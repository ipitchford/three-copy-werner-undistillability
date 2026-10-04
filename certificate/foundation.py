#!/usr/bin/env python3
"""Replay the exact foundations of the portable Werner certificate.

This verifier uses only the adjacent exact_coordinates.py and exported data.
It does not load floating models, historical checkpoints, or solver output.
Full replay is a research job; the importable API does not launch a job.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path
import resource
from resource_report import peak_rss_report
import time

import numpy as np
import sympy as sp
from flint import nmod_mat


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def load_coordinates(bundle):
    path = bundle / 'exact_coordinates.py'
    spec = importlib.util.spec_from_file_location('portable_werner_exact_coordinates', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def rational_matrix(rows):
    return sp.Matrix([[sp.Rational(value) for value in row] for row in rows])


def literal_diagrams(d, cut, permutations):
    """V_p sends input tensor factor i to output factor p[i]; then transpose cut."""
    tuples = tuple(itertools.product(range(d), repeat=4))
    index = {value: i for i, value in enumerate(tuples)}
    matrices = []
    for p in permutations:
        entries = {}
        for column, source in enumerate(tuples):
            target = [0] * 4
            for i in range(4):
                target[p[i]] = source[i]
            row_digits, column_digits = target, list(source)
            for i in cut:
                row_digits[i], column_digits[i] = column_digits[i], row_digits[i]
            key = (index[tuple(row_digits)], index[tuple(column_digits)])
            require(key not in entries, 'Partial transpose did not preserve matrix entries')
            entries[key] = 1
        matrices.append(sp.SparseMatrix(d**4, d**4, entries))
    return matrices


def literal_partial_transpose(matrix, d, cut):
    words = tuple(itertools.product(range(d), repeat=4))
    index = {word: i for i, word in enumerate(words)}
    entries = {}
    for (row, column), value in matrix.todok().items():
        left, right = list(words[row]), list(words[column])
        for slot in cut:
            left[slot], right[slot] = right[slot], left[slot]
        entries[index[tuple(left)], index[tuple(right)]] = value
    return sp.SparseMatrix(d**4, d**4, entries)


def verify_joint_diagram_maps(coordinates):
    """Calibrate both ++ maps from literal tensor indices, not module products."""
    permutations = coordinates.PERMS
    tau, p0 = (0, 1, 3, 2), (2, 3, 1, 0)
    compose = lambda p, q: tuple(p[q[i]] for i in range(4))
    reports = []
    for d in (2, 3):
        ordinary = literal_diagrams(d, (), permutations)
        transposed = literal_diagrams(d, (3,), permutations)
        vtau, vp0 = ordinary[permutations.index(tau)], ordinary[permutations.index(p0)]
        words = tuple(itertools.product(range(d), repeat=4))
        index = {word: i for i, word in enumerate(words)}
        for pi, p in enumerate(permutations):
            diagonal_index = permutations.index(compose(compose(tau, p), tau))
            cross_index = permutations.index(compose(p, p0))
            diagonal = literal_partial_transpose(vtau * ordinary[pi] * vtau, d, (3,))
            cross_product = literal_partial_transpose(ordinary[pi] * vp0, d, (3,))
            # R[i0 i1 i2 i3; j0 j1 j2 j3]
            #   = Y22[i0 i1 i2 j3; i3 j2 j0 j1].
            entries = {}
            for (r, c), value in ordinary[pi].todok().items():
                u, v = words[r], words[c]
                entries[index[(u[0], u[1], u[2], v[0])],
                        index[(v[2], v[3], v[1], u[3])]] = value
            cross_indices = sp.SparseMatrix(d**4, d**4, entries)
            require(diagonal == transposed[diagonal_index], 'Literal joint diagonal map')
            require(cross_product == cross_indices == transposed[cross_index], 'Literal joint cross map')
        reports.append(dict(d=d, diagonal_relations=24, cross_relations=24,
                            direct_cross_index_relations=24))
    return dict(status='EXACT_LITERAL_JOINT_DIAGRAM_MAPS_PASS', tau=tau, p0=p0, reports=reports)


def verify_gaussian_joint_gram():
    """Independent exact 32 by32 ++ Gram calibration using Gaussian integers."""
    words = tuple(itertools.product(range(2), repeat=4))
    index = {word: i for i, word in enumerate(words)}
    mul = lambda x, y: (x[0]*y[0] - x[1]*y[1], x[0]*y[1] + x[1]*y[0])
    conjugate = lambda x: (x[0], -x[1])
    a, b = [(1, 1), (2, -1)], [(2, 3), (-1, 2)]
    def tensor(rays):
        result = []
        for word in words:
            value = (1, 0)
            for ray, i in zip(rays, word):
                value = mul(value, ray[i])
            result.append(value)
        return result
    def outer(left, right=None):
        right = left if right is None else right
        return [[mul(x, conjugate(y)) for y in right] for x in left]
    def pt(matrix):
        return [[matrix[index[i[:3] + (j[3],)]][index[j[:3] + (i[3],)]]
                 for j in words] for i in words]
    def conjugate_swap23(matrix):
        sources = [index[(i[0], i[1], i[3], i[2])] for i in words]
        return [[matrix[i][j] for j in sources] for i in sources]
    y31 = outer(tensor([a, a, a, b]))
    y13 = outer(tensor([b, b, b, a]))
    y22 = outer(tensor([a, a, b, b]))
    first, second = pt(conjugate_swap23(y31)), pt(conjugate_swap23(y13))
    destinations = [index[(s[3], s[2], s[0], s[1])] for s in words]
    cross = pt([[row[j] for j in destinations] for row in y22])
    z1 = tensor([a, a, b, list(map(conjugate, a))])
    z2 = tensor([b, b, a, list(map(conjugate, b))])
    gram = [first[i] + cross[i] for i in range(16)]
    gram += [[conjugate(cross[j][i]) for j in range(16)] + second[i] for i in range(16)]
    require(gram == outer(z1 + z2), 'Gaussian-integer ++ joint Gram identity')
    return dict(status='EXACT_GAUSSIAN_INTEGER_JOINT_GRAM_PASS', order=32,
                a=a, b=b, arithmetic='Pairs of Python integers',
                scope='A finite calibration of the tensor formula; the general formula is established algebraically in the proof.')


def verify_local_modules(module_path, coordinates):
    archive = json.loads(module_path.read_text())
    permutations = tuple(tuple(int(x) for x in p) for p in archive['permutations'])
    expected = tuple(itertools.permutations(range(4)))
    require(permutations == expected == coordinates.PERMS, 'Permutation ordering differs')
    inverse = tuple(permutations.index(tuple(p.index(i) for i in range(4))) for p in permutations)
    seen_cases, reports = set(), []
    for case in archive['representations']:
        d, cut = int(case['d']), tuple(int(x) for x in case['cut'])
        require(d in (2, 3) and len(set(cut)) == len(cut) and all(0 <= x < 4 for x in cut),
                'Invalid local dimension or partial transpose cut')
        require((d, cut) not in seen_cases, 'Duplicate local case')
        seen_cases.add((d, cut))
        diagrams = literal_diagrams(d, cut, permutations)
        labels = set()
        for record in case['modules']:
            label, order = record['label'], int(record['dim'])
            require(label not in labels, 'Duplicate module label')
            labels.add(label)
            embedding = rational_matrix(record['embedding'])
            gram = rational_matrix(record['gram'])
            reps = [rational_matrix(matrix) for matrix in record['basis']]
            require(order > 0 and embedding.shape == (d**4, order), 'Embedding shape')
            require(gram.shape == (order, order) and len(reps) == 24, 'Module matrix shapes')
            require(all(rep.shape == (order, order) for rep in reps), 'Representation shape')
            require(gram == embedding.T * embedding, 'Gram is not the literal embedding Gram')
            pivots = [gram[:k, :k].det() for k in range(1, order + 1)]
            require(all(value > 0 for value in pivots), 'Embedding Gram is not positive definite')
            for p in range(24):
                require(diagrams[p] * embedding == embedding * reps[p],
                        f'Local intertwining failed: d={d}, cut={cut}, label={label}, p={p}')
                require(reps[p].T * gram == gram * reps[inverse[p]],
                        f'Local metric adjoint failed: d={d}, cut={cut}, label={label}, p={p}')
            reports.append(dict(d=d, cut=list(cut), label=label, dimension=order,
                                embedding_relations=24, metric_adjoint_relations=24,
                                positive_leading_minors=[str(x) for x in pivots]))
    require(bool(reports), 'No local modules')
    return dict(status='EXACT_LOCAL_EMBEDDINGS_AND_METRICS_PASS', modules=len(reports),
                cases=len(seen_cases), reports=reports,
                completeness_claim='No irreducible multiplicity completeness is used; each module is a genuine invariant compression.')


def integer_source_column(maps, word, family, coordinates):
    """Independent copy of the exact integer generator map, no old helper imports.

    Output coordinates are normalized physical orbit averages. Canonical sorting
    collects the total coefficient of each average, so no orbit-size multiplier
    appears here. Multiplication of the qubit quotient by 2 makes it integral.
    """
    terms, denominator = coordinates.projected_generator(word, family)
    quotients = []
    for d, multiplier in ((3, 1), (2, 2)):
        q = maps.local[d]['quotient'] * multiplier
        require(all(value.q == 1 for value in q), 'Expected integral scaled local quotient')
        quotients.append([[(i, int(q[i, j])) for i in range(q.rows) if q[i, j]]
                          for j in range(24)])
    q3, q2 = quotients
    result = Counter()
    for term, coefficient in terms.items():
        for choice in itertools.product(q3[term[0]], q3[term[1]], q3[term[2]], q2[term[3]]):
            key = tuple(sorted(item[0] for item in choice[:3])) + (choice[3][0],)
            result[key] += int(coefficient) * math.prod(item[1] for item in choice)
    return {key: value for key, value in result.items() if value}, 2 * denominator


def verify_basis_minors(data, maps, coordinates):
    prime = int(data['prime'])
    require(sp.isprime(prime), 'Modulus is not prime')
    catalog = coordinates.CATALOG
    lookup = {tuple(coordinates.PHYSICAL_KEEP.index(p) for p in word[:3]) +
              (coordinates.AUXILIARY_KEEP.index(word[3]),): i for i, word in enumerate(catalog)}
    require(len(lookup) == len(catalog) == 32200, 'Generator catalog cardinality')
    reports = {}
    for family, size in (('31', 577), ('22', 1220)):
        indices = np.asarray(data['basis' + family + '_indices'])
        buckets = np.asarray(data['basis' + family + '_buckets'])
        signs = np.asarray(data['basis' + family + '_signs'])
        selected_rows = np.asarray(data['basis' + family + '_selected_rows'])
        denominator = int(data['basis' + family + '_denominator'])
        require(indices.shape == (size,) and len(set(map(int, indices))) == size, 'Basis indices')
        require(all(0 <= int(i) < len(catalog) for i in indices), 'Basis catalog index')
        require(buckets.shape == signs.shape == (len(catalog),), 'Sketch arrays')
        require(all(int(s) in (-1, 1) for s in signs), 'Sketch signs')
        require(all(int(b) >= 0 for b in buckets), 'Negative sketch bucket')
        require(selected_rows.shape == (size,) and len(set(map(int, selected_rows))) == size,
                'Selected sketch rows')
        require(all(int(row) >= 0 for row in selected_rows), 'Negative selected row')
        selected = {int(row): i for i, row in enumerate(selected_rows)}
        minor = [[0 for _ in range(size)] for _ in range(size)]
        for j, index in enumerate(indices):
            column, actual_denominator = integer_source_column(maps, catalog[int(index)], family, coordinates)
            require(actual_denominator == denominator, 'Generator denominator differs')
            for key, value in column.items():
                source_row = lookup[key]
                target = selected.get(int(buckets[source_row]))
                if target is not None:
                    minor[target][j] += int(signs[source_row]) * value
        determinant = int(nmod_mat(minor, prime).det())
        require(determinant != 0, f'{family} basis minor is singular modulo {prime}')
        reports[family] = dict(size=size, prime=prime, nonzero_determinant=determinant,
                               integer_column_denominator=denominator)
    require(np.array_equal(data['basis13_indices'], data['basis31_indices']),
            'Family13 must use the same bosonic generator basis as31')
    reports['13'] = dict(reports['31'], same_projector_as31=True)
    return reports


def selected_equality_column(record, family, selected_rows):
    """2912 raw equality rows: three Stiefel groups, AAB31-22, ABB22-13."""
    values = {}
    family_index = {'31': 0, '22': 1, '13': 2}[family]
    for index, value in record['stiefel'].items():
        values[784 * family_index + int(index)] = Fraction(value)
    if family in ('31', '22'):
        sign = 1 if family == '31' else -1
        for index, value in record['aab'].items():
            values[2352 + int(index)] = sign * Fraction(value)
    if family in ('22', '13'):
        sign = 1 if family == '22' else -1
        for index, value in record['abb'].items():
            values[2632 + int(index)] = sign * Fraction(value)
    return [values.get(int(row), Fraction(0)) for row in selected_rows]


def verify_objective_and_equalities(data, maps, coordinates):
    selected_rows = np.asarray(data['equality_basis_rows'])
    row_basis = np.asarray(data['row_basis_integer'])
    matrix_denominator = int(data['matrix_denominator'])
    numerator = np.asarray(data['global_objective_numerator'])
    denominator = int(data['global_objective_denominator'])
    require(selected_rows.ndim == 1 and len(set(map(int, selected_rows))) == len(selected_rows),
            'Equality row indices')
    require(all(0 <= int(row) < 2912 for row in selected_rows), 'Equality row outside raw map')
    require(row_basis.shape == (len(selected_rows), 2374), 'Equality basis shape')
    require(numerator.shape == (2374,) and denominator > 0 and matrix_denominator > 0,
            'Objective/equality rational encodings')
    column = 0
    for family, size in (('31', 577), ('22', 1220), ('13', 577)):
        indices = data['basis' + family + '_indices']
        require(indices.shape == (size,), 'Basis width')
        for index in indices:
            record = maps.evaluate(coordinates.CATALOG[int(index)], family)
            objective = Fraction(record['objective']) if family == '31' else Fraction(0)
            require(objective == Fraction(int(numerator[column]), denominator),
                    f'Objective mismatch in column {column}')
            equalities = selected_equality_column(record, family, selected_rows)
            for row, value in enumerate(equalities):
                require(value == Fraction(int(row_basis[row, column]), matrix_denominator),
                        f'Equality mismatch at row {row}, column {column}')
            column += 1
    require(column == 2374, 'Incomplete source coordinate coverage')
    return dict(status='EXACT_OBJECTIVE_AND_SELECTED_EQUALITIES_PASS',
                columns=column, equality_rows=len(selected_rows), raw_equality_rows=2912,
                objective_normalization='X=54^4 Y; objective includes 1/(8*54^4)',
                equality_rank_claim='No equality rank or completeness assertion is needed: all selected rows are verified valid homogeneous constraints.')


def verify_foundation(bundle, output_path=None):
    bundle = Path(bundle).resolve()
    if output_path is not None:
        output_path = Path(output_path)
        require(not output_path.exists(), 'Refusing to overwrite a foundation receipt')
    started = time.monotonic()
    coordinate_path = bundle / 'exact_coordinates.py'
    module_path = bundle / 'data' / 'local_modules.json'
    basis_path = bundle / 'data' / 'basis_data.npz'
    inputs = {str(path.relative_to(bundle)): sha256(path)
              for path in (coordinate_path, module_path, basis_path)}
    coordinates = load_coordinates(bundle)
    local = verify_local_modules(module_path, coordinates)
    joint_maps = verify_joint_diagram_maps(coordinates)
    gaussian_gram = verify_gaussian_joint_gram()
    print(json.dumps(dict(phase='foundation_local', modules=local['modules'], status='PASS')), flush=True)
    maps = coordinates.ExactMaps()
    dimensions = coordinates.exact_projector_trace_dimensions(maps)
    require({family: value['exact_dimension'] for family, value in dimensions.items()} ==
            {'31': 577, '22': 1220, '13': 577}, 'Exact source dimensions')
    with np.load(basis_path, allow_pickle=False) as data:
        minors = verify_basis_minors(data, maps, coordinates)
        print(json.dumps(dict(phase='foundation_basis', status='PASS', minors=minors)), flush=True)
        features = verify_objective_and_equalities(data, maps, coordinates)
    require(all(sha256(bundle / path) == digest for path, digest in inputs.items()),
            'Foundation inputs changed during replay')
    receipt = dict(status='PORTABLE_EXACT_FOUNDATION_PASS', input_sha256=inputs,
                   source_sha256=sha256(__file__), local=local,
                   joint_diagram_maps=joint_maps, gaussian_joint_gram=gaussian_gram,
                   exact_projector_dimensions=dimensions, basis_minors=minors,
                   objective_and_equalities=features,
                   elapsed_seconds=time.monotonic() - started,
                   **peak_rss_report(),
                   global_certificate_identity_checked=False,
                   multiplier_PSD_checked=False)
    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open('x') as stream:
            json.dump(receipt, stream, indent=2, allow_nan=False)
            stream.write('\n')
    return receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = verify_foundation(args.bundle, args.output)
    print(json.dumps({key: result[key] for key in ('status', 'elapsed_seconds', 'peak_rss_bytes')}), flush=True)
