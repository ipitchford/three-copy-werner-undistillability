#!/usr/bin/env python3
"""Independently replay the portable rational Werner certificate.

Receipts are written outside the immutable input directory. Bounded block
chunks may be resumed, or split into disjoint --start/--stop intervals.
Only phase finalize can report a complete certificate PASS.
"""
import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import resource
from resource_report import peak_rss_report
import signal
import sys
import tempfile
import os
import time
import numpy as np
from flint import fmpz_mat
from foundation import verify_foundation
from exact_tables import block_functionals,clear_vector,require

if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
WIDTHS={'31':577,'22':1220,'13':577}
OFFSETS={'31':0,'22':577,'13':1797}

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda:stream.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()

def dump_new(path,value):
    require(not path.exists(),'Refuse output overwrite: '+str(path))
    with tempfile.NamedTemporaryFile('w',dir=path.parent,prefix=path.name+'.partial-',delete=False) as stream:
        json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n');temporary=Path(stream.name)
    os.link(temporary,path);temporary.unlink()

def npz_new(path,arrays):
    require(not path.exists(),'Refuse output overwrite: '+str(path))
    with tempfile.NamedTemporaryFile(dir=path.parent,prefix=path.name+'.partial-',suffix='.npz',delete=False) as stream:temporary=Path(stream.name)
    np.savez_compressed(temporary,**arrays);os.link(temporary,path);temporary.unlink()

def inventory(bundle):
    manifest_path=bundle/'manifest.json';manifest=json.loads(manifest_path.read_text())
    for relative,digest in manifest['files'].items():
        path=(bundle/relative).resolve();require(path.is_relative_to(bundle.resolve()),'Manifest path escapes bundle')
        require(sha(path)==digest,'Bundle file hash mismatch: '+relative)
    required={'resource_report.py','verify.py','foundation.py','exact_tables.py','exact_coordinates.py','basis_values.py','certificate.json',
        'data/local_modules.json','data/basis_data.npz','data/equality_multiplier.json'}
    require(required<=set(manifest['files']),'Manifest omits a truth dependency')
    certificate=json.loads((bundle/'certificate.json').read_text());retained=certificate['retained_multipliers'];omitted=certificate['omitted_zero_multipliers']
    require(len(retained)==716 and len(omitted)==741,'Certificate multiplier counts')
    require(certificate['source_widths']==[577,1220,577] and certificate['source_family_order']==['31','22','13'],'Source coordinate convention')
    ids=[r['id'] for r in retained]+[r['id'] for r in omitted]
    expected=[f+'_'+str(i).zfill(4) for f,n in [('31',408),('22',541),('13',408)] for i in range(n)]+['joint_'+str(i).zfill(3) for i in range(100)]
    require(sorted(ids)==sorted(expected),'Complete retained/zero-omission partition')
    originals={f:[] for f in WIDTHS};joint=[]
    for rec in retained:
        require(rec['path'] in manifest['files'] and manifest['files'][rec['path']]==rec['sha256'],'Multiplier manifest binding')
        block=json.loads((bundle/rec['path']).read_text());require(block['id']==rec['id'] and block['kind']==rec['kind'],'Multiplier identity')
        origin=block['origin']
        if block['kind']=='source':originals[block['family']].extend(origin['original_source_blocks'])
        else:joint.append(origin['original_source_block'])
    for rec in omitted:
        require(rec['reason']=='exact zero multiplier','Unrecognized omission')
        if rec['kind']=='source':originals[rec['family']].extend(rec['original_source_blocks'])
        else:joint.append(rec['original_source_block'])
    for f,n in [('31',408),('22',585),('13',408)]:require(sorted(originals[f])==list(range(n)),'Original source-cone coverage '+f)
    require(sorted(joint)==list(range(100)),'Original joint-cone coverage')
    return certificate,sha(manifest_path)

def check_commit(output,index,record,manifest_sha):
    receipt_path=output/('block_'+str(index).zfill(4)+'.json')
    if not receipt_path.exists():return False
    receipt=json.loads(receipt_path.read_text());require(receipt['status']=='EXACT_MULTIPLIER_REPLAY_PASS','Non-PASS checkpoint')
    require(receipt['bundle_manifest_sha256']==manifest_sha and receipt['multiplier_sha256']==record['sha256'] and receipt['id']==record['id'],'Checkpoint input drift')
    require(sha(receipt_path.with_suffix('.npz'))==receipt['functional_sha256'],'Checkpoint functional drift')
    require(receipt['PSD']['PSD'] is True,'Checkpoint lacks exact PSD proof')
    return True

def load_modules(archive,block):
    require(len(block['modules'])==4,'Four local modules')
    result=[]
    for ref,d in zip(block['modules'],[3,3,3,2]):
        ci,mi=int(ref['case']),int(ref['module']);require(0<=ci<len(archive['representations']),'Case index')
        case=archive['representations'][ci];require(0<=mi<len(case['modules']),'Module index')
        module=case['modules'][mi]
        require(case['d']==d==ref['d'] and case['cut']==block['cut'],'Literal module dimensions/cut')
        require(module['label']==ref['label'] and module['dim']==ref['dim'],'Exact module identity')
        result.append(module)
    return result

def phase_foundation(bundle,output,manifest_sha):
    path=output/'foundation.json'
    if path.exists():
        saved=json.loads(path.read_text());require(saved['bundle_manifest_sha256']==manifest_sha and saved['status']=='PORTABLE_EXACT_FOUNDATION_PASS','Foundation receipt drift');return saved
    report=verify_foundation(bundle);report['bundle_manifest_sha256']=manifest_sha;dump_new(path,report);return report

def phase_blocks(bundle,output,certificate,manifest_sha,args,started):
    require((output/'foundation.json').exists(),'Run foundation before multiplier replay')
    phase_foundation(bundle,output,manifest_sha)
    archive=json.loads((bundle/'data/local_modules.json').read_text())
    with np.load(bundle/'data/basis_data.npz',allow_pickle=False) as basis:indices={f:basis['basis'+f+'_indices'] for f in WIDTHS}
    records=certificate['retained_multipliers'];stop=min(args.stop,len(records));new=0;existing=0
    def deadline():
        if time.monotonic()-started>args.max_seconds:raise TimeoutError('Bounded replay chunk')
    for index in range(args.start,stop):
        record=records[index]
        if check_commit(output,index,record,manifest_sha):existing+=1;continue
        if new>=args.max_blocks or time.monotonic()-started>args.max_seconds-12:break
        began=time.monotonic();block=json.loads((bundle/record['path']).read_text());modules=load_modules(archive,block)
        try:proof,vectors=block_functionals(block,modules,indices,deadline)
        except TimeoutError:break
        arrays={}
        for f,values in vectors.items():
            require(len(values)==WIDTHS[f],'Complete basis functional')
            numbers,den=clear_vector(values);arrays['functional'+f+'_numerator']=np.asarray(list(map(str,numbers)))
            arrays['functional'+f+'_denominator']=np.asarray(str(den))
        base=output/('block_'+str(index).zfill(4));array_path=base.with_suffix('.npz')
        if array_path.exists():
            with np.load(array_path,allow_pickle=False) as old:require(set(old.files)==set(arrays) and all(np.array_equal(old[k],v) for k,v in arrays.items()),'Interrupted functional artifact mismatch')
        else:npz_new(array_path,arrays)
        receipt=dict(status='EXACT_MULTIPLIER_REPLAY_PASS',id=record['id'],index=index,bundle_manifest_sha256=manifest_sha,
            multiplier_sha256=record['sha256'],functional_sha256=sha(array_path),PSD=proof,functional_families=list(vectors),
            elapsed_seconds=time.monotonic()-began,**peak_rss_report())
        dump_new(base.with_suffix('.json'),receipt);new+=1
        print(json.dumps(dict(phase='block',index=index,id=record['id'],order=proof['order'],status='PASS',seconds=receipt['elapsed_seconds'])),flush=True)
    return dict(new_blocks=new,validated_existing_blocks=existing,range=[args.start,stop])

def finalize(bundle,output,certificate,manifest_sha):
    foundation=phase_foundation(bundle,output,manifest_sha);numbers=np.zeros(2374,dtype=object);denominator=1;receipts=[]
    for index,record in enumerate(certificate['retained_multipliers']):
        require(check_commit(output,index,record,manifest_sha),'Missing replay block '+str(index))
        path=output/('block_'+str(index).zfill(4)+'.npz')
        with np.load(path,allow_pickle=False) as data:
            expected_families=[json.loads((bundle/record['path']).read_text())['family']] if record['kind']=='source' else ['31','22','13']
            require(set(data.files)=={f'functional{f}_{part}' for f in expected_families for part in ('numerator','denominator')},'Unexpected functional fields')
            for f in expected_families:
                ns=np.asarray([int(x) for x in data['functional'+f+'_numerator']],dtype=object);d=int(data['functional'+f+'_denominator'])
                require(len(ns)==WIDTHS[f] and d>0,'Functional vector/denominator')
                common=math.lcm(denominator,d);numbers*=common//denominator;start=OFFSETS[f]
                numbers[start:start+WIDTHS[f]]+=ns*(common//d);denominator=common
        receipts.append(dict(index=index,id=record['id'],receipt_sha256=sha(path.with_suffix('.json')),functional_sha256=sha(path)))
    with np.load(bundle/'data/basis_data.npz',allow_pickle=False) as data:
        B=data['row_basis_integer'];objective=[int(x) for x in data['global_objective_numerator']];objective_den=int(data['global_objective_denominator']);rows=list(map(int,data['equality_basis_rows']))
    equality=json.loads((bundle/'data/equality_multiplier.json').read_text());require(equality['equality_basis_rows']==rows,'Equality multiplier row convention')
    eta=[Fraction(x) for x in equality['eta']];require(len(eta)==272 and B.shape==(272,2374),'Equality dimensions')
    eta_numbers,eta_den=clear_vector(eta);term=fmpz_mat(B.T.tolist())*fmpz_mat([[x] for x in eta_numbers])
    common=math.lcm(denominator,objective_den,eta_den)
    residual=[objective[i]*(common//objective_den)-int(numbers[i])*(common//denominator)-int(term[i,0])*(common//eta_den) for i in range(2374)]
    require(not any(residual),'Exact global functional identity failed')
    receipt=dict(status='COMPLETE_PORTABLE_EXACT_CERTIFICATE_PASS',bundle_manifest_sha256=manifest_sha,
        exact_functional_equations=2374,nonzero_PSD_multipliers=716,zero_omissions=741,
        source_dimensions=[577,1220,577],literal_local_modules_verified=True,rational_basis_completeness_verified=True,
        objective_and_equality_rows_recomputed=True,all_multiplier_PSD_recomputed=True,
        all_multiplier_functionals_recomputed=True,exact_integer_residual_zero=True,
        identity='c = sum(multiplier functionals) + B.T eta',
        foundation_receipt_sha256=sha(output/'foundation.json'),block_receipts=receipts,
        final_identity_common_denominator_bits=common.bit_length())
    path=output/'COMPLETE.json';dump_new(path,receipt);return receipt

def main(args):
    started=time.monotonic();bundle=Path(__file__).resolve().parent;output=args.output.resolve()
    require(not output.is_relative_to(bundle),'Replay receipts must be outside immutable input bundle')
    certificate,manifest_sha=inventory(bundle);output.mkdir(parents=True,exist_ok=True)
    run=dict(schema='portable-exact-replay-v1',bundle_manifest_sha256=manifest_sha)
    run_path=output/'RUN.json'
    if run_path.exists():require(json.loads(run_path.read_text())==run,'Replay directory belongs to another bundle')
    else:dump_new(run_path,run)
    result={}
    if args.phase in ('foundation','all'):result['foundation']=phase_foundation(bundle,output,manifest_sha)['status']
    if args.phase in ('blocks','all'):result.update(phase_blocks(bundle,output,certificate,manifest_sha,args,started))
    if args.phase=='finalize':result=finalize(bundle,output,certificate,manifest_sha)
    elif args.phase=='all' and all((output/('block_'+str(i).zfill(4)+'.json')).exists() for i in range(716)):
        result=finalize(bundle,output,certificate,manifest_sha)
    report={k:v for k,v in result.items() if k!='block_receipts'}
    report.update(phase=args.phase,elapsed_seconds=time.monotonic()-started,
        **peak_rss_report(),bundle_manifest_sha256=manifest_sha)
    dump_new(output/('chunk_'+str(time.time_ns())+'.json'),report);print(json.dumps(report,allow_nan=False),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--phase',choices=['foundation','blocks','finalize','all'],default='all')
    parser.add_argument('--start',type=int,default=0);parser.add_argument('--stop',type=int,default=716)
    parser.add_argument('--max-seconds',type=int,default=150);parser.add_argument('--max-blocks',type=int,default=100000)
    args=parser.parse_args();require(0<=args.start<=args.stop<=716 and args.max_seconds>12 and args.max_blocks>0,'Invalid bounds')
    def alarm(signum,frame):raise TimeoutError('Portable replay hard deadline; committed blocks are preserved')
    signal.signal(signal.SIGALRM,alarm);signal.alarm(args.max_seconds+8);main(args)
