"""Semantic failure controls; run after a complete fresh portable replay."""
import argparse,copy,hashlib,json,math,shutil,sys,tempfile
from fractions import Fraction
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
CERT=ROOT/'certificate'
sys.path.insert(0,str(CERT))
import verify
import foundation
from exact_tables import charpoly_psd,block_functionals,require
from flint import fmpq_mat
from resource_report import peak_rss_report

def rejected(label, expected, operation, reports):
    try: operation()
    except (AssertionError,ValueError,TypeError) as error:
        require(expected in str(error),label+': unexpected failure '+str(error))
        reports.append({'control':label,'rejected':True,'reason':str(error)})
        return
    raise RuntimeError(label+': corrupted input was accepted')

def main(args):
    reports=[]
    certificate,manifest=verify.inventory(CERT)
    complete=json.loads((args.replay/'COMPLETE.json').read_text())
    require(complete['status']=='COMPLETE_PORTABLE_EXACT_CERTIFICATE_PASS' and complete['bundle_manifest_sha256']==manifest,'Controls require this version complete replay')
    rejected('negative Gram eigenvalue','Negative Gram eigenvalue',lambda:charpoly_psd(fmpq_mat([[1,0],[0,-1]])),reports)
    archive=json.loads((CERT/'data/local_modules.json').read_text())
    record=certificate['retained_multipliers'][0]
    block=json.loads((CERT/record['path']).read_text())
    block['Z'][0][0]=str(-1-sum(abs(Fraction(v)) for v in block['Z'][0]))
    with np.load(CERT/'data/basis_data.npz',allow_pickle=False) as z:
        indices={f:z['basis'+f+'_indices'] for f in ('31','22','13')}
        altered={k:z[k].copy() for k in z.files}
    rejected('final multiplier diagonal changed','Negative Gram eigenvalue',lambda:block_functionals(block,verify.load_modules(archive,block),indices),reports)
    altered['global_objective_numerator'][0]=int(altered['global_objective_numerator'][0])+1
    coords=foundation.load_coordinates(CERT);maps=coords.ExactMaps()
    rejected('objective coefficient changed','Objective mismatch',lambda:foundation.verify_objective_and_equalities(altered,maps,coords),reports)
    # Match the forged functional hashes too: this tests the final algebraic
    # identity, not merely the artifact checksum gate.
    with tempfile.TemporaryDirectory(prefix='werner-cross-control-') as temp:
        out=Path(temp)
        for p in args.replay.iterdir():
            if p.is_file() and p.name!='COMPLETE.json':shutil.copyfile(p,out/p.name)
        changed=0
        for i,rec in enumerate(certificate['retained_multipliers']):
            if rec['kind']=='source':continue
            base=out/('block_'+str(i).zfill(4));p=base.with_suffix('.npz')
            with np.load(p,allow_pickle=False) as z:arrays={k:z[k].copy() for k in z.files}
            arrays['functional22_denominator']=np.asarray(str(2*int(arrays['functional22_denominator'])))
            np.savez_compressed(p,**arrays)
            receipt=json.loads(base.with_suffix('.json').read_text())
            receipt['functional_sha256']=verify.sha(p)
            base.with_suffix('.json').write_text(json.dumps(receipt)+'\n');changed+=1
        require(changed==64,'All coupled cross terms altered')
        rejected('coupled cross factor two omitted','Exact global functional identity failed',lambda:verify.finalize(CERT,out,certificate,manifest),reports)
    # Three canonical exact tensor-product boundary examples.
    def product_form(factors):
        value=Fraction(1)
        for m in factors:
            m=np.asarray(m,dtype=object)
            value*=sum(Fraction(int(x*x)) for x in m.ravel())-Fraction(1,2)*int(np.trace(m))**2
        return value
    e00=np.diag([1,0,0]);e01=np.zeros((3,3),dtype=int);e01[0,1]=1
    boundary=[product_form([np.diag([1,1,0]),e00,e00]),product_form([np.diag([1,1,0]),e01,e00]),product_form([np.eye(3,dtype=int),e00,e00])]
    require(boundary==[0,0,Fraction(-3,8)],'Boundary and rank-three negative examples')
    require(peak_rss_report('linux',123)['peak_rss_bytes']==125952,'Linux KiB conversion')
    require(peak_rss_report('darwin',123)['peak_rss_bytes']==123,'macOS byte conversion')
    require(peak_rss_report('unknown',123)['peak_rss_bytes'] is None,'Unknown resource units remain unknown')
    result={'status':'SEMANTIC_NEGATIVE_CONTROLS_PASS','optimized_python':not __debug__,'controls':reports,'boundary_values':list(map(str,boundary)),'resource_unit_checks':'linux, darwin, unknown passed','bundle_manifest_sha256':manifest}
    if args.output:
        require(not args.output.exists(),'Do not overwrite control receipt')
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--replay',type=Path,required=True);p.add_argument('--output',type=Path)
    main(p.parse_args())
