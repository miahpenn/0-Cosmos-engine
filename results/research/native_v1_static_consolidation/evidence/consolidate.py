from pathlib import Path
import ast, re, textwrap, zipfile, hashlib, shutil
root=Path(__file__).parent; w=root/'work'; p=root/'candidate/engine/native_v1';p.mkdir(parents=True,exist_ok=True)
def read(rel):return (w/rel).read_text()
def write(name,s): (p/name).write_text(s.rstrip()+'\n')
def before_main(s):return s.split("if __name__=='__main__':",1)[0]
def replace(s,a,b):
 assert s.count(a)==1,(a,s.count(a));return s.replace(a,b,1)
def defsrc(s,name,cls=None):
 t=ast.parse(s);ns=t.body if cls is None else next(n.body for n in t.body if isinstance(n,ast.ClassDef) and n.name==cls)
 n=next(n for n in ns if isinstance(n,ast.FunctionDef) and n.name==name)
 return textwrap.dedent('\n'.join(s.splitlines()[n.lineno-1:n.end_lineno]))
base='nrc_stage_a11o_a12a_moments/'
for name in ['assembly','geometry_core','native_cmc','matter_source','moments','transport']:
 s=before_main(read(base+name+'.py'));s=re.sub(r'^from (assembly|geometry_core|native_cmc|matter_source|moments|transport) import ',r'from .\1 import ',s,flags=re.M);write(name+'.py',s)
for name in ['basis.json','moments.json']:shutil.copyfile(w/'a12_structural_candidate'/name,p/name)
shutil.copyfile(w/base/"A8_certificate.json",p/"A8_certificate.json")
s=before_main(read('a12_structural_candidate/probe.py'));s=replace(s,'src=pathlib.Path(sys.argv[1]); sys.path.insert(0,str(src)); import geometry_core as gc, native_cmc as nc\nsys.path.insert(0,str(pathlib.Path(__file__).resolve().parent)); import check','from . import geometry_core as gc, native_cmc as nc');write('candidate.py',s)
s=before_main(read('a12_shared_assembly/check.py'));s=replace(s,'src,struct=map(pathlib.Path,sys.argv[1:3]);sys.path.insert(0,str(struct))\nfrom probe import Candidate\nsys.path.insert(0,str(src))','from .candidate import Candidate');s=re.sub(r'^from (assembly|geometry_core|native_cmc|matter_source|moments|transport) import ',r'from .\1 import ',s,flags=re.M);write('shared.py',s)
s=before_main(read('a12_common_residual/common.py'));s=replace(s,'sys.path.insert(0,sys.argv[2]); from probe import Candidate,check','from .candidate import Candidate');write('common.py',s)
c=before_main(read('a12_composed_rate/check.py'));i=before_main(read('a12_structural_scalar_integrated/check.py'));r=before_main(read('a12_rate_interface/check.py'))
# Materialize the exact historical replacements offline, using the literals from the frozen sources.
# No evaluation of the archived replacement scripts or runtime exec.
def literal_replacements(src,fn):
 for n in ast.parse(src).body:
  if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and isinstance(n.value.func,ast.Name) and n.value.func.id==fn:
   yield tuple(ast.literal_eval(x) for x in n.value.args)
split=defsrc(read(base+'geometry_core.py'),'rhs','Assembly');split=replace(split,'def rhs(','def split_rhs(')
vals={n.targets[0].id:ast.literal_eval(n.value) for n in ast.parse(c).body if isinstance(n,ast.Assign) and isinstance(n.targets[0],ast.Name) and n.targets[0].id in ['old','new']}
split=replace(split,vals['old'],vals['new']);split=replace(split,'return dict(p=pt,','diag.update(joint_load_difference=float(np.linalg.norm(joint_delta)),joint_work=float(p@joint_delta),joint_adjusted_energy_closure=float(total-pred-p@joint_delta));return dict(p=pt,')
with (p/'geometry_core.py').open('a') as f:f.write('\n'+split+'\n')
start=c.index('from geometry_core import Assembly');end=c.index('# Preserve all literal nonlinear terms');imports=c[start:end];imports=re.sub(r'^from (assembly|geometry_core|matter_source|moments|transport) import ',r'from .\1 import ',imports,flags=re.M).replace('import geometry_core as gc','from . import geometry_core as gc')
head='"""Static native research rate, restricted V1 witness; not a production backend."""\nimport hashlib\nimport numpy as np, scipy.linalg as la\nfrom numpy.polynomial.legendre import legder,legval\nfrom . import shared as sh, common as cm\n'+imports
classes=c[c.index('class Integrated'):c.index('def run(n):')]
evaluate=defsrc(c,'run')
for a,b in [*literal_replacements(i,'replace'),*literal_replacements(r,'replace')]:evaluate=replace(evaluate,a,b)
# Imports remain local to avoid a circular import at module import time.
evaluate=replace(evaluate,'def evaluate(A,S,raw,sc,mg,g_boundary,scalar_boundary):','def evaluate(A,S,raw,sc,mg,g_boundary,scalar_boundary):\n from .scalar import StructuralRegistry, weak_port')
write('composed.py',head+classes+'\n'+evaluate)
weak=defsrc(i,'weak');weak=replace(weak,'def weak(A,S,raw,sc,mg,reg,al):','def weak_port(A,S,raw,sc,mg,reg,al,scalar_boundary):')
for a,b in [('rf=(Aa@Pc.T).T;rp=','physical_flux=flux.copy();flux=physical_flux if scalar_boundary is None else np.asarray(scalar_boundary,dtype=float);boundary_delta=flux-physical_flux\n rf=(Aa@Pc.T).T;rp='),('wave_rem=rp-strong_load-jump_load','wave_rem=rp-strong_load-jump_load-boundary_delta[:,None]*e'),('boundary_flux=flux.tolist(),','boundary_flux=flux.tolist(),physical_boundary_flux=physical_flux.tolist(),boundary_flux_delta=boundary_delta.tolist(),')]:weak=replace(weak,a,b)
write('scalar.py','"""Factored scalar derivative and explicit wave port."""\nimport numpy as np, scipy.linalg as la\nfrom numpy.polynomial.legendre import leggauss\nfrom . import composed as C\nKAP=C.KAP; mx=C.mx\n'+i[i.index('class ScalarSpace'):i.index('def weak(')]+'\n'+weak)
write('rate.py','"""Immutable rate interface, no temporal stages."""\nimport hashlib\nimport numpy as np\nfrom . import composed as C, scalar as I\n'+r[r.index('def snapshot('):])
d=read('a12_v1_execution_draft/driver.py');d=d[d.index('from matter_source import'):d.index('def advance(')];d=d.replace('from matter_source import','from .matter_source import').replace('from transport import','from .transport import');write('witness.py','"""Homogeneous initial witness and reference RHS only; no evolution entry point."""\nimport numpy as np\nfrom . import rate as R\nC,I=R.C,R.I\n'+d)
with zipfile.ZipFile(next((root/'locked_release').glob('*.zip'))) as z:g=z.read('a12_v1_release_locked/guarded.py').decode()
# Preserve preflight, attempt and ledger functions; do not expose a temporal driver.
write('guards.py','"""Frozen release preflight and diagnostics, no time integrator."""\nimport pathlib, json, os, traceback\nimport numpy as np\nfrom numpy.polynomial.legendre import leggauss\n'+g[g.index('H_LIMIT='):g.index('USED=')]+'\n'+g[g.index('class Stop'):g.index('def sha(')]+'\n'+g[g.index('def clean('):g.index('def validate_reference_solution(')])
write('__init__.py','"""Opt-in V1 research candidate. Not connected to production. No evolution driver."""')
for f in p.glob('*.py'):
 t=ast.parse(f.read_text());assert not re.search(r'\b(?:exec\(|inspect\.getsource|sys\.argv|sys\.path|spec_from_file_location)',f.read_text()),f
print('Static modules:',len(list(p.glob('*.py'))),'data files:',len(list(p.glob('*.json'))))
