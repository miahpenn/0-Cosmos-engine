"""Zero-time comparison only; takes a verified frozen checkpoint ZIP. Run from package root."""
import sys, os, hashlib, json, zipfile, pathlib, io, tempfile, platform
import numpy as np, scipy
from engine.native_v1 import witness as D, guards as G
assert os.environ.get('OPENBLAS_NUM_THREADS')==os.environ.get('OMP_NUM_THREADS')=='1'
assert (platform.python_version(),np.__version__,scipy.__version__)==('3.13.5','2.3.5','1.17.0')
archive=pathlib.Path(sys.argv[1]).read_bytes()
assert hashlib.sha256(archive).hexdigest()=='d83c71a7a4d9b50a65b4682af9f906145370ec35decb4f335afefe991eba446c'
with zipfile.ZipFile(io.BytesIO(archive)) as checkpoint, zipfile.ZipFile(io.BytesIO(checkpoint.read('locked_release/a12_v1_release_locked.zip'))) as z, tempfile.TemporaryDirectory() as t:
    for n in (8,16):
        A,x,g,b=D.initial(n); f=pathlib.Path(t)/f'N{n}.json'
        G.attempt(D,A,x,g,0.,b,f,'dry_initial')
        current=json.loads(f.read_text()); old=json.loads(z.read(f'a12_v1_release_locked/dry_archive/N{n}/initial.json'))
        assert current==old, f'N{n}: frozen initial record differs; no tolerance applied'
        print(f'N{n}: exact JSON value match, zero temporal steps')
