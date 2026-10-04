"""Verify shipped bytes; this is not a mathematical proof checker."""
from pathlib import Path
import hashlib,json
root=Path(__file__).resolve().parent
expected={}
for line in (root/'MANIFEST.sha256').read_text().splitlines():
    digest,name=line.split('  ',1)
    path=(root/name).resolve()
    if not path.is_relative_to(root):raise ValueError('Unsafe manifest path')
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=digest:
        raise ValueError('Manifest mismatch: '+name)
    expected[name]=digest
actual={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file()
        and '.git' not in p.relative_to(root).parts and '__pycache__' not in p.parts
        and p.name!='MANIFEST.sha256'}
if actual!=set(expected):raise ValueError('Manifest coverage mismatch')
print(json.dumps({'status':'PACKAGE_MANIFEST_PASS','files':len(expected)}))
