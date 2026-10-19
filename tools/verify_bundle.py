"""Offline checks only: never import project scripts or connect to Zemax."""
import ast
import hashlib
import json
from pathlib import Path
import re
from urllib.parse import unquote
import zipfile

ROOT = Path(__file__).resolve().parents[1]
manifest = json.loads((ROOT / 'FILE_MANIFEST.json').read_text(encoding='utf-8'))
for item in manifest['files']:
    path = ROOT / item['path']
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != item['sha256']:
        raise RuntimeError(f'Missing or changed archive file: {item["path"]}')
counts = {'python': 0, 'json': 0, 'models': 0}
for item in manifest['files']:
    path = ROOT / item['path']
    if path.suffix == '.py':
        ast.parse(path.read_text(encoding='utf-8-sig'), filename=item['path'])
        counts['python'] += 1
    elif path.suffix == '.json':
        json.loads(path.read_text(encoding='utf-8-sig'))
        counts['json'] += 1
    elif path.suffix.lower() == '.zos':
        counts['models'] += 1
for doc in [ROOT / 'README.md', *sorted((ROOT / 'docs').glob('*.md'))]:
    for target in re.findall(r'\]\(([^)]+)\)', doc.read_text(encoding='utf-8')):
        if '://' in target or target.startswith('#'):
            continue
        path = (doc.parent / unquote(target.split('#')[0])).resolve()
        # The simulation ZIP itself does not recursively contain its own package.
        if target.endswith('codex-4-zemax-simulation-bundle.zip') and not (ROOT / 'packages').exists():
            continue
        if not path.is_relative_to(ROOT) or not path.exists():
            raise RuntimeError(f'Broken local link: {doc.name}: {target}')
package = ROOT / 'packages' / 'codex-4-zemax-simulation-bundle.zip'
if package.exists():
    expected = (package.parent / 'SHA256SUMS.txt').read_text().split()[0]
    if hashlib.sha256(package.read_bytes()).hexdigest() != expected:
        raise RuntimeError('ZIP SHA-256 mismatch')
    with zipfile.ZipFile(package) as bundle:
        if bundle.testzip() is not None:
            raise RuntimeError('ZIP integrity error')
        for item in manifest['files']:
            data = bundle.read('codex-4-zemax/' + item['path'])
            if hashlib.sha256(data).hexdigest() != item['sha256']:
                raise RuntimeError(f'ZIP payload mismatch: {item["path"]}')
print(json.dumps({'status': 'passed', 'files': len(manifest['files']), **counts}, ensure_ascii=False))
