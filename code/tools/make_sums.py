# make_sums.py - write SHA256SUMS.txt for a repo folder (skips .git). Usage: python make_sums.py <folder>
import hashlib, os, sys
R = sys.argv[1]; os.chdir(R); BS = chr(92); lines = []
for root, _, fs in os.walk('.'):
    if '.git' in root.replace(BS, '/').split('/'): continue
    for f in sorted(fs):
        p = os.path.join(root, f).replace(BS, '/')[2:]
        if p == 'SHA256SUMS.txt' or p.endswith('.pyc'): continue
        lines.append(hashlib.sha256(open(p, 'rb').read()).hexdigest() + '  ' + p)
lines.sort(key=lambda x: x[66:])
open('SHA256SUMS.txt', 'w', encoding='utf-8', newline='\n').write('\n'.join(lines) + '\n'); print(len(lines), 'files hashed')
