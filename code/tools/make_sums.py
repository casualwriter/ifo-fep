# make_sums.py - write SHA256SUMS.txt for a repo folder (skips .git). Usage: python make_sums.py <folder>
# Text files are hashed with LF line endings (what git stores), so the sums stay valid after a Windows CRLF edit.
import hashlib, os, sys

BINARY = ('.pdf', '.png', '.zip')
CR, LF, BS = chr(13), chr(10), chr(92)

def main(root):
    os.chdir(root)
    lines = []
    for dirpath, _, files in os.walk('.'):
        if '.git' in dirpath.replace(BS, '/').split('/'):
            continue
        for f in sorted(files):
            p = os.path.join(dirpath, f).replace(BS, '/')[2:]
            if p == 'SHA256SUMS.txt' or p.endswith('.pyc'):
                continue
            data = open(p, 'rb').read()
            if not p.endswith(BINARY):
                data = data.replace((CR + LF).encode(), LF.encode())
            lines.append(hashlib.sha256(data).hexdigest() + '  ' + p)
    lines.sort(key=lambda x: x[66:])
    open('SHA256SUMS.txt', 'w', encoding='utf-8', newline=LF).write(LF.join(lines) + LF)
    print(len(lines), 'files hashed')

if __name__ == '__main__':
    main(sys.argv[1])
