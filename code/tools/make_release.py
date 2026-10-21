"""make_release.py - make a clean release copy of a paper for Zenodo (rules: ifo/CLAUDE.md section 5).

- delete '---' section separators (keep them inside code fences and tables);
- delete the internal checklist block ('> **Internal checklist ...' and the '>' lines after it);
- end each header line before the first '##' with two spaces, so Markdown keeps the line breaks;
- write <stem>-release.md beside the source unless an output path is given. The source file is not changed.

Usage: python make_release.py paper.md [out.md]"""
import os, sys

def clean(text):
    lines = text.split("\n")
    out, in_code, in_check, in_header = [], False, False, True
    for i, ln in enumerate(lines):
        s = ln.strip()
        if s.startswith("```"):
            in_code = not in_code
        if s.startswith("## "):
            in_header = False
        if not in_code and s.startswith("> **Internal checklist"):
            in_check = True; continue
        if in_check:
            if s.startswith(">"): continue
            in_check = False
        if not in_code and s == "---":
            prev = lines[i - 1].strip() if i else ""
            if not prev.startswith("|"):
                continue
        if in_header and s.startswith("**") and not ln.endswith("  "):
            ln = ln + "  "
        out.append(ln)
    res, blank = [], 0                       # squeeze runs of blank lines left by the deletions
    for ln in out:
        blank = blank + 1 if ln.strip() == "" else 0
        if blank <= 1: res.append(ln)
    return "\n".join(res)

if __name__ == "__main__":
    src = sys.argv[1]
    dst = sys.argv[2] if len(sys.argv) > 2 else os.path.splitext(src)[0] + "-release.md"
    t = clean(open(src, encoding="utf-8").read())
    open(dst, "w", encoding="utf-8", newline="\n").write(t)
    print(dst, len(t))
