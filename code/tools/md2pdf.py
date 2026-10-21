"""md2pdf.py - turn a paper in Markdown (with $...$ and $$...$$ math) into a standalone HTML and a PDF.

No Python packages needed. The HTML renders the Markdown with marked.js and the math with KaTeX (both from a CDN),
and Chrome (headless) prints it to PDF. Math is protected from the Markdown parser by placeholders.

Usage:  python md2pdf.py paper.md [out_dir]        -> out_dir/paper.html and out_dir/paper.pdf
Needs:  Google Chrome (path below) and network access for the CDN at print time."""
import html, json, os, subprocess, sys

CHROME = r"C:\Program Files\Google\Chrome\Application\chrome.exe"

PAGE = """<!doctype html>
<html lang="%(lang)s"><head><meta charset="utf-8"><title>%(title)s</title>
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.css">
<script src="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/marked@12.0.2/marked.min.js"></script>
<style>
 @page { size: A4; margin: 18mm 16mm; }
 body { font-family: "Times New Roman", "Noto Serif TC", "PMingLiU", serif; font-size: 10.5pt; line-height: 1.45; color: #111; max-width: 46em; margin: auto; }
 h1 { font-size: 17pt; line-height: 1.25; } h2 { font-size: 13pt; margin-top: 1.4em; } h3 { font-size: 11.5pt; }
 table { border-collapse: collapse; margin: .6em 0; font-size: 9pt; } td, th { border: 1px solid #999; padding: 3px 5px; vertical-align: top; }
 th { background: #eee; } code { font-size: 9pt; } img { max-width: 100%%; }
 blockquote { color: #444; border-left: 3px solid #ccc; margin-left: 0; padding-left: .8em; }
 .katex-display { margin: .5em 0; overflow-x: auto; }
</style></head><body><div id="out"></div>
<script id="src" type="application/json">%(src)s</script>
<script>
const md = JSON.parse(document.getElementById('src').textContent);
const store = [];
// protect math: $$...$$ first, then $...$ (no newline inside inline math)
let t = md.replace(/\\$\\$([\\s\\S]+?)\\$\\$/g, (m, x) => { store.push([x, true]); return '@@M' + (store.length - 1) + '@@'; });
t = t.replace(/\\$([^$\\n]+?)\\$/g, (m, x) => { store.push([x, false]); return '@@M' + (store.length - 1) + '@@'; });
let h = marked.parse(t);
h = h.replace(/@@M(\\d+)@@/g, (m, i) => { const [x, d] = store[+i];
  try { return katex.renderToString(x, {displayMode: d, throwOnError: false}); } catch (e) { return x; } });
document.getElementById('out').innerHTML = h;
document.title = document.title;
</script></body></html>"""

def build(md_path, out_dir=None):
    md_path = os.path.abspath(md_path)
    out_dir = os.path.abspath(out_dir or os.path.dirname(md_path))
    stem = os.path.splitext(os.path.basename(md_path))[0]
    text = open(md_path, encoding="utf-8").read()
    title = next((l[2:].strip() for l in text.splitlines() if l.startswith("# ")), stem)
    lang = "zh-Hant" if stem.endswith("zh-tw") else "en"
    src = json.dumps(text).replace("</", "<\\/")
    page = PAGE % dict(lang=lang, title=html.escape(title), src=src)
    html_path = os.path.join(out_dir, stem + ".html"); pdf_path = os.path.join(out_dir, stem + ".pdf")
    open(html_path, "w", encoding="utf-8").write(page)
    subprocess.run([CHROME, "--headless=new", "--disable-gpu", "--no-pdf-header-footer", "--virtual-time-budget=20000",
                    "--print-to-pdf=" + pdf_path, "file:///" + html_path.replace("\\", "/")], check=True, timeout=180,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    print("html:", html_path); print("pdf :", pdf_path, os.path.getsize(pdf_path), "bytes")

if __name__ == "__main__":
    build(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
