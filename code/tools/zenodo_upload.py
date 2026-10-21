"""zenodo_upload.py - create a Zenodo draft, upload the files and fill the metadata. Publishing is a separate, explicit step.

The access token comes from the environment variable ZENODO_API_KEY (never written to disk, never printed).
By default the script talks to the TEST site (sandbox.zenodo.org, which needs its own token); add --real for zenodo.org.

Usage:
  python zenodo_upload.py <repo_folder> --dry-run              show the plan, send nothing
  python zenodo_upload.py <repo_folder>                        sandbox: create draft, upload files, set metadata
  python zenodo_upload.py <repo_folder> --real                 real Zenodo: the same, still only a DRAFT
  python zenodo_upload.py <repo_folder> --real --publish ID    publish draft ID (permanent; asks you to type PUBLISH)

Metadata is read from <repo_folder>/.zenodo.json. Files uploaded: papers/*/*v1.0*.md and .pdf, the figure, and a zip of code/ + README + LICENSE files.
Only the standard library is used."""
import argparse, io, json, os, sys, urllib.error, urllib.request, zipfile

def api(base, token, method, path, data=None, raw=None, ctype="application/json"):
    url = path if path.startswith("http") else base + path
    body = raw if raw is not None else (json.dumps(data).encode("utf-8") if data is not None else None)
    req = urllib.request.Request(url, data=body, method=method,
                                 headers={"Authorization": "Bearer " + token, "Content-Type": ctype})
    try:
        with urllib.request.urlopen(req, timeout=600) as r:
            t = r.read().decode("utf-8")
            return json.loads(t) if t else {}
    except urllib.error.HTTPError as e:
        sys.exit("Zenodo error %d on %s %s: %s" % (e.code, method, path, e.read().decode("utf-8", "replace")[:600]))

def collect(repo):
    files = []
    for line in sorted(os.listdir(os.path.join(repo, "papers"))):
        d = os.path.join(repo, "papers", line)
        for f in sorted(os.listdir(d)):
            if "v1.0" in f or f.endswith(".png"):
                files.append((f, open(os.path.join(d, f), "rb").read()))
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for top in ["code", "README.md", "README.zh-TW.md", "LICENSE", "LICENSE-papers.md", "CITATION.cff", "requirements.txt", "SHA256SUMS.txt"]:
            p = os.path.join(repo, top)
            if os.path.isdir(p):
                for root, _, fs in os.walk(p):
                    for f in fs:
                        if not f.endswith((".pyc", ".log")):
                            full = os.path.join(root, f); z.write(full, "ifo-fep-code-v1.0/" + os.path.relpath(full, repo).replace("\\", "/"))
            elif os.path.exists(p):
                z.write(p, "ifo-fep-code-v1.0/" + top)
    files.append(("ifo-fep-code-v1.0.zip", buf.getvalue()))
    return files

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("repo"); ap.add_argument("--real", action="store_true"); ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--publish", metavar="ID", help="publish the draft with this deposition id")
    a = ap.parse_args()
    base = "https://zenodo.org/api" if a.real else "https://sandbox.zenodo.org/api"
    site = "REAL zenodo.org" if a.real else "SANDBOX (test; no real DOI)"
    token = os.environ.get("ZENODO_API_KEY", "")
    if not os.path.isdir(a.repo) or not os.path.exists(os.path.join(a.repo, ".zenodo.json")):
        sys.exit("First argument must be the LOCAL repo folder (with .zenodo.json), not a URL. Example: python zenodo_upload.py ..\.. --real")
    meta = json.load(open(os.path.join(a.repo, ".zenodo.json"), encoding="utf-8"))
    files = collect(a.repo)
    print("Site    :", site); print("Token   :", "set" if token else "NOT SET (ZENODO_API_KEY)")
    print("Title   :", meta["title"][:100]); print("License :", meta.get("license"))
    for n, b in files: print("  file  %-48s %9d bytes" % (n, len(b)))
    if a.dry_run: print("Dry run: nothing sent."); return
    if not token: sys.exit("Set ZENODO_API_KEY first (a token made on %s)." % ("zenodo.org" if a.real else "sandbox.zenodo.org"))
    if a.publish:
        print("About to PUBLISH deposition %s on %s. This is permanent: the files and the DOI cannot be removed." % (a.publish, site))
        if input("Type PUBLISH to continue: ").strip() != "PUBLISH": sys.exit("Not published.")
        r = api(base, token, "POST", "/deposit/depositions/%s/actions/publish" % a.publish)
        print("Published. DOI:", r.get("doi"), "URL:", r.get("links", {}).get("html")); return
    dep = api(base, token, "POST", "/deposit/depositions", {})
    bucket = dep["links"]["bucket"]; did = dep["id"]
    for n, b in files:
        api(base, token, "PUT", bucket + "/" + n, raw=b, ctype="application/octet-stream"); print("uploaded", n)
    md = {k: v for k, v in meta.items() if k != "version"}; md["version"] = meta.get("version", "1.0")
    md["publication_date"] = "2026-10-20"
    md["creators"] = meta["creators"]
    md["prereserve_doi"] = True        # reserve the DOI now, so it can be written into the papers before publishing
    r = api(base, token, "PUT", "/deposit/depositions/%s" % did, {"metadata": md})
    print("Draft ready (NOT published). Deposition id:", did)
    print("Reserved DOI:", (r.get("metadata", {}).get("prereserve_doi") or {}).get("doi"))
    print("Check it at:", dep["links"].get("html"))
    print("To publish: python zenodo_upload.py <repo> %s--publish %s" % ("--real " if a.real else "", did))

if __name__ == "__main__":
    main()
