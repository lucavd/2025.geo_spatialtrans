#!/usr/bin/env python3
"""tools/test_okf.py — tests of tools/okf.py on a temporary copy of the bundle (the repo is not touched).
Usage (repo root): python3 -B tools/test_okf.py   -> PASS/FAIL per assertion, exit 1 on any FAIL."""
import pathlib, shutil, subprocess, sys, tempfile
ROOT = pathlib.Path(__file__).resolve().parents[1]
OKF = ROOT / "tools/okf.py"
sys.path.insert(0, str(ROOT / "tools")); import okf  # noqa: E402
res = []


def ck(name, ok, val=""):
    res.append(bool(ok)); print(f"{'PASS' if ok else 'FAIL'} {name} {val}")


def copy_bundle(dst):
    src = ROOT / "docs"
    for d in okf.BUNDLE_DIRS:
        if (src / d).is_dir():
            shutil.copytree(src / d, dst / d)
    for f in okf.BUNDLE_FILES:
        if (src / f).exists():
            shutil.copy(src / f, dst / f)
    for name in okf.EXEMPT:                      # exempt entries: placeholders, never read
        (dst / name).mkdir() if (src / name).is_dir() else (dst / name).write_text("x\n")


def run(b):
    p = subprocess.run([sys.executable, "-B", str(OKF), "check", "--bundle", str(b)], capture_output=True, text=True, cwd="/tmp")
    return p.returncode, p.stdout


with tempfile.TemporaryDirectory() as t:
    b = pathlib.Path(t) / "docs"; b.mkdir(); copy_bundle(b)
    rc, out = run(b)
    ck("T1 check passes on a copy of the repo bundle", rc == 0, out.strip().splitlines()[-1] if out else "")
    ck("T1b root index.md written with okf_version", (b / "index.md").read_text().startswith('---\nokf_version: "0.2"\n---'))
    # mutant M1: bundle file without frontmatter
    m1 = b / "decisions" / "DEC-999.md"; m1.parent.mkdir(exist_ok=True); m1.write_text("# no frontmatter\n")
    rc, out = run(b); ck("T2 mutant M1 (decision without frontmatter) fails", rc == 1 and "DEC-999.md: missing frontmatter" in out)
    m1.unlink()
    # mutant M2: bare decision code in a bundle file
    m2 = b / "steps" / "ZZ.md"; m2.parent.mkdir(exist_ok=True)
    m2.write_text('---\ntype: Step\ntitle: "t"\ndescription: "d"\n---\nVedi DEC-001 per il protocollo.\n')
    rc, out = run(b); ck("T3 mutant M2 (bare DEC-nnn) fails", rc == 1 and "bare decision codes: DEC-001" in out)
    # the same code inside a link or in backticks is fine
    m2.write_text('---\ntype: Step\ntitle: "t"\ndescription: "d"\n---\nVedi [protocollo (DEC-001)](../decisions/DEC-001.md) e `DEC-002`.\n')
    rc, out = run(b); ck("T4 linked / quoted DEC-nnn passes", rc == 0, out.strip().splitlines()[-1])
    m2.unlink()
    # mutant M3: new document outside the bundle and not exempt
    (b / "NEW_NOTE.md").write_text("x\n")
    rc, out = run(b); ck("T5 mutant M3 (document outside the bundle, not exempt) fails", rc == 1 and "NEW_NOTE.md" in out)
    (b / "NEW_NOTE.md").unlink()
    # exempt files are not checked for frontmatter
    rc, out = run(b); ck("T6 exempt documents without frontmatter do not fail", rc == 0)
# link mode: step ids link only where a step file exists
txt = okf.link_text("Come in S1.4c e in S1.4d; vedi DEC-001.", {"DEC-001": "Protocollo"}, {"S1.4d": ("Forma del nucleo", "S1.4d.md")}, "/b", "/b/steps")
ck("T7 link: DEC-001 and S1.4d linked, S1.4c (no step file) left as is",
   txt == "Come in S1.4c e in [Forma del nucleo (S1.4d)](S1.4d.md); vedi [Protocollo (DEC-001)](../decisions/DEC-001.md).", txt)
nf = res.count(False); print(f"RESULT {len(res) - nf}/{len(res)} PASS"); sys.exit(1 if nf else 0)
