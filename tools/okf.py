#!/usr/bin/env python3
"""Open Knowledge Format v0.2 tooling for the geo_spatialtrans documentation bundle (docs/).

Adapted from MMM tools/okf.py (2026-10-20, geo session protocol v2, first v2 session S1.4d).

  python3 -B tools/okf.py check --bundle docs
      Rebuilds every index.md of the bundle from the frontmatter, then checks the bundle.
      Run before each commit; exit status 1 on failure.

  python3 -B tools/okf.py link --bundle docs FILE [FILE ...] [--no-steps]
      In place: bare decision codes (DEC-nnn) become links whose text is the decision title;
      step codes (S1.4d, R2c, F3.1, DOC1 ...) become links to docs/steps/<id>.md, only when
      that step file exists (steps before v2 have no step file).

The bundle is docs/decisions/, docs/steps/, docs/index.md and docs/log.md. Every other
entry of docs/ must be in EXEMPT (documents written before v2, waiting for the migration
step DOC1); an entry that is neither in the bundle nor exempt fails the check, so new
documents cannot silently land outside the bundle.

Checks on bundle files: every non-reserved .md has frontmatter with non-empty type, title
and description; reserved files (index.md, log.md) carry no frontmatter, except okf_version
in the root index.md; no decision code DEC-nnn outside a link; no markdown file of the
bundle ignored by .gitignore. Codes of other registers (BL-nnn, B-nnn, CP-n, P-n) are
warnings until those registers become concepts. Standard library only.
"""
import argparse
import json
import os
import pathlib
import re
import subprocess
import sys

CODE_RE = re.compile(r"(?<![\w\[/-])(DEC-\d{3})(?![\w])")
STEP_RE = re.compile(r"(?<![\w\[/.-])(S\d(?:\.\d)?[a-z]?|R\d[a-z]?|F3\.\d|DOC1)(?![\w.]|-)")
PROTECT_RE = re.compile(r"(`[^`]*`|\[[^\]]*\]\([^)]*\))")
BARE_CODE_RE = re.compile(r"\bDEC-\d{3}\b")
OTHER_CODES_RE = re.compile(r"\b(BL-\d{3}|B-\d{3}|CP-\d+[a-z]?|P-\d+[a-z]?)\b")
RESERVED = {"index.md", "log.md"}
BUNDLE_DIRS = ("decisions", "steps")
BUNDLE_FILES = ("index.md", "log.md")
# documents of docs/ written before protocol v2: exempt until the migration step (DOC1)
EXEMPT = {
    "ROADMAP.md": "Roadmap del progetto: principi, pannello di tessuti, decisioni fino al 2026-10-20, stato delle sessioni.",
    "BACKLOG.md": "Voci aperte e rinviate (BL-nnn), con la sessione di origine e dove affrontarle.",
    "BIO_REFERENCES.md": "Riferimenti biologici (B-nnn) con fonte: DOI o script e dataset che li hanno stimati.",
    "hotspot.md": "Nota storica sugli hotspot della pipeline precedente.",
    "image_gen_explainer.md": "Spiegazione della generazione delle immagini sintetiche della pipeline precedente.",
    "SPATIAL_TRANSCRIPTOMICS_VALIDATION_GUIDE.md": "Guida alla validazione dei dati di trascrittomica spaziale.",
    "sessions": "Note di sessione fino a R2c (protocollo v1).",
    "superpowers": "Design dello Step 1 (specifiche) e materiali collegati.",
}
DIR_HEADINGS = {"steps": "Step", "decisions": "Decisioni"}
DIR_DESCRIPTIONS = {
    "decisions": "Una decisione per file, con titolo e descrizione in parole semplici, chi l'ha presa e se e' in vigore.",
    "steps": "Un file per step (dal protocollo v2): attese pre-registrate, riassunto In breve, approvazione del PI, link al report.",
}


def id_key(code):
    m = re.match(r"DEC-(\d+)", code)
    return (int(m.group(1)), "") if m else (10**6, code)


def read_frontmatter(path):
    """Return (dict of flat keys, body). Quoted values are JSON-decoded; others kept raw."""
    txt = pathlib.Path(path).read_text(encoding="utf-8")
    if not txt.startswith("---\n") or "\n---\n" not in txt[4:]:
        return None, txt
    fm_txt, body = txt[4:].split("\n---\n", 1)
    fm = {}
    for line in fm_txt.splitlines():
        m = re.match(r"^([A-Za-z_]+):\s*(.*)$", line)
        if not m:
            continue
        val = m.group(2).strip()
        if val.startswith('"'):
            try:
                val = json.loads(val)
            except json.JSONDecodeError:
                pass
        fm[m.group(1)] = val
    return fm, body


def decision_titles(bundle):
    titles = {}
    for p in pathlib.Path(bundle, "decisions").glob("DEC-*.md"):
        fm, _ = read_frontmatter(p)
        if fm and fm.get("title"):
            titles[p.stem] = fm["title"]
    return titles


def step_titles(bundle):
    titles = {}
    for p in pathlib.Path(bundle, "steps").glob("*.md"):
        if p.name in RESERVED:
            continue
        fm, _ = read_frontmatter(p)
        if fm and fm.get("title"):
            titles[fm.get("step_id") or p.stem] = (fm["title"], p.name)
    return titles


def link_text(text, titles, steps, bundle, file_dir, link_steps=True):
    """Link bare DEC-nnn (everywhere) and step ids with a step file (outside headings)."""
    bundle, file_dir = os.path.abspath(bundle), os.path.abspath(file_dir)

    def rel(*parts):
        return os.path.relpath(os.path.join(bundle, *parts), file_dir)

    def dec(m):
        code = m.group(1)
        return f"[{titles[code]} ({code})]({rel('decisions', code + '.md')})" if code in titles else code

    def step(m):
        code = m.group(1)
        if code not in steps:
            return code
        title, fname = steps[code]
        return f"[{title} ({code})]({rel('steps', fname)})"

    out, in_fence = [], False
    for line in text.split("\n"):
        if line.lstrip().startswith("```"):
            in_fence = not in_fence
        if in_fence or line.lstrip().startswith("```"):
            out.append(line)
            continue
        heading = line.lstrip().startswith("#")
        parts = PROTECT_RE.split(line)
        for k in range(0, len(parts), 2):
            p = CODE_RE.sub(dec, parts[k])
            if link_steps and not heading:
                p = STEP_RE.sub(step, p)
            parts[k] = p
        out.append("".join(parts))
    return "\n".join(out)


def link_file(path, titles, steps, bundle, link_steps=True):
    path = pathlib.Path(path)
    txt = path.read_text(encoding="utf-8")
    if txt.startswith("---\n") and "\n---\n" in txt[4:]:
        fm_txt, body = txt[4:].split("\n---\n", 1)
        head = "---\n" + fm_txt + "\n---\n"
    else:
        head, body = "", txt
    path.write_text(head + link_text(body, titles, steps, bundle, path.parent, link_steps), encoding="utf-8")


def _strip_protected(body):
    return PROTECT_RE.sub("", body)


def check_file(path, bundle):
    fm, body = read_frontmatter(path)
    problems, warns = [], []
    if path.name in RESERVED:
        is_root_index = path.name == "index.md" and path.parent.resolve() == pathlib.Path(bundle).resolve()
        if fm is not None and not (is_root_index and set(fm) <= {"okf_version"}):
            problems.append("reserved file with frontmatter")
    else:
        if fm is None:
            return ["missing frontmatter"], []
        problems += [f"missing {k}" for k in ("type", "title", "description") if not fm.get(k)]
        for k in ("title", "description"):
            if BARE_CODE_RE.search(str(fm.get(k, ""))):
                warns.append(f"decision code in {k}")
        warns += sorted(set(OTHER_CODES_RE.findall(_strip_protected(body))))
    bare = sorted(set(BARE_CODE_RE.findall(_strip_protected(body))), key=id_key)
    if bare:
        problems.append("bare decision codes: " + ", ".join(bare))
    return problems, warns


def concepts_in(d):
    return [p for p in sorted(pathlib.Path(d).glob("*.md")) if p.name not in RESERVED]


def bundle_files(bundle):
    bundle = pathlib.Path(bundle)
    out = [bundle / f for f in BUNDLE_FILES if (bundle / f).exists()]
    for d in BUNDLE_DIRS:
        if (bundle / d).is_dir():
            out += sorted((bundle / d).rglob("*.md"))
    return out


def stray_entries(bundle):
    """Entries of docs/ that are neither bundle nor exempt."""
    allowed = set(BUNDLE_DIRS) | set(BUNDLE_FILES) | set(EXEMPT)
    return sorted(p.name for p in pathlib.Path(bundle).iterdir() if p.name not in allowed and not p.name.startswith("."))


def build_indexes(bundle):
    bundle = pathlib.Path(bundle)
    for name in BUNDLE_DIRS:
        d = bundle / name
        files = concepts_in(d) if d.is_dir() else []
        if not files:
            continue
        lines = [f"# {DIR_HEADINGS.get(name, name)}\n"]
        key = (lambda p: id_key(p.stem)) if name == "decisions" else (lambda p: p.name)
        for p in sorted(files, key=key):
            fm, _ = read_frontmatter(p)
            fm = fm or {}
            code = f" ({fm.get('step_id') or p.stem})" if name == "steps" else f" ({p.stem})"
            flag = " (superata)" if fm.get("status") == "deprecated" else ""
            lines.append(f"* [{fm.get('title')}{code}]({p.name}){flag} - {fm.get('description')}")
        lines.append("")
        (d / "index.md").write_text("\n".join(lines), encoding="utf-8")
    lines = ['---', 'okf_version: "0.2"', '---', '# Cartelle\n']
    for name in BUNDLE_DIRS:
        if (bundle / name).is_dir() and concepts_in(bundle / name):
            lines.append(f"* [{DIR_HEADINGS[name]}]({name}/) - {DIR_DESCRIPTIONS[name]}")
    lines.append("* [Registro delle modifiche](log.md) - Modifiche al bundle per data, la piu' recente in alto.")
    lines.append("\n# Documenti non ancora in OKF\n")
    lines.append("Scritti prima del protocollo v2; restano come sono fino alla migrazione (step DOC1).\n")
    for name, desc in EXEMPT.items():
        if (bundle / name).exists():
            target = name + ("/" if (bundle / name).is_dir() else "")
            lines.append(f"* [{name}]({target}) - {desc}")
    lines.append("")
    (bundle / "index.md").write_text("\n".join(lines), encoding="utf-8")


def ignored_by_git(bundle):
    """Markdown files of the bundle that git would silently skip (None outside a git repo)."""
    paths = [str(pathlib.Path(bundle, x)) for x in BUNDLE_DIRS + BUNDLE_FILES]
    try:
        out = subprocess.run(["git", "ls-files", "--others", "--ignored", "--exclude-standard", "--"] + paths,
                             capture_output=True, text=True, check=True).stdout.split("\n")
    except Exception:
        return None
    return [p for p in out if p.endswith(".md") or p.endswith("/")]


def check(bundle):
    build_indexes(bundle)
    n_bad, n_files, warn_codes = 0, 0, set()
    stray = stray_entries(bundle)
    if stray:
        n_bad += 1
        print(f"FAIL outside the bundle and not exempt (move to docs/decisions or docs/steps, or add to EXEMPT): {', '.join(stray)}")
    ignored = ignored_by_git(bundle)
    if ignored:
        n_bad += 1
        print(f"FAIL ignored by git (fix .gitignore): {', '.join(ignored)}")
    for p in bundle_files(bundle):
        n_files += 1
        problems, warns = check_file(p, bundle)
        if problems:
            n_bad += 1
            print(f"FAIL {p.relative_to(bundle)}: {'; '.join(problems)}")
        warn_codes.update(w for w in warns if not w.startswith("decision code"))
        for w in warns:
            if w.startswith("decision code"):
                print(f"WARN {p.relative_to(bundle)}: {w}")
    if warn_codes:
        print(f"WARN codes of other registers still unlinked: {', '.join(sorted(warn_codes))}")
    print(f"{n_files} markdown files of the bundle checked, {n_bad} failing")
    return n_bad


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="mode", required=True)
    c = sub.add_parser("check")
    c.add_argument("--bundle", required=True)
    lk = sub.add_parser("link")
    lk.add_argument("--bundle", required=True)
    lk.add_argument("--no-steps", action="store_true")
    lk.add_argument("files", nargs="+")
    a = ap.parse_args()
    if a.mode == "link":
        titles, steps = decision_titles(a.bundle), step_titles(a.bundle)
        for f in a.files:
            link_file(f, titles, steps, a.bundle, link_steps=not a.no_steps)
        return
    sys.exit(1 if check(a.bundle) else 0)


if __name__ == "__main__":
    main()
