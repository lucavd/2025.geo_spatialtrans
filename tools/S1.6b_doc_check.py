#!/usr/bin/env python3
"""tools/S1.6b_doc_check.py - S1.6b, cancello G-1 e controprova (docs/steps/S1.6b.md).

I documenti nominano solo cose che esistono. Da ogni documento si raccolgono:
  * percorsi del repository fra apici inversi (token senza spazi con '/' o con un'estensione nota), anche nei blocchi bash;
  * funzioni R scritte come `nome()` o `pkg::nome()` nel testo e le chiamate `nome(` dentro i blocchi ```r.
Un percorso esiste se e' un file o una cartella dell'albero git del commit (--tree); con segnaposto `<...>` o senza corrispondenza esatta
si accetta un prefisso (`R/01` -> R/01*). Percorsi assoluti o con ~ si cercano sul disco del server. Percorsi ignorati da git (libreria, .venv,
data/real/datasets) si accettano se sono in IGNORED_OK ed esistono nella copia di lavoro. Una funzione esiste se e' definita
(`nome <- function`) in un file .R dell'albero o esportata da un pacchetto di base o della libreria del progetto (controllo in R).

Uso (radice del repository):
  python3 -B tools/S1.6b_doc_check.py [--tree REV] [--text REV|worktree] [--files F ...]   -> PASS/FAIL per riferimento, RESULT
  python3 -B tools/S1.6b_doc_check.py --mutants   -> due mutanti (percorso inventato, funzione inventata) su copie del README: vanno rilevati
  python3 -B tools/S1.6b_doc_check.py --counterproof  -> README.md e CLAUDE.md del 2025 (775eec3) contro l'albero di --tree
"""
import argparse, os, re, subprocess, sys

LIB = "renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"
IGNORED_OK = {LIB, "data/real/datasets", ".venv", "renv/library"}
EXT = (".R", ".r", ".md", ".sh", ".py", ".qmd", ".csv", ".lock", ".json", ".html", ".txt", ".tsv", ".parquet", ".rds", ".js", ".docx")
DEFAULT_FILES = ["README.md", "CLAUDE.md", "reports/S1.6b_guida_revisione.qmd"]
NOT_PATHS = {"n/m"}

def sh(*a):
    return subprocess.run(a, capture_output=True, text=True, check=True).stdout

def tree_of(rev):
    files = [f for f in sh("git", "ls-tree", "-r", "--name-only", rev).split("\n") if f]
    dirs = set()
    for f in files:
        p = f.split("/")
        for i in range(1, len(p)):
            dirs.add("/".join(p[:i]))
    return set(files), dirs

def text_of(rev, f):
    if rev == "worktree":
        return open(f, encoding="utf-8").read() if os.path.exists(f) else None
    try:
        return sh("git", "show", f"{rev}:{f}")
    except subprocess.CalledProcessError:
        return None

def refs(text):
    paths, funs = set(), set()
    blocks = re.findall(r"```(\w*)\n(.*?)```", text, re.S)
    prose = re.sub(r"```.*?```", " ", text, flags=re.S)
    for tok in re.findall(r"`([^`\n]+)`", prose):
        tok = tok.strip()
        m = re.fullmatch(r"((?:[A-Za-z][A-Za-z0-9.]*::)?[A-Za-z_.][A-Za-z0-9_.]*)\(\)", tok)
        if m:
            funs.add(m.group(1)); continue
        if " " in tok or tok.startswith(("http", "-", "$")) or tok in NOT_PATHS:
            continue
        if "/" in tok or tok.endswith(EXT):
            paths.add(tok.rstrip(".,;:)"))
    for lang, body in blocks:
        if lang.lower() == "r":
            body_nc = re.sub(r"#.*", "", body)
            body_nc = re.sub(r'"(?:[^"\\]|\\.)*"', '""', body_nc)
            for m in re.finditer(r"((?:[A-Za-z][A-Za-z0-9.]*::)?[A-Za-z_.][A-Za-z0-9_.]*)\s*\(", body_nc):
                funs.add(m.group(1))
            for m in re.finditer(r'"([A-Za-z0-9_./-]+\.(?:R|csv|json|md))"', body):
                paths.add(m.group(1))
        if lang.lower() in ("bash", "sh"):
            for m in re.finditer(r"(?<![\w/~.-])((?:R|tools|docs|reports|results|data|inst)/[A-Za-z0-9_./<>*-]+)", body):
                paths.add(m.group(1).rstrip(".,;:)"))
    return paths, funs

def path_ok(p, files, dirs):
    if p.startswith(("/", "~")):
        return os.path.exists(os.path.expanduser(p.replace("<", "").replace(">", "")))
    q = p.rstrip("/")
    if q in IGNORED_OK:
        return os.path.exists(q)
    if q in files or q in dirs:
        return True
    pat = re.sub(r"<[^>]*>", "\x00", q)
    pat = re.escape(pat).replace("\x00", ".+").replace(r"\*", ".*")
    rx = re.compile("^" + pat)
    return any(rx.match(f) for f in files) or any(rx.match(d) for d in dirs)

def r_defined(files, rev):
    defs = set()
    for f in files:
        if f.endswith((".R", ".r")):
            try:
                t = sh("git", "show", f"{rev}:{f}")
            except subprocess.CalledProcessError:
                continue
            defs.update(re.findall(r"^\s*([A-Za-z_.][A-Za-z0-9_.]*)\s*(?:<-|=)\s*function", t, re.M))
    return defs

R_CHECK = r"""
args <- commandArgs(TRUE); lib <- args[1]; fs <- args[-1]
if (dir.exists(lib)) .libPaths(lib)
base <- c("base", "stats", "utils", "graphics", "grDevices", "methods", "tools", "parallel")
ex <- unique(unlist(lapply(base, function(p) getNamespaceExports(p))))
for (f in fs) {
  ok <- if (grepl("::", f, fixed = TRUE)) { pp <- strsplit(f, "::", fixed = TRUE)[[1]]
    requireNamespace(pp[1], quietly = TRUE) && pp[2] %in% getNamespaceExports(pp[1]) } else f %in% ex ||
    f %in% c("function", "if", "for", "while", "repeat", "return")
  cat(if (ok) "OK" else "NO", f, "\n")
}
"""

def fun_ok(funs, defs):
    rest = sorted(f for f in funs if f not in defs)
    if not rest:
        return {f: True for f in funs}
    out = subprocess.run(["Rscript", "--vanilla", "-e", R_CHECK, LIB] + rest, capture_output=True, text=True).stdout
    okr = {l.split()[1]: l.split()[0] == "OK" for l in out.splitlines() if l.startswith(("OK", "NO"))}
    return {f: (f in defs) or okr.get(f, False) for f in funs}

def check(texts, files, dirs, defs, verbose=True):
    miss = []; n = 0
    for name, text in texts.items():
        if text is None:
            print(f"FAIL {name}: documento assente"); miss.append((name, "documento")); continue
        paths, funs = refs(text)
        fo = fun_ok(funs, defs)
        for p in sorted(paths):
            n += 1; ok = path_ok(p, files, dirs)
            if verbose or not ok: print(f"{'PASS' if ok else 'FAIL'} {name} percorso {p}")
            if not ok: miss.append((name, p))
        for f in sorted(funs):
            n += 1; ok = fo[f]
            if verbose or not ok: print(f"{'PASS' if ok else 'FAIL'} {name} funzione {f}()")
            if not ok: miss.append((name, f + "()"))
    return n, miss

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", default="HEAD"); ap.add_argument("--text", default="worktree")
    ap.add_argument("--files", nargs="*", default=DEFAULT_FILES)
    ap.add_argument("--mutants", action="store_true"); ap.add_argument("--counterproof", action="store_true")
    a = ap.parse_args()
    files, dirs = tree_of(a.tree); defs = r_defined(files, a.tree)
    if a.mutants:
        base = text_of(a.text, "README.md")
        muts = {"N1 percorso inventato": base + "\nSee `R/04b9_cell_imaginary.R`.\n",
                "N2 funzione inventata": base + "\nCall `simulate_visium_grid()` afterwards.\n"}
        det = 0
        for k, t in muts.items():
            _, miss = check({"README.md": t}, files, dirs, defs, verbose=False)
            hit = any(m[1] in ("R/04b9_cell_imaginary.R", "simulate_visium_grid()") for m in miss)
            det += hit; print(f"{'PASS' if hit else 'FAIL'} mutante {k} rilevato")
        print(f"RESULT {det}/2 mutanti rilevati"); sys.exit(0 if det == 2 else 1)
    if a.counterproof:
        texts = {f + "@775eec3": text_of("775eec3", f) for f in ("README.md", "CLAUDE.md")}
        n, miss = check(texts, files, dirs, defs, verbose=False)
        cell = {k: bool(re.search(r"04b[0-6]?_", t or "")) for k, t in texts.items()}
        print(f"INFO riferimenti controllati {n}, mancanti nel codice consolidato {len(miss)}")
        for k, v in cell.items(): print(f"INFO {k} nomina un file R/04b*.R del livello cellulare: {v}")
        bocciato = len(miss) > 0 or not any(cell.values())
        print(f"{'PASS' if bocciato else 'FAIL'} controprova: il controllo boccia i documenti del 2025"); sys.exit(0 if bocciato else 1)
    texts = {f: text_of(a.text, f) for f in a.files}
    n, miss = check(texts, files, dirs, defs)
    print(f"RESULT {n - len(miss)}/{n} PASS ({len(miss)} mancanti; albero {sh('git', 'rev-parse', '--short', a.tree).strip()}, testo {a.text})")
    sys.exit(0 if not miss else 1)

if __name__ == "__main__":
    main()
