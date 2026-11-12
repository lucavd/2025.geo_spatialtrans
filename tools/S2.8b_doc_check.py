#!/usr/bin/env python3
"""tools/S2.8b_doc_check.py - S2.8b, cancello G-1 e controprova CP-1 (docs/steps/S2.8b.md).

Derivato da tools/S1.6b_doc_check.py. (a) I documenti nominano solo cose che esistono: percorsi del repository fra apici inversi, anche
nei blocchi bash; link markdown relativi (risolti dalla cartella del documento); funzioni R scritte come `nome()` o `pkg::nome()` e chiamate
`nome(` nei blocchi di codice R (```r e ```{r} di Quarto). Un percorso esiste se e' un file o una cartella dell'albero git di --tree
(segnaposto `<...>` e `*` come prefissi); percorsi assoluti o con ~ si cercano sul disco del server; percorsi ignorati da git in IGNORED_OK
se esistono nella copia di lavoro. Una funzione esiste se e' definita (`nome <- function`) in un .R dell'albero o nei documenti controllati,
o esportata da un pacchetto di base o della libreria del progetto. I file inclusi con {{< include ... >}} si controllano anche loro.
(b) Copertura: ogni argomento di simulate_visium_hd_pieces() (dal sorgente di --tree) compare fra apici inversi nella guida d'uso; ogni
funzione marcata @export in R/05b*.R compare nella guida d'uso o nel README.

Uso (radice del repository):
  python3 -B tools/S2.8b_doc_check.py [--tree REV] [--text REV|worktree]   -> PASS/FAIL per riferimento e per voce di copertura, RESULT
  python3 -B tools/S2.8b_doc_check.py --mutants        -> quattro mutanti (percorso inventato, funzione inventata, guida senza `ceiling_q`,
                                                          percorso troncato; il quarto aggiunto dopo la revisione, T01)
  python3 -B tools/S2.8b_doc_check.py --counterproof   -> CP-1: copertura di README.md e CLAUDE.md di 48a5f18 (prima di S2.8b)
"""
import argparse, os, re, subprocess, sys

LIB = "renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"
IGNORED_OK = {LIB, "data/real/datasets", ".venv", "renv/library"}
EXT = (".R", ".r", ".md", ".sh", ".py", ".qmd", ".csv", ".lock", ".json", ".html", ".txt", ".tsv", ".parquet", ".rds", ".js", ".docx")
GUIDE = "vignettes/visium_hd_simulation.qmd"
DEFAULT_FILES = ["README.md", "CLAUDE.md", GUIDE, "reports/S2.8b_guida_revisione.qmd"]
NOT_PATHS = {"n/m"}
ENTRY = "simulate_visium_hd_pieces"
BEFORE = "48a5f18"

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

def with_includes(rev, files):
    out, todo = {}, list(files)
    while todo:
        f = todo.pop(0)
        if f in out:
            continue
        t = text_of(rev, f); out[f] = t
        if t:
            for inc in re.findall(r"\{\{<\s*include\s+([^\s>]+)\s*>\}\}", t):
                todo.append(os.path.normpath(os.path.join(os.path.dirname(f), inc)))
    return out

CODE = re.compile(r"```\{?([A-Za-z]*)[^\n]*\n(.*?)```", re.S)

def refs(text, docname):
    paths, funs, defs, libs = set(), set(), set(), set()
    blocks = CODE.findall(text)
    prose = CODE.sub(" ", text)
    for tok in re.findall(r"`([^`\n]+)`", prose):
        tok = tok.strip()
        m = re.fullmatch(r"((?:[A-Za-z][A-Za-z0-9.]*::)?[A-Za-z_.][A-Za-z0-9_.]*)\(\)", tok)
        if m:
            funs.add(m.group(1)); continue
        if " " in tok and not tok.startswith(("http", "{", "r ")):
            if re.fullmatch(r"(?:R|tools|docs|reports|results|data|inst|vignettes)/[^`]+/", tok):
                paths.add(tok); continue                                    # cartella con uno spazio nel nome (R/test functions/)
            # comandi e frasi fra apici inversi: si controllano le parole che sono percorsi del repository (revisione T02)
            for w in re.findall(r"(?<![\w/~.$-])((?:R|tools|docs|reports|results|data|inst|vignettes)/[A-Za-z0-9_./<>*-]+)", tok):
                paths.add(w.rstrip(".,;:)"))
            continue
        if tok.startswith(("http", "-", "$", "{")) or tok in NOT_PATHS:
            continue
        if "/" in tok or tok.endswith(EXT):
            paths.add(tok.rstrip(".,;:)"))
    for m in re.finditer(r"\]\(([^)\s]+)\)", prose):
        t = m.group(1)
        if t.startswith(("http", "#", "mailto:", "{{")):
            continue
        t = t.split("#")[0]
        if t:
            paths.add("./" + os.path.normpath(os.path.join(os.path.dirname(docname), t)))
    for lang, body in blocks:
        if lang.lower() == "r":
            body_nc = re.sub(r'"(?:[^"\\]|\\.)*"', '""', body)   # prima le stringhe, poi i commenti (revisione T02)
            body_nc = re.sub(r"#.*", "", body_nc)
            body_nc = re.sub(r"`[^`\n]*`", "x_", body_nc)   # nomi fra apici inversi (colonne come `distanza minima (um)`)
            defs.update(re.findall(r"([A-Za-z_.][A-Za-z0-9_.]*)\s*<-\s*function", body_nc))
            libs.update(re.findall(r"\blibrary\(\s*([A-Za-z][A-Za-z0-9.]*)\s*\)", body_nc))
            for m in re.finditer(r"(?<![\w$@.])((?:[A-Za-z][A-Za-z0-9.]*::)?[A-Za-z_.][A-Za-z0-9_.]*)\s*\(", body_nc):
                funs.add(m.group(1))
            for m in re.finditer(r'"([A-Za-z0-9_./-]+\.(?:R|csv|json|md|qmd|txt|rds))"', body):
                p = m.group(1)
                if p.startswith("../") or p.startswith("./"):
                    p = "./" + os.path.normpath(os.path.join(os.path.dirname(docname), p))
                paths.add(p)
        if lang.lower() in ("bash", "sh"):
            for m in re.finditer(r"(?<![\w/~.-])((?:R|tools|docs|reports|results|data|inst|vignettes)/[A-Za-z0-9_./<>*$-]+)", body):
                paths.add(re.sub(r"\$[A-Za-z_]\w*", "*", m.group(1)).rstrip(".,;:)"))   # variabili della shell come *
    return paths, funs, defs | {"lib:" + l for l in libs}

def path_ok(p, files, dirs):
    if p.startswith(("/", "~")):
        return os.path.exists(os.path.expanduser(p.replace("<", "").replace(">", "")))
    q = p[2:] if p.startswith("./") else p
    q = q.rstrip("/")
    if q in IGNORED_OK:
        return os.path.exists(q)
    if q in files or q in dirs:
        return True
    pat = re.sub(r"<[^>]*>", "\x00", q)
    wild = "\x00" in pat or "*" in pat
    pat = re.escape(pat).replace("\x00", ".+").replace(r"\*", ".*")
    # ancorato alla fine (revisione di S2.8b, T01): un prefisso troncato non passa; dopo un segnaposto o un * si accetta il resto del percorso
    rx = re.compile("^" + pat + ("(/.*)?$" if not wild else ".*$"))
    return any(rx.match(f) for f in files) or any(rx.match(d) for d in dirs)

def r_sources(files, rev, pattern=r"\.(R|r)$"):
    out = {}
    for f in files:
        if re.search(pattern, f):
            try:
                out[f] = sh("git", "show", f"{rev}:{f}")
            except subprocess.CalledProcessError:
                pass
    return out

def r_defined(srcs):
    defs = set()
    for t in srcs.values():
        defs.update(re.findall(r"^\s*([A-Za-z_.][A-Za-z0-9_.]*)\s*(?:<-|=)\s*function", t, re.M))
    return defs

def entry_args(srcs):
    for f, t in srcs.items():
        m = re.search(r"^" + ENTRY + r"\s*<-\s*function\(", t, re.M)
        if not m:
            continue
        i, depth, cur, args = m.end(), 1, "", []
        while depth:
            ch = t[i]
            if ch in "([{": depth += 1
            elif ch in ")]}": depth -= 1
            if depth == 1 and ch == ",":
                args.append(cur); cur = ""
            elif depth:
                cur += ch
            i += 1
        args.append(cur)
        return [a.split("=")[0].strip() for a in args if a.strip()]
    return []

def exported(srcs):
    ex = []
    for f, t in sorted(srcs.items()):
        if re.search(r"(^|/)R/05b[^/]*\.R$", f):
            ex += re.findall(r"#'\s*@export\s*\n([A-Za-z_.][A-Za-z0-9_.]*)\s*<-\s*function", t)
    return ex

R_CHECK = r"""
args <- commandArgs(TRUE); lib <- args[1]; pk <- strsplit(args[2], ",")[[1]]; fs <- args[-(1:2)]
if (dir.exists(lib)) .libPaths(lib)
base <- c("base", "stats", "utils", "graphics", "grDevices", "methods", "tools", "parallel", pk[nzchar(pk)])
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
    pk = ",".join(sorted(d[4:] for d in defs if d.startswith("lib:")))
    out = subprocess.run(["Rscript", "--vanilla", "-e", R_CHECK, LIB, pk or ","] + rest, capture_output=True, text=True).stdout
    okr = {l.split()[1]: l.split()[0] == "OK" for l in out.splitlines() if l.startswith(("OK", "NO"))}
    return {f: (f in defs) or okr.get(f, False) for f in funs}

def check_refs(texts, files, dirs, defs, verbose=True):
    miss, n = [], 0
    docdefs = set()
    parsed = {}
    for name, text in texts.items():
        if text is None:
            continue
        parsed[name] = refs(text, name); docdefs |= parsed[name][2]
    for name, text in texts.items():
        if text is None:
            print(f"FAIL {name}: documento assente"); miss.append((name, "documento")); continue
        paths, funs, own = parsed[name]
        fo = fun_ok(funs, defs | {d for d in docdefs if not d.startswith("lib:")} | {d for d in own if d.startswith("lib:")})
        for p in sorted(paths):
            n += 1; ok = path_ok(p, files, dirs)
            if verbose or not ok: print(f"{'PASS' if ok else 'FAIL'} {name} percorso {p}")
            if not ok: miss.append((name, p))
        for f in sorted(funs):
            n += 1; ok = fo[f]
            if verbose or not ok: print(f"{'PASS' if ok else 'FAIL'} {name} funzione {f}()")
            if not ok: miss.append((name, f + "()"))
    return n, miss

def check_coverage(guide, readme, args, ex, verbose=True):
    miss, n = [], 0
    g = guide or ""; r = readme or ""
    for a in args:
        n += 1; ok = re.search(r"`" + re.escape(a) + r"(\s*=[^`]*)?`", g) is not None
        if verbose or not ok: print(f"{'PASS' if ok else 'FAIL'} copertura argomento `{a}` di {ENTRY}() nella guida")
        if not ok: miss.append(("copertura", a))
    for f in ex:
        n += 1; rx = re.compile(r"\b" + re.escape(f) + r"\(")
        ok = rx.search(g) is not None or rx.search(r) is not None
        if verbose or not ok: print(f"{'PASS' if ok else 'FAIL'} copertura funzione esportata {f}() nella guida o nel README")
        if not ok: miss.append(("copertura", f + "()"))
    return n, miss

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", default="HEAD"); ap.add_argument("--text", default="worktree")
    ap.add_argument("--mutants", action="store_true"); ap.add_argument("--counterproof", action="store_true")
    a = ap.parse_args()
    files, dirs = tree_of(a.tree); srcs = r_sources(files, a.tree); defs = r_defined(srcs)
    args = entry_args(srcs); ex = exported(srcs)
    print(f"INFO albero {sh('git', 'rev-parse', '--short', a.tree).strip()}: {len(args)} argomenti di {ENTRY}(), {len(ex)} funzioni esportate in R/05b*.R")
    if not args or not ex:
        print("FAIL argomenti o funzioni esportate non trovati nel sorgente"); sys.exit(1)
    if a.mutants:
        base = with_includes(a.text, DEFAULT_FILES)
        g0 = base[GUIDE]
        muts = {
            "N1 percorso inventato": ({**base, "README.md": base["README.md"] + "\nSee `R/05b9_imaginary_step.R`.\n"}, "R/05b9_imaginary_step.R"),
            "N2 funzione inventata": ({**base, GUIDE: g0 + "\n```{r}\nsimulate_visium_grid_fast(cl, seed = 1)\n```\n"}, "simulate_visium_grid_fast()"),
            "N3 guida senza `ceiling_q`": ({**base, GUIDE: "\n".join(l for l in g0.split("\n") if "ceiling_q" not in l)}, "ceiling_q"),
            "N4 percorso troncato": ({**base, "CLAUDE.md": base["CLAUDE.md"] + "\nVedi `vignettes/visium_hd_simulation.qm`.\n"}, "vignettes/visium_hd_simulation.qm"),
        }
        det = 0
        for k, (texts, target) in muts.items():
            _, m1 = check_refs(texts, files, dirs, defs, verbose=False)
            _, m2 = check_coverage(texts[GUIDE], texts["README.md"], args, ex, verbose=False)
            hit = any(m[1] == target for m in m1 + m2)
            det += hit; print(f"{'PASS' if hit else 'FAIL'} mutante {k} rilevato")
        print(f"RESULT {det}/{len(muts)} mutanti rilevati"); sys.exit(0 if det == len(muts) else 1)
    if a.counterproof:
        r0 = text_of(BEFORE, "README.md"); c0 = text_of(BEFORE, "CLAUDE.md")
        n, miss = check_coverage((r0 or "") + "\n" + (c0 or ""), (r0 or "") + "\n" + (c0 or ""), args, ex, verbose=False)
        print(f"INFO copertura di README.md e CLAUDE.md di {BEFORE}: {n - len(miss)}/{n} voci coperte; scoperte: " + ", ".join(m[1] for m in miss))
        bocciato = len(miss) > 0
        print(f"{'PASS' if bocciato else 'FAIL'} CP-1: il controllo di copertura boccia i documenti di prima di S2.8b"); sys.exit(0 if bocciato else 1)
    texts = with_includes(a.text, DEFAULT_FILES)
    n1, m1 = check_refs(texts, files, dirs, defs)
    n2, m2 = check_coverage(texts.get(GUIDE), texts.get("README.md"), args, ex)
    n, miss = n1 + n2, m1 + m2
    print(f"RESULT {n - len(miss)}/{n} PASS ({len(miss)} mancanti; riferimenti {n1}, copertura {n2}; documenti {len(texts)}; albero {sh('git', 'rev-parse', '--short', a.tree).strip()}, testo {a.text})")
    sys.exit(0 if not miss else 1)

if __name__ == "__main__":
    main()
