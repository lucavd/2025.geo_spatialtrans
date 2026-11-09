"""tools/S2.5_mutants.py — S2.5, G-1: mutanti di R/05b5_simulate_visium_hd.R; ognuno deve far fallire R/testing/test_S2.5.R.
Uso (root del repo): .venv/bin/python tools/S2.5_mutants.py [out.txt]"""
import os, sys, shutil, subprocess, tempfile
F = "R/05b5_simulate_visium_hd.R"
MUT = [
    ("profondita' senza M", "depth_factor = df0 * M, workers = workers)", "depth_factor = df0, workers = workers)"),
    ("semi in progressione", "s <- sample.int(1e8, 4L)", "s <- as.integer(seed) + 0:3"),
    ("blocchi senza M", "cs <- if (is.null(chunk_size)) visium_hd_chunk_size(ntx, chunk_transcripts)", "cs <- if (is.null(chunk_size)) visium_hd_chunk_size(ntx / M, chunk_transcripts)"),
    ("blocco di 8 um sfasato nell'origine", "%/% grid$block + 1L)\n  if (!nrow(df))", "%/% grid$block + 2L)\n  if (!nrow(df))"),
    ("seconda coda ignorata", "p$w2 <- if (ly$second_tail) cp$w2 else 0", "p$w2 <- cp$w2"),
    ("fondo ignorato", "background = ly$background, workers = workers)", "background = TRUE, workers = workers)"),
    ("frame_open non passato", "field = ly$field, frame_open = frame_open, chunk_size = cs", "field = ly$field, frame_open = FALSE, chunk_size = cs"),
    ("griglia non controllata", 'stop("simulate_visium_hd(): `hd` non e\' la griglia di questo livello cellulare.", call. = FALSE)', 'invisible(NULL)'),
    ("generatore non ripristinato (dentro e fuori; il solo interno e' equivalente dopo la correzione G-5)", None, None),
    ("atteso catturato senza M", "cells$n_expected_captured <- cells$n_transcripts_generated / M", "cells$n_expected_captured <- cells$n_transcripts_generated"),
    ("frazioni a 8 um su 4 um2", "ov8$frac_territory <- ov8$area_territory_um2 / (g$bin_um * g$block)^2", "ov8$frac_territory <- ov8$area_territory_um2 / g$bin_um^2"),
    ("aree dei territori non controllate", "anyNA(aH) || max(abs(aH - aT)) > 1e-3)", "FALSE)"),
    ("stato dell'utente senza seme non ripristinato", "            else if (exists(\".Random.seed\", envir = globalenv())) rm(\".Random.seed\", envir = globalenv()) }, add = TRUE)\n  px <-", "            }, add = TRUE)\n  px <-"),
    ("strati senza nome accettati", "(is.null(names(layers)) || any(!nzchar(names(layers))) || anyDuplicated(names(layers)))", "FALSE"),
    ("fattore del preset non usato", "else if (prof_default) visium_hd_depth_factor(archetype) else profiles$depth_factor", "else profiles$depth_factor"),
    ("strisce ignorate", "St <- if (ly$stripes) simulate_stripes(", "St <- if (TRUE) simulate_stripes("),
]
out = sys.argv[1] if len(sys.argv) > 1 else "results/S2.5/S2.5_mutants.txt"
os.makedirs(os.path.dirname(out), exist_ok=True)
src = open(F).read(); lines = []; ok_all = True
env = dict(os.environ, R_LIBS_SITE="/nonexistent", R_LIBS_USER="/nonexistent")
DOUBLE = [("  on.exit({ do.call(RNGkind, as.list(old_kind)); if (had) assign", "  if (FALSE) on.exit({ do.call(RNGkind, as.list(old_kind)); if (had) assign"),
          ("  on.exit({ do.call(RNGkind, as.list(old_kind)); if (had_seed) assign", "  if (FALSE) on.exit({ do.call(RNGkind, as.list(old_kind)); if (had_seed) assign")]
for name, a, b in MUT:
    if a is None:
        ok = all(src.count(x) == 1 for x, _ in DOUBLE)
        if not ok: lines.append(f"ERRORE {name}: frammenti non unici"); ok_all = False; continue
        mutated = src
        for x, y in DOUBLE: mutated = mutated.replace(x, y)
        a, b = "", ""
    else:
        mutated = None
    n = src.count(a) if mutated is None else 1
    if n != 1: lines.append(f"ERRORE {name}: frammento trovato {n} volte"); ok_all = False; continue
    d = tempfile.mkdtemp(prefix="s25mut_"); shutil.copytree("R", f"{d}/R"); open(f"{d}/{F}", "w").write(src.replace(a, b) if mutated is None else mutated)
    p = subprocess.run(["Rscript", "--vanilla", "R/testing/test_S2.5.R"], env=dict(env, S25_RSRC=d), capture_output=True, text=True)
    fails = [l for l in p.stdout.splitlines() if l.startswith("FAIL")] + ([f"ERRORE R: {p.stderr.strip().splitlines()[-1]}"] if p.returncode != 0 and not any(l.startswith("FAIL") for l in p.stdout.splitlines()) and p.stderr.strip() else [])
    det = p.returncode != 0
    ok_all &= det
    lines.append(f"{'RILEVATO' if det else 'NON RILEVATO'} {name}: {len(fails)} asserzioni fallite; {'; '.join(f[:90] for f in fails[:3])}")
    shutil.rmtree(d)
lines.append(f"TOTALE: {sum(l.startswith('RILEVATO') for l in lines)}/{len(MUT)} rilevati")
open(out, "w").write("\n".join(lines) + "\n"); print("\n".join(lines))
sys.exit(0 if ok_all else 1)
