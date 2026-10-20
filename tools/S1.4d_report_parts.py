#!/usr/bin/env python3
"""tools/S1.4d_report_parts.py — estrae dal file dello step (docs/steps/S1.4d.md) le sezioni incluse nel report, cosi' il testo
esiste una volta sola. Uso (root del repo): python3 -B tools/S1.4d_report_parts.py"""
import pathlib, re
t = pathlib.Path("docs/steps/S1.4d.md").read_text(encoding="utf-8").split("\n---\n", 1)[1]
def sec(start, stop_re):
    i = t.index(start); m = re.search(stop_re, t[i + len(start):]); j = i + len(start) + (m.start() if m else len(t))
    return t[i:j].strip() + "\n"
inb = sec("\n# In breve\n", r"\n# ")
rin = sec("\n# Dubbi rinviati\n", r"\n# ").replace("# Dubbi rinviati\n", "", 1)
rev = sec("\n## Revisione avversariale (cancello G-5)\n", r"\n#{1,2} ").replace("## Revisione avversariale (cancello G-5)\n", "", 1)
post = sec("\n## Dopo il risultato e la revisione avversariale", r"\n#{1,2} ")
out = pathlib.Path("reports")
(out / "_S1.4d_inbreve.md").write_text(inb, encoding="utf-8")
ver = sec("\n# Verifiche di Luca\n", r"\n# ").replace("# Verifiche di Luca\n", "", 1)
(out / "_S1.4d_rinviati.md").write_text(rin, encoding="utf-8")
(out / "_S1.4d_verifiche.md").write_text(ver, encoding="utf-8")
(out / "_S1.4d_review.md").write_text(rev + "\n" + post.replace("## Dopo il risultato", "## Note dopo il risultato", 1), encoding="utf-8")
print("ok", len(inb), len(rin), len(rev), len(post))
