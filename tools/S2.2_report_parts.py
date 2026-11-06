"""S2.2: estrae le sezioni del file dello step in frammenti per il report (il testo esiste una volta sola, in docs/steps/S2.2.md).
Uso: python3 tools/S2.2_report_parts.py"""
import re
t = open("docs/steps/S2.2.md").read().split("\n---\n", 1)[1]
t = re.sub(r"\]\(/(decisions|steps)/", r"](../docs/\1/", t)
sec = {}; cur = None
for line in t.split("\n"):
    m = re.match(r"^# (.+)$", line)
    if m: cur = m.group(1); sec[cur] = []; continue
    if cur: sec[cur].append(line)
def w(name, keys, demote=True):
    out = []
    for k in keys:
        body = "\n".join(sec[k])
        if demote: body = re.sub(r"^## ", "### ", body, flags=re.M)
        out.append((f"## {k}\n" if len(keys) > 1 else "") + body)
    open(f"reports/_S2.2_{name}.md", "w").write("\n\n".join(out))
inb = "\n".join(sec["In breve"]).replace("## In breve:", "# In breve:")
open("reports/_S2.2_inbreve.md", "w").write(inb)
w("rinviati", ["Dubbi rinviati"]); w("verifiche", ["Verifiche di Luca"]); w("esiti", ["Esiti"]); w("posthoc", ["Analisi post hoc (dichiarate dopo i risultati; non cambiano gli esiti)"])
w("review", ["Revisione avversariale"]); w("deviazioni", ["Deviazioni dichiarate"]); w("prereg", ["Perimetro", "Dubbi ricevuti", "Attese pre-registrate", "Correzioni prima del calcolo (dichiarate, commit separato)", "Esplorazione dichiarata (prima della pre-registrazione, solo dati veri)"])
print(sorted(sec))
