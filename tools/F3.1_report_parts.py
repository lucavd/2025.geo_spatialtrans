"""tools/F3.1_report_parts.py - F3.1: estrae le sezioni di docs/steps/F3.1.md in reports/F3.1_parts/*.md per il report (il testo esiste una volta sola)."""
import re, os
s = open("docs/steps/F3.1.md").read(); s = s.split("\n---\n", 1)[1] if s.startswith("---") else s
os.makedirs("reports/F3.1_parts", exist_ok=True)
secs = {m.group(1): m.start() for m in re.finditer(r"^# (.+)$", s, flags=re.M)}
order = sorted(secs.items(), key=lambda kv: kv[1])
def body(title):
    i = [k for k, _ in order].index(title); a = order[i][1]; b = order[i + 1][1] if i + 1 < len(order) else len(s)
    t = s[a:b].split("\n", 1)[1]
    return re.sub(r"\(/(decisions|steps)/", r"(../docs/\1/", t)
NAMES = {"In breve": "in_breve", "Dubbi rinviati": "dubbi", "Attese pre-registrate": "attese", "Esiti": "esiti",
         "Correzioni (dopo il calcolo, prima del verdetto)": "correzioni", "Analisi post hoc (dichiarate dopo i risultati; non cambiano gli esiti)": "posthoc",
         "Revisione avversariale": "revisione", "Deviazioni dichiarate": "deviazioni", "Esplorazioni dichiarate (prima della pre-registrazione)": "esplorazioni",
         "Perimetro": "perimetro", "Metodi e misure": "metodi", "Riproduzione": "riproduzione"}
for t, n in NAMES.items():
    # i titoli di livello 2 dentro una sezione diventano di livello 3 nel report
    open(f"reports/F3.1_parts/{n}.md", "w").write(re.sub(r"^## ", "### ", body(t), flags=re.M).replace("\n# ", "\n## "))
print("parts:", len(NAMES))
