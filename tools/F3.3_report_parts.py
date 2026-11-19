"""tools/F3.3_report_parts.py - F3.3: parti del report estratte dal file dello step (il testo esiste una volta sola) e tabelle per tessuto.
Uso: python3 tools/F3.3_report_parts.py"""
import re, pandas as pd
def md(df): return "| " + " | ".join(map(str, df.columns)) + " |\n|" + "---|" * df.shape[1] + "\n" + "\n".join("| " + " | ".join(map(str, r)) + " |" for r in df.itertuples(index=False)) + "\n"
s = open("docs/steps/F3.3.md").read(); O = "reports/F3.3_parts"
def section(title, level=1):
    h = "#" * level + " "; m = re.search(rf"^{re.escape(h + title)}\n(.*?)(?=^#{{1,{level}}} )", s, re.S | re.M); return m.group(1).strip() + "\n"
inb = section("In breve"); open(f"{O}/in_breve.md", "w").write(inb.replace("## In breve", "# In breve", 1))
open(f"{O}/dubbi.md", "w").write(section("Dubbi rinviati"))
open(f"{O}/esiti.md", "w").write(section("Esiti").replace("\n## ", "\n### "))
open(f"{O}/correzioni.md", "w").write(section("Correzioni (dopo il calcolo, prima del verdetto)"))
open(f"{O}/revisione.md", "w").write(section("Revisione avversariale"))
pre = section("Attese pre-registrate").replace("\n## ", "\n### ")
open(f"{O}/prereg.md", "w").write(section("Correzione dopo il cancello G-1, prima di leggere i dati di verifica (2026-11-19)", 2))
# tabella per tessuto e insieme: scarto dentro e fuori, spostamento, esiti
C = pd.read_csv("results/F3.3/F3.3_cells.csv")
TN = {"A1": "intestino", "A2": "tumore del colon", "A3": "stroma", "A4": "linfonodo (500 µm)", "A5": "cervello", "A6": "cuore"}
C["tessuto"] = C.A.map(TN); C["insieme"] = C.set.str.replace(r"_L[12]_A\d$", "", regex=True)
cols = ["tessuto", "measure", "insieme", "level", "n_out", "h", "gap_in_over_h", "verdict_in", "gap_out_over_h", "verdict_out", "shift_over_h", "verdict_shift"]
T = C[cols].copy()
for c in ("gap_in_over_h", "gap_out_over_h", "shift_over_h"): T[c] = T[c].map(lambda v: f"{v:+.2f}")
T["h"] = T.h.map(lambda v: f"{v:.4g}")
T.columns = ["tessuto", "misura", "insieme", "livello", "aree", "margine h", "scarto dentro / h", "esito dentro", "scarto fuori / h", "esito fuori", "spostamento / h", "esito spostamento"]
open(f"{O}/cells_table.md", "w").write(md(T))
V = pd.read_csv("results/F3.3/F3.3_varcomp.csv"); V["tessuto"] = V.A.map(TN)
V = V[["tessuto", "measure", "groups", "n", "sigma_patient", "sigma_area", "rho", "lo", "hi", "verdict"]].round(3)
open(f"{O}/varcomp_table.md", "w").write(md(V))
print("ok")
