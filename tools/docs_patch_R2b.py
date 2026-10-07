"""tools/docs_patch_R2b.py — applica gli aggiornamenti documentali di chiusura R2b (idempotente: salta se il marcatore e' gia' presente).
Uso: python3 tools/docs_patch_R2b.py <bio_rows.md> <sha_chiusura|TBD>
"""
import sys, re
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
bio_rows = Path(sys.argv[1]).read_text().strip(); sha = sys.argv[2] if len(sys.argv) > 2 else "TBD"

def patch(path, marker, fn):
    p = ROOT / path; s = p.read_text()
    if marker in s: print("skip", path); return
    p.write_text(fn(s)); print("patched", path)

# ---------------- BIO_REFERENCES: nuove righe B-043..B-048 + stato + marcatura "sostituito" sulle righe di densita' R2
def bio(s):
    lines = s.split("\n"); i = max(k for k, l in enumerate(lines) if l.startswith("| B-042 |")); lines.insert(i + 1, bio_rows); s = "\n".join(lines)
    for bid in ("B-001", "B-008", "B-015", "B-022", "B-029", "B-036"):
        s = re.sub(rf"(\| {bid} \| A\d \| densità nucleare \(consenso\) \|[^\n]*\| )\*\*provvisorio\*\* \(BL-024\)", r"\1**superato da R2b** (vedi B-043…B-048; consenso R2 fuori dall'intervallo manuale in A2, A3, A4)", s)
    s = s.replace("## Stato\n", "## Stato\n- **R2b (2026-09-21 → 2026-10-07)**: densità nucleare per archetipo da **conteggio manuale** (Luca 24 finestre/2 230 nuclei, F. Pezzuto 12/1 106, alla cieca) riportata come **intervallo** [nuclei evidenti, totale] (B-043…B-048): i segmentatori perdono una classe di nuclei pallidi (OD ematossilina 18–60 % degli appaiati) che gli umani contano; il consenso R2 non è un limite. Le righe di densità R2 (B-001, 008, 015, 022, 029, 036) sono superate; area nucleare (B-002…) ancora provvisoria (BL-025). Dettagli in `reports/R2b.html`.\n", 1)
    return s
patch("docs/BIO_REFERENCES.md", "B-043", bio)

# ---------------- CHANGELOG
chg = f"""## R2b — 2026-09-21 → 2026-10-07 — Verità manuale per la densità nucleare (branch `step1-cell-layer`)

### Codice (tutto in `tools/`; dati derivati fuori git su `/mnt/micron/geo_spatialtrans/R2b/`)
- `R2b_windows.py` (24 finestre seedate a lato adattato, C-R2b.1/4), `R2b_annotator.py` + `R2b_annotator_core.js` (annotatore HTML autonomo alla cieca, v1.1 solo conteggio; core testato con node in `R/testing/test_R2b_core.js`), `R2b_match.py` (appaiamento punto↔oggetto 1:1, inter-rater, sommari; `--selftest` C-R2b.3), `R2b_import_model.py` (VLM → formato annotatore), `R2b_run_vlm.py` (chiamate API OpenAI-compatibili, risposte grezze conservate), `R2b_null_and_figures.py` (CP-R2b.4 nullo dei punti casuali, sovrapposizioni), `R2b_od_strata.py` (CP-R2b.5), `R2b_sheet.py`, `docs_patch_R2b.py`. `R/testing/test_R2b.sh` → `test_R2b.py` (PASS/FAIL dalle tabelle).
- Tracciati in `results/R2b/`: finestre, bersagli, pre-registrazione (con 6 aggiunte datate), annotazioni (Luca, Federica, 3 VLM) e risposte grezze dei modelli, tutte le tabelle, figure, `R2b_test_results.csv`, `R2b_adversarial_review.csv`, checklist della patologa.

### Esiti
- Check C 5 PASS + 1 WARN (C-R2b.1b: maschera bolle R2 parziale, BL-040). Check B: B-R2b.1, 2, 3a, 3b **FAIL**, B-R2b.4 PASS, B-R2b.5 non eseguito. Controprove: CP-R2b.1a WARN (Luca–Federica Δ 2.3 % sul totale, ma 6/12 finestre e 4/6 archetipi entro 10 %; A4 13 %), CP-R2b.1b FAIL (F1 3 µm 0.76), CP-R2b.2 FAIL (13–29 % dei nuclei umani persi da Cellpose **e** StarDist), CP-R2b.3 FAIL (Cellpose in A5/A6: recall basso, non precisione bassa → BL-032/033 riformulati), CP-R2b.4 PASS, **CP-R2b.5 PASS 6/6** (i nuclei persi sono più pallidi ma sopra il fondo: classe ambigua reale, tesi di Luca).
- VLM: GPT-6 Astra localizza (F1 5 µm 0.59–0.92) ma sovraconta 10–47 %; deepseek-flash non è un annotatore (thinking: non ripetibile, ×1.2–3.1; no-thinking: conteggi ricorrenti, posizioni al nullo).
- `docs/BIO_REFERENCES.md` v2: B-043…B-048 (densità come intervallo [evidenti, totale] da due annotatori); B-001/008/015/022/029/036 superate.
- Deviazioni documentate: finestre a lato adattato (non 100 µm); annotatore ridotto a solo conteggio su richiesta di Luca; area nucleare rinviata; doppia corsa DeepSeek-thinking (usata come ripetibilità); A3_w2 fuori archetipo per giudizio di Claude (la patologa non l'ha valutata).
- Chiusura: `{sha}`.

"""
patch("CHANGELOG.md", "## R2b —", lambda s: s.replace("## R2 — 2026-09-20/21", chg + "## R2 — 2026-09-20/21", 1))

# ---------------- BACKLOG
bl = """| BL-034 | R2b | Accordo inter-annotatore umano sulla **posizione**: F1 punto-punto a 3 µm 0.76 (0.84 a 5 µm) fra Luca e Federica a parità di conteggio (Δ 2.3 %). La soglia pre-registrata (0.90 a 3 µm) era pensata per maschere, non per due click su nuclei di 6–11 µm. Per appaiamenti umano↔umano usare 5 µm o una soglia proporzionale al diametro nucleare dell'archetipo. | S1.2 (metriche di confronto) | aperta |
| BL-035 | R2b | **Classe di nuclei pallidi** (13–29 % dei nuclei umani in A1–A3, A6; OD ematossilina 18–60 % degli appaiati) persa da Cellpose e StarDist. Per S1.x la densità del preset va trattata come intervallo [evidenti, totale] (B-043…B-048) e il simulatore deve poter generare la frazione "pallida" come parametro. Per R3 (Voronoi sui nuclei reali) i centroidi dei soli segmentatori sottostimano la densità del 5–25 %: valutare l'aggiunta dei punti manuali o un fattore correttivo per archetipo (totale/evidenti: A1 1.21, A2 1.21, A3 1.34, A4 1.05, A5 1.24, A6 1.09). | R3, S1.2 | aperta |
| BL-036 | R2b | A3_w2 (ROI r3) è epitelio ghiandolare tumorale, non stroma: giudizio di Claude dalla sovrapposizione; la checklist della patologa l'ha lasciata non valutata (12 finestre «sì», 12 vuote, probabile sfasamento di riga). Chiedere a F. Pezzuto il giudizio su A3_w2 e sulle 11 non valutate; A3 riportato con/senza (3 100 / 2 593 /mm²). Eredita BL-031 (ROI misti). | R3 | aperta |
| BL-037 | R2b | Federica ha dichiarato di aver saltato le bolle in A6 ma non ha disegnato esclusioni: la sua densità A6 è su area lorda (−2…−5 % rispetto all'area netta di Luca). In futuri annotatori: esclusioni obbligatorie quando si dichiara di saltare zone. | — | aperta |
| BL-038 | R2b | VLM come annotatori: GPT-6 Astra (eseguito da Luca, versione/data da registrare) sovraconta 10–47 % ma localizza; deepseek-flash non è utilizzabile. Possibile uso di un VLM con visione come **arbitro sì/no** su ritagli 8×8 µm centrati sugli oggetti esclusivi (alternativa all'arbitro OD di R2); Jev (TypeSafe) è solo testo e non è adatto. | R3 (arbitro) | aperta |
| BL-040 | R2b (revisore) | La maschera bolle R2 (`bubble_mask`, BL-026) dichiara 0 % bolle nelle 24 finestre, ma l'annotatore ha escluso bolle sfocate in A6_w2/w3/w4 (1.2–4.9 % dell'area): il rilevatore non vede le bolle senza anello netto. C-R2b.1b → WARN. Migliorare (Hough su gradiente radiale, o soglia su bassa cromaticità + bassa varianza locale) prima di R3 su A6. | R3 | aperta |
| BL-041 | R2b (revisore) | Con 4 finestre per archetipo l'IC t sulla densità è 3–5 volte più largo dell'IC di Poisson (A1 7 100–17 700 vs 11 350–13 560): il consenso R2 cade dentro l'IC t in 5/6 archetipi (fuori solo in A2): B-R2b.2 ha potenza per singolo archetipo solo in A2. Per un riferimento con ±15 % servono ≥ 10–12 finestre per archetipo (stima dalla SD fra finestre), oppure finestre stratificate per sotto-regione (follicolo/paracorticale in A4; villo/cripta in A1). | R2c / R3 | aperta |
| BL-039 | R2b | Area nucleare manuale (B-R2b.5) non eseguita: la fase 2 con calibri (5 nuclei/finestra dai bersagli seedati in `R2b_targets.csv`) resta da fare per chiudere BL-025. | R2c o S1.4 | aperta |
"""
def backlog(s):
    s = s.rstrip("\n") + "\n" + bl
    s = s.replace("| R2b (prima di S1.2) | aperta |", "| R2b (prima di S1.2) | chiusa (R2b: intervallo manuale B-043…B-048; residuo in BL-035) |", 1)   # BL-024
    s = re.sub(r"(\| BL-030 \|[^\n]*\| R2b \| )aperta \|", r"\1chiusa (R2b: A4 28 096, A5 1 185 /mm² manuali; attese R2 smentite confermate) |", s)
    s = re.sub(r"(\| BL-032 \|[^\n]*\| R2b \| )aperta \|", r"\1riformulata (R2b: precisione Cellpose A5 0.89 vs nullo 0.11 → gli esclusivi non sono fantasmi; il problema è il recall 0.77) |", s)
    s = re.sub(r"(\| BL-033 \|[^\n]*\| R2b / R3 \| )aperta \|", r"\1riformulata (R2b: Cellpose A6 recall 0.42, precisione 0.82: perde nuclei, non ne inventa; per A6 usare Space Ranger recall 0.89 o StarDist 0.71) |", s)
    return s
patch("docs/BACKLOG.md", "BL-034", backlog)

# ---------------- ROADMAP: stato sessione, decisioni, prossima
def roadmap(s):
    s = s.replace("| R2b | — | proposta | — | conteggio/contorno manuale (Luca/Federica) su finestre 100×100 µm per archetipo → fissa gli intervalli di BIO_REFERENCES e sceglie il segmentatore per R3; dipende da R2 |",
        "| R2b | 2026-09-21 → 10-07 | **chiusa con riserve** | `reports/R2b.html` | 24 finestre a lato adattato (4/archetipo), Luca 2 230 nuclei + F. Pezzuto 1 106 (alla cieca); 3 segmentatori + 3 VLM appaiati. Check C 5 PASS 1 WARN; B: 1 PASS 4 FAIL; CP: 2 PASS 1 WARN 5 FAIL (attese pre-registrate smentite dai dati, non errori). Risultato: densità = **intervallo** [evidenti, totale] (B-043…B-048); i segmentatori perdono una classe di nuclei pallidi reale (CP-R2b.5); accordo umano 2.3 % sul conteggio, F1 posizione 0.76. Riserve: BL-034…BL-041; revisione avversariale in `results/R2b/R2b_adversarial_review.csv` |")
    dec = """| 2026-10-07 | (R2b) La densità nucleare di riferimento per archetipo è un **intervallo** [nuclei evidenti, nuclei totali] da conteggio manuale (B-043…B-048), non il consenso dei segmentatori; i preset S1.x devono poter generare la frazione di nuclei pallidi come parametro (BL-035) |
| 2026-10-07 | (R2b) Segmentatore per R3: **Cellpose-SAM RGB** in A1, A2, A3, A5; **Space Ranger 4** (StarDist in subordine) in A4 e A6; in ogni caso con fattore correttivo totale/evidenti per archetipo o integrazione dei punti manuali (BL-035) |
| 2026-10-07 | (R2b) Annotazione manuale: un solo compito per volta (solo conteggio), alla cieca, finestre a lato adattato alla densità (70–110 nuclei); giudizi istologici separati in checklist per la patologa |
| 2026-10-07 | (R2b) I VLM generalisti (GPT-6 Astra, deepseek-flash) **non** sono annotatori di riferimento per densità/posizione; al più arbitri sì/no (BL-038) |
| 2026-09-18 | **Una sessione = una chat.**"""
    s = s.replace("| 2026-09-18 | **Una sessione = una chat.**", dec, 1)
    s = re.sub(r"\*\*Prossima sessione\*\*: \*\*R2b\*\*[^\n]*", "**Prossima sessione**: **S1.1** (`extract_regions()`, nessuna dipendenza) oppure **R3** (Voronoi sui nuclei reali: usa le maschere R2 con il fattore correttivo di BL-035; dipende da R2/R2b chiuse). Area nucleare manuale (BL-039) può essere una breve R2c prima di S1.4.", s)
    s = s.replace("Dopo R2 (2026-09-21):", "Dopo R2b (2026-10-07): densità manuale di riferimento per A1–A6 come intervallo (B-043…B-048), annotatore HTML riusabile (`tools/R2b_annotator.py`), 24 finestre con verità a punti su `/mnt/micron/geo_spatialtrans/R2b/`. Dopo R2 (2026-09-21):", 1)
    return s
patch("docs/ROADMAP.md", "| R2b | 2026-09-21", roadmap)
