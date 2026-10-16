#!/usr/bin/env python3
# Aggiornamento dei documenti di coordinamento a chiusura di S1.4c (eseguito una volta; idempotente sui marcatori).
import io, re, sys
def rd(p): return open(p, encoding="utf-8").read()
def wr(p, s): open(p, "w", encoding="utf-8").write(s)
BIO = ['| B-121 | A1 | contatto fra nuclei reali in maschera (Chebyshev ≤ 2 px): frazione di nuclei a contatto cellpose_rgb / stardist_he; solidità mediana; AUC sovrapposizione dei cerchi di pari area → contatto; q99 di N/C (Voronoi dei centroidi) | 0.425 / 0.251 | per ROI cellpose 0.400–0.439, StarDist 0.232–0.274; solidità 0.946 / 0.955; segmentatore di forma (ev) 0.436; AUC CP-0 0.960; q99 N/C 0.757 (n = 33203) | frazione | tools/S1.4c_bio.R + tools/S1.4c_raster.py@S1.4c + maschere native R2, territori S1.4 | S1.4c | dipende dal segmentatore (CP-2 S1.4c: S-1, S-2, S-3 declassate); non usare come bersaglio senza verità di riferimento (BL-120) |', '| B-122 | A2 | contatto fra nuclei reali in maschera (Chebyshev ≤ 2 px): frazione di nuclei a contatto cellpose_rgb / stardist_he; solidità mediana; AUC sovrapposizione dei cerchi di pari area → contatto; q99 di N/C (Voronoi dei centroidi) | 0.615 / 0.320 | per ROI cellpose 0.535–0.762, StarDist 0.271–0.408; solidità 0.954 / 0.961; segmentatore di forma (ev) 0.626; AUC CP-0 0.947; q99 N/C 1.051 (n = 29948) | frazione | tools/S1.4c_bio.R + tools/S1.4c_raster.py@S1.4c + maschere native R2, territori S1.4 | S1.4c | dipende dal segmentatore (CP-2 S1.4c: S-1, S-2, S-3 declassate); non usare come bersaglio senza verità di riferimento (BL-120) |', '| B-123 | A3 | contatto fra nuclei reali in maschera (Chebyshev ≤ 2 px): frazione di nuclei a contatto cellpose_rgb / stardist_he; solidità mediana; AUC sovrapposizione dei cerchi di pari area → contatto; q99 di N/C (Voronoi dei centroidi) | 0.388 / 0.226 | per ROI cellpose 0.207–0.613, StarDist 0.142–0.316; solidità 0.954 / 0.961; segmentatore di forma (ev) 0.397; AUC CP-0 0.967; q99 N/C 0.947 (n = 20018) | frazione | tools/S1.4c_bio.R + tools/S1.4c_raster.py@S1.4c + maschere native R2, territori S1.4 | S1.4c | dipende dal segmentatore (CP-2 S1.4c: S-1, S-2, S-3 declassate); non usare come bersaglio senza verità di riferimento (BL-120) |', '| B-124 | A4 | contatto fra nuclei reali in maschera (Chebyshev ≤ 2 px): frazione di nuclei a contatto cellpose_rgb / stardist_he; solidità mediana; AUC sovrapposizione dei cerchi di pari area → contatto; q99 di N/C (Voronoi dei centroidi) | 0.346 / 0.417 | per ROI cellpose 0.325–0.371, StarDist 0.371–0.461; solidità 0.956 / 0.958; segmentatore di forma (ev) 0.419; AUC CP-0 0.974; q99 N/C 0.816 (n = 99600) | frazione | tools/S1.4c_bio.R + tools/S1.4c_raster.py@S1.4c + maschere native R2, territori S1.4 | S1.4c | dipende dal segmentatore (CP-2 S1.4c: S-1, S-2, S-3 declassate); non usare come bersaglio senza verità di riferimento (BL-120); simulati sui generatori Space Ranger, n_SR/n_StarDist 1.15–1.40 (BL-121) |', '| B-125 | A5 | contatto fra nuclei reali in maschera (Chebyshev ≤ 2 px): frazione di nuclei a contatto cellpose_rgb / stardist_he; solidità mediana; AUC sovrapposizione dei cerchi di pari area → contatto; q99 di N/C (Voronoi dei centroidi) | 0.102 / 0.085 | per ROI cellpose 0.047–0.203, StarDist 0.065–0.108; solidità 0.975 / 0.971; segmentatore di forma (ev) 0.107; AUC CP-0 0.994; q99 N/C 0.474 (n = 1993) | frazione | tools/S1.4c_bio.R + tools/S1.4c_raster.py@S1.4c + maschere native R2, territori S1.4 | S1.4c | dipende dal segmentatore (CP-2 S1.4c: S-1, S-2, S-3 declassate); non usare come bersaglio senza verità di riferimento (BL-120) |', '| B-126 | A6 | contatto fra nuclei reali in maschera (Chebyshev ≤ 2 px): frazione di nuclei a contatto cellpose_rgb / stardist_he; solidità mediana; AUC sovrapposizione dei cerchi di pari area → contatto; q99 di N/C (Voronoi dei centroidi) | 0.060 / 0.063 | per ROI cellpose 0.014–0.166, StarDist 0.044–0.074; solidità 0.946 / 0.953; segmentatore di forma (ev) 0.045; AUC CP-0 0.997; q99 N/C 0.180 (n = 2221) | frazione | tools/S1.4c_bio.R + tools/S1.4c_raster.py@S1.4c + maschere native R2, territori S1.4 | S1.4c | dipende dal segmentatore (CP-2 S1.4c: S-1, S-2, S-3 declassate); non usare come bersaglio senza verità di riferimento (BL-120); simulati sui generatori Space Ranger, n_SR/n_StarDist 1.15–1.40 (BL-121) |']
# ---- ROADMAP
p = "docs/ROADMAP.md"; s = rd(p)
if "(S1.4c, BL-098)" not in s:
    anchor = "| 2026-09-18 | **Una sessione = una chat.**"
    dec = ("| 2026-10-15 | (S1.4c) Perimetro: banco + decisione; implementazione in `04b4` nella sessione dopo (S1.4d) — Luca («Banco + decisione; implementazione in 04b4 nella sessione dopo») |\n"
           "| 2026-10-15 | (S1.4c) Varianti del banco: K-c, K-1, K-2 nel Voronoi (territori invariati, decisione 1) — Luca («K-c, K-1, K-2 nel Voronoi») |\n"
           "| 2026-10-16 | (S1.4c, BL-098) **Forma del nucleo = K-2**: 32-gono di area M2 al generatore, ritagliato nel proprio territorio di Voronoi rientrato di 1e-5 µm e gonfiato per bisezione fino all'area target; tetto φ = min(1, q99 di N/C reale) per archetipo. Contenimento mantenuto (C3, design §7.2 invariato); sostituisce lo scalamento di D-S1.4.2 da S1.4d. Per regola pre-registrata: pannello {B-1, B-2, B-7, B-8} (S-1/S-2/S-3 declassate da CP-2), punteggi K-c 3, K-1 7, K-2 8, P(bootstrap) 1.0; riserve: eleggibilità di K-2 dopo la correzione D-6, 1 punto su K-1 — Luca («Nessun veto: K-2 per S1.4d, fragilità tra le riserve») |\n")
    assert anchor in s; s = s.replace(anchor, dec + anchor, 1)
if "| S1.4c |" not in s:
    m = re.search(r"^\| S1\.4b \|.*$", s, re.M); assert m
    row = ("| S1.4c | 2026-10-15/16 | **chiusa con riserve** | `reports/S1.4c.html` | BL-098 riformulata. Banco K-c (scala) / K-1 (ritaglio) / K-2 (ritaglio + pari area) sui 30 tassellamenti primari × 10 semi, "
           "operatore unico (nuclei simulati rasterizzati sulla griglia nativa, contatto ≤ 2 px, Voronoi dei centroidi osservati). Check C: C-0 330/330, C3 0, C-NO 0, C-area 5.8e-7, C-OP (maschere = R3, scarto 0), test 24 + 17 + 15, mutanti 5/5; "
           "**C-perf WARN** 622.7 s. Check B: K-c non produce contatti (B-7/B-8 FAIL 6/6), K-2 B-4 PASS 6/6 (1.5–4.7 %). CP-0: premessa non smentita (AUC 0.95–0.97 A1–A4); **CP-1 smentita al limite** (3 ≥ 3); **CP-2: S-1/S-2/S-3 dipendono dal segmentatore → descrittive**. "
           "Decisione per regola: **K-2** (8 vs 7 vs 3; dopo la correzione D-6 dell'eleggibilità). Revisione: 2 revisori, 0 discrepanze, 0 bug nella pipeline. Riserve BL-118…BL-127 |")
    s = s[:m.end()] + "\n" + row + s[m.end():]
old = "**Prossima sessione** (aggiornata a S1.4b): "
if old in s:
    new = ("**Prossima sessione** (aggiornata a S1.4c): **S1.4d = implementazione di K-2 in `derive_cell_geometry()`** — modo di forma «ritaglio + pari area» (rientro 1e-5 µm, tetto φ per archetipo nel preset da B-121…B-126, "
           "bisezione accelerata BL-122, trappole GEOS BL-118), test C3/C-NO/C-area, check B di regressione contro i risultati del prototipo S1.4c; dipendenze: S1.4c chiusa con riserve, nessun dato nuovo. "
           "Alternative: BL-120 (verità di riferimento per i contatti), BL-087 (profilo radiale di RNA), poi S1.5. La scelta è di Luca. Testo precedente, superato: ")
    s = s.replace(old, new, 1)
wr(p, s)
# ---- BACKLOG
p = "docs/BACKLOG.md"; s = rd(p)
def set_status(s, bl, status, extra=None):
    m = re.search(r"^\| %s \|.*$" % bl, s, re.M); assert m, bl
    line = m.group(0); cells = line.split(" | ")
    cells[-1] = status + " |"
    if extra: cells[2] = cells[2] + " " + extra
    return s[:m.start()] + " | ".join(cells) + s[m.end():]
if "S1.4c: decisa K-2" not in s:
    s = set_status(s, "BL-098", "decisa in S1.4c (K-2), implementazione in S1.4d", "**S1.4c: decisa K-2** (ritaglio nel territorio + pari area; B-4 1.5–4.7 %, contatti riprodotti in B-7/B-8; regola pre-registrata, nessun veto di Luca). Resta aperta la forma fine (BL-119, BL-123).")
    s = set_status(s, "BL-111", "superata in S1.4c", "S1.4c: con il ritaglio nel territorio i nuclei sono disgiunti per costruzione (nessun vincolo di scala fra nuclei).")
    s = set_status(s, "BL-112", "risolta nel banco S1.4c; da fissare nel design in S1.4d", "S1.4c: operatore unico = nuclei rasterizzati, centroidi dei pixel, Voronoi dei centroidi osservati con `tessellate_voronoi()` (C-OP: sui reali riproduce S1.4 a 1e-9).")
    new = [
    "| BL-118 | S1.4c (D-1, D-2) | Trappole GEOS (pacchetto `geos` 0.2.5 / libgeos 3.11): (a) il ritaglio con `geos_intersection_prec` su griglia 1e-6 fa uscire i vertici nuovi dal territorio fino a ~1e-6 µm (C3 36/~1 000) → ritagliare contro il territorio rientrato di 1e-5 µm; (b) `geos_buffer(-ε)` su geometrie con modello di precisione restituisce **in silenzio** poligoni vuoti (23/674) → buffer con precisione 0 e `stop()` sui vuoti. | S1.4d / S1.5 (helper F1) | aperta |",
    "| BL-119 | S1.4c (D-1, fig. 3) | La faccia di contatto reale non cade sulla bisettrice dei centroidi: nuclei con > 5 % di pixel fuori dal Voronoi dei centroidi 0.15–0.28 nei reali contro 0.03–0.09 con K-2 (A1–A3). Candidato: K-2 nel diagramma di potenza (BL-117) o piano di contatto pesato dalla taglia. | dopo S1.4d | aperta |",
    "| BL-120 | S1.4c (CP-2) | I riferimenti di contatto e solidità dipendono dal segmentatore (S-1 cellpose / StarDist: A1 0.42 / 0.25, A2 0.61 / 0.32, A3 0.39 / 0.23; StarDist più convesso in 5/6): S-1, S-2, S-3 declassate. Serve una verità di riferimento (coppie a contatto annotate a mano sulle finestre di R2b). | R2c | aperta |",
    "| BL-121 | S1.4c (CODE-7, fig. 3) | A4/A6: simulati sui generatori Space Ranger (n_SR / n_StarDist 1.15–1.40), vicini reali StarDist; S-1 A4 simulato 0.75 contro 0.42. Con l'operatore sulla sola maschera StarDist B-8 reale A4 0.335 → 0.148 e B-7 0.397 → 0.188 (punteggi invariati, B-7 K-2 A4 WARN → FAIL). Fissare un'unica maschera per posizioni, aree e vicini in A4/A6. Legata a BL-100, BL-115. | R3b / R2c | aperta |",
    "| BL-122 | S1.4c (C-perf) | K-2 su S (449 536 celle): 622.7 s, 5.6 GB (WARN), fino a 53 intersezioni GEOS per cella ritagliata. Prima dell'integrazione: punto di partenza analitico (area disco ∩ poligono), bisezione solo sulle celle ritagliate, tolleranza relativa. | S1.4d | aperta |",
    "| BL-123 | S1.4c (fig. 1, CODE-6) | Forma: i nuclei simulati sono dischi ritagliati, i reali sono allungati (ecc. mediana 0.56–0.81, B-115…B-120) e meno convessi (solidità 0.946–0.975 contro 0.960–0.976); la solidità di skimage a ~200 px non distingue un disco (0.961 ± 0.011). Valutare ellissi ritagliate (BL-063) e una metrica di forma sub-pixel. | dopo S1.4d | aperta |",
    "| BL-124 | S1.4c (RC-10, CODE-8) | Lezione di metodo: lo score di CP-0 è stato definito nelle dev notes come (r_i + r_j) / d invece della differenza pre-registrata, senza marcarlo come deviazione (esito invariato: AUC 0.942–0.975). Ogni ridefinizione rispetto alla pre-registrazione va dichiarata come deviazione. | prossime pre-registrazioni | aperta |",
    "| BL-125 | S1.4c (D-6, CODE-9/10/11, RC-05/12/20) | Lezione di metodo: (a) l'eleggibilità va testata prima del run anche su un caso con C-perf WARN (D-6 corretta dopo i risultati ha capovolto K-1 → K-2); (b) pre-registrare un criterio fra le varianti candidate (K-2 su K-1 di 1 punto, P 0.93); (c) una controprova che può essere smentita all'uguaglianza (CP-1, 3 ≥ 3) e su metriche poi declassate va formulata sul pannello effettivo. | prossime pre-registrazioni | aperta |",
    "| BL-126 | S1.4c (CODE-13) | I test di `tools/S1.4c_decision.R` non rilevano 3 mutanti (bootstrap senza reinserimento, bootstrap non appaiato, metro con la media): aggiungere i test prima di riusare le funzioni (effetto sui dati di S1.4c nullo). | prossima decisione per regola | aperta |",
    "| BL-127 | S1.4c (RC-07, φ) | B-4 di K-2 in A2 0.047 (5 fold singoli > 0.05; PASS dipende dalla media sui fold). In A2 φ = 1 (q99 di N/C reale 1.05): K-2 può riempire il territorio (citoplasma nullo). Fissare φ per archetipo nel preset e riportare la frazione con N/C ≥ 0.95 in S1.4d. | S1.4d | aperta |"]
    s = s.rstrip("\n") + "\n" + "\n".join(new) + "\n"
wr(p, s)
# ---- BIO_REFERENCES
p = "docs/BIO_REFERENCES.md"; s = rd(p)
if "| B-121 |" not in s:
    m = re.search(r"^\| B-120 \|.*$", s, re.M); assert m
    s = s[:m.end()] + "\n" + "\n".join(BIO) + s[m.end():]
    st = "## Stato\n"; assert st in s
    s = s.replace(st, st + "- **S1.4c (2026-10-15/16)**: B-121…B-126 = contatto fra nuclei reali misurato sulle maschere native con l'operatore di S1.4c (frazione a contatto con due segmentatori, solidità, AUC di CP-0, q99 di N/C usato come tetto di K-2). **Dipendono dal segmentatore** (CP-2): valori descrittivi finché non c'è una verità di riferimento (BL-120).\n", 1)
wr(p, s)
# ---- CHANGELOG
p = "CHANGELOG.md"; s = rd(p)
if "## S1.4c" not in s:
    entry = open("/tmp/s14c_changelog.md", encoding="utf-8").read()
    i = s.index("## S1.4b"); s = s[:i] + entry + "\n" + s[i:]
wr(p, s)
print("docs aggiornati")
