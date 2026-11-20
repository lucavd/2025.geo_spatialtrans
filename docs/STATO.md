---
type: State
title: "Stato di geo_spatialtrans"
description: "Memoria fra le sessioni: fase, dati, numeri, bilancio dell'errore, minacce, strade fallite, prossima mossa."
status: draft
generated: { by: geo-agent/progetto-master, at: 2026-11-20T21:30:00+01:00 }
---
# Fase

Sonda, budget 5, firmata a G0 il 2026-11-20. Scommessa chiusa: SC1 (adozione, contratto sul colon), SC2 (scelta del tessuto: polmone 10x dopo Xenium). Contratto del colon in `docs/archivio/PAPER_colon_SC1.md`: riserva se la registrazione fallisce.

# Dati (polmone 10x, nota tecnica CG000709, 1 blocco FFPE, 1 donatore)

- H1 = Visium_HD_Human_Lung_Cancer_HD_Only_Experiment1; X1 = Visium_HD_Human_Lung_Cancer_post_Xenium_v1_Experiment1 + Xenium_V1_Human_Lung_Cancer_FFPE (289 geni); H2 = Visium_HD_Human_Lung_Cancer_HD_Only_Experiment2; X2 = Visium_HD_Human_Lung_Cancer_post_Xenium_Prime_5K_Experiment2 + Xenium_Prime_Human_Lung_Cancer_FFPE (5001 geni).
- CDN: `https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/<s>/<s>_<file>`; Xenium `.../samples/xenium/3.0.0/<s>/<s>_<file>`. curl funziona; urllib di Python riceve 403. Pagine 10x bloccate da un controllo anti-bot.
- UMI medi per bin da 8 um (metrics_summary, SR 3.0.0; mediane assenti): H1 211,7; X1 129,3; H2 274,9; X2 70,4. Colon P2 488,4.
- Registrazione non fornita. Catena: `he_imagealignment.csv` (H&E -> Xenium) + allineamento Space Ranger (H&E -> quadratini). L'H&E dello Xenium ha dimensioni e pixel (0,2738 um) del tissue_image del Visium HD (sotto-agente SC2, da verificare in S1).
- Nessun tipo cellulare fornito; riferimento RCTD: HLCA (CELLxGENE).
- Metadati in `/mnt/micron/geo_spatialtrans/sc2_meta/`; inventario SC2 in ESISTENTE.

# Metro

D = radice con segno della media su 12 voci (tumore, stroma x M-1...M-6) di (Delta^2 - Var Delta)/h^2; Delta = braccio meno riferimento; h = precisione della media delle 5 aree della sezione di taratura (H2). `tools/riesame_D.py`; controlli in `docs/riesame/FATTI_B.md`. Aree di tumore e stroma del polmone da definire con una regola scritta prima.

# Numeri (h)

Polmone: pavimento, tetto, soglia non misurati. Colon (riserva): pavimento 7,77; tetto 2,23; simulatore ereditato 7,83.

# Bilancio dell'errore

Polmone: non ancora. Riferimento dal colon (`results/riesame/B_per_measure_common.csv`): tumore, concordanza con il nucleo 43 %; tumore, vicini dello stesso tipo 18 %; singoletti RCTD stroma 15 % e tumore 15 %; marcatori estranei 9 %.

# Minacce alla validita'

1. Registrazione composta non misurata: cancello di S1.
2. Un donatore: il tetto e' fra sezioni seriali, probabilmente piu' stretto che fra pazienti, quindi soglia piu' difficile.
3. Verita' solo su sezioni meno sensibili (1,64x v1, 3,90x 5K): A5 vale a quella profondita'.
4. Circolarita': il metro usa bin2cell + RCTD (A4).
5. h dipende dal disegno a 5 aree; la meta' tenuta da parte riduce l'area per tetto e taratura.

# Strade fallite

- Generatore dello strato cellulare (colon): 7,83 contro 7,77 del banale (`docs/riesame/CHIUSURA_GENERATORE.md`).
- Taratura dello spostamento (colon): marcatori estranei del tumore 8,1 con, 3,3 senza.
- Cattura contro Xenium in quadratini sul colon: indecidibile (XD-26). Sul polmone diventa A3.
- Colon con verita' sulla stessa sezione: non esiste (Long 2025 = P2 su sezione seriale; ricerca SC2).

# Prossima mossa

S1, sonda 1 di 5: cancello della registrazione.

1. Prima di scaricare, scrivere e committare: (a) la regola della meta' spaziale tenuta da parte per ogni sezione; (b) il metodo del residuo: nuclei Xenium mappati sull'H&E contro nuclei StarDist sull'H&E, e scarto H&E -> quadratini da 2 um con i picchi di UMI nucleari (MALAT1, NEAT1); 3 regioni disgiunte nella meta' di sviluppo; p95 <= 2 um, affine locale ammessa (XD-4).
2. Scaricare per X1 e X2: binned_outputs, tissue_image.btf, he_imagealignment.csv, cells, cell_boundaries, transcripts.parquet (circa 15 GB) in `/mnt/micron/geo_spatialtrans/data_lung/`.
3. Esito: passa in X1 e X2, si prosegue; passa in una sola, si prosegue con quella e il tetto di A5 si riscrive prima di misurarlo; fallisce in entrambe, ritorno al colon.

Poi S2 metro e check smoke sul polmone; S3 pavimento e tetto; S4 componenti veri e bilancio; S5 A3 e A5, G1. Pulizie: `CLAUDE.md` aggiornato al polmone e sotto 100 righe; storici in `docs/archivio/`.

# Aperti

- `results/F3.2/sheets/` non tracciato (PNG del 18-19 novembre, di F3.2): lasciato com'e'.
- Riserva per un tetto fra pazienti con verita': Vannan 2025 TMA5 (4 donatori, fibrosi); qualita' del Visium HD non letta.
