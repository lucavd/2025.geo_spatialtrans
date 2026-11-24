---
type: State
title: "Stato di geo_spatialtrans"
description: "Memoria fra le sessioni: fase, dati, numeri, bilancio dell'errore, minacce, strade fallite, prossima mossa."
status: draft
generated: { by: geo-agent/progetto-master, at: 2026-11-24T01:30:00+01:00 }
---
# Fase

Sonda, budget 7 (proroga del 2026-11-23). Registro: S1, S2 (`infra`), S3b, S3c, S4 (`misura`). Usate 5 sessioni, ne restano 2.

# Dati (polmone 10x, CG000709, 1 blocco FFPE, 1 donatore)

- H1 (T) e H2 (V) in `/mnt/micron/geo_spatialtrans/data_lung/H1`, `/H2`; md5 in `results/sonda/S3/MD5_download_H.txt`. HLCA core in `/mnt/micron/geo_spatialtrans/ref_lung/`.
- MALAT1, NEAT1, KRT18 non sono nel pannello di sonde.
- Metà tenuta da parte: solo H2 (x < m − 50 µm, m = 3372 µm), X1, X2. H1 letta intera, mai in taratura (griglie di H1 e H2 ruotate di 180°).

# Metro (S3c, `docs/sonda/REGOLA_S3c.md`)

- Tabella e h di S3b (`results/sonda/S3b/area_values.csv`, `h.csv`). Stimatore appaiato `tools/sonda/S3c_D.py`: d_k = braccio − V per area, q = (media(d)² − var(d)/5)/h², D = radice con segno della media di q sulle 12 celle.
- Bootstrap per coppie stratificato per tessuto, B = 2000, seme 20261123, semi-ampiezze ed se × f = 1,584 (scelto da COV). Cancelli G-ID, C-1p, C-2, C-4, C-5 e CP-perm passati. Revisione indipendente di S3c e S4 fatta (S4: 235/235).

# Numeri (h)

Tumore + stroma, `results/sonda/S3c/D.json`: pavimento D(B) = 3,151 (2,908; 3,440), tetto D(T) = 0,336 (−0,263; 0,821), soglia 1,336, pavimento − δ = 2,151. A1 oggi = D(B) = 3,151 ± 0,141: non regge (il simulatore è il banale).

Componenti veri, S4 (`docs/sonda/REGOLA_S4.md`, `results/sonda/S4/D.json`; banale con un componente preso dall'area omologa di H2, un seme):

| braccio | componente vero | D | intervallo | quota del divario chiusa |
|---|---|---|---|---|
| Y | mappa dei tipi | 1,771 | 1,514; 2,200 | 0,69 |
| G | nuclei | 3,014 | 2,766; 3,351 | 0,09 |
| U | mappa degli UMI | 2,636 | 1,565; 3,582 | 0,30 (0,02; 0,60) |
| A | nuclei + tipi + UMI | 0,961 | 0,691; 1,315 | 0,92 |
| AE | A + espressione per cellula | 0,509 | 0,359; 0,652 | 0,99 |

Confronti: A e AE contro B «diverso»; AE contro T «uguale»; G contro B «uguale»; Y e U contro B indeterminati. Attese: P1, P2, P4, P5 confermate; P3 (U non aiuta) indeterminata. Cancelli G-T 16/16, G-ID 50/50, G-S3c scarto 0, G-MD5 uguali, G-REV 235/235. Controprova CP-AE passata.

# Check

`check smoke` = `tools/sonda/S3c_smoke.sh` (circa 5 min). Non contiene S4: i bracci S4 sono diagnostica, non catena. `check full` da scrivere.

# Bilancio dell'errore (per componente, S4; `results/sonda/S4/per_measure.csv`)

1. Disposizione dei tipi nello spazio: 69 % del divario. Sistema singoletti RCTD e vicini dello stesso tipo in entrambi i tessuti.
2. UMI concentrati sui nuclei (quota assegnata alle cellule, M-1, 19 % del divario per voce): nessuna sostituzione singola la muove (Y, G, U); torna solo con nuclei e UMI veri insieme (A). In A coincide con V per costruzione: la quota è un limite superiore.
3. Espressione della singola cellula: AE − A = −0,45 h, «uguale»; 7 % del divario, soprattutto la concordanza col nucleo nello stroma. Sotto il 10 %: nessuna sessione.
4. Geometria dei nuclei da sola: 9 %. Nessuna sessione.

Con i componenti veri migliori D = 0,509 (AE) e 0,961 (A), sotto la soglia: lo spazio per l'articolo c'è, se un generatore impara 1 e 2.

# Minacce alla validità

1. HLCA è tessuto sano: le cellule tumorali non hanno un tipo (pesa su M-4).
2. I componenti veri vengono dall'area stessa di V: Y porta anche la composizione locale dell'area, non solo il raggruppamento. Un generatore impara le statistiche, non la mappa. Prima prova della costruzione: mappa generata contro mappa vera ruotata o di un'altra area.
3. Un donatore: tetto fra sezioni seriali.
4. A5 (X1 contro X2): controllare l'orientamento delle griglie prima di usarle.
5. U chiude il 30 % del divario, quasi tutto sui singoletti, senza tipi veri; intervallo 2-60 %: non spiegato, non cambia l'ordine delle voci.

# Strade fallite

- Generatore dello strato cellulare (colon): 7,83 contro 7,77 del banale.
- Tetto di S3 su aree non omologhe e κ del colon: C-4 rosso.
- Stimatore indipendente su aree omologhe (S3b): tetto −0,77, distorto.

# Prossima mossa

Restano 2 sessioni di sonda: (1) A3 e A5 sullo Xenium, prima l'orientamento delle griglie di X1 e X2; (2) G1. A G1: con i componenti veri migliori D ≤ soglia, la regola del contratto indica costruzione; budget proposto da rifare sul bilancio per componente (2 voci grandi × 2 colpi, più conferma).

# Aperti

- `results/F3.2/sheets/` non tracciato: lasciato com'è.
- `tools/okf.py check` rosso da SC1 (esenzioni non modificate).
- Revisione indipendente dei numeri di S1, S2, S3b non fatta (prima di G1).
- X2 senza registrazione misurata.
