---
type: State
title: "Stato di geo_spatialtrans"
description: "Memoria fra le sessioni: fase, dati, numeri, bilancio dell'errore, minacce, strade fallite, prossima mossa."
status: draft
generated: { by: geo-agent/progetto-master, at: 2026-11-21T12:30:00+01:00 }
---
# Fase

Sonda, budget 5, firmata a G0 il 2026-11-20. Nel registro S1 e S2, entrambe `infra`. Prossima: sessione 3 di 5. Allarme aperto: più di una sessione su tre senza misura. Riserva: colon (`docs/archivio/PAPER_colon_SC1.md`).

# Dati (polmone 10x, CG000709, 1 blocco FFPE, 1 donatore)

- H1, H2 solo Visium HD (non ancora scaricati). X1 = Visium HD dopo Xenium v1 (286 geni in comune con le sonde); X2 = dopo Xenium 5K (4892 in comune).
- X1, X2 in `/mnt/micron/geo_spatialtrans/data_lung/` (md5 in `MD5_download.txt`); quadratini da 2 e 8 µm in `*/visium/binned/`.
- MALAT1 e NEAT1 non sono nel pannello di sonde Visium HD.

# Registrazione

- S1 (`docs/sonda/REGOLA_S1.md`): catena pubblicata spostata di circa 10 µm; X1 passa con affine locale (p95 0,84 / 0,68 / 0,39 µm, controlli 9/9 e 0-4 %). X2 non giudicabile (errore standard 1,1-3,4 µm).
- S2 (`docs/sonda/REGOLA_S2.md`, `tools/sonda/S2_registration.py`, `results/sonda/S2/`): metri che pesano i geni per i conteggi, A (correlazione) e B (Poisson), σ 2 e 4 µm. Sintetico: tutti passano; sintetico tipo X2 a 100 µm, errore 1,37 µm (S1) contro 0,22-0,56. X1: controlli 9/9 per tutti; ridanno S1 entro 0,3 µm A2 e B2 (ammessi); A4 e B4 no (0,31-0,33 µm in x, R3), A4 cambia l'esito di R2. B2 su X1: p95 locale 0,54 / 0,85 / 0,56 µm. Il percorso S1 dentro S2 ridà S1 esattamente.
- X2, errore standard a 250 µm: A2 1,09 / 0,89 / 0,53; B2 1,19 / 1,21 / 0,43 (S1: 1,88 / 3,44 / 1,14). Nessun L* (anche: R2 ha 8 riquadri usati su 16-25). Esito: X2 non giudicabile, primo colpo dell'ipotesi «metro pesato».
- Diagnosi (`S2_diag_counts.py`, solo conteggi): il filtro per gene ereditato (≥ 10 conteggi in entrambe le mappe) tiene in X2 il 0-4 % delle UMI a 100 µm e il 6-34 % a 200-250 µm (X1: 42-90 %). Senza filtro, X2 avrebbe 25/25 e 16/16 riquadri validi a 200 e 250 µm, 60-98/100 a 100 µm. Il metro pesato non ha avuto una prova piena.

# Metro

D = radice con segno della media su 12 voci (tumore, stroma x M-1...M-6) di (Delta^2 - Var Delta)/h^2; `tools/riesame_D.py`; controlli in `docs/riesame/FATTI_B.md`. Aree di tumore e stroma del polmone da definire con una regola scritta prima.

# Check

`check smoke` parziale: solo il test sintetico del metro della registrazione in uso (S1, `tools/sonda/S1_test_cc.py`). Catena e metro D da aggiungere. `check full` da scrivere.

# Numeri (h)

Polmone: pavimento, tetto, soglia non misurati. Colon (riserva): pavimento 7,77; tetto 2,23; simulatore ereditato 7,83.

# Bilancio dell'errore

Polmone: non ancora. Riferimento dal colon (`results/riesame/B_per_measure_common.csv`): tumore, concordanza con il nucleo 43 %; tumore, vicini dello stesso tipo 18 %; singoletti RCTD stroma 15 % e tumore 15 %; marcatori estranei 9 %.

# Minacce alla validita'

1. X2 senza registrazione misurata: A5 (tetto = concordanza X1 contro X2) e la parte X2 di A3 dipendono da X2.
2. La catena pubblicata sbaglia di circa 10 µm: la verità Xenium va corretta con un'affine stimata anche dall'espressione Visium HD (lieve circolarità per A3).
3. Un donatore: tetto fra sezioni seriali.
4. Verita' solo su sezioni meno sensibili.
5. Circolarita' del metro (bin2cell + RCTD, A4).

# Strade fallite

- Generatore dello strato cellulare (colon): 7,83 contro 7,77 del banale (`docs/riesame/CHIUSURA_GENERATORE.md`).
- Taratura dello spostamento (colon): marcatori estranei del tumore 8,1 con, 3,3 senza.
- Cattura contro Xenium in quadratini sul colon: indecidibile (XD-26).
- Residuo della registrazione con MALAT1/NEAT1: geni assenti dal pannello.

# Prossima mossa

Proposta a Luca (2026-11-21), in attesa: sessione 3 = metro D e check smoke sul polmone con X1 (aree di tumore e stroma con regola scritta prima, affine di X1 stimata sulla metà di sviluppo), per avere una misura e chiudere l'allarme. La seconda e ultima prova su X2 (metri A2 e B2 senza filtro per gene, validità sui conteggi totali, stessa soglia 0,5 µm, regola scritta prima) diventa il primo passo della sessione su A3 e A5; se fallisce, il tetto di A5 si riscrive prima di misurarlo.

# Aperti

- `results/F3.2/sheets/` non tracciato (PNG del 18-19 novembre, di F3.2): lasciato com'e'.
- `tools/okf.py check` fallisce da SC1: `PAPER.md`, `STATO.md`, `archivio`, `sonda` fuori dalle esenzioni. Elenco non modificato.
- Aggregazione dei quadratini da 8 µm sfasata di circa 2 µm rispetto a floor(i/4) dei 2 µm.
- Revisione indipendente dei numeri di S1 e S2 non fatta (prevista prima di G1).
- X1, R3: metri pesati e S1 differiscono di circa 0,3 µm in x.
- Riserva per un tetto fra pazienti: Vannan 2025 TMA5.
