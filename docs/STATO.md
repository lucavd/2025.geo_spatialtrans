---
type: State
title: "Stato di geo_spatialtrans"
description: "Memoria fra le sessioni: fase, dati, numeri, bilancio dell'errore, minacce, strade fallite, prossima mossa."
status: draft
generated: { by: geo-agent/progetto-master, at: 2026-11-21T09:30:00+01:00 }
---
# Fase

Sonda, budget 5, firmata a G0 il 2026-11-20. S1 (cancello della registrazione) chiusa il 2026-11-21: X1 passa, X2 non giudicabile. Riserva: colon (`docs/archivio/PAPER_colon_SC1.md`).

# Dati (polmone 10x, CG000709, 1 blocco FFPE, 1 donatore)

- H1, H2 solo Visium HD (non ancora scaricati). X1 = Visium HD dopo Xenium v1 (289 geni, 286 in comune con le sonde Visium HD); X2 = Visium HD dopo Xenium 5K (5001 geni, 4892 in comune).
- X1, X2 in `/mnt/micron/geo_spatialtrans/data_lung/` (22 file, dimensioni = CDN, md5 in `MD5_download.txt`); quadratini da 2 e 8 µm estratti in `*/visium/binned/`.
- MALAT1 e NEAT1 non sono nel pannello di sonde Visium HD.

# Registrazione (S1, `docs/sonda/REGOLA_S1.md`, `tools/sonda/S1_registration.py`)

- Catena Xenium µm -> H&E: inversa di `he_imagealignment.csv` su pixel Xenium da 0,2125 µm. Regola della variante fermata (62 % e 59 % su tessuto, soglia 0,8) perché lo Xenium copre più dell'area di cattura (38-41 % delle cellule fuori griglia; 100 % e 99,9 % su tessuto fra quelle dentro): errore di pre-registrazione, Luca: proseguire.
- H&E dello Xenium = tissue_image del Visium HD (stesse dimensioni, spostamento nullo). H&E -> quadratini: proiettiva esatta (2e-11 px).
- X1, riquadri da 200 µm (errore standard 0,29-0,41 µm): catena pubblicata spostata di circa (3,5; 9,2) µm, residuo p95 9,0-11,2 µm nelle 3 regioni. Affine locale (riquadri non usati nella stima): p95 0,84 / 0,68 / 0,39 µm. Controlli positivi 9/9 entro 0,3 µm, negativi 0-4 % entro 2 µm. Passa. Una sola traslazione globale (non pre-registrata) darebbe p95 0,88 / 1,79 / 1,53 µm.
- X2: regola di precisione non soddisfatta (errore standard 1,1-3,4 µm a 200-250 µm; a 100 µm nessun riquadro valido): poche UMI (70 per quadratino da 8 µm) e pochi geni con almeno 10 conteggi. Cancello non giudicato.
- Risultati: `results/sonda/S1/` (copie da `/mnt/micron/geo_spatialtrans/S1_sonda/`).

# Metro

D = radice con segno della media su 12 voci (tumore, stroma x M-1...M-6) di (Delta^2 - Var Delta)/h^2; `tools/riesame_D.py`; controlli in `docs/riesame/FATTI_B.md`. Aree di tumore e stroma del polmone da definire con una regola scritta prima.

# Numeri (h)

Polmone: pavimento, tetto, soglia non misurati. Colon (riserva): pavimento 7,77; tetto 2,23; simulatore ereditato 7,83.

# Bilancio dell'errore

Polmone: non ancora. Riferimento dal colon (`results/riesame/B_per_measure_common.csv`): tumore, concordanza con il nucleo 43 %; tumore, vicini dello stesso tipo 18 %; singoletti RCTD stroma 15 % e tumore 15 %; marcatori estranei 9 %.

# Minacce alla validita'

1. X2 senza registrazione misurata: A5 (tetto = concordanza X1 contro X2) dipende da X2.
2. La catena pubblicata sbaglia di circa 10 µm: la verità Xenium va corretta con una trasformazione stimata anche dall'espressione Visium HD (lieve circolarità per A3).
3. Un donatore: tetto fra sezioni seriali.
4. Verita' solo su sezioni meno sensibili (X2: 70 UMI per quadratino da 8 µm).
5. Circolarita' del metro (bin2cell + RCTD, A4).

# Strade fallite

- Generatore dello strato cellulare (colon): 7,83 contro 7,77 del banale (`docs/riesame/CHIUSURA_GENERATORE.md`).
- Taratura dello spostamento (colon): marcatori estranei del tumore 8,1 con, 3,3 senza.
- Cattura contro Xenium in quadratini sul colon: indecidibile (XD-26).
- Residuo della registrazione con i picchi nucleari MALAT1/NEAT1: geni assenti dal pannello.

# Prossima mossa

Decisione di Luca attesa: misurare X2 con un metro più preciso (una sessione, regola scritta prima: pesi dei geni per conteggi invece della standardizzazione per gene, oppure gruppi di geni) oppure proseguire solo con X1 riscrivendo prima il tetto di A5. Poi metro e check smoke sul polmone, con la correzione della registrazione di X1 (affine stimata sulla metà di sviluppo).

# Aperti

- `results/F3.2/sheets/` non tracciato (PNG del 18-19 novembre, di F3.2): lasciato com'e'.
- `tools/okf.py check` fallisce da SC1: `PAPER.md`, `STATO.md`, `archivio`, `sonda` fuori dalle esenzioni. Elenco non modificato.
- Aggregazione dei quadratini da 8 µm sfasata di circa 2 µm rispetto a floor(i/4) dei 2 µm: da capire prima di usare i due livelli insieme.
- Revisione indipendente dei numeri di S1 non fatta (prevista prima di G1).
- Riserva per un tetto fra pazienti: Vannan 2025 TMA5.
- S1 non registrata in results/misure.tsv: progetto.py chiudi rifiuta senza check (check smoke da scrivere in S2); registrarla appena il check esiste.
