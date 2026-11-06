Tabelle: `results/S2.4/S2.4_verdicts.csv` (una riga per giudizio), `S2.4_verdict_counts.csv`, `S2.4_predictions.csv`, `S2.4_P5_glass.csv`, `S2.4_CP12.csv`, `S2.4_CP3.csv`, `S2.4_measures.csv`, `S2.4_G2.csv`.

### Regole

- **R-1, passo 1 delle strisce.** Curva monotona; le righe di S2.0 (ρ 0, −0.2, −0.3) ritrovate identiche (scarto massimo 9e-17). Passo 1 per il generatore (righe; colonne): intestino −0.305; −0.275, tumore e stroma −0.313; −0.212, linfonodo −0.390; −0.461, cervello −0.310; −0.296, cuore −0.174; −0.208. Nessun valore fuori dalla curva; rispetto ai valori estrapolati della sezione 9 cambiano al più di 0.01 (`results/S2.4/S2.4_R1_generator.csv`). Decisione tecnica [DEC-055](../docs/decisions/DEC-055.md).
- **R-2, parte lenta dal vetro.** Identificata in intestino (righe, colonne) e in tumore e stroma (righe, colonne); non nel cervello (correlazione fra metà 0.19 e 0.30), nel linfonodo (righe: 7 righe con entrambe le metà) e nel cuore (vetro insufficiente). 2 sezioni su 5 per asse: **la parte lenta non entra nel simulatore** ([DEC-056](../docs/decisions/DEC-056.md)). Correlazione della parte fine dal vetro con quella dal tessuto 0.34-0.90: i fattori del vetro sono in parte cattura per riga e colonna (`S2.4_R2_slow.csv`).
- **R-3, fondo e seconda coda insieme** (verosimiglianza profilo su ℓ₂, intervalli bootstrap delle tessere, B = 200, `S2.4_R3_fits.csv`, `S2.4_R3_boot.csv`): intestino b 0 [0, 0.079], w₂ 0.036 [0.017, 0.039], ℓ₂ 1 826 µm [1 030, 2 047]; tumore e stroma 0.174 [0.022, 0.259], 0.023 [0.014, 0.041], 429 µm [172, 761]; linfonodo 0 [0, 0], 0.036 [0.033, 0.041], 1 575 µm [1 405, 1 981]; cervello 0.147 [0.117, 0.178], 0.041 [0.036, 0.048], 305 µm [272, 342]; cuore 0 [0, 0.137], 0.077 [0.012, 0.128], 113 µm [71, 179]. Tessuto osservato / previsto 1.010-1.047. In tre sezioni su cinque il fondo costante va a zero: il vetro è spiegato da una coda lunga. Decisione di Luca [DEC-053](../docs/decisions/DEC-053.md).

### Cancelli

| id | esito | dettaglio |
|---|---|---|
| G-1 | superato | R: 37/37 (`S2.4_G1_tests.txt`); Python: 6/6 (`S2.4_G1_pytests.txt`); mutanti 12/12 rilevati (`S2.4_mutants.txt`). Prima del verdetto: due test riscritti e due difetti corretti (correzione C-1) |
| G-2 | **fallito come scritto** | 120 esecuzioni: conservazione 120/120, nessun NA, fondo entro 3 errori Monte Carlo (massimo 2.39); catturati: una esecuzione su 120 a 3.01 contro la soglia scritta di 3 (A4_f1, variante senza seconda coda). Sulle 90 esecuzioni con strisce z medio −0.22, deviazione standard 0.86; con la soglia per i confronti multipli scritta in testa ai cancelli (Φ⁻¹(1 − 0.0005/120) = 4.46) sarebbe superato. Decisione di Luca |
| G-3 | superato | R-1 (3 tabelle), R-3 (2 tabelle senza la colonna dei secondi) e pannello (30 aree, seme k = 0) rifatti in un'altra cartella con 3 e 2 processi invece di 6 e 5: tutto identico; .json identici senza tempi e memoria (`S2.4_G3.txt`) |
| G-4 | superato | (a) 362 file letti, md5 identici prima e dopo (208 in comune con il controllo finale di S2.3, identici); (b) batterie di S2.1 21/21, S2.2 36/36, S2.3 39/39 e 9/9; (c) due aree del pannello di S2.3 rifatte con il codice di S2.4: 4/4 file identici |
| G-5 | superato, con rilievi | tre revisori (sezione «Revisione avversariale»): con le definizioni dello step ogni ricalcolo torna entro la tolleranza; le differenze oltre la tolleranza vengono da definizioni che il testo non fissava (spiegate e qui precisate) |

### Previsioni

Conteggi per esito (giudizi per asse, fascia o quantità):

| previsione | confermata | indeterminata | smentita | altro |
|---|---|---|---|---|
| P-1 parità | 5 | 7 | 0 | |
| P-2 deviazione standard della parte fine | 8 | 4 | 0 | |
| P-3 passo 1 della parte fine | 11 | 1 | 0 | |
| P-4 fasce profonde dei lumi | 0 | 7 | 1 | 9 descrittivi, 1 non valutabile |
| P-5 vetro delle sezioni dopo R-3 | 18 | 9 | 1 | 1 non valutabile (una sola tessera) |
| P-6 quadratini vuoti | 0 | 3 | 3 | |
| P-7 tempo e memoria | 8 | 1 | 3 | |

- **P-1...P-3, strisce**: nessuna smentita su 36 giudizi. Attesa smentita nello stroma: indeterminata (righe −0.030 ± 0.019 contro ±0.020; colonne −0.023 ± 0.020 contro ±0.015). Righe del tumore indeterminate come atteso (−0.066 ± 0.050 contro ±0.073). Righe dell'intestino in P-2 al limite (−0.013 ± 0.005 contro ±0.013), come atteso.
- **P-4, lumi profondi**: tumore 8-16 µm 0.72 ± 0.30, 16-32 µm 0.76 ± 0.40, 32-64 µm 0.85 ± 0.48 (indeterminate; in E-2, senza seconda coda: 0.50, 0.41, 0.39); intestino 1.20 ± 0.25 e 1.05 ± 0.64 (indeterminate); cuore 8-16 µm 0.95 (indeterminata), 16-32 µm 2.27 (indeterminata, 2 aree), 32-64 µm 3.48 (smentita per eccesso, 1 area): il cuore è il controllo negativo.
- **P-5, vetro**: confermate tutte le fasce dell'intestino e del linfonodo salvo una indeterminata (linfonodo 256-512 µm); cervello confermato salvo 32-64 µm (0.86, indeterminata); tumore e stroma confermate fino a 64 µm, indeterminate a 64-512 µm (256-512 µm: 0.70 [0.59, 0.82]), **smentita oltre 512 µm (0.42 [0.33, 0.47], 8 tessere)**, come atteso; cuore indeterminato (poche tessere), 256-512 µm non valutabile.
- **P-6, quadratini vuoti**: smentita per difetto nell'intestino (0.11 ± 0.03), nel tumore (0.10 ± 0.02) e nello stroma (0.40 ± 0.13); indeterminata nel linfonodo (0.71), nel cervello (0.96) e nel cuore (1.23). Le strisce non creano i vuoti veri.
- **P-7, tempo e memoria per area, un processo, server condiviso** (carico medio fino a 740 per un'altra utente): tempo confermato entro 300 s in cinque tessuti (cuore 24 s ... linfonodo 133 s), indeterminato nel tumore (280 ± 44 s); memoria confermata nell'intestino (5 317 MB), nel linfonodo (3 075 MB) e nel cuore (3 429 MB), **smentita nel tumore (11 533 MB), nello stroma (10 497 MB) e nel cervello (15 888 MB)**.

### Controprove

- **CP-1**: con le strisce spente lo stimatore trova deviazione standard della parte fine 0.013-0.031 e parità entro 1 ± 0.002 in 6 tessuti su 6: superata, P-1...P-3 hanno potere.
- **CP-2**: nel tumore, con gli stessi semi, la seconda coda alza la fascia 16-32 µm da 0.024 a 0.046 (differenza 0.022 ± 0.002) e la 32-64 µm da 0.016 a 0.038 (0.022 ± 0.001): superata.
- **CP-3**: negativo non identificato in 5 sezioni su 5 (specificità); positivo identificato in intestino, tumore e stroma, cervello (entrambi gli assi, deviazione standard ritrovata 0.066-0.082 contro 0.08 più la parte fine mediata) e nelle colonne del linfonodo; **non** nelle righe del linfonodo e nel cuore (vetro insufficiente). Fallita come scritta (atteso 5 su 5): in quelle due sezioni R-2 non ha potere, e lì la parte lenta è «non valutabile», non assente.

### Misure

- **Sovra-generazione** (BL-209): M medio per area 3.20 (cuore) ... 5.00 (stroma), intervallo 2.91-5.78; catturati / attesi della scala di R4 1.000 in tutti i tessuti; il tempo della posizione e dello spostamento cresce di 2.0-5.5 volte rispetto alle strisce spente. A tessere di 500 µm il massimo locale risparmierebbe il 18-22 % dei trascritti generati (E-4) ma l'assottigliamento non sarebbe più esatto: con i valori di R-3 la seconda coda porta oltre 280 µm il 2-4 % dei trascritti (`S2.4_measures.csv`; calcolo analitico in «Analisi post hoc»). Il simulatore usa il massimo dell'area intera.
- **UMI per quadratino in tessuto**, simulato di S2.4 / simulato di S2.3: 0.948 (cuore, dove il fondo è andato a 0) ... 1.003; simulato / vero 0.85 (intestino) ... 1.01 (stroma).
