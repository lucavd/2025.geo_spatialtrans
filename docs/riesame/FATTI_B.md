---
type: Facts
title: "Riesame, Parte B: i fatti sul generatore"
description: "Numeri riletti da results/riesame (B_D_common.json, B_per_measure_common.csv, B_G4.csv, xen_match.csv); dati ai difensori senza la storia del progetto."
status: draft
generated: { by: geo-agent/skill-v2, at: 2026-11-20T17:30:00+01:00 }
sources:
  - id: regola
    resource: "/riesame/REGOLA_B.md"
    title: "Riesame, Parte B: regola di esito scritta prima del calcolo"
---
# Fatti della Parte B

Domanda. Un simulatore di dati Visium HD (spatial transcriptomics a quadratini da 2 µm) con uno strato cellulare esplicito (generatore di cellule, nuclei, tipi; profili di espressione; spostamento dei trascritti; modello di cattura con strisce e fondo). Serve a testare metodi di analisi. La domanda del riesame: il simulatore da' qualcosa che non si ottiene piu' semplicemente?

Bracci, tumore del colon e suo stroma, 18 024 geni comuni:
- V: sezione vera di taratura (5 aree per tessuto).
- T: sezioni vere di altri due pazienti (P1, P5): il tetto, cioe' quanto due dati veri distano fra loro.
- S: il simulatore di oggi (5 aree x 2 semi).
- B: braccio banale scritto prima: centri da Poisson omogeneo con la densita' del preset, territori di Voronoi, nucleo = cerchio con l'area media ritagliato, tipi estratti a caso dalla composizione, trascritti dal profilo medio del tipo messi a caso nella cellula, stessa cattura del simulatore (strisce, fondo, profondita') con lo spostamento dei trascritti spento (5 aree x 2 semi; un'area con un solo seme per un errore geometrico del pacchetto).
Stessa catena di analisi per tutti (bin2cell + RCTD + etichette trasferite).

Misura D: radice con segno della media, sulle 12 celle tessuto x misura (M-1 quota UMI assegnati alle cellule, M-2 UMI per cellula, M-3 marcatori estranei, M-4 quota singoletti RCTD, M-5 concordanza tipo RCTD con il nucleo, M-6 vicini dello stesso tipo), di (Delta^2 - Var Delta)/h^2, Delta = media del braccio - media di V; h = precisione della media di V (semi-ampiezza dell'intervallo bootstrap). Piu' basso e' meglio. Incertezza: bootstrap parametrico per aree, B = 2000.
Soglie scritte prima: margine delta = 1,0 h. Uguale se |Dab| + w <= 1; diverso se |Dab| - w > 1; altrimenti indeterminato (w = semi-ampiezza al 95 %). Onere della prova sul componente nostro: vale solo con un 'diverso' a suo favore.

Regola di esito scritta prima (sul generatore): D(S) uguale a D(XC) -> si tiene; D(S) peggio del riferimento e meglio di D(B) (entrambi diversi) -> si tiene solo per prove controllate; D(S) uguale, indeterminato o peggio rispetto a D(B) -> il generatore non batte il banale: si toglie. Verdetto unico fra: continuare, restringere, cambiare strada, chiudere.

## Risultati (results/riesame/B_D_common.json)

| gruppo | D(T) tetto | D(S) simulatore | D(B) banale | S - T (w) | S - B (w) | B - T (w) |
|---|---|---|---|---|---|---|
| tumore + stroma | 2.23 (1.85; 2.71) | 7.83 (7.44; 8.24) | 7.77 (7.41; 8.15) | +5.60 (0.50) different | +0.06 (0.29) equal | +5.54 (0.51) different |
| tumore | 2.10 (1.39; 3.01) | 9.89 (9.27; 10.49) | 9.83 (9.30; 10.42) | +7.78 (0.83) different | +0.05 (0.44) equal | +7.73 (0.83) different |
| stroma | 2.35 (2.00; 2.72) | 4.99 (4.78; 5.23) | 4.91 (4.70; 5.12) | +2.65 (0.38) different | +0.09 (0.22) equal | +2.56 (0.41) different |

Per misura (distanza con segno in h, results/riesame/B_per_measure_common.csv):

| tessuto | misura | T | S | B |
|---|---|---|---|---|
| tumore | UMI assegnati | 2.3 | 2.3 | 0.5 |
| tumore | UMI per cellula | 3.3 | -0.4 | -0.3 |
| tumore | marcatori estranei | 2.3 | 8.1 | 3.3 |
| tumore | singoletti | 1.3 | 10.2 | 10.2 |
| tumore | concordanza col nucleo | 1.9 | 17.1 | 17.7 |
| tumore | vicini stesso tipo | 0.2 | 10.9 | 12.3 |
| stroma | UMI assegnati | -0.1 | 1.0 | 1.0 |
| stroma | UMI per cellula | 1.5 | -0.2 | -0.2 |
| stroma | marcatori estranei | 1.3 | 1.9 | 1.1 |
| stroma | singoletti | 1.9 | 10.4 | 10.6 |
| stroma | concordanza col nucleo | 3.6 | 2.5 | 0.5 |
| stroma | vicini stesso tipo | 3.6 | 5.6 | 5.4 |

Controlli: ricalcolo della catena contro i valori precedenti 7/7 entro 1e-9 (massimo 1,7e-15); due pazienti veri fra loro non diversi (P1 - P5 = -0,42, w 0,68); stimatore senza distorsione su spostamenti costruiti 0, 1, 2 h (stimati: quadrato medio -0,03; 0,99; 2,01); md5 dei 550 file letti uguali prima e dopo.

Domanda sulla cattura (il modello di strisce, fondo e spostamento batte il dato Xenium vero diviso in quadratini e diradato alla stessa profondita'?): indecidibile con i dati disponibili. Le registrazioni Xenium -> Visium HD pubblicate da 10x appartengono a pacchetti Xenium da circa 1,1 milioni di cellule che non sono su GEO; i pacchetti GEO (P1 307 762, P2 340 837, P5 275 998 cellule) non corrispondono a nessuna registrazione (identificativi comuni 0,02-0,03 %, livello del caso). Vincolo che la riapre: quei pacchetti.

Fatti dal riesame precedente (Parte A, file FATTI.md): prima prova d'uso: 36 giudizi, 14 confermati, 4 indeterminati, 18 smentiti; valori veri dentro la variabilita' simulata: 18 su 228; nessun confronto pubblicato con simulatori esistenti (uno solo, interno, del luglio 2025).

# Revisione indipendente (G-5) e difensori

- G-5, sotto-agente con codice proprio (`results/riesame/G5/`): 18 D e 12 scarti coincidono con i dichiarati entro 1,8e-15; i 12 w entro l'errore Monte Carlo (20 semi, deviazione standard di w 0,006-0,032 h; chi-quadro 8,15 su 12, p = 0,77); esiti identici con tutti i 20 semi. Rilievi senza effetto sui numeri: s_p dello stroma include le 3 aree da 6,5 mm (come il codice); A2_r2 di B con un solo seme (B-13); intervalli di D ricentrati sulla mediana, solo descrittivi. Esito: passa.
- Difensore di «continuare»: il difetto sta nei componenti comuni a S e B (profili, cattura, catena), non nella geometria; propone una prova B' con un solo componente cambiato. Ammette che, con la regola scritta prima, i fatti portano a «chiudere» e che B' sarebbe una scelta fatta dopo aver visto i dati.
- Difensore di «chiudere»: S uguale a B in tutti e tre i gruppi, quindi il generatore si toglie; la cattura e' senza risposta. Ammette che la chiusura poggia sulla clausola delle due sessioni, non su una prova contro la cattura.
