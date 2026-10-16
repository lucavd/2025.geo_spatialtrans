# R2c — Pre-registrazione: verità manuale sul contatto fra nuclei (BL-120)

**Approvata da Luca il 2026-10-16** («si»), prima di qualsiasi annotazione, calcolo di candidati sui segmentatori o risultato.
Scelte di disegno di Luca (2026-10-16): compito «sì/no su coppie proposte»; copertura «nuclei focali + Federica a metà».

## Obiettivo unico
Verità di riferimento manuale per il contatto fra nuclei sulle 24 finestre di R2b, per decidere per archetipo quale segmentatore
(cellpose_rgb, stardist_he) dà una frazione di nuclei a contatto (S-1) credibile, e dove cade K-2 (S1.4c). Dipendenze: S1.4c chiusa
con riserve, finestre e punti di R2b presenti — soddisfatte.

## Disegno
- **Punti**: punti manuali di Luca R2b (`results/R2b/R2b_points.csv`, rater = Luca), in µm del ROI.
- **Candidati**: lati della triangolazione di Delaunay dei punti di Luca per finestra con distanza d ≤ d_max(archetipo).
  Regola per d_max (fissata ora): q99.5 della distanza fra i due punti manuali appaiati alle due etichette di ogni coppia a contatto
  (Chebyshev ≤ δ = 2 px, operatore di S1.4c) in cellpose_rgb ∪ stardist_he, sulle coppie con entrambi i punti appaiati, + 1 µm.
- **Nuclei focali**: 15 per finestra (60 per archetipo), campionati con seme fisso 20261016 fra i punti di Luca a distanza ≥ d_max
  dal bordo della finestra e dalle zone escluse. Si giudicano tutte le coppie candidate che contengono almeno un nucleo focale.
- **Compito**: ritaglio ingrandito centrato sulla coppia, due marcatori piccoli, nessun contorno; risposte «si toccano» /
  «separati» / «non giudicabile». Ordine casuale (seme), alla cieca rispetto ai segmentatori. ~10 % delle coppie ripresentate
  (intra-rater); coppie-trappola con d > 20 µm (fra punti della stessa finestra), ~5 % del totale.
- **Annotatori**: Luca tutte le coppie; Federica le coppie delle sue 2 finestre per archetipo di R2b (flag `rater2`), HTML autonomo.
- **S-1 manuale** per nucleo focale = almeno una coppia «si toccano». Primaria: «non giudicabile» escluse (nucleo escluso se
  tutte le sue coppie sono non giudicabili e nessuna è «si toccano»); intervallo con le varianti «non giudicabile = contatto» /
  «= separati». IC 95 %: bootstrap a grappolo per finestra (B = 2000, finestre ricampionate con reinserimento entro archetipo).
- **S-1 dei segmentatori** sugli stessi nuclei focali appaiati in quel segmentatore, sulle maschere native del ROI con il codice
  raster di S1.4c (contatto con qualunque etichetta, δ = 2 px). Space Ranger 4 come terzo confronto descrittivo in A4/A6 (BL-121).

## Check C — correttezza software
| id | attesa | soglia |
|---|---|---|
| C-1 copertura | frazione delle coppie a contatto di ciascun segmentatore (entrambi i punti appaiati) che cade fra i candidati | PASS ≥ 0.99, WARN 0.97–0.99, FAIL < 0.97, per archetipo × segmentatore |
| C-2 annotatore | ritagli pixel-identici alla finestra R2b; ordine riproducibile dal seme; HTML senza informazioni dei segmentatori (grep); test del core JS | tutti PASS |
| C-3 operatore | S-1 dei segmentatori calcolato dal codice R2c = valori per nucleo di S1.4c (`cp2`) sugli stessi nuclei | identico |
| C-4 mutanti | bootstrap senza reinserimento; bootstrap non a grappolo; κ con marginali scambiate | 3/3 rilevati |
| C-5 eleggibilità (BL-125a) | regole di candidati/focali eseguite e ispezionate su una finestra per archetipo prima di generare l'annotatore | eseguito |

## Check B — previsioni (di Claude, non misure)
- **P-1**: S-1 manuale ≥ 0.30 in A1–A4 e ≤ 0.15 in A5, A6 (soglie di P1 di S1.4c).
- **P-2**: S-1 manuale strettamente fra StarDist e cellpose in ≥ 3 di A1–A4 (cellpose più largo, StarDist più convesso; S1.4c).
- **P-3**: sulle coppie con entrambi i punti appaiati in entrambi i segmentatori: sensibilità cellpose > StarDist, specificità StarDist > cellpose.
- **B-dec** (per archetipo): un segmentatore è «credibile» se il suo S-1 cade nell'IC 95 % manuale; se entrambi, quello con |Δ|
  minore; se nessuno, riferimento = intervallo manuale. Esiti: cellpose / StarDist / nessuno / indeciso (entrambi credibili e
  |Δ| che differiscono < 0.02).
- **B-K2**: S-1 di K-2 dai raster di S1.4c (10 semi, mediana), nuclei simulati con centroide nelle finestre (stessa regola di
  bordo dei focali). PASS dentro l'IC manuale; WARN dentro l'IC allargato di |S-1 Luca − S-1 Federica| dell'archetipo; FAIL
  altrimenti. Limite dichiarato: insieme di nuclei diverso (generatori dai segmentatori).

## Controprove
| id | test | soglia |
|---|---|---|
| CP-1 | κ di Cohen Luca–Federica sulle coppie comuni (non giudicabile esclusa) | ≥ 0.60 PASS; 0.40–0.60 WARN; < 0.40 FAIL → nessuna verità, nessuna decisione. Convenzione Landis & Koch 1977 (doi:10.2307/2529310) |
| CP-2 | κ intra-annotatore sulle ripetute | ≥ 0.70 |
| CP-3 | coppie-trappola giudicate «separati» | ≥ 95 % |
| CP-4 | AUC di −d per «si toccano» in A1–A4 | ≥ 0.80 per archetipo |
| CP-5 | B-dec con le sole risposte di Federica sul suo sottoinsieme | stesso esito in ≥ 3/4 di A1–A4 |

## Potenza (dichiarata prima)
60 focali, S-1 ≈ 0.3: semi-ampiezza IC ~0.12–0.14. Divario cellpose–StarDist ~0.18 (A1), ~0.28 (A2): risolvibile; A3 ~0.13:
rischio «indeciso». S-2/S-3 restano declassate (non misurate da un giudizio sì/no).

## Perimetro
Solo `tools/R2c_*`, `results/R2c/`, `R/testing/test_R2c.sh`; pacchetto invariato. Ogni ridefinizione rispetto a questo testo
va dichiarata come deviazione (BL-124).
