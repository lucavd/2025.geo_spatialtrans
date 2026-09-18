# BIO_REFERENCES — valori di riferimento biologici

Regola: una riga entra qui solo con **fonte tracciabile**: DOI di un articolo, oppure `script` + `dataset` (in `data/real/README.md`) che ha prodotto la stima. Niente valori a memoria.

| id | Archetipo | Metrica | Valore | Intervallo / distribuzione | Unità | Fonte (DOI o script@commit + dataset) | Sessione | Note |
|---|---|---|---|---|---|---|---|---|
| B-001 | A1 | densità nucleare | | | cell/mm² | | | |

## Stato
- **S0 (2026-09-18)**: schema creato, **nessun valore inserito**. Le ipotesi numeriche presenti nel design doc (densità del catalogo `cell_types`, `nucleus_to_eq_ratio`, intervallo C8 [5, 25] µm) e in `full_test.R` (sparsità, UMI) **non** sono riferimenti: restano ipotesi finché una riga di questa tabella non le copre.

## Metriche previste per archetipo (colonne della tabella, da riempire in R2–R5)
Geometria (R2–R3): densità nucleare (cell/mm²), area nucleare (µm²), distanza al primo vicino (µm), g(r) / Ripley K, area ed eq_radius del territorio Voronoi (µm², µm), rapporto N/C, eccentricità del territorio.
Espressione (R4–R5): library size a 2 µm e 8 µm (UMI/bin), frazione di bin vuoti, sparsità, UMI nucleo vs citoplasma, `diffusion_sigma_um`.

## Convenzioni
- `id`: `B-###` progressivo. `Fonte`: DOI **oppure** `script@commit` + nome del dataset come in `data/real/README.md`.
- Una riga stimata da dati reali riporta sempre n (cellule/bin) e l'intervallo (IQR o IC 95 %).
