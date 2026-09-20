# data/real — dataset reali di riferimento (sessione R1, 2026-09-18)

I dati stanno su `/mnt/micron/geo_spatialtrans/data_real/` (volume `micron`, XFS, 3.5 TB) e sono raggiunti dal repo tramite il symlink `data/real/datasets` (fuori git). In git restano solo questo README, `checksums.md5` e `inventory/`.

Fonte: dataset pubblici 10x Genomics Visium HD, pagine lette il 2026-09-18 (inventario completo: `inventory/datasets_10x_pages_2026-09-18.csv`, 46 pagine = 32 dataset_id, 14 con due elaborazioni Space Ranger; `inventory/files_10x_pages_2026-09-18.csv`, 582 URL con md5 pubblicato, 526 contenuti distinti perché alcuni input sono condivisi fra pagine — confrontare i md5 per URL, non per nome file). Licenza dichiarata in pagina: CC BY 4.0 per tutti i dataset scaricati.

## Mappa dataset → archetipi (decisione R1, approvata da Luca il 2026-09-18)

| Arch. | Dataset (Space Ranger) | Tessuto | Specie | Preserv. | Pubbl. | File | GB | Immagine microscopica (px) | µm/px (`scalefactors_json`) | Passo griglia 2 µm ricavato | Bin 2 µm in tessuto | Pagina |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A1 | `Visium_HD_Mouse_Small_Intestine` (3.0.0) | Intestine | Mouse | FFPE | 2024-03-25 | 6 | 7.8 | `Visium_HD_Mouse_Small_Intestine_tissue_image.btf` (21943, 23618, 3) | 0.27376 | 1.9995 µm | 5,479,660 | [pagina](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-mouse-intestine) |
| A2 | `Visium_HD_Human_Colon_Cancer` (4.0.1) | Colon | Human | FFPE | 2025-07-03 | 8 | 36.7 | `Visium_HD_Human_Colon_Cancer_tissue_image.btf` (48740, 75250, 3) | 0.27370 | 1.9991 µm | 8,176,117 | [pagina](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-human-crc-v4) |
| A3 | `Visium_HD_Human_Colon_Cancer` (4.0.1) | Colon | Human | FFPE | 2025-07-03 | 8 | 36.7 | `Visium_HD_Human_Colon_Cancer_tissue_image.btf` (48740, 75250, 3) | 0.27370 | 1.9991 µm | 8,176,117 | [pagina](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-human-crc-v4) |
| A4 | `Visium_HD_Human_Lymph_Node_FFPE` (4.0.1) | Lymph Node | Human | FFPE | 2025-07-03 | 8 | 7.2 | `Visium_HD_Human_Lymph_Node_FFPE_tissue_image.tif` (26128, 28248, 3) | 0.27378 | 2.0000 µm | 5,482,441 | [pagina](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-human-lymph-node-v4) |
| A5 | `Visium_HD_6p5mm_Mouse_Brain` (4.1.0) | Brain | Mouse | FFPE | 2026-07-14 | 8 | 12.1 | `Visium_HD_6p5mm_Mouse_Brain_tissue_image.btf` (22074, 22354, 3) | 0.27402 | 1.9998 µm | 6,343,031 | [pagina](https://www.10xgenomics.com/datasets/visium-hd-cytassist-6p5mm-mouse-brain) |
| A6 | `Visium_HD_6p5mm_Human_Heart` (4.1.0) | Heart | Human | FFPE | 2026-04-28 | 8 | 16.1 | `Visium_HD_6p5mm_Human_Heart_tissue_image.btf` (49099, 39139, 3) | 0.27381 | 2.0004 µm | 10,615,406 | [pagina](https://www.10xgenomics.com/datasets/visium-hd-cytassist-6p5mm-human-heart) |

Note per archetipo:

- **A1** — epitelio semplice denso (mucosa intestino tenue, cripte/villi); muscularis = regione candidata A6 secondaria
- **A2** — CRC FFPE (campione Nature Genetics 2025); contiene mucosa normale adiacente (A1 umano), stroma (A3) e aggregati linfoidi (A4 secondario)
- **A3** — stroma/desmoplasia come regione interna di A2 — NESSUN dataset dedicato
- **A4** — linfonodo: follicoli e paracorticale; riserva: Visium_HD_Human_Tonsil_Fresh_Frozen (3.1.1, 8 GB tif + 17 GB binned)
- **A5** — corteccia (neuroni+glia); preferito a Visium_HD_Mouse_Brain la cui immagine microscopica è 145 MB
- **A6** — cardiomiociti allungati = controllo negativo anisotropo; nota pagina: microbolle nell H&E

Copertura: A1, A2, A4, A5, A6 con dataset dedicato; **A3 (stroma) non ha dataset dedicato** e sarà estratto come regione interna del CRC (A2). La pre-registrazione R1 prevedeva l'assenza di dataset dedicati per A4 e A6: previsione **smentita** (esistono linfonodo FFPE e cuore FFPE).

## Verifica istologica dei ritagli (controprova 2, 2026-09-19)

Ritagli H&E a risoluzione nativa (400 µm, `results/R1/crops/`) rivisti da **Federica Pezzuto** (anatomopatologa) tramite Luca Vedovelli:

- A1 intestino tenue murino: confermato. - A2 CRC umano: confermato (ritagli: stroma con vaso, mucosa normale con cellule caliciformi, ghiandole tumorali).
- A4 linfonodo: confermato; il ritaglio 1 **non è rappresentativo** («qualcosa di coagulato, grasso o capsula, comunque sicuramente non linfociti») pur essendo dentro la maschera `in_tissue` di Space Ranger → la maschera `in_tissue` non equivale a "tessuto dell'archetipo" (BL-014).
- A5 corteccia murina: confermato. - A6 cuore umano: confermato; i cerchi grigi sono **artefatti di montaggio** (bolle sotto il coprivetrino), da mascherare in R2 (BL-017).

## File scaricati per dataset

Immagine microscopica H&E (`*_tissue_image.btf|tif`; Olympus VS200 20×/0.8 NA secondo la pagina 10x per A1, A4, A5, A6 — la pagina del CRC (A2/A3) non ha il blocco Imaging: microscopio non dichiarato), immagine CytAssist (`*_image.tif`, 3000×3200), `binned_outputs.tar.gz` (2/8/16 µm), `spatial.tar.gz`, `segmented_outputs.tar.gz` e `barcode_mappings.parquet` (solo Space Ranger ≥ 4), `metrics_summary.csv`, `web_summary.html`. **Non** scaricati: FASTQ, `.cloupe`, `molecule_info.h5`. Totale 38 file, 79.84 GB. I tar sono estratti accanto all'originale in `spatial_toplevel/`, `binned_outputs_x/`, `segmented_outputs_x/` (marker `.extracted`).

## Integrità e scala

- `checksums.md5`: md5 di ogni file **ricalcolato dopo il download** e uguale al md5 pubblicato da 10x (38/38, `results/R1/R1_md5.csv`). Verifica: `cd data/real && md5sum -c checksums.md5` (~15 min).
- Scala: `microns_per_pixel` da `square_002um/spatial/scalefactors_json.json` (0.2737–0.2740 µm/px). Controprova indipendente: passo della griglia 2 µm stimato per minimi quadrati da `tissue_positions.parquet` × µm/px = 1.9991–2.0004 µm (anisotropia ≤ 1.0001).
- Metadati immagine: il tag TIFF `XResolution` è 96 dpi (= 264.6 µm/px) in tutte le immagini — placeholder, non usabile; l'OME-XML `PhysicalSizeX` è invece presente e corretto (0.27377–0.27381 µm, entro 0.08 % da `scalefactors_json`; rilievo del revisore R1). `estimate_pixel_size()` (C-A) può leggere l'OME-XML ma deve verificarlo contro la griglia dei bin, mai fidarsi del tag di risoluzione (BL-016).

## Revisione avversariale (2026-09-20)

Sotto-agente indipendente in sola lettura: 30/30 affermazioni confermate (md5 completo, µm/px, passo griglia con metodo alternativo, dimensioni immagini, conteggi bin, inventario) — `results/R1/R1_adversarial_review.csv`. Otto rilievi non bloccanti: accolti 1, 3, 4, 6, 7; in BACKLOG BL-020 (orientamento array→pixel: A6 ruotato di 180°, rotazione ≤ 0.75° in A2), BL-021 (0.30–0.36 % dei bin `in_tissue` fuori immagine in A2/A4), BL-022 (A2 con seconda pagina TIFF illeggibile; A4 TIFF classico a strisce), BL-023 (versione Space Ranger dal web_summary via euristica).

## Riproduzione

```bash
bash tools/R1_download.sh data/real/inventory/R1_download_manifest.csv /mnt/micron/geo_spatialtrans/data_real /mnt/micron/geo_spatialtrans/data_real/_logs/R1_download_status.tsv   # riprendibile
bash tools/setup_python_env.sh   # .venv (uv) con tifffile, pyarrow, h5py, scikit-image
bash R/testing/test_R1.sh        # verifiche C/B/controprova; R1_FAST=1 salta il ricalcolo md5
```
