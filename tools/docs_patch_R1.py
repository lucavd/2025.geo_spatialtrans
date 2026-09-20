#!/usr/bin/env python3
"""Applica le modifiche documentali della sessione R1 (idempotente: salta se il marker R1 è già presente)."""
import re, sys, pathlib
root = pathlib.Path.home() / "2025.geo_spatialtrans"
def patch(path, fn):
    p = root / path; s = p.read_text(); n = fn(s)
    if n != s: p.write_text(n); print("patched", path)
    else: print("unchanged", path)

# ---------- CHANGELOG ----------
CH = """## R1 — 2026-09-18/19 — Inventario e download dei dataset Visium HD di riferimento (branch `step1-cell-layer`)

### Dati
- Inventario di 46 pagine dataset Visium HD 10x (582 file con md5 pubblicato), letto con Claude in Chrome il 2026-09-18: `data/real/inventory/datasets_10x_pages_2026-09-18.csv`, `files_10x_pages_2026-09-18.csv`.
- Scaricati 5 dataset (38 file, 79.84 GB) su `/mnt/micron/geo_spatialtrans/data_real` (symlink `data/real/datasets`), md5 = pubblicato 38/38 (download + ricalcolo indipendente): A1 `Visium_HD_Mouse_Small_Intestine` (SR 3.0.0), A2/A3 `Visium_HD_Human_Colon_Cancer` (4.0.1), A4 `Visium_HD_Human_Lymph_Node_FFPE` (4.0.1), A5 `Visium_HD_6p5mm_Mouse_Brain` (4.1.0), A6 `Visium_HD_6p5mm_Human_Heart` (4.1.0). Niente FASTQ/cloupe/molecule_info.
- `data/real/README.md` (mappa dataset → archetipi, scala, verifica istologica), `data/real/checksums.md5`, `data/real/inventory/R1_download_manifest.csv`, `R1_archetype_map.csv`, `R1_download_status.tsv`.

### Codice / ambiente
- `tools/R1_download.sh` (download riprendibile con verifica md5), `tools/R1_verify.py` (metadati BigTIFF, scala, controprova griglia, ritagli), `R/testing/test_R1.sh` (un comando, PASS/FAIL; `R1_FAST=1` salta il md5).
- Ambiente Python di progetto: `tools/setup_python_env.sh` (uv, Python 3.12, `.venv/` ignorato) + `tools/requirements-r1.txt` / `tools/requirements.lock` (tifffile, imagecodecs, zarr, pyarrow, h5py, scikit-image).
- Tabelle di verifica tracciate in `results/R1/*.csv`; ritagli in `results/R1/crops/`.

### Esiti
- Check C 30/30 PASS; controprova 1 (scala dalla griglia dei bin) 5/5 entro 0.05 %; controprova 2 (istologia) confermata da Federica Pezzuto (anatomopatologa) per A1, A2, A4, A5, A6 con due rilievi (A4 ritaglio non rappresentativo; A6 bolle di montaggio).
- Deviazioni registrate: budget disco ridefinito (79.8 GB su micron); previsione B-R1.2 smentita (esistono dataset A4 e A6); soglia B-R1.3 senza fonte (il paper di Cellpose non fissa un minimo di µm/px).
- Reperto: i tag TIFF di risoluzione sono placeholder (96 dpi) in tutte le immagini → `estimate_pixel_size()` non può usarli.

### Documenti
- `reports/R1.qmd` + `R1.html`; `docs/sessions/2026-09-18_R1.md`; ROADMAP §1, §3, §6, §7; BACKLOG BL-014…BL-019; BIO_REFERENCES (stato).

"""
def ch(s):
    if "## R1 —" in s: return s
    return s.replace("## S0 — 2026-09-18", CH + "## S0 — 2026-09-18", 1)
patch("CHANGELOG.md", ch)

# ---------- ROADMAP ----------
def rm(s):
    if "| R1 | 2026-09-18" in s: return s
    s = s.replace("| R1 | — | prossima | — | inventario e download dataset Visium HD pubblici; nessuna dipendenza |",
        "| R1 | 2026-09-18/19 | **chiusa con riserve** | `reports/R1.html` | 46 pagine inventariate; 5 dataset / 79.8 GB su micron, md5 38/38; check C 30/30; controprova scala 5/5; istologia confermata (F. Pezzuto). Riserve: BL-014 (in_tissue ≠ archetipo), BL-015 (soglia µm/px senza fonte), BL-018 (A3 senza dataset dedicato), BL-019 (A1 solo SR 3.0.0) |\n| R2 | — | prossima | — | segmentazione nucleare su H&E (Cellpose/StarDist, GPU) per archetipo; dipende da R1 (chiusa) |")
    s = s.replace("**Prossima sessione: R1** (inventario dataset reali).", "**Prossima sessione: R2** (segmentazione nucleare; ambiente Python GPU da creare in `.venv` con `tools/setup_python_env.sh`) oppure **S1.1** (`extract_regions()`), senza dipendenze reciproche.")
    s = s.replace("Nessun dato reale ancora sul server (R1).", "Dopo R1 (2026-09-19): 5 dataset Visium HD reali su `/mnt/micron/geo_spatialtrans/data_real` (A1, A2/A3, A4, A5, A6; 79.8 GB; scala 0.274 µm/px verificata), inventario in `data/real/`.")
    # decisioni
    dec = ("| 2026-09-18 | (R1) Dati reali sul volume `/mnt/micron/geo_spatialtrans/data_real` (symlink `data/real/datasets`); budget R1 ridefinito a 79.8 GB prima del download |\n"
           "| 2026-09-18 | (R1) Mappa dataset → archetipi: A1 Mouse_Small_Intestine (SR 3.0.0), A2/A3 Human_Colon_Cancer (4.0.1), A4 Human_Lymph_Node_FFPE (4.0.1), A5 6p5mm_Mouse_Brain (4.1.0), A6 6p5mm_Human_Heart (4.1.0); A3 = regione stromale interna al CRC |\n"
           "| 2026-09-18 | (R1) Per ogni dataset si scaricano: immagine microscopica, immagine CytAssist, binned_outputs, spatial, segmented_outputs + barcode_mappings (SR ≥ 4), metrics, web_summary. Mai FASTQ/cloupe/molecule_info |\n"
           "| 2026-09-19 | (R1) La scala µm/px di un dataset reale viene da `scalefactors_json.json` verificata con la griglia dei bin; i tag TIFF non sono una fonte ammissibile |\n"
           "| 2026-09-19 | (R1) I riferimenti biologici R2–R5 si stimano su regioni annotate per archetipo, non sull'intera maschera `in_tissue` |\n")
    s = re.sub(r"(\| 2026-09-18 \| \(S0\) Il check S0[^\n]*\n)", r"\1" + dec, s, count=1)
    # tabella archetipi: colonna dataset reale
    s = s.replace("Le densità e i raggi qui sotto sono **ipotesi di partenza** da confermare in `BIO_REFERENCES.md` (binario R, sessione R1-R2), non valori accertati.",
                  "Le densità e i raggi qui sotto sono **ipotesi di partenza** da confermare in `BIO_REFERENCES.md` (binario R, sessione R2-R3), non valori accertati. Dataset reali assegnati in R1 (2026-09-18): A1 `Visium_HD_Mouse_Small_Intestine`; A2 e A3 `Visium_HD_Human_Colon_Cancer` (A3 = regioni stromali); A4 `Visium_HD_Human_Lymph_Node_FFPE`; A5 `Visium_HD_6p5mm_Mouse_Brain`; A6 `Visium_HD_6p5mm_Human_Heart`. Dettagli in `data/real/README.md`.")
    return s
patch("docs/ROADMAP.md", rm)

# ---------- BACKLOG ----------
BL = """| BL-014 | R1 | La maschera `in_tissue` di Space Ranger include materiale non parenchimale (linfonodo: ritaglio su grasso/capsula «sicuramente non linfociti», F. Pezzuto). I riferimenti R2–R3 per archetipo vanno stimati su regioni annotate (manuali o da `extract_regions()`), non sull'intera maschera. | R2 (annotazione regioni) ↔ S1.1 | aperta |
| BL-015 | R1 | Soglia B-R1.3 (µm/px ≤ 0.5 per la segmentazione nucleare) senza fonte: il paper di Cellpose (10.1038/s41592-020-01018-x) riporta solo il diametro mediano di training (30 px) e il resize a taglia comune. Sostituire con test empirico: segmentazione a risoluzione nativa (0.274 µm/px) vs sottocampionata 2×/4× sullo stesso ritaglio. | R2 | aperta |
| BL-016 | R1 | I tag TIFF di risoluzione (XResolution/ResolutionUnit) sono 96 dpi in tutte e 5 le immagini 10x: placeholder. `estimate_pixel_size()` (design C-A) non deve leggerli; fonti ammissibili: `scalefactors_json.json` + griglia dei bin, o input esplicito. Aggiornare il design §C-A. | S1.5 (design v1.2) | aperta |
| BL-017 | R1 | Cuore (A6): bolle di montaggio nell'H&E (confermate dalla patologa come artefatto di montaggio). Vanno rilevate e mascherate prima della segmentazione (cerchi grigi a bordo netto, diametro ~50–200 µm). | R2 | aperta |
| BL-018 | R1 | A3 (stroma lasso) non ha dataset 10x dedicato: coperto solo da regioni stromali del CRC (A2). Valutare se lo stroma desmoplastico tumorale è un surrogato accettabile dello stroma lasso, o cercare un secondo dataset (es. `Visium_HD_11mm_Human_TA` con polmone/mammella/colon). | R2–R3 | aperta |
| BL-019 | R1 | A1 (`Visium_HD_Mouse_Small_Intestine`) esiste solo con Space Ranger 3.0.0: niente `segmented_outputs`/`barcode_mappings`. Il confronto R2 «Cellpose vs segmentazione Space Ranger» sarà possibile per A2, A4, A5, A6 ma non per A1; in alternativa rielaborare A1 con Space Ranger 4 dai FASTQ (23.9 GB) o dai bin 2 µm. | R2 | aperta |
"""
def bl(s):
    if "| BL-014 |" in s: return s
    return re.sub(r"(\| BL-013 \|[^\n]*\n)", r"\1" + BL, s, count=1)
patch("docs/BACKLOG.md", bl)

# ---------- BIO_REFERENCES ----------
def br(s):
    if "R1 (2026-09-19)" in s: return s
    return s.replace("## Metriche previste per archetipo", "- **R1 (2026-09-19)**: dataset reali disponibili (vedi `data/real/README.md`), **ancora nessun valore biologico**. Unico dato accertato: scala 0.2737–0.2740 µm/px per le 5 immagini (da `scalefactors_json.json`, verificata con la griglia dei bin, `results/R1/R1_scale.csv`).\n\n## Metriche previste per archetipo", 1)
patch("docs/BIO_REFERENCES.md", br)
