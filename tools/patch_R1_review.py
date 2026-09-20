#!/usr/bin/env python3
"""Correzioni accolte dalla revisione avversariale R1 (idempotente). Eseguire dalla root del repo."""
import re, pathlib
def patch(path, pairs, must=True):
    p = pathlib.Path(path); s = p.read_text(); n = s
    for old, new in pairs:
        if old in n: n = n.replace(old, new, 1)
        elif must and new not in n: raise SystemExit(f"{path}: pattern non trovato: {old[:70]!r}")
    if n != s: p.write_text(n); print("patched", path)
    else: print("unchanged", path)

# --- tools/R1_verify.py: docstring, fit su tutti i bin, passo per asse, rotazione ---
patch("tools/R1_verify.py", [
 ("Scrive results/R1/*.csv, results/R1/crops/*.png, data/real/README.md, data/real/checksums.md5.",
  "Scrive results/R1/*.csv e results/R1/crops/*.png (README e checksums.md5 sono generati a parte)."),
 ('''    p = pos.dropna(subset=["pxl_col_in_fullres", "pxl_row_in_fullres"])
    if len(p) > 300_000:
        p = p.sample(300_000, random_state=42)
''', '''    p = pos.dropna(subset=["pxl_col_in_fullres", "pxl_row_in_fullres"])   # tutti i bin (rilievo revisore R1)
'''),
 ('''    return dict(pitch_px_sv1=sv[0], pitch_px_sv2=sv[1], pitch_px=float(np.sqrt(abs(np.linalg.det(A)))),
                rms_resid_px=float(np.sqrt((resid ** 2).mean())), n_fit=len(p))''',
  '''    ax1 = float(np.hypot(*A[:, 0])); ax2 = float(np.hypot(*A[:, 1]))     # passo lungo array_col e array_row
    rot = float(np.degrees(np.arctan2(A[1, 0], A[0, 0])))                  # orientamento dell'asse array_col nell'immagine
    return dict(pitch_px_sv1=sv[0], pitch_px_sv2=sv[1], pitch_px=float(np.sqrt(abs(np.linalg.det(A)))),
                pitch_px_axis1=ax1, pitch_px_axis2=ax2, grid_rotation_deg=rot,
                rms_resid_px=float(np.sqrt((resid ** 2).mean())), n_fit=len(p))'''),
 ('''            pitch_um = gp["pitch_px"] * upp; sc["pitch_um_from_grid"] = pitch_um''',
  '''            pitch_um = gp["pitch_px"] * upp; sc["pitch_um_from_grid"] = pitch_um
            sc["pitch_um_axis1"] = gp["pitch_px_axis1"] * upp; sc["pitch_um_axis2"] = gp["pitch_px_axis2"] * upp'''),
])

# --- README dati ---
patch("data/real/README.md", [
 ("- **I tag TIFF di risoluzione sono 96 dpi (= 264.6 µm/px) in tutte le immagini: placeholder, non usabili.** `estimate_pixel_size()` (C-A) non può basarsi sui metadati dell'immagine (BL-016).",
  "- Metadati immagine: il tag TIFF `XResolution` è 96 dpi (= 264.6 µm/px) in tutte le immagini — placeholder, non usabile; l'OME-XML `PhysicalSizeX` è invece presente e corretto (0.27377–0.27381 µm, entro 0.08 % da `scalefactors_json`; rilievo del revisore R1). `estimate_pixel_size()` (C-A) può leggere l'OME-XML ma deve verificarlo contro la griglia dei bin, mai fidarsi del tag di risoluzione (BL-016)."),
 ("Fonte: dataset pubblici 10x Genomics Visium HD, pagine lette il 2026-09-18 (inventario completo: `inventory/datasets_10x_pages_2026-09-18.csv`, 46 pagine; `inventory/files_10x_pages_2026-09-18.csv`, 582 file con md5 pubblicato).",
  "Fonte: dataset pubblici 10x Genomics Visium HD, pagine lette il 2026-09-18 (inventario completo: `inventory/datasets_10x_pages_2026-09-18.csv`, 46 pagine = 32 dataset_id, 14 con due elaborazioni Space Ranger; `inventory/files_10x_pages_2026-09-18.csv`, 582 URL con md5 pubblicato, 526 contenuti distinti perché alcuni input sono condivisi fra pagine — confrontare i md5 per URL, non per nome file)."),
 ("## Riproduzione", "## Revisione avversariale (2026-09-20)\n\nSotto-agente indipendente in sola lettura: 30/30 affermazioni confermate (md5 completo, µm/px, passo griglia con metodo alternativo, dimensioni immagini, conteggi bin, inventario) — `results/R1/R1_adversarial_review.csv`. Otto rilievi non bloccanti: accolti 1, 3, 4, 6, 7; in BACKLOG BL-020 (orientamento array→pixel: A6 ruotato di 180°, rotazione ≤ 0.75° in A2), BL-021 (0.30–0.36 % dei bin `in_tissue` fuori immagine in A2/A4), BL-022 (A2 con seconda pagina TIFF illeggibile; A4 TIFF classico a strisce), BL-023 (versione Space Ranger dal web_summary via euristica).\n\n## Riproduzione"),
])

# --- CHANGELOG ---
patch("CHANGELOG.md", [
 ("- Reperto: i tag TIFF di risoluzione sono placeholder (96 dpi) in tutte le immagini → `estimate_pixel_size()` non può usarli.",
  "- Reperto: il tag TIFF di risoluzione è placeholder (96 dpi) in tutte le immagini, mentre l'OME-XML `PhysicalSizeX` è corretto (entro 0.08 % da `scalefactors_json`) → `estimate_pixel_size()` può usare l'OME-XML, mai il tag, e deve verificare contro la griglia dei bin.\n- Revisione avversariale (2026-09-20, sotto-agente indipendente): 30/30 affermazioni confermate, 0 discrepanze, 8 rilievi non bloccanti (5 accolti, 4 in BACKLOG BL-020…BL-023)."),
 ("- `reports/R1.qmd` + `R1.html`; `docs/sessions/2026-09-18_R1.md`; ROADMAP §1, §3, §6, §7; BACKLOG BL-014…BL-019; BIO_REFERENCES (stato).",
  "- `reports/R1.qmd` + `R1.html`; `docs/sessions/2026-09-18_R1.md`; ROADMAP §1, §3, §6, §7; BACKLOG BL-014…BL-023; BIO_REFERENCES (stato); `results/R1/R1_cp2_histology.csv` (verifica istologica), `results/R1/R1_adversarial_review.csv`."),
])

# --- BACKLOG ---
BL = """| BL-020 | R1 (revisore) | Orientamento array→pixel non uniforme fra dataset: in A1–A5 `pxl_col` cresce con `array_col` e `pxl_row` decresce con `array_row`; in A6 i segni sono invertiti (rotazione di 180°); rotazione residua fino a ≈ 0.75° (A2). Documentare e gestire nella mappa bin↔immagine prima di R2 (maschere, regioni). | R2 | aperta |
| BL-021 | R1 (revisore) | 0.36 % (A2, ≈ 29 000) e 0.30 % (A4) dei bin `in_tissue` cadono fuori dall'immagine microscopica: l'area di cattura sporge dal campo acquisito. Da escludere esplicitamente nell'allineamento immagine↔bin. | R2 | aperta |
| BL-022 | R1 (revisore) | Formati immagine eterogenei: A2 `tissue_image.btf` ha una seconda pagina TIFF non leggibile da tifffile; A4 è TIFF classico (non BigTIFF) a strisce. Il lettore immagini di R2 deve gestire entrambi e leggere sempre la serie 0. | R2 | aperta |
| BL-023 | R1 (revisore) | La versione di Space Ranger è letta dal `web_summary.html` con un'euristica (stringa 3.x.y/4.x.y più frequente) quando manca un'etichetta: fragile. Leggere la versione da una fonte strutturata (JSON del web_summary o `metrics_summary`) e confrontarla con la pagina 10x. | R2 | aperta |
"""
patch("docs/BACKLOG.md", [
 ("| BL-016 | R1 | I tag TIFF di risoluzione (XResolution/ResolutionUnit) sono 96 dpi in tutte e 5 le immagini 10x: placeholder. `estimate_pixel_size()` (design C-A) non deve leggerli; fonti ammissibili: `scalefactors_json.json` + griglia dei bin, o input esplicito. Aggiornare il design §C-A. | S1.5 (design v1.2) | aperta |",
  "| BL-016 | R1 | Il tag TIFF di risoluzione (XResolution/ResolutionUnit) è 96 dpi in tutte e 5 le immagini 10x: placeholder. L'OME-XML `PhysicalSizeX` è invece corretto (entro 0.08 % da `scalefactors_json`, rilievo del revisore R1). `estimate_pixel_size()` (design C-A): fonti ammissibili in ordine = input esplicito, `scalefactors_json.json` + griglia dei bin, OME-XML verificato contro la griglia; mai il tag di risoluzione. Aggiornare il design §C-A. | S1.5 (design v1.2) | aperta |"),
])
s = pathlib.Path("docs/BACKLOG.md").read_text()
if "| BL-020 |" not in s:
    s = re.sub(r"(\| BL-019 \|[^\n]*\n)", r"\1" + BL, s, count=1); pathlib.Path("docs/BACKLOG.md").write_text(s); print("patched docs/BACKLOG.md (BL-020..023)")

# --- ROADMAP: decisione sulla scala ---
patch("docs/ROADMAP.md", [
 ("| 2026-09-19 | (R1) La scala µm/px di un dataset reale viene da `scalefactors_json.json` verificata con la griglia dei bin; i tag TIFF non sono una fonte ammissibile |",
  "| 2026-09-19 | (R1) La scala µm/px di un dataset reale viene da `scalefactors_json.json` verificata con la griglia dei bin; l'OME-XML `PhysicalSizeX` è ammesso solo se verificato contro la griglia; il tag TIFF di risoluzione non è mai una fonte |"),
 ("Riserve: BL-014 (in_tissue ≠ archetipo), BL-015 (soglia µm/px senza fonte), BL-018 (A3 senza dataset dedicato), BL-019 (A1 solo SR 3.0.0) |",
  "Riserve: BL-014 (in_tissue ≠ archetipo), BL-015 (soglia µm/px senza fonte), BL-018 (A3 senza dataset dedicato), BL-019 (A1 solo SR 3.0.0); revisione avversariale 30/30 confermate, BL-020…BL-023 |"),
])

# --- nota di sessione ---
patch("docs/sessions/2026-09-18_R1.md", [
 ("| Revisione avversariale | vedi `results/R1/R1_adversarial_review.csv` e report | |",
  "| Revisione avversariale | **PASS** | sotto-agente indipendente (2026-09-20), sola lettura: 30/30 affermazioni confermate con strumenti diversi (md5sum completo 79.84 GB, gdalinfo, mediana Δpxl per il passo), 0 discrepanze; 8 rilievi non bloccanti → 5 accolti prima del commit (OME-XML PhysicalSizeX corretto; fit CP-1 su tutti i bin e per asse; CP-2 in `R1_cp2_histology.csv`; docstring; conteggi inventario 46 pagine / 582 URL / 526 contenuti), 4 in BACKLOG (BL-020…BL-023) |"),
 ("- BL-018 A3 senza dataset dedicato. BL-019 A1 solo SR 3.0.0 (niente segmentazione 10x).",
  "- BL-018 A3 senza dataset dedicato. BL-019 A1 solo SR 3.0.0 (niente segmentazione 10x).\n- Dal revisore: BL-020 orientamento array→pixel (A6 ruotato 180°, rotazione ≤ 0.75° in A2); BL-021 0.30–0.36 % bin in_tissue fuori immagine (A2, A4); BL-022 formati TIFF eterogenei; BL-023 versione Space Ranger da fonte strutturata."),
])
print("done")
