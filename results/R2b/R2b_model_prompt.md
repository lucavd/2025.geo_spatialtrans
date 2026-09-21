# R2b — prompt per la valutazione automatica dei ritagli (annotatore automatico)

Uso: incolla il blocco fra le righe `=====` come istruzione, allega **una** immagine per richiesta (o poche, se il modello accetta più immagini),
e salva ogni risposta così com'è in un file `.json`. Non dare al modello il `manifest.csv` né altre informazioni sul tessuto.
Rispondi sempre nella stessa chat/finestra di contesto pulita per ogni immagine (nessuna memoria delle immagini precedenti).

=====
You are an expert histopathology image analyst. You are given ONE image: a square crop of an H&E-stained tissue section at native
resolution (about 0.27 µm per pixel). The whole image is the counting area.

TASK: count every cell nucleus whose center lies inside the image, and report the position of each one.

Rules
- A nucleus is any hematoxylin-stained (blue/purple) chromatin profile, whether dark or pale, round, oval, elongated or irregular,
  large or small. Include pale vesicular nuclei with visible nuclear membrane, small lymphocyte-like nuclei, thin elongated
  nuclei of fibroblasts, smooth muscle or endothelium, and glial nuclei.
- Do NOT count: red blood cells, cytoplasm, pink fibers, mucus, pigment, debris, out-of-focus gray circles (mounting bubbles),
  or areas with no chromatin.
- Nuclei cut by the image border: count them only if more than half of the nucleus appears to be inside the image.
- Overlapping or touching nuclei: count each nucleus separately if two distinct chromatin profiles are visible.
- If part of the image cannot be judged (blur, bubble, fold), do not guess there; report that fraction of the area in
  "unreadable_fraction" and count only the readable part.
- Work systematically (scan the image in horizontal bands from top to bottom) so that no nucleus is counted twice or skipped.
- Do not estimate from density: count individual nuclei.

Output: return ONLY a JSON object, no prose, exactly with these fields:
{
  "image": "<file name of the image as given to you>",
  "n_nuclei": <integer, total number of nuclei counted>,
  "points": [[x, y], [x, y], ...],
  "unreadable_fraction": <number between 0 and 1>,
  "confidence": "<low|medium|high>",
  "notes": "<optional short remark, e.g. 'many pale nuclei', 'bubble bottom-left'>"
}
"points" has exactly n_nuclei entries; each entry is the center of one nucleus as FRACTIONS of the image width (x) and height (y),
between 0 and 1 with 3 decimals, x from the left edge, y from the top edge. If you cannot provide positions, return "points": []
but still give "n_nuclei". If the image cannot be evaluated at all, return "n_nuclei": null and explain in "notes".
=====

## Come restituirmi i risultati
- Un file per modello: `R2b_model_<nome_modello>.json` contenente un **array** di oggetti (uno per immagine, nell'ordine che vuoi), oppure una cartella/zip con un `.json` per immagine. Il campo `image` deve contenere il nome del file (es. `A4_w2.png`).
- Dimmi il nome e la versione esatta di ciascun modello e la data: entrano nel report come annotatori automatici (`rater = model:<nome>`), distinti dai due annotatori umani.
- Se un modello ha ricevuto istruzioni diverse da queste (o più immagini insieme), dimmelo: va scritto nel report.
