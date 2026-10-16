# R2c — note di sviluppo (prima delle annotazioni)

Ordine: pre-registrazione `1aa0eed` → primo run regola pre-registrata + addendum D-1..D-4 `cbd60c2` → questa fase (candidati v2,
annotatore, check C-2/C-3/C-5). Nessuna risposta manuale esiste a questo punto.

## Difetti trovati dai check prima delle annotazioni
- **D-1** (C-1 FAIL, A1 cellpose 3/99, A2 StarDist 1/30): i lati di Delaunay non coprono tutte le coppie a contatto → tolta la
  restrizione (addendum).
- **D-2** (C-5): d_max di A6 da 1 sola coppia (5.09 µm) → 0 coppie giudicabili → d_max combinato (addendum).
- **D-3** fascia d_max < d ≤ d_geo verificata a campione, **D-4** trappole oltre d_geo (addendum).
- **D-5** (C-2e FAIL al primo run): il grep di cecità trovava «catch» = parola chiave JavaScript `try { } catch (e)` (2
  occorrenze per file, verificate a mano); il controllo esclude ora solo il costrutto `catch (`, tutto il resto invariato.

## Dettagli di implementazione dichiarati
- Convenzione delle coordinate: i punti R2b sono continui con il pixel k = [k, k+1) (annotatore R2b: drawImage a 0,0); i
  centroidi delle maschere di S1.4c mettono il centro del pixel a k (scarto 0.5 px = 0.14 µm, preesistente, trascurabile per
  l'appaiamento ≤ 3 µm). I marcatori R2c usano la convenzione dei clic (C-2c: errore 1e-13 µm).
- Lato del ritaglio per archetipo = max(2.5 d_max, 30 µm): A1 30.0, A2 40.8, A3 46.6, A4 30.0, A5 38.2, A6 67.8 µm.
- A6: 59 focali invece di 60 (A6_w4 con 14 eleggibili: d_max 27.1 µm da bordo e da bolle escluse).
- Cecità: l'HTML non contiene tipo di item, finestra, archetipo o segmentatori (C-2e). Limiti noti: le ripetute usano la stessa
  immagine dell'originale (visibile solo leggendo il sorgente) e upp/lato del ritaglio differiscono per dataset; l'archetipo è
  comunque riconoscibile dall'istologia e non è un obiettivo di cecità.
- Ritagli PNG senza perdita (C-2a/b: identici al ROI e alla finestra R2b); zip 54.2 MB (Luca, 1 475 immagini per 1 607 item),
  28.4 MB (Federica, 763 per 834).
