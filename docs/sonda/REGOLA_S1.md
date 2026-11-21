---
type: Rule
title: "S1, cancello della registrazione: metà tenuta da parte e residuo"
description: "Regole scritte e committate prima di scaricare i dati del polmone: metà spaziale tenuta da parte per ogni sezione, regioni, misura del residuo della registrazione, controlli, esito."
status: stable
generated: { by: geo-agent/progetto-master, at: 2026-11-21T08:00:00+01:00 }
sources:
  - id: paper
    path: docs/PAPER.md
  - id: xd4
    path: docs/riesame/DECISIONI_B_XENIUM.md
---
# Perché

Il contratto (`docs/PAPER.md`[^paper], «Cancelli scritti ora») chiede, prima di ogni altro lavoro della sonda, che la catena Xenium -> H&E -> quadratini da 2 µm abbia un residuo al 95° percentile di al massimo 2 µm in 3 regioni disgiunte, con l'affine locale ammessa come in XD-4[^xd4]. Se fallisce in X1 e in X2 si torna al colon.

Deviazione da `STATO.md` (prossima mossa, punto b): MALAT1 e NEAT1 non sono nel pannello di sonde del Visium HD (Human Transcriptome v2.0, 0 righe); le sonde «unspliced» sono esoniche, non intronici. Il residuo della seconda tratta non si misura quindi con i picchi nucleari. Lo sostituisce una misura da capo a fondo sulla stessa sezione: trascritti Xenium portati sulla griglia contro UMI del Visium HD sugli stessi geni (286 geni in X1, 4892 in X2). Misura la catena intera, che è ciò che il cancello giudica.

# Metà tenuta da parte (vale per H1, H2, X1, X2)

- Quadratini da 8 µm `in_tissue` di Space Ranger (solo la maschera del tessuto, nessuna espressione). x = (array_col + 0,5) × 8 µm. m = mediana di x sui quadratini in tissue.
- Sviluppo: x < m − 50 µm. Tenuta da parte: x > m + 50 µm. La striscia centrale di 100 µm non appartiene a nessuna delle due.
- Ai quadratini da 2 µm si applica la stessa soglia in µm sul loro centro. Nessuno script legge l'espressione della metà tenuta da parte fino alla conferma.

# Catena

1. Xenium (µm) -> H&E a piena risoluzione: matrice 3×3 di `he_imagealignment.csv`. Varianti: matrice diretta o inversa, unità Xenium in pixel da 0,2125 µm o in µm (4 varianti). Si sceglie quella che porta la quota più alta di centroidi delle cellule Xenium su quadratini da 8 µm in tissue; deve essere almeno 0,8 e superare la seconda di almeno 0,2, altrimenti ci si ferma.
2. H&E -> quadratini da 2 µm: `tissue_positions.parquet` di `square_002um` (centro di ogni quadratino in pixel a piena risoluzione); trasformazione proiettiva stimata su tutti i quadratini, usata se il residuo massimo è sotto 0,01 pixel, altrimenti interpolazione.
3. Le due immagini H&E (`he_image.ome.tif` dello Xenium e `tissue_image.btf` del Visium HD) sono la stessa immagine se hanno le stesse dimensioni e la correlazione di fase su 9 ritagli da 2048 pixel dà spostamento nullo (≤ 1 pixel). Se non lo sono, si aggiunge un'affine stimata con la correlazione di fase sui 9 ritagli; il cancello si giudica comunque sulla catena intera.

# Regioni

Quadrati da 1 mm × 1 mm interamente nella metà di sviluppo, con angolo su un reticolo da 100 µm. Coperti: quadratini da 8 µm in tissue con almeno un centroide Xenium (catena del punto 1-2) nel loro intorno 3 × 3. Candidati: quota coperta ≥ 0,9 (se non bastano, 0,8, poi lato 750 µm; ogni rilassamento si annota). R1 = candidato più vicino al baricentro dei quadratini coperti di sviluppo; R2 = disgiunto da R1 e più lontano da R1; R3 = disgiunto da entrambi e con la distanza minima da R1 e R2 più grande. Parità: x minore, poi y minore.

# Misura del residuo

- Trascritti Xenium con `qv` ≥ 20 e identificativo Ensembl nel pannello di sonde Visium HD; UMI del Visium HD da `square_002um/filtered_feature_bc_matrix.h5` sugli stessi geni.
- Spostamento grossolano per regione: correlazione incrociata delle mappe totali a 8 µm, finestra ± 200 µm.
- Riquadri: la regione si divide in riquadri di lato L. In ogni riquadro, geni con almeno 10 conteggi in entrambe le mappe; ogni mappa per gene standardizzata; CC(s) = somma su geni e quadratini di zX_g(x + s) · zV_g(x), s intero in quadratini nella finestra ± 12 µm attorno allo spostamento grossolano; picco raffinato con una quadrica sul 3 × 3 (se il vertice esce dal 3 × 3, il massimo intero). Riquadro valido se entrambe le mappe hanno almeno 500 conteggi sui geni usati. Picco sul bordo della finestra: residuo > 2 µm.
- Residuo del riquadro = |s| in µm (quadratini da 2 µm).
- Precisione: per ogni riquadro, trascritti Xenium divisi a caso in due metà e UMI Visium assottigliate con una binomiale 1/2 (seme 20261121); errore standard radiale ≈ deviazione standard delle differenze fra le due metà, diviso 2. L = il più piccolo fra 100, 200 e 250 µm per cui la media quadratica dell'errore standard radiale è ≤ 0,5 µm in tutte e 3 le regioni della sezione. Si sceglie prima di guardare i residui.

# Affine locale (XD-4)

Se il 95° percentile globale supera 2 µm in una regione: affine a 6 parametri dei centri dei riquadri verso centri + s, stimata sui riquadri pari della scacchiera e valutata sui dispari, poi il contrario. Il residuo è quello dei riquadri non usati nella stima.

# Controlli del metro

- Positivo: spostamenti noti aggiunti ai trascritti Xenium, (1,0; 0), (0; 2,5), (1,5; −1,5) µm; la mediana fra riquadri della differenza di s deve ritrovarli entro 0,3 µm per asse in ogni regione.
- Negativo: riquadro Xenium della regione r contro riquadro Visium con lo stesso indice della regione r + 1 (ciclico); al più il 10% di questi riquadri con |s − grossolano| ≤ 2 µm.
- Se un controllo fallisce il metro non vale: il cancello non si giudica finché non lo si ripara.

# Esito

- Una sezione passa se in tutte e 3 le regioni il 95° percentile dei residui è ≤ 2 µm, globale oppure dopo l'affine locale (riquadri non usati nella stima). Si riportano anche le quote di riquadri entro 1, 2 e 5 µm.
- X1 e X2 passano: si prosegue. Una sola: si prosegue con quella e il tetto di A5 si riscrive prima di misurarlo. Nessuna: ritorno al colon.
- Solo se una sezione fallisce: nuclei StarDist sull'H&E contro nuclei Xenium (tratta 1) per dire dove si perde. Diagnosi, nessuna scelta dipende da essa.

# Limiti dichiarati ora

- Il residuo è la media sul riquadro: deformazioni più piccole del riquadro non si vedono.
- Se serve l'affine locale, è stimata anche con l'espressione del Visium HD: lieve circolarità per A3, da dichiarare.

[^paper]: docs/PAPER.md
[^xd4]: docs/riesame/DECISIONI_B_XENIUM.md, XD-4

# Errore di pre-registrazione (2026-11-21, prima di misurare il residuo)

La regola della variante (Catena, punto 1) si è fermata: la variante migliore, matrice inversa con pixel Xenium da 0,2125 µm, porta su quadratini in tissue il 62,2% delle cellule Xenium in X1 e il 59,0% in X2 (soglia 0,8); distacco dalla seconda 0,36 e 0,39 (soglia 0,2). Causa: lo Xenium copre più tessuto dell'area di cattura del Visium HD; il denominatore contava anche le cellule fuori dalla griglia.

| variante | dentro la griglia X1 / X2 | su tessuto fra quelle dentro | tessuto Visium con una cellula Xenium vicina |
|---|---|---|---|
| inversa, px 0,2125 | 0,622 / 0,590 | 1,000 / 0,999 | 0,871 / 0,874 |
| diretta, µm | 0,429 / 0,228 | 0,601 / 0,884 | 0,094 / 0,082 |
| diretta, px 0,2125 | 0,178 / 0,094 | 0,609 / 0,456 | 0,340 / 0,139 |
| inversa, µm | 0,376 / 0,000 | 0,004 / - | 0,000 / 0,000 |

Decisione di Luca, in chat il 2026-11-21: «Proseguire con la variante scelta». Si usa la variante che la regola mette prima; il resto della regola resta com'è scritto. Nessun residuo era stato calcolato al momento della decisione.

Altri fatti della preparazione: le due immagini H&E sono la stessa immagine (stesse dimensioni, spostamento nullo nei 9 ritagli) in X1 e X2; la proiettiva quadratini -> pixel riproduce Space Ranger entro 2e-11 pixel; il centro di un quadratino da 8 µm dista circa 7,4 pixel (2 µm) dalla media dei suoi 2 µm presunti (aggregazione non allineata a floor(i/4)): irrilevante per la divisione in metà (striscia di 100 µm) e per la copertura.
