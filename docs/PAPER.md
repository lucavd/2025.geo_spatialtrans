---
type: Contract
title: "Contratto dell'articolo geo_spatialtrans"
description: "Contratto sul polmone 10x, scelto in SC2."
status: signed
generated: { by: geo-agent/progetto-master, at: 2026-11-20T21:00:00+01:00 }
---
# Messaggio

Un simulatore di Visium HD con verità per cellula riproduce un tumore polmonare quanto una seconda sezione e ordina i metodi come la verità.

Lettore: chi sviluppa o confronta metodi per Visium HD. Rivista: Genome Biology (A), Bioinformatics (B). Simulatore = catena banale più i pezzi che la battono.

# Dati

Polmone 10x (CG000709), un blocco FFPE, sezioni seriali: H1, H2 solo Visium HD; X1 (Xenium v1, 289 geni), X2 (Xenium 5K) con Visium HD sulla stessa sezione. Riserva: colon, `docs/archivio/PAPER_colon_SC1.md`.

# Misura principale

D (h), meglio verso il basso; metro di `docs/riesame/FATTI_B.md`: 12 voci, tumore e stroma × 6 misure dopo bin2cell + RCTD (HLCA). Unità: l'area (5); bootstrap B = 2000. Taratura su H2. Pavimento (banale) e tetto (H1 contro H2) da misurare. Soglia = tetto + δ, δ = 1,0 h (`docs/riesame/REGOLA_B.md`); sensibilità con δ = 0,5 e 2.

# Affermazioni

- A1* (Fig. 1): D ≤ soglia, tarato su H2, H1 come riferimento.
- A2 (Fig. 2): nessuna voce oltre il tetto più largo per voce.
- A3 (Fig. 3): batte lo Xenium in quadratini sulla stessa sezione (X1, X2), geni del pannello; scarto ≥ δ.
- A4 (Fig. 4): regge con una catena mai usata in taratura; ≤ suo tetto + δ.
- A5 (Fig. 5): classifica dei metodi sul simulato concorde con quella contro lo Xenium; ≥ concordanza X1 contro X2 meno un margine fissato prima.

# Cancelli scritti ora

Registrazione, prima cosa della sonda: Xenium, H&E, quadratini da 2 µm; residuo al 95° percentile ≤ 2 µm in 3 regioni disgiunte (affine locale come XD-4). Se fallisce, si torna al colon.

G1, con i componenti veri o spenti migliori:

- pavimento e tetto distanti meno di δ: chiudere;
- D ≤ soglia: costruzione;
- soglia < D ≤ pavimento − δ: riposizionare sull'esito B;
- D > pavimento − δ: chiudere.

# Esiti

- A: A1-A5 reggono.
- B: A1 non regge, D ≤ pavimento − δ; si dichiara lo scarto.
- C: G1 dice chiudere, oppure a costruzione finita D > pavimento − δ.

# Limiti accettati

- Un solo donatore: il tetto è fra sezioni seriali, non fra pazienti.
- Visium HD dopo Xenium meno sensibile (UMI per bin da 8 µm: 129 e 70 contro 212 e 275).
- La classifica dei metodi è misurata su un solo blocco.

# Dati tenuti da parte

Una metà spaziale di ogni sezione, scelta con una regola scritta prima di leggere i dati.

# Budget proposto

Scommessa 2 (SC1, SC2); sonda 5 (registrazione, check, pavimento e tetto, componenti veri, A3 e A5); costruzione 10 (5 voci × 2 colpi, confermato a G1); conferma 2; scrittura firmata a G3. Dati circa 25 GB.

# Decisioni del PI

- 2026-11-20 (SC2): polmone 10x dopo Xenium, cancello sulla registrazione, riserva il colon; G0 firmato, sonda 5.
