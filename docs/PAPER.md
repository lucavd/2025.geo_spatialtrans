---
type: Contract
title: "Contratto dell'articolo geo_spatialtrans"
description: "Contratto sul polmone 10x."
status: signed
generated: { by: geo-agent/progetto-master, at: 2026-11-20T21:00:00+01:00 }
---
# Messaggio

Un simulatore di Visium HD con verità per cellula riproduce un tumore polmonare quanto una seconda sezione e ordina i metodi come la verità.

Lettore: chi sviluppa o confronta metodi Visium HD. Rivista: Genome Biology (A), Bioinformatics (B). Simulatore = banale più i pezzi che lo battono.

# Dati

Polmone 10x (CG000709), un blocco FFPE, sezioni seriali: H1, H2 solo Visium HD; X1 (Xenium v1, 289 geni), X2 (Xenium 5K) con Visium HD sulla stessa sezione. Riserva: colon, `docs/archivio/PAPER_colon_SC1.md`.

# Misura principale

D (h), meglio verso il basso; metro di `docs/riesame/FATTI_B.md`: 12 voci, tumore e stroma × 6 misure dopo bin2cell + RCTD (HLCA). Unità: l'area (5); bootstrap B = 2000. Soglia = tetto + δ, δ = 1,0 h (`docs/riesame/RIESAME_C.md` §1.4); sensibilità δ = 0,5.

# Affermazioni

- A1* (Fig. 1): D ≤ soglia, tarato su H2, H1 come riferimento.
- A2 (Fig. 2): nessuna voce oltre il tetto più largo per voce.
- A4 (Fig. 4): regge con una catena mai usata in taratura; ≤ suo tetto + δ.
- A5 (Fig. 5): classifica di 5 regole di assegnazione (ENACT, bin2cell) sul simulato concorde con quella contro lo Xenium di X1: ≥ concordanza di X1 fra metà disgiunte − 0,2.

# Cancelli scritti ora

Registrazione, prima cosa della sonda: Xenium, H&E, quadratini da 2 µm; residuo al 95° percentile ≤ 2 µm in 3 regioni disgiunte (affine locale come XD-4). Se fallisce, si torna al colon.

G1, con i componenti veri o spenti migliori:

- pavimento e tetto distanti meno di δ: chiudere;
- D ≤ soglia: costruzione;
- soglia < D ≤ pavimento − δ: riposizionare sull'esito B;
- D > pavimento − δ: chiudere.

# Esiti

- A: A1, A2, A4, A5 reggono.
- B: A1 non regge, D ≤ pavimento − δ; si dichiara lo scarto.
- C: G1 dice chiudere, oppure a costruzione finita D > pavimento − δ.

# Limiti accettati

- Un solo donatore: il tetto è fra sezioni seriali, non fra pazienti.
- Visium HD dopo Xenium meno sensibile (UMI per bin da 8 µm: 129 e 70 contro 212 e 275).
- La classifica dei metodi è misurata su un solo blocco.
- Nessun confronto con lo Xenium in quadratini: 284 geni, 2 % degli UMI, la catena non ci gira (S5).

# Dati tenuti da parte

Metà di H2, X1, X2 (S1); H1 intera (S3b).

# Budget proposto

Scommessa 2; sonda 5; costruzione 10 (5 voci × 2 colpi, confermato a G1); conferma 2; scrittura firmata a G3.

# Decisioni del PI

- 2026-11-20 (SC2): polmone 10x dopo Xenium, cancello sulla registrazione, riserva il colon; G0: sonda 5.
- 2026-11-21 (S1): X2 con un metro più preciso, una sessione.
- 2026-11-23 (S3, S3b, S3c): h senza κ del colon, aree omologhe H1-H2; metà nascosta solo in H2; sonda prorogata a 7 sessioni.
- 2026-11-25 (S6): δ = 1 h (RIESAME_C §1.4); sonda a 9.
- 2026-11-25 (S8): A3 tolta, fra i limiti; A5 con ENACT più bin2cell, margine 0,2.
