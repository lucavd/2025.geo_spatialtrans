---
type: Contract
title: "Contratto dell'articolo geo_spatialtrans"
description: "Contratto in bozza, non firmato."
status: draft
generated: { by: geo-agent/skill-v2, at: 2026-11-20T19:00:00+01:00 }
---
# Messaggio

Un simulatore di Visium HD con verità per cellula produce, nel tumore del colon, risposte dei metodi d'analisi vicine al vero quanto un secondo paziente.

Lettore: chi sviluppa o confronta metodi per Visium HD. Rivista: Genome Biology (A), Bioinformatics (B).

Simulatore = catena banale più i pezzi che la battono.

# Misura principale

D (h), più basso è meglio; metro in `docs/riesame/FATTI_B.md`: 12 voci, tumore e stroma × 6 misure dopo bin2cell + RCTD. Unità: l'area (5); bootstrap per aree, B = 2000. Pavimento 7,77 (banale, fisso). Tetto 2,23 (1,85; 2,71). Soglia 3,23 = tetto + δ, δ = 1,0 h (Luca, `docs/riesame/REGOLA_B.md`); sensibilità con δ = 0,5 e 2.

# Affermazioni

| id | testo | figura | serve |
|---|---|---|---|
| A1* | D ≤ soglia, tarato su P2 e misurato con P1 e P5 come riferimento | Fig. 1 | ≤ 3,23 |
| A2 | Nessuna voce oltre il tetto più largo per voce | Fig. 2 | 0 voci oltre 3,6 h |
| A3 | Batte lo Xenium in quadratini (ENACT), geni del pannello | Fig. 3 | scarto ≥ δ |
| A4 | Regge con una catena mai usata in taratura | Fig. 4 | ≤ suo tetto + δ |

# Cancello G1, regola scritta ora

Nella sonda ogni componente comune si sostituisce con la versione vera (profili: cellule vere di P2; geometria: nuclei veri) o si spegne (cattura, spostamento). Con la combinazione migliore:

- D ≤ 3,23: costruzione.
- 3,23 < D ≤ 6,77: riposizionare sull'esito B.
- D > 6,77: chiudere.

# Esiti

- A: A1-A4 reggono.
- B: D ≤ 6,77, A1 non regge. Messaggio: scarto dal vero misurato e dichiarato; non validato per classificare metodi.
- C: G1 dice chiudere, oppure a costruzione finita D > 6,77.
- A3 indecidibile (pannello senza spazio): si scrive, non sposta l'esito.

# Limiti accettati

- Solo tumore del colon e stroma: gli altri tessuti hanno una sola sezione.
- Nessuna classifica contro una verità vera (esiste solo nel polmone).
- Una sezione di taratura; h dipende dal disegno a 5 aree.

# Dati tenuti da parte

P1 e P5: una lettura residua (DEC-093), per la conferma di A1. Rene umano e fegato di ratto mai letti.

# Budget proposto

| fase | sessioni | motivo |
|---|---|---|
| scommessa | 2 | questa e la firma |
| sonda | 3 | check; componenti veri; ENACT e seconda catena |
| costruzione | 10 | 5 voci × 2 colpi (4 voci portano D solo a 3,28); confermato a G1 |
| conferma | 2 | |
| scrittura | firmata a G3 | |

Calcolo: 2-3 ore-core per braccio. PI: un'ora per cancello.

# Decisioni del PI

(nessuna ancora)
