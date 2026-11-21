---
type: Rule
title: "S2, cancello della registrazione su X2 con un metro più preciso"
description: "Regola scritta e committata prima di calcolare: candidati che pesano i geni per i conteggi, ammissione su sintetico e X1, scelta per precisione su X2, giudizio di X2."
status: stable
generated: { by: geo-agent/progetto-master, at: 2026-11-21T12:00:00+01:00 }
sources:
  - id: s1
    path: docs/sonda/REGOLA_S1.md
  - id: stato
    path: docs/STATO.md
---
# Perché

In S1 X2 non è giudicabile: errore standard radiale 1,1-3,4 µm a 200-250 µm, nessun riquadro valido a 100 µm (`results/sonda/S1/X2_se_by_L.csv`). Decisione di Luca, 2026-11-21: «proviamo x2 con altro metodo più preciso». Ipotesi: il metro di S1 standardizza ogni gene, quindi in X2 (fino a 4892 geni, quasi tutti con pochi conteggi per riquadro) i geni rari pesano quanto quelli abbondanti e il loro rumore domina. Un metro che pesa i geni per i conteggi deve dare errore standard ≤ 0,5 µm a 200 µm o meno in tutte e 3 le regioni di X2.

# Cosa resta come in S1[^s1]

Catena, variante, regioni, metà tenuta da parte, spostamento grossolano (`L_choice.json` di S1, mappe totali), finestra ± 12 µm, filtro dei geni (≥ 10 conteggi in entrambe le mappe nel riquadro), validità (≥ 500 conteggi), raffinamento a quadrica, regola del bordo, stima della precisione (metà Xenium a caso, Visium binomiale 1/2, seme 20261121), L in {100, 200, 250} µm con soglia 0,5 µm, affine locale a scacchiera, controlli positivi e negativi con le stesse tolleranze, regola di esito.

# Candidati (cambia solo la funzione di somiglianza nel riquadro)

Mappa Xenium per gene su una finestra allargata di 3σ, lisciata con un nucleo gaussiano di deviazione σ, poi ritagliata. Due funzioni, ciascuna con σ = 2 µm e σ = 4 µm (quattro candidati):

- A, correlazione sui conteggi: CC(s) = Σ_g Σ_b X̃_g(b + s) (V_g(b) − V̄_g). Nessuna standardizzazione: ogni gene pesa per i suoi conteggi.
- B, verosimiglianza di Poisson: Λ_g = X̃_g + 0,1 · media di X̃_g sulla finestra; ℓ(s) = Σ_g [Σ_b V_g(b) log Λ_g(b + s) − N_g log Σ_b Λ_g(b + s)], N_g = UMI del gene nel riquadro (scala per gene profilata).

Il metro di S1 gira accanto, come riferimento.

# Ammissione (mai sui residui di X2)

Un candidato è ammesso se:

1. passa il test sintetico di S1 (`tools/sonda/S1_test_cc.py`: spostamento noto e spostamento aggiunto entro 0,3 µm, negativo ≤ 10 %);
2. su X1 con L = 200 µm: controlli positivi 9/9 entro 0,3 µm per asse e negativi ≤ 10 % in ogni regione;
3. su X1 ridà S1: per ogni regione la mediana fra riquadri della differenza (candidato − S1) è entro 0,3 µm per asse, e l'esito di X1 non cambia.

# Scelta (sulla sola precisione)

Per ogni candidato ammesso e ogni L: errore standard radiale (media quadratica fra riquadri) e riquadri usati in ciascuna regione di X2. L* = il più piccolo L con errore standard ≤ 0,5 µm e almeno metà dei riquadri usati in tutte e 3 le regioni (la seconda condizione è nuova: impedisce di arrivare alla precisione scartando riquadri). Vince il candidato con L* più piccolo; parità: il massimo fra regioni dell'errore standard più basso, poi B prima di A, poi σ minore. Sono permessi solo errori standard e conteggi di riquadri: nessuno spostamento di X2 si stampa prima della scelta.

# Giudizio di X2

Solo il candidato scelto, con L*: residui, affine locale, controlli e esito come in S1. Se i controlli di X2 falliscono il metro non vale e X2 non è giudicato.

- Nessun candidato ammesso o nessun L* su X2: X2 resta non giudicabile; è il primo colpo di questa ipotesi.
- X2 passa: si prosegue con X1 e X2.
- X2 non passa o non è giudicabile: si prosegue con X1, e il tetto di A5 (concordanza X1 contro X2) si riscrive prima di misurarlo.

# Diagnosi (nessuna scelta ne dipende)

Sintetico simile a X2 (2000 geni a bassa espressione, circa 70 UMI per quadratino da 8 µm): deviazione standard fra riquadri dei cinque metri.

# Limiti dichiarati ora

- σ e la pseudoconta 0,1 sono fissati qui, senza taratura.
- La scelta fra quattro candidati sulla precisione di X2 è una scelta sui dati di X2, ma non sui residui che il cancello giudica.

[^s1]: docs/sonda/REGOLA_S1.md
