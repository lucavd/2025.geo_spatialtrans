
Misura: sulle sei aree di 1 mm² di P-1 di S2.7 (prima area di ogni tessuto: intestino A1_r1, tumore del colon A2_r1, stroma A3_r1, linfonodo A4_f1, cervello A5_r1, cuore A6_r1), 3 semi, quattro configurazioni: motore di S2.5 che calcola da sé griglia e verità geometrica (E0); funzione a pezzi con pezzi di 504 µm (E504), 1 008 µm (E1008), un pezzo solo (EInf). Un processo per esecuzione, un thread, al più 6 esecuzioni insieme, nessun altro utente sulla macchina.

- t(a, E) = media dei 3 tempi di orologio; picco(a, E) = massimo sui semi del picco di memoria residente (processo principale e processi dei pezzi).
- migliore(a) = la configurazione a pezzi più veloce fra E504, E1008, EInf con picco ≤ 8 GB; rapporto(a) = t(a, migliore) / t(a, E0).
- **Se rapporto ≤ 1.5 in almeno 5 tessuti su 6**: `simulate_visium_hd_pieces()` diventa l'unico ingresso documentato (la guida di S2.8b dice quale lato usare per aree piccole e per sezioni intere). **Altrimenti**: due ingressi documentati, `simulate_visium_hd()` per aree singole fino a circa 1 mm² e `simulate_visium_hd_pieces()` per aree più grandi e aree di cattura intere, con la nota che a parità di seme i due motori danno realizzazioni diverse della stessa legge (S2.7, P-1).
- Se migliore(a) non è E504 in almeno 5 tessuti su 6, la scelta automatica del lato dei pezzi diventa una voce di BACKLOG per l'ottimizzazione (S1.7), con il guadagno misurato.

Perché 1.5 (scelta tecnica dell'agente, da registrare come decisione con il risultato): un confronto di metodi ripete la simulazione da decine a centinaia di volte per tessuto. Oltre una volta e mezza il costo del motore di S2.5, un secondo ingresso per le aree piccole ripaga la sua complessità; sotto, conta di più avere un solo motore, il cui risultato non dipende da come si divide l'area (S2.7, G-2). Le sezioni intere non entrano nella regola: le fa solo la funzione a pezzi (S2.6: il motore di S2.5 non regge in memoria il tumore e il cervello, gli altri tre con 26-107 GB).

