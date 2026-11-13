
Codice: `tools/F3.1_sim.R`, `tools/F3.1_export_sim.R`, `tools/F3.1_reg_sim.py`, `tools/F3.1_export_real.py`, `tools/F3.1_reference.py`, `tools/F3.1_b2c.py`, `tools/F3.1_rctd.R`. Ambiente di bin2cell: `.venv-b2c` (bin2cell 0.3.5, scanpy 1.12.5). Stesse funzioni sui due lati, a partire dallo stesso formato (matrice per quadratino, quadratini con nucleo e in tessuto, nuclei con tipo).

Misure con previsione (per area; valore del tessuto = media delle 5 aree, sul simulato prima la media dei 2 semi):

| id | misura | definizione |
|---|---|---|
| M-1 | UMI assegnati | quota degli UMI dei quadratini in tessuto che finiscono in una cellula dopo l'allargamento |
| M-2 | UMI per cellula | media degli UMI per cellula dopo l'allargamento (cellule con almeno metà del nucleo in tessuto) |
| M-3 | marcatori estranei | per cellula con tipo e almeno 5 UMI sui marcatori: quota degli UMI sui 50 marcatori degli altri tipi (marcatori dal riferimento); media sulle cellule |
| M-4 | singoletti | quota dei blocchi di 8 µm non scartati che RCTD chiama singoletto |
| M-5 | concordanza con il nucleo | quota dei blocchi (non scartati, con un nucleo tipizzato che copre almeno 4 dei 16 quadratini) in cui il primo tipo di RCTD è il tipo del nucleo |
| M-6 | coerenza spaziale | quota delle coppie di blocchi adiacenti (non scartati) con lo stesso primo tipo |

Misure senza previsione (descrittive o esplicative): numero di cellule per area; blocchi sotto 20 UMI e scartati; M-3 sui profili veri generati (solo simulato); M-1...M-3 con la registrazione corretta dello scarto di R5 (solo vero); verità del simulato: purezza (quota degli UMI assegnati a una cellula che vengono da lei), recupero (quota degli UMI catturati di una cellula finiti nella sua cellula), quota da cellule di altro tipo, accuratezza del primo tipo di RCTD contro il tipo che domina gli UMI del blocco.

