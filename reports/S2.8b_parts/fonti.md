| frase del riassunto | fonte (file, riga o colonna) | valore letto |
|---|---|---|
| il README di prima nominava una sola delle 15 opzioni e delle 20 funzioni | `/mnt/micron/geo_spatialtrans/S2.8b/official3/CP1.log` | «1/35 voci coperte» |
| due esecuzioni della guida danno conteggi identici | `official3/gate_guide.log` | «PASS impronte e controlli uguali nelle due rese» |
| una copia con un esempio sbagliato non si produce | `official3/gate_cp2.log` | «PASS CP-2 … multiplo di 56» |
| una quarantina di imprecisioni | `results/S2.8b/S2.8b_review_discrepancies.csv` | 42 righe, una delle quali (N-04) errore del revisore |
| una frase di dettaglio messa in coda | stesso file, riga N3 | «rinviata» (BL-255) |
| il cervello impiega 335 secondi per mm² con i pezzi predefiniti | `results/S2.7/S2.7_P3_P4.csv`, riga `sec_A5_p504w8`, colonna `s_per_mm2_tissue` | 334.64 |
| 300 del bilancio | documento di progetto dello Step 2, sezione 12 («meno di 5 minuti» per mm²) | 300 s |
| circa un decimo in più | 334.64 / 300 | 1.115 |
| con pezzi più grandi il cervello restava sotto | `results/S2.7/S2.7_P3_P4.csv`, riga `sec_A5_p1008w8` | 263.32 |
| aggiornamento di R il 5 novembre | `/var/log/apt/history.log` del server | r-base-core 4.6.1-9 → 4.6.2-1, Start-Date 2026-11-05 09:45:40 |
| le cellule dell'esempio sono le stesse di fine ottobre | `official3/gate_quickstart.log` | «cell_df identico a S1.6b: TRUE» (file di S1.6b del 2026-10-30) |
| cambia un solo byte | `official3/gate_quickstart.log` | «byte diversi dal file di S1.6b: 1 (posizioni: 10)» |
| anche quelle di una prova del livello cellulare | `/mnt/micron/geo_spatialtrans/S2.8b/probe/T7b.log` | `equal_461=TRUE` in 3 livelli su 3 |
| le prove della griglia passano tutte | `official3/gate_battery.log` | 21, 36, 39, 37, 40, 33, 78, 46 PASS, 0 FAIL |
| quelle delle cellule con i soliti errori noti più il confronto | `official3/battery1.log` | «FAIL trovati [T13 T7b] ammessi [T13]», unica regressione |
| il ramo di sviluppo contiene lo Step 2; il principale non è stato toccato | `official3/G4_merge.log` | vedi G-4 |
| una produzione della guida richiede circa 3 minuti | `/mnt/micron/geo_spatialtrans/S2.8b/probe/render.log` | 2 min 53 s |
| da 40 minuti a un'ora | `official2/start.txt`, `end.txt`; `official3/start.txt`, `end.txt` | 10:52-11:32; 11:36-12:39 |
| tre giri | cartelle `official`, `official2`, `official3` | commit `038f150`, `9102831`, `fea12bd` |
| la stima era un'ora | file dello step, «Stima dei tempi» | «circa un'ora di orologio» |
