| cosa | perché | cosa serve per chiuderlo | step destinatario | voce di BACKLOG |
|---|---|---|---|---|
| `main` fermo al 2025 | scelta di Luca: non toccarlo adesso | decidere come aggiornarlo | alla prima versione usabile per i metodi (dopo F3.1) | BL-253 |
| cartella di ripresa senza `run.rds` riusata senza controlli | il codice di `R/` è fermo in questo step | rifiuto e test | primo step che tocca il codice dello Step 2 | BL-254 |
| documentazione nel codice della funzione a pezzi (`keep_profiles` con l'uscita leggera, fondo nell'uscita leggera, versione di R nelle impronte) e frase della guida sulle funzioni dei passi (N3) | codice fermo; la frase della guida richiederebbe un nuovo commit provato | roxygen e frase corretti | primo step che tocca il codice dello Step 2 | BL-255 |
| `test_S2.3.R` esce con 0 anche con FAIL | codice dei test fermo | `quit(status = )` | primo step che tocca il codice dello Step 2 | BL-256 |
| confronti di md5 di file RDS fra versioni di R (T7b di `test_S1.5d.R`, avvio rapido) | codice dei test fermo | confronto del contenuto | primo step che tocca il codice (urgente per chi rifà la batteria dello Step 1) | BL-257 |
| strumenti dei documenti (portata, file non tracciati), link rotto nella guida per Daniele dello Step 1, durata delle batterie dello Step 2 nel README (da 2 a 16 minuti nel terzo giro, non fino a 25) | documenti fuori dal bundle dopo il commit provato; conversione dei documenti | correzioni | DOC1 | BL-258 |
