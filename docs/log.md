# Registro delle modifiche

## 2026-10-23 (risultati e chiusura di S1.5b)

* **Chiusura**: [Con i tipi cellulari nei preset, il simulatore dà a ogni tipo la taglia di nucleo che ha nei tessuti veri? (S1.5b)](steps/S1.5b.md) chiuso con l'approvazione di Luca («va bene», lettura del riassunto «In breve» senza dubbi); nessuna decisione per Luca. ROADMAP §7: S1.5b chiusa, prossima sessione proposta Step 1.5 (semina eterogenea e dipendente dal tipo).
* **Aggiornamento**: [Con i tipi cellulari nei preset, il simulatore dà a ogni tipo la taglia di nucleo che ha nei tessuti veri? (S1.5b)](steps/S1.5b.md): cancelli, regola del braccio, verdetti (banco v2 dopo la correzione della semina), deviazioni dichiarate, post hoc, revisione avversariale (due revisori, due giri), riassunto «In breve» senza decisioni per Luca, rinvii BL-153…BL-156. Decisione tecnica [Tipi nei preset: braccio ricentrato, semina cieca al tipo, tipi predefiniti (DEC-019)](decisions/DEC-019.md). Report `reports/S1.5b.html`. In attesa della lettura di Luca.

## 2026-10-22 (apertura di S1.5b)

* **Aggiornamento**: [S1.5b](steps/S1.5b.md): ok di Luca su previsioni e regola dei bracci («Approvo così»); decisioni [Cervello: i nuclei persi da Cellpose ricevono un tipo dall'espressione (DEC-017)](decisions/DEC-017.md) e [Stroma del tumore: preset misto come misurato (DEC-018)](decisions/DEC-018.md).
* **Creazione**: file dello step [Con i tipi cellulari nei preset, il simulatore dà a ogni tipo la taglia di nucleo che ha nei tessuti veri? (S1.5b)](steps/S1.5b.md) con dubbi ricevuti (BL-099, BL-138, BL-055, BL-149, BL-150), definizioni operative, cancelli G-1…G-5, previsioni P-1…P-8, controprova CP-1, due decisioni per Luca (composizione del cervello, stroma misto), scritti prima di ogni simulazione; riferimenti per tipo, valori reali per ROI e stime a mano degli attesi in `results/S1.5b/`. Previsioni in attesa dell'ok di Luca.
* **Creazione**: decisione tecnica [Tipi nei preset: gruppi di taglia, densità del tipo puro, distanza minima unica (DEC-016)](decisions/DEC-016.md).
* **Nota**: i documenti di R6 datano la chiusura al 2026-10-23, ma i commit di chiusura sono del 2026-10-22 (orologio del server); le date di R6 non si correggono qui.

## 2026-10-23 (chiusura di R6)

* **Chiusura**: [Di che tipo sono le cellule vere nei sei tessuti, e i tipi spiegano le taglie dei nuclei? (R6)](steps/R6.md) chiuso con l'approvazione di Luca («ok per tutte e 3 le questioni aperte, chiudi»); decisioni [I riferimenti per tipo cellulare vengono dalle etichette del criterio corretto (DEC-013)](decisions/DEC-013.md), [Nel colon l'epitelio e' tumorale o normale per l'istologia della ROI (DEC-014)](decisions/DEC-014.md), [Nel linfonodo i riferimenti per B e T restano provvisori e i linfociti hanno una sola taglia (DEC-015)](decisions/DEC-015.md). BIO_REFERENCES B-161…B-185 adottati; ROADMAP §7: R6 chiusa, prossima sessione proposta S1.5b (tipi nei preset).

## 2026-10-23 (risultati di R6)

* **Aggiornamento**: [Di che tipo sono le cellule vere nei sei tessuti, e i tipi spiegano le taglie dei nuclei? (R6)](steps/R6.md): risultati di cancelli, previsioni e controprove, correzione post hoc del criterio di etichettatura (PH-1), revisione avversariale con tre correzioni ed esecuzione ripetuta, riassunto «In breve» con tre decisioni per Luca, dubbi rinviati (BL-147…BL-152). Decisione tecnica [Come si attribuisce l'espressione a 2 µm ai nuclei veri e come si classificano (DEC-012)](decisions/DEC-012.md). Report `reports/R6.html`.

## 2026-10-23 (apertura di R6)

* **Creazione**: file dello step [Di che tipo sono le cellule vere nei sei tessuti, e i tipi spiegano le taglie dei nuclei? (R6)](steps/R6.md) con dubbi ricevuti, perimetro, cancelli G-1…G-5, previsioni P-1…P-6, controprove; ok di Luca («si va bene 1 e 2») su previsioni, catalogo dei tipi e soglie.

## 2026-10-21 (chiusura di S1.5)

* **Chiusura**: [Il simulatore del livello cellulare funziona tutto insieme, con i sei preset dei tessuti, e sa riconoscere la propria scala? (S1.5)](steps/S1.5.md) chiuso con l'approvazione di Luca; decisione [Il caso limite «scala sbagliata di 10 volte, scatta C8» fallito in S1.5 e' un errore della pre-registrazione: la scala viene dall'utente o dai metadati (DEC-011)](decisions/DEC-011.md). Design v1.6 aggiornato di conseguenza; nella skill di sessione la regola sulle scelte decise dai dati.

## 2026-10-21 (risultati di S1.5)

* **Aggiornamento**: [Il simulatore del livello cellulare funziona tutto insieme, con i sei preset dei tessuti, e sa riconoscere la propria scala? (S1.5)](steps/S1.5.md): risultati dei cancelli, delle previsioni e della controprova; analisi post hoc PH-1...PH-3; revisione avversariale (due revisori) con le correzioni del codice (commit `ab03c27`) e la controprova delle correzioni; note D-11...D-14; riassunto «In breve» con una decisione per Luca (il caso limite «scala sbagliata di 10 volte → C8» non e' realizzabile); dubbi rinviati in BACKLOG (BL-143...BL-146). Design v1.6, BIO_REFERENCES B-155...B-160, report `reports/S1.5.html`.

## 2026-10-21 (apertura di S1.5)

* **Aggiornamento**: [Il simulatore del livello cellulare funziona tutto insieme, con i sei preset dei tessuti, e sa riconoscere la propria scala? (S1.5)](steps/S1.5.md): ok di Luca sulle previsioni («Approvo così»); le tre scelte del modello, chieste per errore a Luca, riscritte come decise dai dati su sua correzione: distanza minima ([Nei sei preset la distanza minima fra i centri delle cellule e' quella stimata dai tessuti veri (DEC-009)](decisions/DEC-009.md)), conversione del modello di taglia con una regola sul pannello, nuclei in piu' parti ([I nuclei simulati in piu' pezzi restano come sono: sono piu' rari che nelle maschere reali (DEC-010)](decisions/DEC-010.md), `results/S1.5/S1.5_real_multipart.csv`).
* **Creazione**: file dello step [Il simulatore del livello cellulare funziona tutto insieme, con i sei preset dei tessuti, e sa riconoscere la propria scala? (S1.5)](steps/S1.5.md) con dubbi ricevuti e loro ripartizione, definizioni operative, cancelli G-1…G-5, previsioni P-1…P-8, controprova CP-1 e tre decisioni chieste a Luca, scritti prima di ogni calcolo; riferimenti reali per ROI in `results/S1.5/S1.5_real_refs_roi.csv`. Previsioni in attesa dell'ok di Luca.
* **Creazione**: decisione [S1.5 costruisce orchestratore, controlli, preset e stima della scala; la composizione per tipi cellulari passa a due step successivi (DEC-008)](decisions/DEC-008.md), scelta di Luca in chat.

## 2026-10-21 (chiusura di R3b)

* **Chiusura**: [Nel linfonodo e nel cuore, da quale insieme di nuclei veri prendiamo posizioni, aree e vicini? (R3b)](steps/R3b.md) chiuso con l'approvazione di Luca; decisioni [Nel linfonodo e nel cuore una sola maschera nucleare: i nuclei di Space Ranger con la forma di StarDist, con la taglia stimata sui soli nuclei di forma StarDist (DEC-006)](decisions/DEC-006.md) e [Il controllo di integrita' di R3b fallito sul centro dei nuclei nel tessuto valido e' un errore di scrittura del controllo, non della maschera (DEC-007)](decisions/DEC-007.md). BIO_REFERENCES B-133…B-154; BACKLOG BL-100, BL-115, BL-121, BL-130 chiuse; ROADMAP §7: R3b chiusa, prossima sessione proposta S1.5.

## 2026-10-20 (apertura di R3b)

* **Aggiornamento**: [Nel linfonodo e nel cuore, da quale insieme di nuclei veri prendiamo posizioni, aree e vicini? (R3b)](steps/R3b.md): risultati, analisi post hoc, revisione avversariale (due revisori, numeri concordanti), note D-9...D-11, riassunto «In breve» con due decisioni per Luca; report `reports/R3b.html`. BACKLOG BL-141, BL-142; BL-064 rinviata a S1.5.
* **Creazione**: file dello step [Nel linfonodo e nel cuore, da quale insieme di nuclei veri prendiamo posizioni, aree e vicini? (R3b)](steps/R3b.md) con dubbi ricevuti, maschere a confronto, cancelli G-1…G-5, previsioni P-1…P-4, controprova CP-1 e regola di decisione, scritti prima di ogni calcolo; previsioni approvate da Luca in chat il 2026-10-20 («Approvo così»). Forma delle cellule del cuore (BL-064) rinviata a S1.5 da Luca.

## 2026-10-20 (apertura di S1.4d, primo step con il protocollo v2)

* **Chiusura**: [La forma del nucleo ritagliato nel territorio e riportato all'area giusta entra nel pacchetto e riproduce i contatti contati a mano? (S1.4d)](steps/S1.4d.md) chiuso con l'approvazione di Luca; decisioni [Il controllo di precisione sulle aree fallito in S1.4d è un errore della soglia, non del programma (DEC-004)](decisions/DEC-004.md) e [Nel tumore del colon nessun limite all'area del nucleo; la coda del modello di taglia si corregge con i tipi cellulari (DEC-005)](decisions/DEC-005.md). ROADMAP §7: S1.4d chiusa, prossima sessione proposta R3b.
* **Aggiornamento**: [La forma del nucleo ritagliato nel territorio e riportato all'area giusta entra nel pacchetto e riproduce i contatti contati a mano? (S1.4d)](steps/S1.4d.md): risultati, revisione avversariale (due revisori, numeri concordanti), correzioni D-8…D-14, riassunto «In breve» con due decisioni per Luca. Decisione tecnica [Come il nucleo ritagliato e riportato all'area giusta è implementato nel pacchetto (DEC-003)](decisions/DEC-003.md).
* **Creazione**: file dello step [La forma del nucleo ritagliato nel territorio e riportato all'area giusta entra nel pacchetto e riproduce i contatti contati a mano? (S1.4d)](steps/S1.4d.md) con dubbi ricevuti, cancelli G-1…G-5, previsioni P-1…P-5 e controprova, scritti prima di ogni calcolo; previsioni approvate da Luca in chat il 2026-10-20 («Approvo così»).
* **Creazione**: decisione [Protocollo di sessione versione 2 (DEC-002)](decisions/DEC-002.md), testo preparato con Luca in chat il 2026-10-20; non ancora riletto da Luca nella forma registrata.
* **Creazione**: decisione tecnica [Strumento di controllo della documentazione in formato OKF (DEC-001)](decisions/DEC-001.md) e lo strumento `tools/okf.py` con i suoi test (`tools/test_okf.py`, 8/8).
* **Creazione**: bundle OKF (`docs/decisions/`, `docs/steps/`, `docs/index.md`, questo registro).
