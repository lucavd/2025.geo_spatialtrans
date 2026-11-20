---
type: Review
title: "Riesame di geo_spatialtrans, Parte A: criterio, candidati di articolo, verdetti per componente"
description: "Riesame su carta del progetto dopo 43 sessioni: il tetto misurato dai dati di F3.3, le affermazioni possibili, tre candidati di articolo, un verdetto provvisorio per componente e le opzioni per il PI."
status: draft
generated: { by: geo-agent/skill-v2, at: 2026-11-20T13:30:00+01:00 }
pi_approval: { by: human:luca-vedovelli, at: 2026-11-20, evidence: "in chat: «1. Sì: Parte B come scritta, con il braccio X1. 2. Sì: margine di equivalenza 1,0 h e onere della prova sui nostri componenti.»" }
sources:
  - id: consegna
    resource: "docs/RIESAME_geo_spatialtrans.md (consegna del 2026-11-20, fuori dal repository)"
    title: "Riesame di geo_spatialtrans: consegna per la sessione"
  - id: fatti
    resource: "/riesame/FATTI.md"
    title: "Riesame, Parte A: i fatti riletti dai file"
  - id: esistente
    resource: "/ESISTENTE.md"
    title: "Che cosa esiste già per simulare Visium HD a livello di cellula"
  - id: regola
    resource: "/riesame/REGOLA_B.md"
    title: "Riesame, Parte B: regola di esito scritta prima del calcolo"
  - id: f33
    resource: "/steps/F3.3.md"
    title: "Il simulatore di oggi regge fuori dalle aree su cui è stato tarato?"
---
# Che cosa è

Riesame, non step: decide se il progetto continua, si restringe, cambia strada o chiude. Nessuna modifica al simulatore, nessuna voce di backlog. Le due decisioni aperte di F3.3 e la prova della semina a compartimenti restano sospese fino al verdetto (conferma di Luca in chat, 2026-11-20). Numeri: [FATTI](/riesame/FATTI.md)[^fatti]. Esistente: [ESISTENTE](/ESISTENTE.md)[^esistente].

# 1. Il criterio

## Il tetto

Oggi «come il vero» vuol dire «entro h dalla sezione di taratura». h è la semi-ampiezza dell'intervallo della media di 5 aree della stessa sezione: misura la precisione di una media, non la variabilità biologica.

Il tetto misurato (stimatore quadratico corretto, unità h, dalle tabelle di F3.3)[^f33]:

| tessuto | vero T da V (tetto) | vero T da vero T' | simulato da V |
|---|---|---|---|
| tumore del colon | 2,25 (1,63; 2,95) | 2,76 (2,26; 3,32) | 10,27 (9,70; 10,83) |
| stroma | 1,98 (1,67; 2,31) | 2,06 (1,69; 2,49) | 4,79 (4,62; 4,98) |
| cervello di topo | 3,53 (2,93; 4,20) | non misurabile | 1,32 (0,98; 1,72) |
| tumore e stroma | 2,12 (1,78; 2,55) | 2,44 (2,11; 2,80) | 8,01 (7,64; 8,38) |

Lettura:

- Due sezioni vere dello stesso tessuto distano circa 2-3,5 h. Un giudizio «smentito» con uno scarto sotto questo livello non è un difetto rispetto a un secondo paziente: il criterio di F3.1 è troppo severo di circa 2 h.
- Il simulatore però è a 8 h dalla sua sezione di taratura nel tumore e nello stroma: quasi 4 volte il tetto. Allentare il criterio al tetto non lo salva.
- Per misura: 7 su 18 combinazioni restano oltre il tetto con il loro intervallo, tutte misure di mescolanza (marcatori estranei, singoletti, concordanza con il nucleo, vicini dello stesso tipo) nel tumore e nello stroma. Dentro il tetto: trascritti per cellula nei tre tessuti, due misure del cervello. Indeterminate: 6.
- Nel cervello il simulato è vicino alla sezione di taratura (1,3 h) e lontano dalla seconda sezione (4,3 h), che dista 3,5 h dalla prima: lì il tetto è grande e il simulato è tarato sulla prima sezione.

Regola nuova del criterio (proposta, vale per la Parte B): un braccio si giudica contro il tetto D(T), mai contro h da solo.

## Tessuti con una sola sezione

Intestino di topo, linfonodo, cuore: il tetto non si misura. Affermazioni possibili:

- «difetto più grande di ogni tetto osservato altrove» (oltre 3,5 h): intestino, UMI assegnati (8,0 h), marcatori estranei (5,3 h), vicini dello stesso tipo (−9,0 h). Linfonodo e cuore: nessuna.
- controlli con una ragione strutturale: l'anisotropia del cuore, che il simulatore isotropo non può riprodurre (controprova di F3.3, 0,11 contro 0,10 dentro campione).
- non possibili: «come il vero», «dentro il tetto», difetti piccoli. I difetti del linfonodo (marcatori estranei 3,4 h, vicini −2,4 h) stanno dentro la fascia dei tetti degli altri tessuti.
- per un articolo: questi tre tessuti entrano come esempi qualitativi o si tolgono; servirebbe una seconda sezione pubblica (non cercata in questo riesame).

## Candidati di articolo

**A. Il modello della cattura su disposizioni vere** (strada «cambiare strada»).

- Messaggio: un modello della cattura di Visium HD applicato a cellule Xenium vere produce quadratini da 2 µm vicini al vero quanto un secondo paziente, con la verità cellula per cellula.
- Lettore: chi sviluppa o sceglie metodi di assegnazione dei quadratini alle cellule e di deconvoluzione per Visium HD. Rivista: Bioinformatics o Genome Biology.
- Misura principale: D, scarto da V in unità di h sulle misure dei metodi, tumore e stroma; pavimento D(X1), Xenium in quadratini alla profondità di V; tetto D(T).
- Affermazioni e figure: (1) lo Xenium in quadratini, come in ENACT, resta oltre il tetto (figura dei bracci con la fascia del tetto); (2) con il modello della cattura si arriva al tetto (stessa figura, per misura); (3) la classifica dei metodi sui dati costruiti coincide con quella sul dato vero con la verità Xenium (figura di concordanza dei ranghi); (4) l'errore di bin2cell cresce con lo spostamento dei trascritti, e i valori misurati nei sei tessuti cadono dove previsto (curva di sensibilità); (5) le misure di riferimento di Visium HD in sei tessuti: densità, aree nucleari, spostamento (tabella).
- Stato: nessuna affermazione misurata. Esistente: lo strumento completo non trovato; ENACT ha già fatto la strada corta senza cattura.

**B. Il generatore controllabile con lo scarto dichiarato** (strada «continuare» o «restringere»).

- Messaggio: un simulatore parametrico di Visium HD con strato cellulare esplicito permette prove controllate dei metodi, con lo scarto dal vero misurato e dichiarato.
- Lettore: chi sviluppa metodi e vuole manopole (densità, taglia dei nuclei, spostamento). Rivista: Bioinformatics.
- Misura principale: D(S) fra pavimento D(B) e tetto D(T); risposta dei metodi ai parametri.
- Affermazioni e figure: (1) il simulatore batte la soluzione banale (figura dei bracci); (2) lo scarto dal tetto per misura è dichiarato (figura per misura); (3) l'errore dei metodi risponde a densità, taglia dei nuclei e spostamento come previsto (curve); (4) il cuore come controllo negativo dell'isotropia (figura).
- Stato: l'affermazione (2) oggi vale 8,0 h contro un tetto di 2,1 nel tumore e nello stroma. Obiezione nota (review arXiv 2606.09675): un simulatore favorisce i metodi che ne condividono le ipotesi.

**C. Non c'è articolo.** Si chiude con `docs/CHIUSURA.md`; si riusano le misure sui Visium HD veri, il banco di prova e i dati di verifica sigillati.

# 2. Esistente

Registro: [ESISTENTE](/ESISTENTE.md). Risposte brevi:

1. L'articolo completo: non trovato. ENACT (Bioinformatics 2025) ha già Xenium del colon diviso in quadratini da 2 µm con verità per cellula, senza spostamento né cattura, con i confini veri dati al metodo.
2. Lo strumento con griglia 2 µm, verità per cellula, spostamento e cattura: non trovato. Il più vicino è ALBIS (2026): molecole con cellula d'origine, cellule sferiche, nessuna cattura.
3. Metro e verità pronti: metri sì (precisione e richiamo sui conteggi di ENACT; abbinamento ai nuclei Xenium di STCS). Verità sulla stessa sezione: Visium HD dopo Xenium, solo polmone 10x. Dataset sintetici Visium HD depositati: non trovati.
4. Come altri hanno superato l'ostacolo: molecole vere divise in quadratini; Visium HD dopo Xenium sulla stessa sezione; sezioni adiacenti con nuclei annotati a mano; simulazioni parametriche semplici.

# 3. Dati per la Parte B

Xenium degli stessi blocchi per P1, P2 (taratura) e P5, sezioni seriali a 5-10 µm, pannello di 422 geni, 408 nei 18 024 comuni, contorni dei nuclei e delle cellule e trascritti con coordinate. Sul pannello: tumore, epiteliale 11 e fibroblasti 13 marcatori su 50; stroma, endoteliale 7, epiteliale 12, fibroblasti 8, immunitario 18. RCTD con i soli geni del pannello (un'area per tessuto, conteggi Visium HD ridotti): tipo principale uguale nel 98 % e nel 96 % dei blocchi, ma i singoletti salgono a 99,9 % e 99,98 %. I geni del pannello portano il 3,5-4,1 % degli UMI. Quindi RCTD regge per il tipo principale, ma sul pannello quasi non vede i miscugli: la misura dei singoletti perde il suo spazio, e la regola di esito ne tiene conto.

# 4. Verdetti provvisori per componente

| componente | sessioni | che cosa dicono F3.1 e F3.3 | verdetto | che cosa lo ribalta |
|---|---|---|---|---|
| Generatore dello strato cellulare (regioni, semina, tessellazione, nuclei) | 14 | I trascritti per cellula tornano in 6 tessuti su 6 e restano dentro il tetto fuori campione; i difetti grandi non sono geometrici | sostituire per i dati realistici; tenere per prove controllate solo se batte il banale | Parte B: D(S) uguale a D(XC), oppure D(S) migliore di D(B) di più del margine di equivalenza |
| Disposizione dei tipi cellulari | 3 | Causa principale in F3.1; vicini dello stesso tipo 11,0 h nel tumore (tetto 0,5) e 5,6 h nello stroma (tetto 2,9); lo stroma di P1 e P5 è più mescolato di P2 | sostituire con disposizioni vere | Parte B: lo Xenium in quadratini non migliora le misure di mescolanza rispetto al simulato |
| Profili di espressione | 1 | Profili dei linfociti B giusti; nel cervello l'eccesso di estranei viene dai profili | tenere | Parte B: il braccio banale, senza spostamento, ha lo stesso eccesso di marcatori estranei del simulato |
| Spostamento dei trascritti | 1 | Seconda causa in F3.1: troppi trascritti dei vicini nei nuclei; marcatori estranei oltre il tetto nel tumore (8,1 h contro 1,2) e nello stroma (1,9 contro 0,9) | correggere, per tipo cellulare | Parte B: il braccio banale ha gli stessi marcatori estranei del simulato |
| Cattura e griglia (strisce, fondo, profondità) | 4 | Trascritti per cellula costruiti dalla taratura; profondità fino a 7 volte diversa fra sezioni dello stesso tumore | tenere; profondità come parametro per sezione | Parte B: D(XC) non migliore di D(X1) |
| Riferimenti misurati sui Visium HD veri | 10 | Base di tutte le misure; lo spostamento laterale a 2 µm misurato non trovato in letteratura | tenere | Un lavoro pubblicato con le stesse misure |
| Banco di prova dei metodi e dati di verifica | 3 | Ha trovato i difetti e il sovradattamento alla sezione; il margine h da solo è troppo severo | tenere, con il tetto come criterio | Il controllo di rumore della Parte B fallisce (due pazienti veri distinti fra loro oltre il margine) |
| Infrastruttura, scala, consolidamento | 7 | non giudicata | resta com'è | non pertinente |

Scommessa della consegna (tenere la cattura, togliere il generatore): coerente con questi verdetti, non ancora misurata.

# 5. Regola di esito della Parte B

[Riesame, Parte B: regola di esito scritta prima del calcolo](/riesame/REGOLA_B.md)[^regola]. Committata prima di ogni calcolo della Parte B.

# 6. Opzioni per il PI

| opzione | costo | guadagno atteso | che cosa si perde |
|---|---|---|---|
| 1. Parte B come scritta (consigliata) | 1-2 sessioni; circa 4 ore di calcolo con 4 processi; 0 ore del PI | decide con i numeri se il modello della cattura vale un articolo e se il generatore resta | niente; due sessioni senza risposta chiudono |
| 2. Chiudere ora | 0 | si smette di spendere; si riusano riferimenti e banco di prova | il candidato A, l'unico con uno strumento che non esiste già |
| 3. Parte B ridotta: solo generatore, geni completi, senza Xenium | meno di 1 sessione | decide solo sul generatore | la domanda sulla cattura, cioè la scommessa |

# In breve

## In breve: il simulatore dà qualcosa che non si ottiene dividendo in quadratini un dato Xenium vero?

Non c'è ancora un contratto, quindi nessuna affermazione dell'articolo è misurata. Oggi abbiamo misurato il tetto, cioè quanto distano fra loro due pazienti veri, e cercato chi ha già fatto lo stesso lavoro.

**A cosa serviva.** Dopo 43 sessioni il progetto non ha un articolo definito. Prima di un altro miglioramento dovevamo sapere due cose. Il criterio con cui giudichiamo il simulatore è giusto? E il simulatore dà qualcosa che una strada più corta non dà?

**Cosa abbiamo trovato.**

1. **Il criterio era troppo severo, ma non è questo che fa fallire il simulatore.**
   - Il margine usato finora misura quanto è precisa la media di 5 aree di una sezione, non quanto variano i tessuti veri. Due pazienti veri dello stesso tumore distano fra loro circa due margini.
   - Il simulatore dista dal suo paziente di taratura circa otto margini nel tumore del colon e nel suo stroma: quasi quattro volte quanto distano due pazienti veri.
   - Restano oltre il tetto le misure di mescolanza: trascritti di altri tipi nelle cellule ricostruite, miscugli di tipi nei quadrati di 8 µm, vicini dello stesso tipo. Il numero di trascritti per cellula invece è dentro il tetto in tutti i tessuti dove si misura.
2. **Per l'intestino di topo, il linfonodo e il cuore abbiamo una sola sezione: lì si possono dire solo le cose grosse.**
   - Nell'intestino tre difetti sono più grandi di ogni distanza fra pazienti vista altrove; nel linfonodo e nel cuore nessuno lo è.
   - Il cuore resta utile come controllo, perché le cellule vere sono allungate e quelle simulate no.
3. **Nessuno ha già pubblicato il simulatore completo, ma la strada corta esiste già.**
   - Un lavoro del 2025 (ENACT) ha diviso in quadratini da 2 µm le molecole Xenium dello stesso tumore del colon, con la verità cellula per cellula. Non ha però né lo spostamento dei trascritti né un modello di come il vetrino cattura i trascritti.
   - Uno strumento con quadratini da 2 µm, verità per cellula, spostamento e cattura insieme non l'abbiamo trovato. Il più vicino (ALBIS, 2026) ha cellule sferiche e nessun modello della cattura.
   - Nella tua libreria Zotero non ci sono lavori su Visium HD.
4. **I dati per la prova decisiva ci sono, con un limite.**
   - Il paziente di taratura e altri due hanno lo Xenium degli stessi blocchi, con i contorni dei nuclei e la posizione di ogni trascritto.
   - Il pannello Xenium ha però solo 408 geni utili, che portano il 4 % dei trascritti. Su quei geni il programma che assegna i tipi ai quadrati riconosce bene il tipo principale, ma quasi non vede i miscugli di due tipi. La prova lo tiene in conto: se sul pannello non si vedono più i difetti noti del simulatore, il generatore si giudica sui geni completi.
5. **Verdetto provvisorio: si tiene il modello della cattura e si sostituisce il generatore delle cellule con disposizioni vere.** È la scommessa della consegna. I numeri di oggi la sostengono per il generatore. Per la cattura non dicono ancora nulla.

**Cosa devi decidere.**

1. **Facciamo la prova decisiva come scritta? Consiglio di sì.**
   - Confronta sulle stesse misure il paziente di taratura, un altro paziente, il simulatore di oggi, una soluzione banale, lo Xenium in quadratini e lo Xenium con il nostro modello della cattura.
   - Costa una o due sessioni e circa 4 ore di calcolo sul server con 4 processi. Non chiede lavoro manuale.
   - Alternativa: chiudere ora. Si perde l'unico candidato di articolo con uno strumento che non esiste già.
   - Alternativa: una prova ridotta, solo sul generatore e senza Xenium. Costa meno di una sessione, ma lascia aperta proprio la domanda sulla cattura.
2. **Due bracci si dicono uguali se differiscono meno di un margine, e un dubbio non si conta a favore del simulatore? Consiglio di sì.**
   - Un margine è circa metà della distanza fra due pazienti veri. Il rumore di una differenza fra due bracci, misurato oggi, vale circa mezzo margine, quindi sia «uguale» sia «diverso» sono esiti possibili.
   - Ogni nostro pezzo deve dimostrare di servire: se il confronto resta incerto, quel pezzo conta come non utile.
   - Aggiungo un braccio che costa poco: lo Xenium in quadratini portato alla stessa profondità del Visium HD. Il modello della cattura deve battere questo, non lo Xenium grezzo, che perderebbe già solo per la profondità.
   - Alternativa: una soglia di mezzo margine. È più severa, ma con 5 aree per sezione molti confronti uscirebbero incerti per costruzione.

**Comunque**

- Le misure sui Visium HD veri (densità, aree dei nuclei, spostamento dei trascritti), il banco di prova dei metodi e i dati di verifica sigillati si riusano in ogni caso.
- Il rene e il fegato restano sigillati. La prova decisiva legge per la seconda volta i pazienti 1 e 5: dopo, a ognuno resta una lettura.
- La semina a compartimenti e le due decisioni rimaste aperte nella misura fuori campione restano sospese fino al verdetto.

Appena decidi questi 2 punti, il riesame su carta è chiuso.

**Sui tempi di calcolo.** La stima è mia e viene dai tempi misurati quando abbiamo misurato il simulatore fuori campione: 48 esecuzioni dei due metodi in circa 2 ore con 4 processi. La prova decisiva ne chiede circa 100, su meno geni, più lo scaricamento dello Xenium del paziente di taratura (16 GB, circa 12 minuti).

[^fatti]: [Riesame, Parte A: i fatti riletti dai file](/riesame/FATTI.md).
[^esistente]: [Che cosa esiste già per simulare Visium HD a livello di cellula](/ESISTENTE.md).
[^f33]: [Il simulatore di oggi regge fuori dalle aree su cui è stato tarato? (F3.3)](/steps/F3.3.md).
[^regola]: [Riesame, Parte B: regola di esito scritta prima del calcolo](/riesame/REGOLA_B.md).

# Decisioni di Luca (2026-11-20)

1. «Sì: Parte B come scritta, con il braccio X1.»
2. «Sì: margine di equivalenza 1,0 h e onere della prova sui nostri componenti.» Con un controllo positivo in più per la cattura (G-6 nella [regola di esito](/riesame/REGOLA_B.md)).
3. Protocollo: la skill `progetto-master`, copiata in `.claude/skills/progetto-master/SKILL.md` (commit `3f00f5e`).
