# Riesame di geo_spatialtrans: consegna per la sessione

Scritto il 2026-11-20 fuori dal repository. Viene da una lettura del branch `f3-methods` (commit `40d3b68`) e da una ricerca dell'esistente di circa un'ora. Ogni numero ha la sua fonte: rileggilo dal file prima di usarlo. È allineata alla skill `progetto-master` salvata lo stesso giorno.

## Che cosa è questa sessione

Questa sessione non è uno step. È un riesame: decide se il progetto continua, si restringe, cambia strada o si chiude. Non modifica il simulatore e non apre voci di backlog.

Carica la skill `progetto-master`: questo è il suo riesame, e i suoi sei passi valgono tutti. Dove questa consegna è più specifica, vale la consegna. Tre precisazioni:

- Lo strumento `tools/progetto.py` non è ancora nel repository e questa sessione non lo usa: niente `apri`, `check`, `chiudi`. L'adozione (init, contratto, stato) si fa dopo il verdetto, se il verdetto non è «chiudere».
- Della skill `geo-spatialtrans-session` valgono solo le regole che proteggono i dati: lavoro solo sul server, dati reali in sola lettura, sigillo dei dati di verifica (DEC-093, regola R-2), rene e fegato sigillati. Il suo rituale di step non si applica.
- Si lavora su un branch nuovo, `riesame`, creato da `f3-methods`.

Le due decisioni aperte di F3.3 e la prova della semina a compartimenti restano sospese fino al verdetto. Chiedi al PI di confermarlo in apertura. È l'unica domanda di apertura.

## Perché adesso

- Dopo 43 sessioni il progetto non ha un articolo definito: mancano messaggio, lettore e rivista. L'obiettivo scritto è «dati simulati controllabili e biologicamente solidi».
- La prima prova d'uso ha smentito 18 giudizi su 36. Fuori campione, sugli altri pazienti, il simulato peggiora di circa tre margini.
- Non esiste un pavimento: nessuno ha misurato, sulle stesse misure, quanto vale una soluzione banale o una già pubblicata.
- In nessun branch c'è un documento sullo stato dell'arte. Una ricerca ha trovato una strada più corta verso lo stesso scopo.
- Il protocollo attuale propone un altro step di miglioramento. Prima serve sapere se il simulatore ha un lettore.

## La domanda

Che cosa dà questo simulatore che non si ottiene dividendo in quadratini di 2 µm un dato Xenium vero?

Se la risposta è «niente che si misuri», il progetto si chiude. Se è «qualcosa», quel qualcosa è l'articolo.

## Fatti da rileggere

| fatto | valore | fonte |
|---|---|---|
| Sessioni dal 2026-09-18 | 43 | `docs/ROADMAP.md`, sezione 7 |
| Prima prova d'uso, 6 misure per 6 tessuti | 14 giudizi confermati, 4 indeterminati, 18 smentiti | `docs/steps/F3.1.md` |
| Valori del simulatore stimati sulle stesse 30 aree usate per verificarlo | 380 su 478 | `docs/steps/F3.2.md` |
| Fuori campione, stesse sezioni | nessun sovradattamento visibile; uno scarto sotto un margine non si vedrebbe | `docs/steps/F3.3.md` |
| Fuori campione, altri pazienti e sezioni | il simulato peggiora di circa tre margini | `docs/steps/F3.3.md` |
| Trascritti per cellula nel tumore del colon | circa 1.100 nel paziente di taratura; 640, 2.000 e 270 nelle altre sezioni | `docs/steps/F3.3.md` |
| Valori veri dentro la variabilità delle aree simulate | 18 su 228 | `docs/steps/F3.3.md` |
| Vicini dello stesso tipo nello stroma | 89% nel paziente di taratura, 72-75% negli altri, 65% nel simulato | `docs/steps/F3.3.md` |
| Backlog | 273 voci, 155 aperte | `docs/BACKLOG.md` |
| Sessioni su geometria di nuclei e territori | 15 (S1.3 e seguenti, S1.4 e seguenti, R2, R2b, R2c, R3, R3b): lettura di chi scrive | `docs/ROADMAP.md`, sezione 7 |
| Confronti con simulatori esistenti | uno solo, con SimSpace, del 2025-07-28 | branch `simspace`, file `SIMPSPACE.md` |

## Che cosa esiste

La colonna «letto» dice che cosa ha aperto chi scrive. Rifai la ricerca: non fidarti di questa tabella.

| lavoro | che cosa fa | letto | dove |
|---|---|---|---|
| ENACT (Kamel et al., Bioinformatics 2025; preprint ottobre 2024) | Valuta l'assegnazione dei quadratini alle cellule su due dataset sintetici simili a Visium HD, ricavati da Xenium e da seqFISH+. I notebook stanno in `src/synthetic_data`. Il testo principale non tratta lo spostamento dei trascritti. | testo principale; supplemento no | https://doi.org/10.1093/bioinformatics/btaf094 e https://github.com/Sanofi-Public/enact-pipeline |
| Dati di Oliveira et al. sul tumore del colon | Visium HD, Xenium e cellula singola da sezioni seriali degli stessi blocchi. Sono gli stessi campioni del nostro tumore del colon. | solo i titoli delle serie GEO | GSE280315 (Visium HD), GSE280314 (Xenium), GSE280311 (Flex) |
| SpatialSimBench (Genome Biology 2025) e Open Problems | Confronto di 13 simulatori spaziali su 10 dataset con 35 metriche, con metodi di controllo. | abstract | https://doi.org/10.1186/s13059-025-03505-w e https://openproblems.bio/benchmarks/spatial_simulators |
| SRTsim, scCube, Spider, SimSpace, STADiffuser | Simulatori di conteggi a livello di cellula o di spot. | abstract o estratti | https://doi.org/10.1186/s13059-023-02879-z, https://doi.org/10.1038/s41467-024-49445-0, https://doi.org/10.1093/bioinformatics/btaf562, bioRxiv 2025.07.18.665587, https://www.nature.com/articles/s41467-026-76829-1 |
| Review sulla segmentazione (giugno 2026) | Cita solo scCube e Spider come simulatori. Avverte che un simulatore favorisce i metodi che ne condividono le ipotesi. | sì | https://arxiv.org/abs/2606.09675 |
| STCS (preprint 2026) | Valida su un Visium HD con segmentazione ricavata da uno Xenium accoppiato. | estratto | https://www.biorxiv.org/content/10.64898/2026.02.26.708370v2 |
| Benchmark delle piattaforme subcellulari (preprint 2024) | Misura lo spostamento dei trascritti in Visium HD e Stereo-seq. | estratto | https://www.biorxiv.org/content/10.1101/2024.12.23.630033v1.full |

Non è stato trovato un simulatore parametrico di Visium HD con uno strato cellulare esplicito e un modello della cattura. La ricerca è stata di una decina di interrogazioni su web e PubMed: «non trovato» non vuol dire «non esiste».

## Parte A: riesame su carta

Una sessione. Nessun calcolo nuovo sul simulatore.

1. Attacca il criterio prima del metodo.
   - Oggi «come il vero» vuol dire «come una sezione di un paziente». F3.3 mostra che un altro paziente vero è lontano dalla sezione di taratura. Calcola dai risultati di F3.3 il tetto: lo scarto, in margini, fra una sezione vera e un'altra sezione vera dello stesso tessuto. Un simulatore non si giudica più severamente di un secondo paziente.
   - Intestino, linfonodo e cuore hanno una sola sezione: lì il tetto non si misura. Scrivi quali affermazioni restano possibili per questi tre tessuti.
   - Scrivi due o tre candidati di articolo, ognuno con messaggio in una frase, lettore, misura principale e da 3 a 5 affermazioni con la loro figura. Un candidato è «non c'è articolo».
2. Rifai la ricerca dell'esistente e lascia una riga per voce in `docs/ESISTENTE.md`. Rispondi a quattro domande: l'articolo è già scritto? Esiste lo strumento? Esistono metro e verità già pronti? Come hanno superato altri lo stesso ostacolo? Verifica in particolare:
   - come ENACT costruisce i suoi dataset sintetici (supplemento e notebook);
   - se un simulatore pubblicato produce quadratini di 2 µm con l'assegnazione vera alle cellule e con lo spostamento dei trascritti;
   - la libreria Zotero del PI.
3. Verifica i dati per la Parte B: quali pazienti hanno lo Xenium degli stessi blocchi, quanti geni ha il pannello, quanti sono in comune con i 18.024 geni comuni di F3.3, se ci sono i contorni dei nuclei. Controlla che RCTD e i 50 marcatori per tipo reggano sul pannello ridotto.
4. Dai un verdetto provvisorio per componente, con l'osservazione che lo ribalterebbe.

| componente | sessioni | che cosa dicono F3.1 e F3.3 | verdetto | che cosa lo ribalta |
|---|---|---|---|---|
| Generatore dello strato cellulare (regioni, semina, tessellazione, nuclei) | | | | |
| Disposizione dei tipi cellulari | | | | |
| Profili di espressione | | | | |
| Spostamento dei trascritti | | | | |
| Cattura e griglia (strisce, fondo, profondità) | | | | |
| Riferimenti misurati sui Visium HD veri | | | | |
| Banco di prova dei metodi e dati di verifica | | | | |

5. Scrivi la regola di esito della Parte B, con le soglie, prima di ogni calcolo. Fai il commit.
6. Fermati. Porta al PI al massimo 3 opzioni, ognuna con costo e guadagno atteso. «Chiudere» è una delle tre. Scrivi il messaggio secondo la sezione «Scrivere per il PI» della skill.

## Parte B: la prova che decide

Una sessione, due al massimo. Tessuto: il tumore del colon e il suo stroma, perché hanno lo Xenium degli stessi blocchi e più pazienti.

Misura principale: lo scarto dal vero in unità di margine, sulle misure di F3.3, con la regola V, il margine h e lo stimatore quadratico corretto (DEC-098). Più basso è meglio. Tutti i bracci usano gli stessi geni: il pannello Xenium dentro i 18.024 geni comuni.

| braccio | che cos'è | costo |
|---|---|---|
| V | Visium HD vero del paziente di taratura. È il riferimento. | ricalcolo sui geni del pannello |
| T | Visium HD vero di un altro paziente. Il suo scarto da V è il tetto. | ricalcolo sui geni del pannello |
| S | Il simulatore di oggi. | ricalcolo sui geni del pannello |
| B | Soluzione banale: cellule a caso con la densità del preset, profilo medio per tipo, cattura senza spostamento. È il pavimento. | nuovo, piccolo |
| X | Xenium in quadratini: i trascritti Xenium sommati in quadratini di 2 µm, con i nuclei Xenium dati a bin2cell come in F3.1. | nuovo |
| XC | Xenium più cattura: cellule e trascritti Xenium, poi il modello della cattura del nostro Step 2. | nuovo; si fa solo se costa meno di una sessione |

Regola di esito. D è lo scarto dal vero. Le soglie di «meglio» e «uguale» si fissano sul rumore fra semi e fra aree, prima del calcolo.

| domanda | risultato | conseguenza |
|---|---|---|
| Il modello della cattura serve? | D(X) al livello di D(T) | La strada corta basta. Il simulatore non ha un lettore come sostituto del dato vero. |
| | D(X) sopra D(T), e D(XC) meglio di D(X) | Il modello della cattura vale: è il candidato articolo. |
| | D(X) sopra D(T), e D(XC) uguale a D(X) | Il modello della cattura non aggiunge. |
| Il generatore dello strato cellulare serve? | D(S) uguale a D(XC) | Il generatore regge il confronto con la disposizione vera: si tiene. |
| | D(S) fra D(B) e D(XC) | Si tiene solo per prove controllate, con lo scarto dichiarato. |
| | D(S) uguale a D(B) | Il generatore non batte il banale: si toglie. |

Prima del verdetto, i due sotto-agenti del passo 4 della skill difendono «continuare» e «chiudere» sui fatti completi. Nella Parte A non servono.

Verdetto, uno solo:

- chiudere, se la strada corta basta oppure se il modello della cattura non aggiunge;
- cambiare strada, se il modello della cattura vale e il generatore no: l'articolo usa disposizioni vere;
- continuare, se valgono tutti e due: allora si scrive il contratto in `docs/PAPER.md`, con pavimento, tetto e soglia;
- restringere, se il risultato vale solo per il tumore del colon.

Se lo Xenium non è usabile, la Parte B si fa con V, T, S e B. Decide solo sul generatore, e il verdetto lo dice.

## Che cosa non si fa

- Non si prova la semina a compartimenti e non si tocca il simulatore.
- Non si chiede al PI lavoro manuale: né annotazioni, né conteggi, né letture a vista.
- Non si aprono rene e fegato.
- Non si propone una terza sessione per la Parte B. Con due sessioni senza risposta il verdetto è «chiudere», con il motivo.
- Non si contano i costi già sostenuti.

## Che cosa consegna

- Parte A: `docs/riesame/FATTI.md` con i soli numeri riletti, `docs/ESISTENTE.md`, la tabella dei verdetti, la regola di esito con il suo commit, un riassunto «In breve» per il PI.
- Parte B: la tabella dei bracci con D e la sua incertezza, il verdetto, poi `docs/PAPER.md` oppure `docs/CHIUSURA.md`.

`docs/CHIUSURA.md`, se serve, sta in una pagina: che cosa abbiamo misurato, perché non basta, che cosa si riusa, quale vincolo riapre il progetto se cade.

## Che cosa si riusa in ogni caso

- Le misure sui Visium HD veri: densità, aree nucleari, spostamento dei trascritti.
- Il banco di prova dei metodi (`tools/F3.1_*` e `tools/F3.3_*`).
- I dati di verifica sigillati e le regole di DEC-093.

## Scommessa di chi scrive

È un'ipotesi da smentire, non un fatto: si tiene il modello della cattura e si lascia il generatore dello strato cellulare, sostituito da disposizioni vere.

## Limiti di questa consegna

- Chi scrive non ha visto il server né le trascrizioni delle sessioni.
- Dei dati Xenium degli stessi blocchi ha visto solo i titoli delle serie GEO.
- Di ENACT non ha letto il supplemento. Di scCube e degli altri simulatori ha letto abstract o estratti.
- Il paziente di taratura e il numero di geni del pannello Xenium sono da verificare.
