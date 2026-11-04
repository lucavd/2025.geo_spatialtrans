# Step 2 — Dalla mappa delle cellule ai conteggi di Visium HD: documento di progetto

**Versione**: 0.1 (2026-11-04, sessione S2.0; scritta con le misure di S2.0 in `docs/steps/S2.0.md`). **Stato**: proposta, da mandare a Daniele prima del codice. **Ramo**: `step2-visium-grid`.
**Scopo**: dire che cosa fa lo Step 2 del simulatore, con quali regole e con quali numeri, prima di scriverne il codice. Lo Step 1 (design `2026-04-29-step1-cell-layer-design.md`, v1.9) produce le cellule: posizione, tipo, territorio, nucleo. Lo Step 2 le trasforma nei conteggi che produrrebbe un vetrino Visium HD: quadratini (bin) di 2 µm, somma a 8 µm, e la verità di riferimento con cui si mettono alla prova i metodi di analisi.

**Come leggerlo.** La sezione 1 dice che cosa esce. La sezione 2 elenca le decisioni con la loro fonte. Le sezioni 3-8 descrivono i sei passi del modello. La sezione 9 è la tabella dei parametri per tessuto. Le sezioni 10-15 trattano lumi, verità, prestazioni, tolleranze, vecchia pipeline e piano di implementazione. La sezione 16 è la mappa della revisione di Daniele (2025, sezioni 3 e 4) e dei dubbi ricevuti.

## 1. Che cosa produce lo Step 2

Una lista `hd` con:

| oggetto | contenuto | a che cosa serve |
|---|---|---|
| `counts_2um` | matrice sparsa geni × bin di 2 µm, tutti i bin dell'area simulata (tessuto, lumi e vetro) | i dati «osservati» |
| `counts_8um` | somma 4 × 4 dei bin di 2 µm | i dati a 8 µm, come Space Ranger |
| `bins` | una riga per bin: riga e colonna della griglia, centro in µm, frazione di tessuto, `in_tissue` (frazione ≥ 0.5) | geometria dei dati |
| `truth$overlap` | per ogni coppia (bin, cellula): frazione dell'area del bin coperta dal territorio della cellula e dal suo nucleo | verità geometrica della revisione di Daniele (3.2.A) |
| `truth$origin` | per ogni coppia (bin, cellula): numero di trascritti catturati nel bin che vengono da quella cellula, più la quota del fondo | verità dopo lo spostamento: assegnazione dei bin alle cellule, deconvoluzione |
| `truth$cell_counts` | matrice geni × cellule dei trascritti generati per cellula (prima dello spostamento e della cattura tecnica) | profili veri |
| `truth$layers` | strati tecnici attivi e loro valori (fattori delle strisce per riga e colonna, fondo, parametri dello spostamento) | per mettere alla prova i metodi che correggono quegli strati |
| `metadata` | `pixel_size_um` ereditato dallo Step 1, seme, preset, versione | riproducibilità |

## 2. Decisioni

| n. | decisione | fonte |
|---|---|---|
| D2-1 | La griglia di 2 µm si costruisce dalla scala dello Step 1 (`metadata$pixel_size_um`, senza valore di default); l'area simulata è il riquadro delle cellule più un margine di vetro. Nessuna sotto-campionatura delle cellule. | decisione 10 del design dello Step 1 (C-A); chiude BL-005 e BL-010 |
| D2-2 | I trascritti si generano per cellula, si posizionano dentro la cellula, si spostano, poi si contano nei bin. Il bin non ha un profilo proprio. | revisione di Daniele 3.2.B, 3.2.C |
| D2-3 | Ogni tipo cellulare ha il suo numero di trascritti per cellula in tutti i tessuti (mediana e distribuzione di R4). | [DEC-035](/decisions/DEC-035.md) punto 1 |
| D2-4 | Il nucleo ha densità di trascritti r volte quella del citoplasma; r per tessuto invertito in S2.0 con la regola di spostamento scelta (sezione 9), sulla statistica pesata per cellula. Nel cuore il controllo interno è smentito: valore con riserva. | [DEC-034](/decisions/DEC-034.md), R5, S2.0 |
| D2-5 | Lo spostamento ha due livelli: una sfocatura gaussiana stretta per tutti i trascritti (circa 0.5-1 µm) e, per i trascritti che finiscono fuori dal tessuto, l'alone misurato ai bordi (gaussiana di 3-5 µm più coda esponenziale). La regola con cui un trascritto passa all'alone è quella scelta dai dati in S2.0 (sezione 6). Il cuore usa i valori mediani delle altre sezioni. | [DEC-037](/decisions/DEC-037.md); S2.0, M2: regola R-A, vince la zona di bordo ([DEC-039](/decisions/DEC-039.md)) |
| D2-6 | Le strisce della griglia a 2 µm sono uno strato tecnico che si può spegnere: alternanza fra righe e colonne pari e dispari più un fattore casuale per ogni riga e ogni colonna (sezione 7). | [DEC-035](/decisions/DEC-035.md) punto 2; S2.0, M1: P-1 smentita in 5 sezioni su 5, regola R-B ([DEC-039](/decisions/DEC-039.md)) |
| D2-7 | Dentro la cellula la densità dei trascritti varia su scala di 1-2 µm più del caso: uno strato moltiplicativo (gamma, media 1, varianza v per tessuto) è acceso nei preset e si può spegnere. | S2.0, M3: P-4 smentita per eccesso in 6 tessuti su 6, regola R-C, post hoc PH-2 e PH-3 ([DEC-039](/decisions/DEC-039.md)) |
| D2-8 | Nessun «dropout» come processo separato: gli zeri in eccesso osservati vengono da sottotipi o programmi dentro i tipi, quindi i profili per cellula hanno sottotipi o programmi. | R4 (PH-2), BL-194 |
| D2-9 | La profondità (trascritti per cellula) dei preset viene dalle aree misurate di R4; un fattore di profondità permette di passare alla media della sezione. | BL-195 (scelta da confermare da Luca) |
| D2-10 | Fondo sul vetro: un numero costante di trascritti per bin (misurato sul vetro lontano in R5), con la composizione media del tessuto. | R5 (D-3) |
| D2-11 | Lo spazio acellulare dentro il tessuto resta assente per costruzione (decisione 1 del design dello Step 1): lumi e vetro ricevono solo trascritti spostati e fondo. Limite dichiarato. | BL-077, BL-086, BL-201 |
| D2-12 | Tutte le sovrapposizioni geometriche fra bin e poligoni usano l'helper con arrotondamento a 1e-6 µm dello Step 1 (`overlay_prec`), mai `sf::st_intersection` sui territori. | S1.3b, S1.3d, S1.5 |

## 3. Passo 1: la griglia (revisione 3.2.A, 3.2.E; C-A)

La griglia è quella di Visium HD: quadrati di 2 µm, righe e colonne numerate dall'origine dell'area simulata, nessuno spazio fra i bin. Il centro del bin (i, j) è ((i + 0.5) · 2, (j + 0.5) · 2) µm nel riferimento delle cellule dello Step 1, già in µm grazie a `pixel_size_um`. L'area simulata è il riquadro dei territori più un margine di vetro di 100 µm, arrotondato a multipli di 8 µm così che la somma a 8 µm sia completa. Un'area di cattura intera (6.5 mm) ha 3 250 × 3 250 = 10.6 milioni di bin; il simulatore lavora a tessere (sezione 12).

La parità di riga e di colonna e i fattori delle strisce si legano all'indice del bin nella griglia simulata. La vecchia funzione `create_sampling_grid()` (`R/05_grid_sampling.R`) non si usa: la sua scala non coincideva con quella del resto della pipeline (`pixel_size_um = 1` contro 10, BL-005) e il numero di cellule non dipendeva dal parametro `n_cells` (BL-010). Nello Step 2 il numero di cellule è quello che lo Step 1 ha piazzato, e la scala è una sola.

Verità geometrica (revisione 3.2.A): per ogni bin, la frazione della sua area coperta da ciascun territorio e da ciascun nucleo, calcolata sulla griglia di 2 µm con l'helper di sovrapposizione dello Step 1 (D2-12). Somma delle frazioni di un bin ≤ 1; la differenza è vetro, lume o regione esclusa.

## 4. Passo 2: i profili per cellula (revisione 3.2.B, 4.1, 4.4, 4.5)

Ogni cellula riceve:

1. **un totale di trascritti** N_c, estratto dalla distribuzione per tipo e tessuto misurata in R4 (mediana e dispersione dei trascritti per cellula, righe B-275...B-335 di `BIO_REFERENCES.md`; [DEC-035](/decisions/DEC-035.md)). Il totale è sulla scala dei trascritti catturati (R4 misura la cattura, non le molecole vere), quindi non serve una seconda efficienza di cattura media. Profondità: D2-9 e la scelta di Luca su BL-195;
2. **un profilo di espressione** p_c (frazioni per gene) dal suo tipo, con un sottotipo o un programma (BL-194: nell'intestino gli zeri in eccesso vengono da marcatori di sottotipi dentro i tipi grossolani);
3. **i conteggi per gene** x_c ~ Multinomiale(N_c, p_c) o binomiale negativa con dispersione per gene (la vecchia `R/06d_dispersion_params.R` resta riusabile).

La fonte dei profili per tipo è una scelta di Luca (riassunto di S2.0). Il resto del documento non dipende dalla scelta.

Il numero di zeri a 2 µm e a 8 µm (revisione 4.4 e 4.5) non è un parametro: è una conseguenza di trascritti per cellula, densità delle cellule, spostamento e somma; il passo 6 lo controlla contro R4 (B-270...B-329). La variabilità fra bin della stessa cellula (revisione 4.1) è il passo 4 più lo strato del passo 3.

## 5. Passo 3: dove stanno i trascritti dentro la cellula (revisione 3.2.D, 4.2)

I trascritti di una cellula si posizionano nel suo territorio con densità uniforme nel citoplasma e r volte più alta nel nucleo (D2-4). Valori di r per tessuto nella sezione 9. La posizione non dipende dal gene: l'arricchimento nucleare per gene (revisione 4.2) resta un'estensione, non è nella prima versione (nessun riferimento misurato).

**Variabilità dentro la cellula oltre il caso (revisione 4.1; S2.0 M3).** Con trascritti posizionati a densità uniforme (r volte nel nucleo), spostati e catturati, i conteggi dei bin dentro lo stesso nucleo variano da 2.9 a 10.9 volte più del simulato (indice di dispersione osservato/simulato, 30 ROI: P-4 smentita per eccesso in 6 tessuti su 6). Circa metà dell'eccesso viene dalle strisce non periodiche (sezione 7), che il prototipo pre-registrato non conteneva: con i fattori di riga e colonna misurati il rapporto scende a 1.8-4.1 (post hoc PH-2). Il resto è una variazione moltiplicativa della densità con varianza v da 0.033 (cervello) a 0.095 (stroma); 0.153 nel cuore, dove il controllo interno è smentito. È correlata fra bin adiacenti (correlazione in più 0.09-0.38 a 2 µm) e quasi assente a 4 µm (post hoc PH-3): una struttura di 1-2 µm. Non è un gradiente lungo, come darebbe un nucleo segmentato male, e non è un rumore indipendente per bin. Se sia biologica (granuli, nucleoli, strutture del nucleo) o tecnica alla stessa scala questi dati non lo dicono.

Il simulatore la riproduce con un campo moltiplicativo gamma di media 1 e varianza v per tessuto sulla griglia fine dei trascritti, con correlazione di circa 1 µm, prima dello spostamento. v si calibra nello step di implementazione (S2.3) sull'indice di dispersione e sulla correlazione a 2 µm, con le strisce accese, perché il valore efficace dipende dal resto del modello; i valori della sezione 9 sono il punto di partenza. La misura è dentro il nucleo; il campo si applica a tutta la cellula (assunzione dichiarata). Strato spegnibile (`intracell_noise = FALSE`).

Opzione (spenta per difetto): gradiente radiale dal nucleo alla periferia nel cervello e nel cuore (BL-087); la misura del profilo per tessuto resta allo step che scriverà questa funzione.

## 6. Passo 4: lo spostamento dei trascritti (revisione 4.3; C-B; DEC-037)

Ogni trascritto si sposta dal punto in cui nasce al punto in cui viene catturato. R5 ha misurato due cose che una sola campana non concilia: ai bordi del tessuto verso il vetro i conteggi si allargano di 3-5 µm con una coda del 6-15 % fino a decine di µm; dentro il tessuto il contrasto fra nucleo e citoplasma richiede una sfocatura sotto 1 µm. Luca ha deciso due livelli ([DEC-037](/decisions/DEC-037.md)); S2.0 ha scelto dai dati la regola che li lega, fra tre candidate legate ai tre meccanismi possibili:

- **due popolazioni**: ovunque una parte dei trascritti si sposta poco (0.7 µm) e una parte di più; nessuna dipendenza dal luogo;
- **zona di bordo**: i trascritti nati vicino al bordo del tessuto (peso exp(−d/D) alla distanza d dal bordo) si spostano con il nucleo largo, gli altri con quello stretto (bordo sottile o sfrangiato);
- **trasporto sul vetro**: tutti si spostano poco; quelli che atterrano fuori dal tessuto si spostano ancora (vetro libero).

In tutte la coda esponenziale vale per tutti i trascritti. S2.0 (M2) ha adattato le tre regole e il nucleo di R5 alle 314 finestre di bordo di R5 (nessun errore; il nucleo di R5 ricalcolato coincide con R5 entro 0.008 unità di log-verosimiglianza). Vince la **zona di bordo** in tutte e quattro le sezioni giudicate (somma degli AIC). Mediana per finestra contro il nucleo di R5: −116 nell'intestino, −691 nel tumore del colon, −49 nel cervello; nel linfonodo +1.8, cioè lì le regole non si distinguono. Le due popolazioni si riducono quasi al nucleo di R5 (frazione stretta 0.66-0.92); il trasporto sul vetro è peggiore del nucleo di R5 in tutte le sezioni (sfocatura al limite inferiore). La zona di bordo è profonda 7-25 µm (D) con una sfocatura larga di 4.5-6.2 µm; oltre qualche D dal bordo resta la sola sfocatura stretta di 0.7 µm.

Il controllo interno (contrasto fra nucleo e citoplasma senza i bin a ridosso del bordo del nucleo, con r invertito con la zona di bordo) è confermato nel tumore, nello stroma e nel cervello, indeterminato nell'intestino (0.90) e nel linfonodo (0.88), smentito nel cuore (0.66: il modello non descrive l'interno dei nuclei del cuore). La controprova sui lumi interni (CP-1), mai usati nell'adattamento, dà con la zona di bordo nessuna delle 8 combinazioni sezione × fascia smentita (6 confermate, 2 indeterminate); il nucleo di R5 ricalcolato ridà esattamente stime e intervalli di R5.

In pratica: un trascritto nato a distanza d dal bordo del tessuto o di un lume usa la sfocatura larga σ_e con probabilità exp(−d/D), altrimenti quella stretta σ_in; poi la coda esponenziale (massa w, lunghezza ℓ) vale per tutti. La regola è misurata sul bordo esterno del tessuto; nei lumi la si applica con la distanza dal bordo del lume, e CP-1 dice che funziona entro 8 µm. Valori per sezione nella sezione 9; il cuore usa la mediana delle altre quattro.

## 7. Passo 5: cattura e strati tecnici (revisione 3.2.C)

1. **Strisce** (D2-6). Ogni bin ha un fattore di cattura c_ij = p(parità di i, parità di j) · exp(α_i + β_j): la parte periodica p ha quattro classi normalizzate a media 1, con i valori della sezione; α_i e β_j sono fattori casuali per riga e per colonna. S2.0 (M1) ha misurato le due parti nelle cinque sezioni. La parte periodica (rapporto dispari/pari) va da 0.95 (righe del cuore) a 1.75 (colonne dell'intestino) e, salvo nel linfonodo, è la stessa in tutta la sezione (deviazione standard fra tessere di 300 µm ≤ 0.02). La parte non periodica non è trascurabile: deviazione standard 0.18-0.28 in scala logaritmica, riproducibile fra le due metà di ogni riga e colonna (correlazione 0.73-0.96; 0.39 nelle colonne del linfonodo), con correlazione negativa fra righe o colonne adiacenti (da −0.07 a −0.29). Gli stessi fattori si ritrovano sul vetro lontano dal tessuto (analisi post hoc PH-1), quindi sono cattura tecnica per riga e colonna, non tessuto. Il simulatore genera α_i e β_j gaussiani con la deviazione standard e la correlazione a passo 1 della sezione (media mobile di ordine 1). Nel linfonodo la parità cambia nello spazio (deviazione standard fra tessere 0.09 nelle righe, 0.05 nelle colonne): il preset la fa variare a tessere di 300 µm con quella deviazione standard. Valori nella sezione 9. Lo strato si spegne con un argomento (`stripes = FALSE`): i metodi che correggono le strisce si provano con lo strato acceso, gli altri si possono provare anche senza.
2. **Fondo sul vetro** (D2-10): b trascritti per bin in media (Poisson), con la composizione per gene media del tessuto; sul tessuto il fondo si somma ai trascritti delle cellule.
3. **Campionamento**: i trascritti spostati e il fondo si contano nel bin in cui cadono; i fattori delle strisce agiscono come assottigliamento binomiale (si generano N_c · max(c) / media(c) trascritti e si tiene ciascuno con probabilità c_bin / max(c)), così i conteggi restano interi e la verità di origine resta esatta.

## 8. Passo 6: somma a 8 µm e uscite (revisione 3.2.E)

`counts_8um` somma i 16 bin di 2 µm di ogni bin di 8 µm (righe e colonne a blocchi di 4, come Space Ranger). La verità di origine e quella geometrica si sommano allo stesso modo. Controllo: ogni bin di 8 µm è uguale alla somma esatta dei suoi 16 (identità su interi; R5 l'ha verificata sui dati veri).

## Tabella dei parametri (sezione 9)

Valori per tessuto; la colonna «fonte» rimanda alla riga di `docs/BIO_REFERENCES.md`, alla decisione o al file del repository che li ha prodotti.

| parametro | intestino (A1) | tumore (A2) | stroma (A3) | linfonodo (A4) | cervello (A5) | cuore (A6) | fonte |
|---|---|---|---|---|---|---|---|
| trascritti per cellula | per tipo, mediana e dispersione di R4 | per tipo, mediana e dispersione di R4 | per tipo, mediana e dispersione di R4 | per tipo, mediana e dispersione di R4 | per tipo, mediana e dispersione di R4 | per tipo, mediana e dispersione di R4 | B-275...B-335; [DEC-035](/decisions/DEC-035.md) |
| rapporto nucleo/citoplasma r (regola C2) | 0.75 | 0.76 | 0.82 | 0.55 | 1.48 | 2.61* | B-384, B-386, B-388, B-390, B-392, B-394; *controllo interno smentito |
| sfocatura stretta per tutti σ_in (µm) | 0.7 | 0.7 | 0.7 | 0.7 | 0.7 | 0.7 | valore fisso di S2.0 dentro 0.5-1 µm di [DEC-037](/decisions/DEC-037.md); B-367...B-372 |
| zona di bordo: sfocatura larga σ_e (µm) | 6.23 | 4.72 | 4.72 | 4.49 | 4.99 | 4.86† | B-379, B-380, B-381, B-382, B-383; †mediana delle altre 4 sezioni ([DEC-037](/decisions/DEC-037.md)) |
| zona di bordo: profondità D (µm) | 24.75 | 8.94 | 8.94 | 22.04 | 6.86 | 15.49† | B-379, B-380, B-381, B-382, B-383; †mediana delle altre 4 sezioni ([DEC-037](/decisions/DEC-037.md)) |
| coda esponenziale: massa w | 0.087 | 0.035 | 0.035 | 0.109 | 0.112 | 0.098† | B-379, B-380, B-381, B-382, B-383; †mediana delle altre 4 sezioni ([DEC-037](/decisions/DEC-037.md)) |
| coda esponenziale: lunghezza ℓ (µm) | 30.28 | 19.28 | 19.28 | 38.08 | 33.07 | 31.67† | B-379, B-380, B-381, B-382, B-383; †mediana delle altre 4 sezioni ([DEC-037](/decisions/DEC-037.md)) |
| fondo sul vetro b (UMI per bin di 2 µm) | 0.161 | 0.197 | 0.197 | 0.042 | 0.309 | 0.117 | B-342, B-348, B-354, B-360, B-366 |
| strisce: parità righe dispari/pari | 1.417 | 1.578 | 1.578 | 0.980‡ | 1.178 | 0.947 | B-374, B-375, B-376, B-377, B-378; ‡variabile nello spazio |
| strisce: parità colonne dispari/pari | 1.750 | 1.724 | 1.724 | 1.125‡ | 1.231 | 0.967 | B-374, B-375, B-376, B-377, B-378 |
| strisce: SD dei fattori per riga α (log) | 0.205 | 0.225 | 0.225 | 0.225 | 0.248 | 0.195 | B-374, B-375, B-376, B-377, B-378 |
| strisce: SD dei fattori per colonna β (log) | 0.188 | 0.196 | 0.196 | 0.282 | 0.184 | 0.190 | B-374, B-375, B-376, B-377, B-378 |
| strisce: correlazione a passo 1 di α; β | -0.19; -0.25 | -0.29; -0.18 | -0.29; -0.18 | -0.21; -0.12 | -0.26; -0.27 | -0.07; -0.22 | B-374, B-375, B-376, B-377, B-378 |
| variabilità dentro il nucleo v (gamma moltiplicativa) | 0.044 | 0.069 | 0.095 | 0.063 | 0.033 | 0.153* | B-385, B-387, B-389, B-391, B-393, B-395; post hoc PH-2 |
| scala della variabilità dentro il nucleo | 1-2 µm | 1-2 µm | 1-2 µm | 1-2 µm | 1-2 µm | 1-2 µm | B-385, B-387, B-389, B-391, B-393, B-395; post hoc PH-3 (correlazione in più a 2 µm, quasi nulla a 4 µm) |

Il tumore (A2) e lo stroma (A3) stanno sulla stessa sezione: spostamento, fondo e strisce sono quelli della sezione; r e v sono per tessuto.

## 10. Lumi, vetro e spazio acellulare (BL-077, BL-086, BL-201)

Nello Step 1 ogni punto del tessuto appartiene a un territorio cellulare (decisione 1 del design dello Step 1): non c'è matrice senza cellule dentro le regioni (BL-077). Lo Step 2 non la aggiunge: lumi e vetro ricevono solo i trascritti spostati dalle cellule vicine e il fondo. Entro 8 µm dal bordo dei lumi questo riproduce i lumi veri entro un quarto (R5 CP-2: 9 combinazioni su 10; controprova di S2.0 con la regola scelta: con la zona di bordo nessuna delle 8 combinazioni sezione × fascia smentita, 6 confermate e 2 indeterminate). Oltre 8 µm i lumi veri del tumore del colon hanno 1.5 volte i trascritti previsti (materiale nelle ghiandole, BL-201) e quelli del cuore meno: limite dichiarato; un'opzione «materiale nei lumi» (sorgente debole uniforme nei lumi) resta per quando avremo lumi annotati. I conteggi a 2 µm non distinguono un vuoto piccolo dentro il tessuto (BL-086): i metodi che cercano lo spazio acellulare vanno provati con vuoti grandi.

## 11. Verità di riferimento e metodi che servono (revisione 3.2.A, 4.3)

| metodo da mettere alla prova | verità |
|---|---|
| assegnazione dei bin alle cellule, segmentazione a partire dai conteggi | `truth$origin` (quota di ogni cellula in ogni bin), `truth$overlap` |
| deconvoluzione a 8 µm | `truth$origin` sommata a 8 µm e aggregata per tipo |
| correzione delle strisce | `truth$layers$stripes` (parità, fattori per riga e colonna) |
| correzione del fondo e dell'RNA ambientale | `truth$layers$background`, quota del fondo in `truth$origin` |
| correzione dello spostamento (diffusione) | parametri in `truth$layers$diffusion`; differenza fra `truth$overlap` e `truth$origin` |
| analisi di espressione per tipo e sottotipo | `truth$cell_counts`, tipo e sottotipo per cellula |

La revisione di Daniele (3.2.A) chiamava «verità» la sola frazione geometrica. Con lo spostamento la geometria non basta: un bin riceve trascritti da cellule che non tocca. Per questo lo Step 2 dà entrambe.

## 12. Prestazioni (BL-045, BL-057, BL-072, BL-144)

Ordini di grandezza (da misurare negli step d'implementazione): un'area di 1 mm² ha 250 000 bin di 2 µm e da circa 3 000 (cuore) a 12 000 cellule (intestino, linfonodo); i trascritti vanno da circa 1 milione (cuore, linfonodo) a circa 30 milioni (cervello: 6 742 per cellula, B-320). Un'area di cattura intera è 42 volte più grande. Bilancio proposto, come per lo Step 1: un'area di 1 mm² in meno di 5 minuti e 8 GB; l'area intera a tessere di 500 µm con un bordo di 200 µm (oltre 5 lunghezze della coda più lunga, 38 µm), un processo per tessera. Lo Step 1 sul linfonodo 4 000 × 4 000 px supera già il bilancio (310 s, 8.5 GB, BL-144): il ritaglio dei territori (BL-072) e le regioni a 2 µm (BL-045) vanno ottimizzati prima dell'area intera; la griglia dei vicini su estensioni enormi (BL-057) non serve allo Step 2, che lavora a tessere.

## 13. Tolleranze dei controlli (BL-146)

- Conteggi: identità esatte su interi (trascritti generati = trascritti spostati = catturati + persi per assottigliamento; somma a 8 µm).
- Frazioni geometriche: somma delle frazioni di un bin entro g · P_bin + g² con g = 1e-6 µm e P_bin = 8 µm (perimetro), come la regola di S1.5 per i pezzi piccoli (BL-093); la somma per regione si lega alla griglia, non all'area della regione (risposta a BL-146).
- Distribuzioni simulate contro riferimenti: errore Monte Carlo fra semi e intervalli di rilevanza pre-registrati negli step.

## 14. Che cosa si riusa della vecchia pipeline (BL-001, BL-005, BL-010)

| file | che cosa fa | nello Step 2 |
|---|---|---|
| `R/05_grid_sampling.R` | griglia di campionamento con scala propria | sostituito dal passo 1 (BL-005, BL-010) |
| `R/06a_expression_params.R`, `R/06k_expression_profiles_wrapper.R` | parametri e involucro della generazione per spot | sostituiti dall'orchestratore dello Step 2 |
| `R/06b_expression_baseline.R` | profili sintetici per tipo con geni marcatori | riusato solo se Luca sceglie i profili sintetici |
| `R/06d_dispersion_params.R` | dispersione binomiale negativa | riusabile per la dispersione per gene |
| `R/06e_library_size.R` | profondità per spot (8 000 UMI, Visium classico) | sostituito dai trascritti per cellula di R4 (B-275...B-335) |
| `R/06f_dropout_models.R` | dropout e RNA ambientale | non usato: niente dropout separato (D2-8); il fondo sostituisce l'RNA ambientale (D2-10). La divergenza fra le due configurazioni di prova (`dropout_range` 0.45-0.65 contro 0.25-0.45, BL-001) non riguarda più il simulatore |
| `R/06g_gene_modules.R`, `R/06h_spatial_correlation.R` | moduli di geni e campi gaussiani | opzione per i programmi dentro i tipi (BL-194) |
| `R/06i_hybrid_cells.R` | cellule ibride | non usato: le miscele nascono dalla geometria e dallo spostamento |
| `R/06j_expression_generation.R` | campionamento per spot | il campionamento binomiale negativo si riusa per cellula |

## 15. Piano di implementazione (una funzione per sessione)

| sessione | funzione | controlli | riferimenti | controprova |
|---|---|---|---|---|
| S2.1 | `build_hd_grid()` e `truth$overlap` | somma delle frazioni, identità della somma a 8 µm, geometrie valide | — | griglia spostata di mezzo bin: le frazioni si ricalcolano identiche a meno della traslazione |
| S2.2 | `simulate_cell_profiles()` | totali per cellula, multinomiale, sottotipi | trascritti per cellula per tipo (R4); profili scelti da Luca | profili veri recuperati dalla media delle cellule simulate |
| S2.3 | `place_and_move_transcripts()` | conservazione dei trascritti, densità nucleo/citoplasma, nucleo di spostamento | R5 (bordi, contrasto interno), regola di S2.0 | controllo interno e lumi di S2.0 rifatti sulle cellule simulate |
| S2.4 | `capture_and_bin()` | strisce, fondo, assottigliamento, somma a 8 µm | R4 (strisce, vuoti, UMI per bin), M1 di S2.0 | strisce stimate sui dati simulati con lo stimatore di S2.0 |
| S2.5 | `simulate_visium_hd()` e pannello A1-A6 | tutti i controlli dello Step 2 | R4 per tessuto (UMI per bin a 2 e 8 µm, frazione di bin vuoti, geni per bin, rapporto nucleo/corona) | il cuore come controllo negativo dove il Voronoi isotropo fallisce |

## 16. Mappa della revisione di Daniele e dei dubbi ricevuti

| elemento | dove |
|---|---|
| revisione 3.2.A (mappatura bin → cellula) | sezioni 3, 11 |
| revisione 3.2.B (profili per cellula) | sezione 4 |
| revisione 3.2.C (aggregazione cellula → bin, fondo, rumore di cattura) | sezioni 6, 7 |
| revisione 3.2.D (variazione dentro la cellula, nucleo e citoplasma) | sezione 5 |
| revisione 3.2.E (somma 2 → 8 µm) | sezione 8 |
| revisione 4.1 (variazione dentro la cellula) | sezioni 4, 5 |
| revisione 4.2 (arricchimento nucleare per gene) | sezione 5 (estensione, non nella prima versione) |
| revisione 4.3 (bin di confine, diffusione) | sezioni 6, 11 |
| revisione 4.4 (trascritti per bin a 8 µm) | sezioni 4, 9 |
| revisione 4.5 (sparsità a 8 µm) | sezioni 4, 9 |
| BL-197 regola d'uscita | sezione 6 |
| BL-192 strisce | sezione 7 |
| BL-194 sottotipi | sezione 4 |
| BL-195 profondità | sezioni 2 (D2-9), 4 |
| BL-001, BL-005, BL-010 vecchia pipeline | sezioni 3, 14 |
| BL-146 tolleranze | sezione 13 |
| BL-201, BL-077, BL-086 lumi e spazio acellulare | sezione 10 |
| BL-087 gradiente radiale | sezione 5 |
| BL-045, BL-057, BL-072, BL-144 prestazioni | sezione 12 |
| [DEC-034](/decisions/DEC-034.md), [DEC-035](/decisions/DEC-035.md), [DEC-037](/decisions/DEC-037.md) | sezione 2 |
| [DEC-038](/decisions/DEC-038.md) scelte tecniche di S2.0 | sezioni 6, 7, 13 |
| [DEC-039](/decisions/DEC-039.md) regole decise dai dati in S2.0 (R-A, R-B, R-C) | sezioni 5, 6, 7, 9 |

## 17. Limiti dichiarati

- Spazio acellulare dentro il tessuto assente (sezione 10); lumi oltre 8 µm del tumore del colon e del cuore non riprodotti (BL-201).
- Arricchimento nucleare per gene assente (revisione 4.2).
- Variabilità dentro la cellula misurata solo dentro il nucleo e applicata a tutta la cellula; origine (biologica o tecnica) non nota (sezione 5).
- Cuore: il modello non descrive l'interno dei nuclei (controllo interno smentito); r e v del cuore con riserva.
- Cuore: spostamento con i valori delle altre sezioni ([DEC-037](/decisions/DEC-037.md)); cellule isotrope (limite dello Step 1).
- I limiti dello Step 1 passano ai conteggi: contatti fra nuclei in eccesso nei tessuti densi, territori troppo uniformi con la semina uniforme (design dello Step 1, 7.4).

## 18. Storico

| data | versione | che cosa | chi |
|---|---|---|---|
| 2026-11-04 | 0.1 | prima bozza con le misure di S2.0 | geo-agent/skill-v2 (S2.0) |