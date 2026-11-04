# Step 2 — Dalla mappa delle cellule ai conteggi di Visium HD: documento di progetto

**Versione**: 0.2 (2026-11-04, sessione S2.0, dopo la revisione avversariale; scritta con le misure di S2.0 in `docs/steps/S2.0.md`). **Stato**: proposta, da mandare a Daniele prima del codice. **Ramo**: `step2-visium-grid`.
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
| D2-5 | Lo spostamento ha due livelli, legati dalla regola scelta dai dati in S2.0 (zona di bordo, sezione 6): tutti i trascritti hanno la sfocatura gaussiana stretta σ_in = 0.7 µm, salvo quelli nati vicino al bordo del tessuto o di un lume (peso exp(−d/D) alla distanza d dal bordo), che hanno la sfocatura larga σ_e (4.5-6.2 µm); una coda esponenziale (massa w, lunghezza ℓ) vale per tutti. La regola precisa la formulazione di DEC-037 («per i trascritti che finiscono fuori dal tessuto»): l'alone dipende da dove il trascritto nasce, non da dove arriva; DEC-037 affidava la regola a questo documento (BL-197). Sostituisce il nucleo gaussiano unico `diffusion_sigma_um` della decisione 11 del design dello Step 1. Il cuore usa i valori mediani delle altre sezioni. | [DEC-037](/decisions/DEC-037.md); S2.0, M2: regola R-A ([DEC-039](/decisions/DEC-039.md)) |
| D2-6 | Le strisce della griglia a 2 µm sono uno strato tecnico che si può spegnere: alternanza fra righe e colonne pari e dispari più un fattore casuale per ogni riga e ogni colonna (sezione 7). | [DEC-035](/decisions/DEC-035.md) punto 2; S2.0, M1: P-1 smentita in 5 sezioni su 5, regola R-B ([DEC-039](/decisions/DEC-039.md)) |
| D2-7 | La densità dei trascritti varia su scala di 1-2 µm più del caso, dentro i nuclei (misura pre-registrata) e probabilmente in tutto il tessuto (revisione G-5): uno strato moltiplicativo (gamma, media 1, varianza per tessuto) è acceso nei preset e si può spegnere; valori al netto delle strisce e forma di campo correlato, tarati in S2.3 ([DEC-040](/decisions/DEC-040.md)). | S2.0, M3: P-4 smentita per eccesso in 6 tessuti su 6, regola R-C, post hoc PH-2 e PH-3 ([DEC-039](/decisions/DEC-039.md)) |
| D2-8 | Nessun «dropout» come processo separato: gli zeri in eccesso osservati vengono da sottotipi o programmi dentro i tipi, quindi i profili per cellula hanno sottotipi o programmi. | R4 (PH-2), BL-194 |
| D2-9 | I trascritti per cellula restano per tipo (R4); un fattore di profondità per tessuto li porta a riprodurre gli UMI veri per bin di 2 µm in tessuto (sezione 9, righe di riferimento). Il fattore si tara in S2.2 e si verifica in S2.5; senza, il simulatore produrrebbe meno RNA del vero (stima provvisoria della revisione G-5: 67-82 %, un terzo nel cuore). | BL-195, BL-167; [DEC-040](/decisions/DEC-040.md) |
| D2-10 | Fondo sul vetro: un numero costante di trascritti per bin (misurato sul vetro lontano in R5), con la composizione media del tessuto. | R5 (D-3) |
| D2-11 | Lo spazio acellulare dentro il tessuto resta assente per costruzione (decisione 1 del design dello Step 1): lumi e vetro ricevono solo trascritti spostati e fondo. Limite dichiarato. | BL-077, BL-086, BL-201 |
| D2-12 | Tutte le sovrapposizioni geometriche fra bin e poligoni usano l'helper con arrotondamento a 1e-6 µm dello Step 1 (`overlay_prec`), mai `sf::st_intersection` sui territori. | S1.3b, S1.3d, S1.5 |

## 3. Passo 1: la griglia (revisione 3.2.A, 3.2.E; C-A)

La griglia è quella di Visium HD: quadrati di 2 µm, nessuno spazio fra i bin; i è la riga (asse y), j la colonna (asse x), numerate da 0 dall'origine (x₀, y₀) dell'area simulata. Il centro del bin (i, j) è (x₀ + (j + 0.5) · 2, y₀ + (i + 0.5) · 2) µm nel riferimento delle cellule dello Step 1, già in µm grazie a `pixel_size_um`. L'area simulata è il riquadro dei territori più un margine di vetro di 280 µm (multiplo di 8 µm, così la somma a 8 µm è completa). La coda in due dimensioni lascia oltre la distanza r la frazione (1 + r/ℓ) exp(−r/ℓ) della sua massa: con ℓ = 41.5 µm (la più lunga, linfonodo) il 31 % oltre 100 µm, il 4.7 % oltre 200 µm, lo 0.9 % oltre 280 µm (6.7 ℓ); i trascritti che escono dall'area si contano come persi (sezione 13) così che la somma a 8 µm sia completa. Un'area di cattura intera ha 3 350 × 3 350 = 11.2 milioni di bin (la matrice di Space Ranger letta in S2.0, 6.7 mm di lato; l'area nominale è 6.5 mm); il simulatore lavora a tessere (sezione 12).

La parità di riga e di colonna e i fattori delle strisce si legano all'indice del bin nella griglia simulata. La vecchia funzione `create_sampling_grid()` (`R/05_grid_sampling.R`) non si usa: la sua scala non coincideva con quella del resto della pipeline (`pixel_size_um = 1` contro 10, BL-005) e il numero di cellule non dipendeva dal parametro `n_cells` (BL-010). Nello Step 2 il numero di cellule è quello che lo Step 1 ha piazzato, e la scala è una sola.

Verità geometrica (revisione 3.2.A): per ogni bin, la frazione della sua area coperta da ciascun territorio e da ciascun nucleo, calcolata sulla griglia di 2 µm con l'helper di sovrapposizione dello Step 1 (D2-12). Somma delle frazioni di un bin ≤ 1; la differenza è vetro, lume o regione esclusa.

## 4. Passo 2: i profili per cellula (revisione 3.2.B, 4.1, 4.4, 4.5)

Ogni cellula riceve:

1. **un totale di trascritti** N_c, estratto dalla distribuzione per tipo e tessuto misurata in R4 (mediana e dispersione dei trascritti per cellula: per tessuto B-275, B-285, B-295, B-307, B-320, B-332, per tipo B-276...B-278, B-286...B-288, B-296...B-300, B-308...B-313, B-321...B-325, B-333...B-335 di `BIO_REFERENCES.md`; [DEC-035](/decisions/DEC-035.md)). Il totale è sulla scala dei trascritti catturati (R4 misura la cattura, non le molecole vere), quindi non serve una seconda efficienza di cattura media. Profondità: D2-9;
2. **un profilo di espressione** p_c (frazioni per gene) dal suo tipo, con un sottotipo o un programma (BL-194: nell'intestino gli zeri in eccesso vengono da marcatori di sottotipi dentro i tipi grossolani);
3. **i conteggi per gene** x_c ~ Multinomiale(N_c, p_c) o binomiale negativa con dispersione per gene (la vecchia `R/06d_dispersion_params.R` resta riusabile).

La fonte dei profili per tipo la decidono i dati, in S2.2, con una regola scritta prima del calcolo ([DEC-040](/decisions/DEC-040.md)): candidate i profili medi dei nuclei veri di R6 (impronta intera, oppure solo l'interno del nucleo, meno contaminato dai vicini) e un atlante a cellula singola del tessuto; vince quella che, passata nel modello del simulatore, riproduce meglio i valori veri di R4 per bin (geni rilevati, sparsità a 2 e 8 µm, somiglianza fra bin vicini). Il resto del documento non dipende dalla scelta. Il confronto del modello di rumore per gene (binomiale negativa contro scDesign3) è BL-203.

Il numero di zeri a 2 µm e a 8 µm (revisione 4.4 e 4.5) non è un parametro: è una conseguenza di trascritti per cellula, densità delle cellule, spostamento e somma; il passo 6 lo controlla contro R4 (UMI per bin a 2 e 8 µm e sparsità: B-270...B-272, B-280...B-282, B-290...B-292, B-302...B-304, B-315...B-317, B-327...B-329). La variabilità fra bin della stessa cellula (revisione 4.1) è il passo 4 più lo strato del passo 3.

## 5. Passo 3: dove stanno i trascritti dentro la cellula (revisione 3.2.D, 4.2)

I trascritti di una cellula si posizionano nel suo territorio con densità uniforme nel citoplasma e r volte quella del citoplasma nel nucleo (D2-4): r < 1 (nucleo meno denso) in intestino, tumore, stroma e linfonodo, r > 1 in cervello e cuore. Valori di r per tessuto nella sezione 9. La posizione non dipende dal gene: l'arricchimento nucleare per gene (revisione 4.2) resta un'estensione, non è nella prima versione (nessun riferimento misurato).

**Variabilità dentro la cellula oltre il caso (revisione 4.1; S2.0 M3).** Con trascritti posizionati a densità uniforme (r volte nel nucleo), spostati e catturati, i conteggi dei bin dentro lo stesso nucleo variano da 2.9 a 10.9 volte più del simulato (indice di dispersione osservato/simulato, 30 ROI: P-4 smentita per eccesso in 6 tessuti su 6). Una parte dell'eccesso viene dalle strisce non periodiche (sezione 7): circa metà in intestino, tumore, stroma e linfonodo (0.44-0.60), sette decimi nel cervello, un quarto nel cuore, che il prototipo pre-registrato non conteneva: con i fattori di riga e colonna misurati il rapporto scende a 1.8-4.1 (post hoc PH-2). Il resto è una variazione moltiplicativa della densità con varianza v da 0.033 (cervello) a 0.095 (stroma); 0.153 nel cuore, dove il controllo interno è smentito. È correlata fra bin adiacenti (correlazione in più 0.09-0.38 a 2 µm) e molto meno a 4 µm (al più 0.09, nel cervello; post hoc PH-3): una struttura di 1-2 µm. Non è un gradiente lungo, come darebbe un nucleo segmentato male, e non è un rumore indipendente per bin. Se sia biologica (granuli, nucleoli, strutture del nucleo) o tecnica alla stessa scala questi dati non lo dicono.

Il simulatore la riproduce con un campo moltiplicativo gamma di media 1 sulla griglia fine dei trascritti, con correlazione spaziale corta (nucleo gaussiano di circa 1 µm, da calibrare), prima dello spostamento. v è la varianza in più per bin di 2 µm, non la varianza puntuale del campo, che è più grande; la varianza del campo e la sua lunghezza di correlazione si calibrano nello step di implementazione (S2.3) sull'indice di dispersione e sulla correlazione a 2 µm, con le strisce accese, perché il valore efficace dipende dal resto del modello; i valori della sezione 9 sono il punto di partenza. La misura pre-registrata è dentro il nucleo. La revisione G-5 trova un eccesso simile o più grande anche nella corona di 8 µm (v 0.12-0.21 con i fattori di riga e colonna) e nello spazio acellulare (0.04-0.23), con correlazione positiva fra bin adiacenti anche lì (7 ROI; provvisorio, perché per la corona manca un nullo simulato): il campo si applica quindi a tutte le sorgenti di trascritti, cioè a tutto il tessuto, non al solo nucleo. Nell'intestino v è probabilmente sottostimato di circa il 20 %: i bin profondi sono scelti sulla maschera dei nuclei non spostata dello scarto della sezione (−2.2 µm), che gonfia l'indice simulato (revisione G-5). Strato spegnibile (`intracell_noise = FALSE`).

Opzione (spenta per difetto): gradiente radiale dal nucleo alla periferia nel cervello e nel cuore (BL-087); la misura del profilo per tessuto resta allo step che scriverà questa funzione.

## 6. Passo 4: lo spostamento dei trascritti (revisione 4.3; C-B; DEC-037)

Ogni trascritto si sposta dal punto in cui nasce al punto in cui viene catturato. R5 ha misurato due cose che una sola campana non concilia: ai bordi del tessuto verso il vetro i conteggi si allargano di 3-5 µm con una coda del 6-15 % (0.38 nel cuore, poco determinata) fino a decine di µm; dentro il tessuto il contrasto fra nucleo e citoplasma richiede una sfocatura sotto 1 µm. Luca ha deciso due livelli ([DEC-037](/decisions/DEC-037.md)); S2.0 ha scelto dai dati la regola che li lega, fra tre candidate legate ai tre meccanismi possibili:

- **due popolazioni**: ovunque una parte dei trascritti si sposta poco (0.7 µm) e una parte di più; nessuna dipendenza dal luogo;
- **zona di bordo**: i trascritti nati vicino al bordo del tessuto (peso exp(−d/D) alla distanza d dal bordo) si spostano con il nucleo largo, gli altri con quello stretto (bordo sottile o sfrangiato);
- **trasporto sul vetro**: tutti si spostano poco; quelli che atterrano fuori dal tessuto si spostano ancora (vetro libero).

In tutte la coda esponenziale vale per tutti i trascritti. S2.0 (M2) ha adattato le tre regole e il nucleo di R5 alle 314 finestre di bordo di R5 (nessun errore; il nucleo di R5, che la pre-registrazione non chiedeva di rifare, è stato ricalcolato come controllo e coincide con R5 entro 0.008 unità di log-verosimiglianza). Vince la **zona di bordo** in tutte e quattro le sezioni giudicate (somma degli AIC). Mediana per finestra contro il nucleo di R5: −116 nell'intestino, −691 nel tumore del colon, −49 nel cervello; nel linfonodo +1.8, cioè lì le regole non si distinguono. Le due popolazioni si riducono quasi al nucleo di R5: la frazione con la sfocatura larga vale 0.66-0.92 e la parte stretta è nulla in mediana in intestino, tumore e linfonodo (0.14 nel cervello). Il trasporto sul vetro è peggiore del nucleo di R5 in tutte le sezioni e degenera (sfocatura sul vetro e lunghezza della coda ai limiti inferiori, coda di massa 0.46-0.86). In 5 finestre su 314 la zona di bordo è peggiore del nucleo di R5 di 0.6-1.3 unità, sempre con D al limite di 10⁴ µm: è l'effetto atteso della correzione AM-1 (con D grande le due regole non coincidono esattamente), non un difetto dell'ottimizzatore. La zona di bordo è profonda 7-25 µm (D) con una sfocatura larga di 4.5-6.2 µm; oltre qualche D dal bordo resta la sola sfocatura stretta di 0.7 µm.

La revisione avversariale ha trovato soluzioni ferme con la coda spenta (w = 0): tutte le partenze ereditavano w dal nucleo di R5 della finestra, e con w = 0 l'ottimizzatore non si muove. Il nuovo adattamento con partenze a coda accesa (post hoc PH-6) migliora la verosimiglianza di oltre 1 unità in 13 coppie finestra × regola su 107 con w = 0 (fino a 174 unità in una finestra del tumore, anche con il nucleo di R5). Non cambia nessun esito pre-registrato né la regola vincente; cambiano la lunghezza della coda (intestino 30.3 → 32.9 µm, tumore 19.3 → 24.7, linfonodo 38.1 → 41.5) e la sua massa nel tumore (0.035 → 0.040). I valori della sezione 9 e i controlli che seguono sono quelli corretti.

Il controllo interno (contrasto fra nucleo e citoplasma senza i bin a ridosso del bordo del nucleo, con r invertito con la regola) era pre-registrato nel tumore e nel cervello, come condizione per ammettere una regola. Con la zona di bordo è confermato in entrambi (0.999 e 0.976). Con le due popolazioni (P-3) è confermato nel cervello (1.010) e «non calcolabile come scritto» nel tumore (una ROI su 5 non si inverte; sulle altre 4, 0.867 [0.841, 0.886], che sarebbe smentita): R-A ha letto «non calcolabile» come non smentita e ha ammesso le due popolazioni, ma la vincente non cambia. Negli altri tessuti il controllo è descrittivo, giudicato con lo stesso intervallo: con la zona di bordo stroma 0.972, intestino 0.90 e linfonodo 0.88 (indeterminati), cuore 0.66 (il modello non descrive l'interno dei nuclei del cuore). La controprova sui lumi interni (CP-1), mai usati nell'adattamento, dà con la zona di bordo nessuna delle 8 combinazioni sezione × fascia smentita (6 confermate, 2 indeterminate); il nucleo di R5 ricalcolato ridà esattamente stime e intervalli di R5. Nel cervello (2 lumi, descrittivo) la zona di bordo va peggio del nucleo di R5 nella fascia 4-8 µm (0.71 contro 0.81).

Il modello adattato (C-2 della pre-registrazione) è una miscela: per una sorgente di trascritti T sulla griglia fine, l'intensità attesa è

λ ∝ (1 − w) [G(σ_e) ∗ (T π) + G(σ_in) ∗ (T (1 − π))] + w E(ℓ) ∗ T, con π(d) = exp(−d/D),

dove G(σ) è la gaussiana isotropa in due dimensioni con deviazione standard σ per asse, E(ℓ) il nucleo E(r) ∝ exp(−r/ℓ) normalizzato nel piano (K2 di R5), ∗ la convoluzione e d la distanza della sorgente dal bordo del tessuto. Per un singolo trascritto: con probabilità w si sposta secondo E(ℓ); altrimenti secondo G(σ_e) con probabilità π(d), secondo G(σ_in) con probabilità 1 − π(d). La regola è misurata sul bordo esterno del tessuto; nei lumi la si applica con la distanza dal bordo del lume, e CP-1 dice che funziona entro 8 µm. Valori per sezione nella sezione 9; il cuore usa la mediana delle altre quattro.

## 7. Passo 5: cattura e strati tecnici (revisione 3.2.C)

1. **Strisce** (D2-6). Ogni bin ha un fattore di cattura c_ij = p(parità di i, parità di j) · exp(α_i + β_j), poi normalizzato a media 1 sui bin di tessuto dell'area simulata. La parte periodica p è separabile, p = p_r(i mod 2) · p_c(j mod 2), con i due rapporti dispari/pari della sezione (R4 misura le quattro classi: se non fossero separabili si usano le quattro classi misurate); α_i e β_j sono fattori casuali per riga e per colonna. La media di exp(α + β) prima della normalizzazione è exp((σ_α² + σ_β²)/2), 1.04-1.06 con i valori della sezione 9: per questo si normalizza. S2.0 (M1) ha misurato le due parti nelle cinque sezioni. La parte periodica (rapporto dispari/pari) va da 0.95 (righe del cuore) a 1.75 (colonne dell'intestino) e, salvo nel linfonodo, è la stessa in tutta la sezione (deviazione standard fra tessere di 300 µm ≤ 0.021). La parte non periodica non è trascurabile: deviazione standard 0.18-0.28 in scala logaritmica (P-1, come pre-registrata), riproducibile fra le due metà di ogni riga e colonna. La revisione avversariale ha mostrato che la stima alternata dei fattori non converge nella parte lenta: una tendenza comune a molte righe passa indifferentemente nei fattori o nel campo del tessuto lisciato a 10 µm, e deriva con il numero di iterazioni. P-1 resta smentita in ogni variante, ma per il simulatore servono quantità ben definite. Per questo si usa la sola parte fine (post hoc PH-5): lo scarto di ogni fattore dalla media mobile su 21 righe o colonne, stabile fra 10 e 20 iterazioni entro 0.003. La parte fine ha deviazione standard 0.18-0.25, correlazione fra le due metà 0.85-0.97 (0.39 nelle colonne del linfonodo) e passo 1 da −0.22 a −0.46. Su sezioni sintetiche con fattori noti lo stesso procedimento ritrova la deviazione standard e il passo 1 della parte fine vera, e il filtro sposta il passo 1 di circa −0.05 per fattori indipendenti; con questa calibrazione il passo 1 da dare al generatore va da −0.17 (cuore, righe) a −0.46 (linfonodo, colonne), in parte estrapolato oltre l'ultimo punto calibrato (−0.3). Le righe e le colonne con tessuto ma zero UMI escono dalla maschera prima della stima: il loro fattore sentinella abbassava la parità del linfonodo (righe 0.980 → 0.996, colonne 1.125 → 1.105) e del cervello (righe 1.178 → 1.201). I fattori sono correlati con quelli del vetro lontano dal tessuto (post hoc PH-1, correlazione 0.30-0.83): sono almeno in parte cattura tecnica per riga e colonna. Il simulatore genera α_i e β_j gaussiani, centrati, con la deviazione standard della parte fine e il passo 1 calibrato ρ, come media mobile di ordine 1: x_i = (e_i + θ e_{i−1}) / √(1 + θ²), con e_i normali indipendenti e θ = (1 − √(1 − 4ρ²)) / (2ρ). Le tendenze lente della cattura lungo righe e colonne non sono identificate e il simulatore non le riproduce (limite). Nel linfonodo la parità cambia nello spazio (deviazione standard fra tessere del logaritmo del rapporto dispari/pari 0.09 nelle righe, 0.05 nelle colonne; nelle colonne P-2 è indeterminata ma l'intervallo esclude 0): il preset la fa variare con un campo gaussiano liscio del logaritmo del rapporto, lunghezza di correlazione 300 µm (la scala delle tessere della misura), con quelle deviazioni standard. Valori nella sezione 9. Lo strato si spegne con un argomento (`stripes = FALSE`): i metodi che correggono le strisce si provano con lo strato acceso, gli altri si possono provare anche senza.
2. **Fondo sul vetro** (D2-10): b trascritti per bin in media (Poisson), con la composizione per gene media del tessuto; sul tessuto il fondo si somma ai trascritti delle cellule.
3. **Campionamento**: i trascritti spostati e il fondo si contano nel bin in cui cadono; i fattori delle strisce agiscono come assottigliamento binomiale (si generano N_c · max(c) / media(c) trascritti e si tiene ciascuno con probabilità c_bin / max(c)), così i conteggi restano interi e la verità di origine resta esatta. Il massimo di c su un'area intera, con 3 350 fattori per asse, può valere diverse volte la media: il fattore di sovra-generazione max(c)/media(c) si calcola per tessera (non sull'area intera) ed entra nel bilancio delle prestazioni (sezione 12).

## 8. Passo 6: somma a 8 µm e uscite (revisione 3.2.E)

`counts_8um` somma i 16 bin di 2 µm di ogni bin di 8 µm (righe e colonne a blocchi di 4, come Space Ranger). La verità di origine e quella geometrica si sommano allo stesso modo. Controllo: ogni bin di 8 µm è uguale alla somma esatta dei suoi 16 (identità su interi; R5 l'ha verificata sui dati veri).

## Tabella dei parametri (sezione 9)

Valori per tessuto; la colonna «fonte» rimanda alla riga di `docs/BIO_REFERENCES.md`, alla decisione o al file del repository che li ha prodotti.

| parametro | intestino (A1) | tumore (A2) | stroma (A3) | linfonodo (A4) | cervello (A5) | cuore (A6) | fonte |
|---|---|---|---|---|---|---|---|
| trascritti per cellula | per tipo, mediana e dispersione di R4 | per tipo, mediana e dispersione di R4 | per tipo, mediana e dispersione di R4 | per tipo, mediana e dispersione di R4 | per tipo, mediana e dispersione di R4 | per tipo, mediana e dispersione di R4 | per tessuto B-275, B-285, B-295, B-307, B-320, B-332; per tipo B-276...B-278, B-286...B-288, B-296...B-300, B-308...B-313, B-321...B-325, B-333...B-335; [DEC-035](/decisions/DEC-035.md) |
| rapporto nucleo/citoplasma r (regola C2) | 0.75 | 0.75 | 0.82 | 0.55 | 1.48 | 2.61* | B-384, B-386, B-388, B-390, B-392, B-394; *controllo interno smentito |
| sfocatura stretta per tutti σ_in (µm) | 0.7 | 0.7 | 0.7 | 0.7 | 0.7 | 0.7 | valore unico per scelta ([DEC-038](/decisions/DEC-038.md)): R5 PH-5 `results/R5/R5_PH5_narrowcore.csv` (sigma_fixed 0.7) e pre-registrazione di M2, dentro 0.5-1 µm di [DEC-037](/decisions/DEC-037.md); valori per tessuto in B-367...B-372 |
| zona di bordo: sfocatura larga σ_e (µm) | 6.23 | 4.72 | 4.72 | 4.49 | 4.99 | 4.86† | B-379, B-380, B-381, B-382, B-383; †mediana delle altre 4 sezioni ([DEC-037](/decisions/DEC-037.md)), calcolata da `results/S2.0/wfix/S2.0_edges_sections.csv` |
| zona di bordo: profondità D (µm) | 24.75 | 8.94 | 8.94 | 22.04 | 6.86 | 15.49† | B-379, B-380, B-381, B-382, B-383; †mediana delle altre 4 sezioni ([DEC-037](/decisions/DEC-037.md)), calcolata da `results/S2.0/wfix/S2.0_edges_sections.csv` |
| coda esponenziale: massa w | 0.087 | 0.040 | 0.040 | 0.109 | 0.112 | 0.098† | B-379, B-380, B-381, B-382, B-383; †mediana delle altre 4 sezioni ([DEC-037](/decisions/DEC-037.md)), calcolata da `results/S2.0/wfix/S2.0_edges_sections.csv` |
| coda esponenziale: lunghezza ℓ (µm) | 32.85 | 24.73 | 24.73 | 41.50 | 33.07 | 32.96† | B-379, B-380, B-381, B-382, B-383; †mediana delle altre 4 sezioni ([DEC-037](/decisions/DEC-037.md)), calcolata da `results/S2.0/wfix/S2.0_edges_sections.csv` |
| fondo sul vetro b (UMI per bin di 2 µm) | 0.161 | 0.197 | 0.197 | 0.042 | 0.309 | 0.117 | B-342, B-348, B-354, B-360, B-366 |
| strisce: parità righe dispari/pari | 1.417 | 1.578 | 1.578 | 0.996‡ | 1.201 | 0.947 | B-397, B-398, B-399, B-400, B-401; ‡variabile nello spazio; righe e colonne con 0 UMI tolte (PH-5) |
| strisce: parità colonne dispari/pari | 1.750 | 1.724 | 1.724 | 1.105‡ | 1.231 | 0.967 | B-397, B-398, B-399, B-400, B-401 |
| strisce: SD della parte fine dei fattori per riga α (log) | 0.193 | 0.221 | 0.221 | 0.211 | 0.242 | 0.183 | B-397, B-398, B-399, B-400, B-401; parte lenta non identificata (PH-5) |
| strisce: SD della parte fine dei fattori per colonna β (log) | 0.185 | 0.191 | 0.191 | 0.251 | 0.180 | 0.188 | B-397, B-398, B-399, B-400, B-401 |
| strisce: passo 1 per il generatore, α; β (dalla parte fine misurata, calibrata su sezioni sintetiche) | -0.30; -0.27§ | -0.31; -0.21§ | -0.31; -0.21§ | -0.38; -0.46§ | -0.31; -0.29§ | -0.17; -0.21 | B-397, B-398, B-399, B-400, B-401; `results/S2.0/S2.0_PH5_calib.csv`; §almeno un valore oltre l'ultimo punto calibrato (passo 1 vero −0.3) |
| variabilità dentro il nucleo v, come scritta da R-C (prototipo senza fattori di riga e colonna) | 0.102 | 0.144 | 0.170 | 0.159 | 0.115 | 0.201 | B-385, B-387, B-389, B-391, B-393, B-395; `results/S2.0/S2.0_RC.json` |
| variabilità dentro il nucleo v, al netto delle strisce (post hoc PH-2; punto di partenza della taratura, DEC-040) | 0.044 | 0.069 | 0.095 | 0.063 | 0.033 | 0.153* | B-385, B-387, B-389, B-391, B-393, B-395; post hoc PH-2 |
| scala della variabilità dentro il nucleo | 1-2 µm | 1-2 µm | 1-2 µm | 1-2 µm | 1-2 µm | 1-2 µm | B-385, B-387, B-389, B-391, B-393, B-395; post hoc PH-3 (correlazione in più a 2 µm, al più 0.09 a 4 µm) |
| riferimento (controllo, non parametro): UMI per bin di 2 µm in tessuto, media | 27.14 | 45.46 | 19.53 | 11.40 | 41.79 | 3.42 | B-270, B-280, B-290, B-302, B-315, B-327 |
| riferimento (controllo, non parametro): UMI per bin di 8 µm in tessuto, media | 424.81 | 719.50 | 312.39 | 182.42 | 667.44 | 54.78 | B-271, B-281, B-291, B-303, B-316, B-328 |
| riferimento (controllo, non parametro): sparsità bin × gene a 2 µm; 8 µm | 0.9987; 0.9830 | 0.9976; 0.9679 | 0.9990; 0.9857 | 0.9994; 0.9910 | 0.9979; 0.9716 | 0.9998; 0.9975 | B-272, B-282, B-292, B-304, B-317, B-329 |

Il tumore (A2) e lo stroma (A3) stanno sulla stessa sezione: spostamento, fondo e strisce sono quelli della sezione; r e v sono per tessuto. Le ultime tre righe non sono parametri: sono i valori veri che il simulatore deve riprodurre (passo 6 e S2.5).

## 10. Lumi, vetro e spazio acellulare (BL-077, BL-086, BL-201)

Nello Step 1 ogni punto del tessuto appartiene a un territorio cellulare (decisione 1 del design dello Step 1): non c'è matrice senza cellule dentro le regioni (BL-077). Lo Step 2 non la aggiunge: lumi e vetro ricevono solo i trascritti spostati dalle cellule vicine e il fondo. Entro 8 µm dal bordo dei lumi questo riproduce i lumi veri entro un quarto (R5 CP-2: 9 combinazioni su 10; controprova di S2.0 con la regola scelta: con la zona di bordo nessuna delle 8 combinazioni sezione × fascia smentita, 6 confermate e 2 indeterminate). Oltre 8 µm i lumi veri del tumore del colon hanno 1.4-1.5 volte i trascritti previsti con la zona di bordo (materiale nelle ghiandole, BL-201), quelli dell'intestino e del linfonodo 1.1-1.25 volte a 8-32 µm, quelli del cuore meno (0.49-0.74): limite dichiarato; un'opzione «materiale nei lumi» (sorgente debole uniforme nei lumi) resta per quando avremo lumi annotati. I conteggi a 2 µm non distinguono un vuoto piccolo dentro il tessuto (BL-086): i metodi che cercano lo spazio acellulare vanno provati con vuoti grandi.

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

Ordini di grandezza dalle righe di `docs/BIO_REFERENCES.md` (da misurare negli step d'implementazione). Un'area di 1 mm² ha 250 000 bin di 2 µm. Le cellule per mm² vanno da circa 1 000 (cuore, B-036 e B-048; cervello, B-029 e B-047) a 24 000-28 000 (linfonodo, B-022 e B-046), con 11 000-12 000 nell'intestino (B-001, B-043). I trascritti catturati per mm² (UMI per bin di 2 µm × 250 000) vanno da 0.86 milioni (cuore, B-327) a 11.4 milioni (tumore del colon, B-280): 10.4 milioni nel cervello (B-315), 6.8 nell'intestino (B-270), 2.85 nel linfonodo (B-302). Con l'assottigliamento delle strisce i trascritti generati sono max(c)/media(c) volte di più (sezione 7). Un'area di cattura intera (11.2 milioni di bin) è 45 volte più grande.

Bilancio proposto, come per lo Step 1: un'area di 1 mm² in meno di 5 minuti e 8 GB; l'area intera a tessere di 500 µm con un bordo di 280 µm (meno dell'1 % della coda più lunga oltre il bordo, sezione 3), un processo per tessera. Lo Step 1 sul linfonodo 4 000 × 4 000 px supera già il bilancio (310 s, 8.5 GB, BL-144): il ritaglio dei territori (BL-072) e le regioni a 2 µm (BL-045) vanno ottimizzati prima dell'area intera; la griglia dei vicini su estensioni enormi (BL-057) non serve allo Step 2, che lavora a tessere.

## 13. Tolleranze dei controlli (BL-146)

Scelte tecniche registrate in [DEC-038](/decisions/DEC-038.md).

- Conteggi: identità esatte su interi. Trascritti generati = catturati + persi per assottigliamento + usciti dall'area simulata (quelli che escono da una tessera si contano nella tessera vicina); somma a 8 µm uguale alla somma dei 16 bin di 2 µm. `truth$cell_counts` contiene i trascritti generati per cellula, prima dello spostamento e dell'assottigliamento: il totale atteso catturato di una cellula è N_c (scala di R4, D2-3), quello generato N_c · max(c)/media(c) (sezione 7).
- Frazioni geometriche: l'area coperta in un bin (somma delle aree dei pezzi) si confronta con l'area attesa entro g · Σ P + n · g², con g = 1e-6 µm, Σ P la somma dei perimetri dei pezzi e n il loro numero (la regola di S1.5 per due pezzi, τ = g · (P_i + P_j) + g², BL-093, estesa a n pezzi); per le frazioni la tolleranza si divide per l'area del bin (4 µm²). La somma per regione si lega alla griglia, non all'area della regione (risposta a BL-146). Su un'area di cattura intera le coordinate assolute arrivano a 1e4 µm e l'errore relativo dell'area dei poligoni a 1.4e-7 (BL-104): le sovrapposizioni si calcolano in coordinate locali della tessera.
- Distribuzioni simulate contro riferimenti: errore Monte Carlo fra semi e intervalli di rilevanza pre-registrati negli step.

## 14. Che cosa si riusa della vecchia pipeline (BL-001, BL-005, BL-010)

| file | che cosa fa | nello Step 2 |
|---|---|---|
| `R/05_grid_sampling.R` | griglia di campionamento con scala propria | sostituito dal passo 1 (BL-005, BL-010) |
| `R/06a_expression_params.R`, `R/06k_expression_profiles_wrapper.R` | parametri e involucro della generazione per spot | sostituiti dall'orchestratore dello Step 2 |
| `R/06b_expression_baseline.R` | profili sintetici per tipo con geni marcatori | riusabile per i controlli (profili noti), non come fonte dei preset |
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
| S2.2 | `simulate_cell_profiles()` | totali per cellula, multinomiale, sottotipi | trascritti per cellula per tipo (R4); profili scelti dalla regola pre-registrata di S2.2 (DEC-040) | profili veri recuperati dalla media delle cellule simulate |
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
| revisione 4.4 (trascritti per bin a 8 µm) | sezioni 4, 9 (righe dei riferimenti a 2 e 8 µm) |
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
| BL-104 aree in coordinate assolute | sezione 13 |
| BL-167 densità dei preset | sezione 17 |
| [DEC-034](/decisions/DEC-034.md), [DEC-035](/decisions/DEC-035.md), [DEC-037](/decisions/DEC-037.md) | sezione 2 |
| [DEC-038](/decisions/DEC-038.md) scelte tecniche di S2.0 (σ_in unico, margini, tolleranze) | sezioni 9, 13 |
| [DEC-039](/decisions/DEC-039.md) regole decise dai dati in S2.0 (R-A, R-B, R-C) | sezioni 5, 6, 7, 9 |
| [DEC-040](/decisions/DEC-040.md) profili, profondità e variabilità a 1-2 µm decisi con i dati | sezioni 2, 4, 5, 9 |

## 17. Limiti dichiarati

- Spazio acellulare dentro il tessuto assente (sezione 10); lumi oltre 8 µm del tumore del colon e del cuore non riprodotti (BL-201).
- Arricchimento nucleare per gene assente (revisione 4.2).
- Strisce: solo la parte fine dei fattori di riga e colonna (periodi sotto circa 40 µm); le tendenze lente della cattura non sono identificate dallo stimatore e mancano.
- Variabilità a 1-2 µm misurata (pre-registrata) solo dentro il nucleo, applicata a tutto il tessuto sulla base di una verifica provvisoria della revisione; origine (biologica o tecnica) non nota (sezione 5).
- Cuore: il modello non descrive l'interno dei nuclei (controllo interno smentito); r e v del cuore con riserva.
- Fondo vicino al bordo: il simulatore usa il fondo del vetro lontano (sezione 9), ma nelle finestre di bordo il fondo adattato è più alto in intestino, tumore e cervello (0.22, 0.46, 0.55 UMI per bin contro 0.16, 0.20, 0.31): vicino al tessuto il vetro simulato avrà meno conteggi del vero (forse una coda non modellata).
- Bilancio dei trascritti per area: densità delle cellule (Step 1) per trascritti per cellula (R4, nucleo più corona di 8 µm) può dare meno UMI per area dei dati veri, fino a un terzo nel cuore, dove la corona di 8 µm è piccola per i cardiomiociti (stima provvisoria della revisione G-5); da controllare come previsione in S2.2 e S2.5 con BL-195 e BL-167 (densità dei preset riferite al tessuto pieno).
- Revisione 4.1: la metrica di Daniele (correlazione fra i profili per gene dei bin della stessa cellula, contro cellule diverse dello stesso tipo) non è misurata in S2.0, che misura la dispersione dei totali per bin; il campo moltiplicativo della sezione 5 non crea variazione di composizione. Si misura in S2.2 con i profili scelti.
- Revisione 3.2.B e 3.2.D: gli effetti spaziali dentro un tipo (campi gaussiani, `R/06h_spatial_correlation.R`) sono solo un'opzione per i programmi; la localizzazione periferica o di membrana non è trattata, salvo il gradiente radiale (BL-087).
- Cuore: spostamento con i valori delle altre sezioni ([DEC-037](/decisions/DEC-037.md)); cellule isotrope (limite dello Step 1).
- I limiti dello Step 1 passano ai conteggi: contatti fra nuclei in eccesso nei tessuti densi, territori troppo uniformi con la semina uniforme (design dello Step 1, 7.4).

## 18. Storico

| data | versione | che cosa | chi |
|---|---|---|---|
| 2026-11-04 | 0.1 | prima bozza con le misure di S2.0 | geo-agent/skill-v2 (S2.0) |
| 2026-11-04 | 0.2 | correzioni dalla revisione avversariale G-5 (documento, nucleo, strisce, bordi), PH-4, PH-5, PH-6, DEC-038, DEC-039, DEC-040 | geo-agent/skill-v2 (S2.0) |