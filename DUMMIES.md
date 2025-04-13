# Spiegazione semplice del progetto di simulazione di trascrittomica spaziale

## Cos'è questo progetto?

Questo progetto crea dati artificiali che imitano come i geni si attivano in cellule di diverso tipo all'interno di un tessuto. È come creare una mappa dettagliata che mostra quali geni sono "accesi" o "spenti" in ogni singola cellula di un tessuto, tenendo conto della loro posizione.

## Perché è utile?

1. I ricercatori possono testare nuovi metodi di analisi su questi dati simulati
2. Si possono creare situazioni difficili da trovare nei dati reali
3. Si può risparmiare tempo e denaro rispetto a esperimenti di laboratorio costosi

## Come funziona in pochi passi:

1. **Carica un'immagine**: Usiamo un'immagine (come una foto di tessuto) come punto di partenza

2. **Crea una griglia**: Posizioniamo una griglia sopra l'immagine per rappresentare dove si troverebbero le cellule

3. **Identifica i tipi di cellule**: Raggruppiamo le cellule in diversi "tipi" in base al loro colore/intensità nell'immagine

4. **Genera profili di espressione genica**: Per ogni tipo di cellula, creiamo un modello di quali geni sono attivi

5. **Aggiungi effetti biologici realistici**:
   - Alcune cellule vicine si influenzano a vicenda
   - Alcuni geni tendono ad attivarsi insieme
   - A volte l'espressione di un gene non viene rilevata (dropout)

6. **Crea risultati**: Produciamo dati che sembrano quelli che otterresti da un vero esperimento di laboratorio

7. **Genera report**: Alla fine, creiamo un rapporto con grafici e statistiche su ciò che è stato simulato

## Un esempio concreto:

Immagina che vogliamo simulare un tessuto con tre tipi di cellule (come un tessuto con cellule immunitarie, strutturali e nervose):

1. Usiamo un'immagine che rappresenta il tessuto
2. Identifichiamo le diverse aree e le etichettiamo come "Tipo 1", "Tipo 2" e "Tipo 3"
3. Per il Tipo 1 diciamo che i geni A, B e C sono molto attivi
4. Per il Tipo 2 diciamo che i geni D, E e F sono molto attivi
5. Per il Tipo 3 diciamo che i geni G, H e I sono molto attivi
6. Aggiungiamo variabilità in modo che non tutte le cellule dello stesso tipo siano identiche
7. Creiamo una tabella finale dove ogni riga è una cellula (con la sua posizione) e ogni colonna è un gene, con numeri che indicano quanto è attivo quel gene in quella cellula

Il risultato è come una grande tabella di Excel dove puoi vedere l'attività di ogni gene in ogni cellula, organizzata in base alla posizione delle cellule nel tessuto.

## Cosa rende speciale questo progetto?

1. **È modulare**: Ogni parte del processo è un pezzo separato che può essere migliorato indipendentemente
2. **È personalizzabile**: Puoi cambiare facilmente parametri come il numero di cellule, geni o tipi cellulari
3. **È riproducibile**: Usando gli stessi parametri otterrai sempre gli stessi risultati
4. **È automatizzato**: Una volta configurato, il processo funziona da solo dall'inizio alla fine

In sintesi, questo progetto è come un laboratorio virtuale per creare dati di espressione genica che sembrano reali ma che puoi controllare completamente.

# Dal tessuto all'espressione genica: il processo dettagliato

## Come trasformiamo un'immagine in espressione genica

Quando partiamo da un'immagine reale di un tessuto, seguiamo questi passaggi:

### 1. Preparazione dell'immagine

- **Caricamento**: L'immagine viene caricata in formato digitale (PNG, JPEG)
- **Conversione**: Se a colori, viene convertita in scala di grigi
- **Normalizzazione**: I valori dell'intensità vengono standardizzati tra 0 e 1
- **Soglia**: Applichiamo una soglia per distinguere il tessuto dallo sfondo vuoto

### 2. Identificazione dei tipi cellulari tramite clustering

- L'algoritmo analizza i diversi livelli di grigio nell'immagine
- Raggruppa pixel simili in "cluster" che rappresentano diversi tipi di cellule
- Ogni pixel viene etichettato con un numero (1, 2, 3, ecc.) che indica il tipo cellulare

### 3. Creazione della griglia di campionamento (sampling)

Qui entra in gioco il "sampling" che hai notato. Ci sono due approcci principali:

#### Approccio "grid_mode" (simulazione di Visium HD)
- **Cosa fa**: Crea una griglia regolare, come una scacchiera, sull'immagine
- **Perché**: Imita le tecnologie reali come Visium HD, dove i sensori sono posizionati in punti fissi
- **Come funziona**: 
  - La griglia ha una risoluzione specifica (es. 2μm tra un punto e l'altro)
  - Ogni punto della griglia diventa una "cella" nella nostra simulazione
  - Il tipo cellulare di quel punto viene determinato dal cluster sottostante

#### Approccio casuale (senza grid_mode)
- **Cosa fa**: Sceglie punti casuali all'interno dell'area del tessuto
- **Perché**: Simula tecnologie che possono campionare cellule in posizioni casuali
- **Come funziona**: 
  - Si selezionano casualmente N punti (dove N = n_cells) all'interno dell'area del tessuto
  - Per ogni punto, si determina il tipo cellulare dal pixel sottostante

### 4. Relazione tra immagine e griglia

L'immagine del tessuto e la griglia sono allineate così:
- Ogni pixel dell'immagine ha coordinate (x, y)
- Ogni punto della griglia ha coordinate (x, y)
- Il tipo cellulare di un punto della griglia è determinato dal cluster del pixel in quella posizione
- Se un punto della griglia cade fuori dal tessuto, viene ignorato

### 5. Generazione dell'espressione genica

Per ogni punto della griglia (che ora rappresenta una cellula):
- **Profilo di base**: Determiniamo un profilo di espressione base in base al tipo cellulare
- **Geni marcatori**: Ogni tipo cellulare ha geni specifici che sono molto espressi
- **Variabilità**: Aggiungiamo variazioni casuali per riflettere le differenze tra cellule dello stesso tipo
- **Effetti spaziali**: Le cellule vicine tendono ad avere profili di espressione più simili

### 6. Dalla griglia alla matrice di espressione

La matrice finale di espressione rappresenta:
- **Righe**: Ogni riga è una "cella" (un punto della griglia)
- **Colonne**: Ogni colonna è un gene
- **Valori**: Quanto è espresso quel gene in quella cella (numeri interi positivi)
- **Coordinate spaziali**: Manteniamo le coordinate (x, y) di ogni cella

## I geni corrispondono a box sulla griglia?

No, i geni non corrispondono ai box della griglia. La relazione è diversa:

- **Box della griglia** = cellule/spot di tessuto (la posizione fisica)
- **Geni** = caratteristiche misurate per ogni box (sono le colonne della matrice finale)

Per ogni punto della griglia (cellula), misuriamo l'espressione di tutti i geni. È come se stessimo compilando un questionario per ogni cellula, dove le domande sono "Quanto è espresso il gene A?", "Quanto è espresso il gene B?", e così via.

## Esempio pratico

Immagina un'immagine di un tessuto dove:
- Le aree scure sono cellule di tipo 1
- Le aree grigie sono cellule di tipo 2
- Le aree chiare sono cellule di tipo 3

Quando mettiamo la griglia sopra:
1. Un punto che cade su un'area scura diventa una "cellula di tipo 1"
2. Per questa cellula, generiamo valori di espressione per tutti i geni
3. I geni marcatori del tipo 1 avranno valori più alti
4. Se due punti della griglia sono vicini, i loro profili di espressione saranno simili

Alla fine, otteniamo una tabella dove:
- Ogni riga è un punto della griglia (una cellula)
- Le prime colonne indicano le coordinate x, y e il tipo cellulare
- Le restanti colonne mostrano il livello di espressione di ogni gene

Questa tabella, insieme alle informazioni spaziali, è il prodotto finale della simulazione, pronto per essere analizzato con gli stessi strumenti usati per i dati reali di trascrittomica spaziale.