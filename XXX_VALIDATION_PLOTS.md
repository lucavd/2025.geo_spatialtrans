# Validazione Visiva della Coerenza dei Dati Spaziali Simulati (Senza Dati Reali)

Questo documento descrive quattro tipi di plot chiave che possono essere generati (ad esempio in R) per valutare visivamente la coerenza interna e la verosimiglianza dei dati di trascrittomica spaziale simulati, senza la necessità di confrontarli direttamente con un dataset sperimentale reale.

## 1. Plot Spaziale dell'Espressione dei Geni Marker

* **Scopo:** Verificare se i geni designati come specifici per certi tipi cellulari/regioni ("marker") durante la simulazione sono effettivamente espressi nelle corrette posizioni spaziali.
* **Come si presenta il plot:** È una mappa che mostra la posizione (coordinate x, y) di tutti gli spot/punti simulati. Ogni punto è colorato in base al livello di espressione (es. conteggi normalizzati e log-trasformati) di *un singolo gene marker*. È necessario creare un grafico separato per ogni marker che si desidera ispezionare. Una scala di colori consistente (es. `viridis` o `magma` in R, dal blu/basso al giallo/alto) è essenziale.
* **Cosa Osservare (Segno di Coerenza):** Per un gene simulato come marcatore del "Tipo A", i punti nelle regioni spaziali definite come "Tipo A" (note dalla ground truth della simulazione) dovrebbero avere colori che indicano alta espressione. Le altre regioni dovrebbero mostrare bassa espressione. Il pattern visivo dei colori deve corrispondere alla mappa spaziale dei tipi cellulari simulati. Una chiara localizzazione è un buon segno; una distribuzione casuale o in regioni errate indica un problema.

## 2. Plot della Relazione Media-Varianza dei Geni

* **Scopo:** Controllare se la variabilità statistica dei conteggi genici simulati è realistica, mostrando la caratteristica sovradispersione (varianza > media) tipica dei dati RNA-seq e della modellazione Negative Binomial.
* **Come si presenta il plot:** Un grafico a dispersione (scatter plot) dove ogni punto rappresenta un gene. L'asse X mostra l'espressione media del gene (su tutti gli spot), l'asse Y mostra la varianza dell'espressione dello stesso gene. Entrambi gli assi sono tipicamente in scala logaritmica. Viene spesso aggiunta una linea di riferimento diagonale che indica Varianza = Media (comportamento Poissoniano).
* **Cosa Osservare (Segno di Coerenza):** La maggior parte dei punti (geni) deve trovarsi significativamente *al di sopra* della linea Varianza = Media. Questo conferma la sovradispersione attesa. Il grafico dovrebbe mostrare un trend generale in cui la varianza aumenta all'aumentare della media. Se i punti si raggruppano attorno o sotto la linea di riferimento, la simulazione della varianza non è biologicamente plausibile per dati trascrittomici.

## 3. Plot della Relazione Dropout-Espressione Media

* **Scopo:** Verificare se l'artefatto tecnico del dropout (mancata rilevazione di un gene, risultante in un conteggio zero) è simulato in modo plausibile, essendo più probabile per geni meno espressi.
* **Come si presenta il plot:** Un grafico a dispersione dove ogni punto è un gene. L'asse X mostra l'espressione media del gene (tipicamente log-trasformata, es. log(media + 1)). L'asse Y mostra la frazione (o percentuale) di spot in cui quel gene ha un conteggio pari a zero. Data la potenziale sovrapposizione di molti punti, è utile usare trasparenza (alpha) o una visualizzazione basata sulla densità dei punti.
* **Cosa Osservare (Segno di Coerenza):** Si deve osservare una chiara tendenza *negativa*. I geni molto espressi (a destra sull'asse X) devono avere una frazione di zeri vicina a zero (in basso sull'asse Y). I geni poco espressi (a sinistra sull'asse X) devono mostrare una frazione di zeri più alta e più variabile. La forma complessiva di questa relazione è un indicatore chiave della plausibilità del modello di dropout.

## 4. Plot UMAP/t-SNE Colorato per Tipo Cellulare Simulata

* **Scopo:** Valutare se i profili di espressione globali generati per i diversi tipi cellulari/cluster *simulati* sono abbastanza distinti da poter essere separati da metodi di analisi standard. È un controllo funzionale della simulazione.
* **Come si presenta il plot:** Si applica un algoritmo di riduzione dimensionale (es. UMAP o t-SNE) alla matrice di espressione (spot x geni) per ottenere coordinate 2D per ogni spot. Si crea un grafico a dispersione di queste coordinate (es. UMAP1 vs UMAP2). Fondamentale: **i punti vengono colorati in base all'etichetta del tipo cellulare/cluster assegnata durante la simulazione (la ground truth)**.
* **Cosa Osservare (Segno di Coerenza):** I punti dello stesso colore (stesso tipo simulato) dovrebbero tendere a raggrupparsi nello spazio UMAP/t-SNE. La nettezza della separazione tra i gruppi di colori diversi dipenderà dalla "difficoltà" della simulazione (es. quanto erano distinguibili i geni marker). Anche se i cluster non sono perfettamente isolati, una chiara tendenza al raggruppamento per tipo simulato indica che il framework sta generando profili di espressione coerenti con le etichette definite. Se i colori appaiono completamente mescolati, la simulazione potrebbe non star creando differenze significative tra i tipi.