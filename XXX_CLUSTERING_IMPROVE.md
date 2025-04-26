# Proposta di Miglioramento per il Clustering Iniziale nel Framework di Simulazione Spaziale

## Introduzione

Questo documento riassume le proposte di miglioramento discusse per la fase iniziale di clustering basata su immagine nel framework di simulazione spaziale `2025.geo_spatialtrans`. L'obiettivo è aumentare il realismo della segmentazione spaziale iniziale, affrontando le limitazioni del K-means++ basato unicamente sull'intensità, pur mantenendo un buon compromesso con la velocità computazionale e riconoscendo che questo è solo il primo passo di una simulazione più complessa.

## 1. Opzioni Aggiuntive per l'Algoritmo di Clustering

Si propone di rendere selezionabile l'algoritmo di clustering iniziale, offrendo alternative che incorporano informazioni spaziali.

* **Metodo Attuale (Default): `kmeans++`**
    * Utilizza K-means++ basandosi unicamente sulla similarità di intensità dei pixel in scala di grigi.
    * *Limitazione:* Può raggruppare pixel spazialmente distanti se hanno intensità simili.

* **Opzione Proposta 1: `spatial_kmeans` (K-means Spazialmente Vincolato)**
    * Utilizza K-means modificato con una metrica di distanza combinata:
        $$D_{combinata}^2 = D_{intensità}^2 + w^2 \cdot D_{spaziale}^2$$
    * `w` è un parametro di peso che bilancia l'importanza della vicinanza spaziale rispetto alla similarità di intensità.
    * *Beneficio:* Produce cluster più compatti e spazialmente coesi.
    * *Velocità:* Moderato aumento del tempo di calcolo rispetto al K-means standard.

* **Opzione Proposta 2: `slic` (Simple Linear Iterative Clustering)**
    * Algoritmo che genera "superpixel" (piccole regioni compatte e omogenee).
    * Utilizza una distanza combinata (intensità + spazio) ma limita la ricerca dei cluster a regioni locali.
    * *Beneficio:* Garantisce intrinsecamente la coesione spaziale ed è efficiente.
    * *Velocità:* Generalmente veloce, comparabile o superiore a K-means su immagini grandi.

* **Implementazione Suggerita:**
    * Introdurre un nuovo parametro di configurazione, ad esempio `clustering_method`, con valori possibili `"kmeans++"` (default), `"spatial_kmeans"`, `"slic"`.
    * Modificare/estendere la funzione `R/functions/04_clustering.R` per includere queste opzioni.

## 2. Opzione per la Determinazione Automatica del Numero di Cluster (`k`)

Si propone di aggiungere una funzionalità *opzionale* per stimare il numero di cluster (`k`), alternativa all'impostazione manuale.

* **Metodo Attuale:** Richiede l'impostazione manuale del numero di cluster tramite il parametro `k_cell_types`. Questo rimane l'approccio **principale e raccomandato** per scenari di simulazione controllata e benchmarking, dove la complessità deve essere nota a priori.

* **Opzione Proposta: Stima Automatica di `k`**
    * Introdurre un parametro booleano, ad esempio `estimate_k = FALSE` (default), che se impostato a `TRUE` attiva la stima automatica.
    * *Metodi di Stima:* Utilizzare tecniche standard come:
        * Metodo Elbow (basato su WCSS con metrica appropriata)
        * Punteggio Silhouette (basato su metrica appropriata)
        * Gap Statistic
    * *Workflow:* Quando `estimate_k = TRUE`, il framework esegue la stima usando il metodo di clustering selezionato e può:
        1.  Utilizzare direttamente il `k` stimato per il clustering finale.
        2.  Suggerire il `k` stimato all'utente, richiedendo comunque una conferma o permettendo una sovrascrittura manuale.
    * *Beneficio:* Offre una guida basata sui dati per scegliere `k` in scenari esplorativi, pur mantenendo il controllo manuale come standard per la riproducibilità.

## Conclusione

L'introduzione di queste opzioni aumenterebbe la flessibilità del framework, permettendo agli utenti di generare basi spaziali per le simulazioni con un maggiore grado di coerenza spaziale e offrendo assistenza nella scelta della complessità dello scenario (`k`), adattandosi meglio a diverse esigenze e preferenze.