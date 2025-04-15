# Documentazione: Miglioramento dei Moduli di Co-espressione Genica

Questa documentazione descrive l'implementazione migliorata dei moduli di co-espressione genica nel pacchetto 2025.geo_spatialtrans, che permette una simulazione più realistica delle reti geniche e delle interazioni tra geni.

## Funzionalità introdotte

### 1. Moduli gerarchici
- Possibilità di creare strutture gerarchiche di moduli (moduli principali e sotto-moduli)
- Riflette la struttura reale dei pathway biologici, dove i geni sono organizzati in gruppi funzionali annidati

### 2. Interazioni tra moduli (crosstalk)
- Matrice di rete che modella le interazioni tra i diversi moduli
- Permette la propagazione degli effetti tra moduli correlati
- Densità di rete personalizzabile

### 3. Distribuzione realistica delle dimensioni dei moduli
- Distribuzione esponenziale: modelli più piccoli e numerosi rispetto a pochi moduli grandi
- Simile alla distribuzione osservata nei dati reali di espressione genica
- Possibilità di utilizzare anche distribuzioni uniformi

### 4. Fattori latenti
- Modella componenti latenti (non osservate) che influenzano l'espressione di gruppi di geni
- Ogni fattore latente può influenzare più moduli con intensità differenti
- Consente di simulare processi biologici come vie di segnalazione o stati cellulari

### 5. Sovrapposizione tra moduli
- Possibilità per i geni di appartenere a più moduli contemporaneamente
- Riflette la realtà biologica in cui lo stesso gene può partecipare a diversi pathway

## Parametri

I nuovi parametri introdotti per controllare queste funzionalità sono:

```r
cell_specific_params = list(
  use_gene_modules = TRUE,          # Attiva/disattiva moduli di co-espressione
  n_gene_modules = 5,               # Numero di moduli di base
  module_correlation = 0.7,         # Intensità della correlazione tra geni dello stesso modulo
  module_hierarchical = FALSE,      # Attiva/disattiva gerarchia dei moduli
  module_overlap = 0.1,             # Proporzione di geni che possono appartenere a più moduli
  module_size_distribution = "exponential", # Tipo di distribuzione delle dimensioni ("exponential" o "uniform")
  n_latent_factors = 3,             # Numero di fattori latenti
  module_network_density = 0.2,     # Densità della rete di interazione tra moduli
  latent_factor_strength = 0.8      # Intensità dell'effetto dei fattori latenti
)
```

## Come funziona

1. **Creazione dei moduli di base**: I geni vengono assegnati ai moduli di base con dimensioni che seguono la distribuzione specificata.

2. **Creazione della gerarchia** (se abilitata): Alcuni moduli principali vengono divisi in sotto-moduli, mantenendo la relazione gerarchica.

3. **Generazione della rete di interazione**: Viene creata una matrice di interazione tra i moduli, dove ogni connessione rappresenta un "crosstalk" biologico.

4. **Generazione dei fattori latenti**: Vengono creati fattori latenti che influenzano gruppi di moduli correlati.

5. **Propagazione degli effetti**: Gli effetti dei fattori latenti vengono propagati attraverso la rete di moduli, simulando complesse interazioni regolative.

## Esempi di utilizzo

### Simulazione con molti piccoli moduli interagenti
```r
params <- list(
  use_gene_modules = TRUE,
  n_gene_modules = 10,
  module_hierarchical = TRUE,
  module_size_distribution = "exponential",
  n_latent_factors = 5,
  module_network_density = 0.3
)
```

### Simulazione con pochi grandi moduli indipendenti
```r
params <- list(
  use_gene_modules = TRUE,
  n_gene_modules = 3,
  module_hierarchical = FALSE,
  module_size_distribution = "uniform",
  n_latent_factors = 2,
  module_network_density = 0.1
)
```

## Impatto sulla simulazione

Questa implementazione migliorata consente di:

1. Generare dati di espressione genica che riflettono la co-regolazione osservata in dati reali
2. Simulare stati cellulari specifici attraverso l'attivazione di particolari fattori latenti
3. Modellare le reti di regolazione genica con diversi livelli di complessità
4. Simulare esperimenti con perturbazioni di specifici pathway
5. Generare dati di riferimento per testare metodi di identificazione di moduli
