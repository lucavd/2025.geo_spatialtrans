# Strategia di Implementazione: Modulo Hotspot Cellulari

## Panoramica

Questo documento delinea una strategia dettagliata per l'implementazione di un modulo dedicato alla simulazione di hotspot cellulari all'interno del framework di simulazione di dati spaziali trascrittomici. Gli hotspot rappresentano nicchie cellulari spazialmente definite con proprietà trascrittomiche distinte, come microambienti tumorali, siti di infiammazione, o regioni specializzate di un tessuto.

## Architettura del Modulo

### File di Implementazione
```
R/functions/06q_cellular_hotspots.R
```

Il modulo sarà implementato come componente indipendente ma integrabile con i moduli esistenti, seguendo le convenzioni di nomenclatura e struttura del progetto.

### Flusso di Integrazione

1. Il modulo verrà chiamato dopo la generazione dei profili di espressione base
2. Modificherà l'espressione di geni selezionati (per indice) nelle regioni degli hotspot
3. Integrerà informazioni spaziali per creare pattern di espressione realistici
4. Manterrà la compatibilità con i moduli esistenti e i formati di output, lavorando esclusivamente con indici di geni anziché nomi biologici

## Funzioni Core del Modulo

### 1. Definizione degli Hotspot Spaziali

```r
define_spatial_hotspots <- function(
  coordinates,                # DataFrame con coordinate x, y delle celle
  n_hotspots = 3,             # Numero di hotspot da generare
  hotspot_config = list(
    size_range = c(50, 200),  # Range dimensioni hotspot (raggio in μm)
    shape = "circular",       # Forma: "circular", "elliptical", "irregular"
    intensity_profile = "gaussian", # Profilo intensità: "gaussian", "step", "linear"
    overlap_allowed = FALSE,  # Permette sovrapposizione tra hotspot
    border_gradient = TRUE,   # Gradiente all'interfaccia
    custom_centers = NULL     # Coordinate opzionali per i centri
  ),
  spatial_constraint = NULL,  # Vincolo spaziale opzionale (es. regioni permesse)
  random_seed = 42            # Seed per riproducibilità
) {
  # Implementazione:
  # 1. Posizionamento dei centri degli hotspot (casuale o specificato)
  # 2. Definizione dei confini e della forma
  # 3. Calcolo delle appartenenze e intensità per ogni cella
  # 4. Gestione dei gradienti ai bordi
  
  # Restituisce una lista con:
  # - Matrice di appartenenza delle celle agli hotspot (valori 0-1)
  # - Metadati degli hotspot (centri, dimensioni, tipi)
  # - Dati di distanza (distanza di ogni cella dal centro dell'hotspot più vicino)
}
```

### 2. Creazione di Moduli Genici Funzionali

```r
create_gene_modules <- function(
  n_genes,                  # Numero totale di geni nella simulazione
  module_config = list(
    proliferation = list(   # Modulo di geni relativi alla proliferazione
      size = 20,            # Dimensione del modulo (numero di geni)
      selection = "random"  # Metodo di selezione: "random", "sequential", "custom"
    ),
    hypoxia = list(         # Modulo di geni relativi all'ipossia
      size = 15, 
      selection = "random"
    ),
    # Altri moduli funzionali
    immune = list(size = 25, selection = "random"),
    ecm = list(size = 15, selection = "random"),
    stemness = list(size = 10, selection = "random")
  ),
  custom_indices = list(),  # Indici personalizzati per i moduli (se selection = "custom")
  random_seed = 42          # Seed per riproducibilità
) {
  # Implementazione:
  # 1. Inizializzazione della lista di moduli
  # 2. Per ogni modulo funzionale:
  #    - Se selection = "random", seleziona casualmente gli indici
  #    - Se selection = "sequential", usa blocchi sequenziali
  #    - Se selection = "custom", utilizza gli indici forniti
  # 3. Verifica che non ci siano sovrapposizioni tra moduli se richiesto
  
  # Restituisce:
  # - Lista di moduli, ciascuno contenente un vettore di indici di geni
}
```

### 3. Modificazione dell'Espressione Basata sugli Hotspot

```r
modify_expression_in_hotspots <- function(
  expression_matrix,       # Matrice di espressione (geni x celle)
  hotspot_membership,      # Matrice di appartenenza agli hotspot
  gene_modules,            # Lista di moduli genici (per indice)
  expression_config = list(
    # Configurazione per ciascun modulo e tipo di hotspot
    hotspot_1 = list(
      upregulated = c("proliferation", "stemness"),  # Moduli da sovraesprimere
      downregulated = c("immune"),                   # Moduli da sottoesprimere
      up_factor = 2.5,      # Fattore di sovraespressione
      down_factor = 0.4     # Fattore di sottoespressione
    ),
    hotspot_2 = list(
      # Configurazione per hotspot 2
    )
  ),
  gradient_effect = TRUE,  # Applicare effetto gradiente basato su distanza dal centro
  random_seed = 42         # Seed per riproducibilità
) {
  # Implementazione:
  # 1. Per ogni hotspot e cella associata:
  #    - Calcolo del fattore di appartenenza (considerando gradiente)
  #    - Identificazione dei geni target (dai moduli definiti)
  #    - Applicazione dei fattori di modificazione
  # 2. Gestione delle interferenze tra hotspot sovrapposti
  # 3. Normalizzazione delle modifiche per mantenere proprietà statistiche
  
  # Restituisce:
  # - Matrice di espressione modificata
  # - Metadati delle modifiche applicate (per validazione)
}
```

### 4. Simulazione del Microambiente degli Hotspot

```r
simulate_hotspot_microenvironment <- function(
  expression_matrix,        # Matrice espressione (geni x celle)
  coordinates,              # Coordinate spaziali delle celle
  hotspot_membership,       # Appartenenza delle celle agli hotspot
  gene_modules,             # Moduli genici funzionali
  microenvironment_config = list(
    hypoxia = list(
      enabled = TRUE,       # Attiva simulazione ipossia
      core_intensity = 0.8, # Intensità al centro
      gradient_range = 0.6, # Estensione del gradiente
      target_modules = c("hypoxia")  # Moduli target
    ),
    immune_infiltration = list(
      enabled = TRUE,       # Attiva infiltrazione immunitaria
      intensity = 0.5,      # Intensità infiltrazione
      location = "border",  # Localizzazione: "border", "core", "random"
      target_modules = c("immune")  # Moduli target
    ),
    ecm_remodeling = list(
      enabled = TRUE,       # Attiva rimodellamento ECM
      intensity = 0.7,      # Intensità rimodellamento
      target_modules = c("ecm")  # Moduli target
    )
  ),
  random_seed = 42           # Seed per riproducibilità
) {
  # Implementazione:
  # 1. Applicazione di pattern di espressione per ciascun processo biologico
  #    usando gli indici dei geni dai moduli specificati
  # 2. Creazione di gradienti biologicamente plausibili
  # 3. Integrazione delle diverse alterazioni in modo coerente
  # 4. Conservazione delle proprietà statistiche della distribuzione
  
  # Restituisce:
  # - Matrice di espressione con microambiente simulato
  # - Metadati delle alterazioni applicate
}
```

### 5. Simulazione delle Interazioni Cellula-Cellula

```r
simulate_cell_interactions <- function(
  expression_matrix,        # Matrice espressione (geni x celle)
  coordinates,              # Coordinate spaziali
  hotspot_membership,       # Appartenenza agli hotspot
  gene_modules,             # Moduli genici funzionali
  interaction_config = list(
    interaction_range = 50, # Raggio di interazione (μm)
    interaction_pairs = list(
      pair_1 = list(
        source_module = "stemness",  # Modulo dei geni "mittente"
        target_module = "proliferation", # Modulo dei geni "ricevente"
        effect = "up",    # Effetto: "up" o "down"
        strength = 0.7    # Forza dell'effetto
      ),
      # Altre coppie di interazione
    ),
    decay_with_distance = TRUE,  # L'effetto diminuisce con la distanza
    decay_factor = 0.1          # Fattore di diminuzione per unità di distanza
  ),
  random_seed = 42           # Seed per riproducibilità
) {
  # Implementazione:
  # 1. Costruzione della rete di interazioni tra celle vicine
  # 2. Per ogni coppia di celle entro il raggio di interazione:
  #    - Calcolo della forza dell'interazione basata sulla distanza
  #    - Per ogni coppia di interazione nei parametri:
  #      - Modifica dell'espressione dei geni target in base all'espressione dei geni sorgente
  # 3. Gestione delle interazioni cumulative
  
  # Restituisce:
  # - Matrice di espressione modificata dalle interazioni
  # - Rete di interazioni cellula-cellula (opzionale)
}
```

### 6. Funzione Wrapper Principale

```r
generate_cellular_hotspots <- function(
  expression_data,          # Lista con expression_matrix e coordinates
  n_genes,                  # Numero totale di geni nella simulazione
  hotspot_params = list(    # Parametri configurazione completa
    enabled = TRUE,         # Attiva/disattiva l'intero modulo
    n_hotspots = 3,         # Numero di hotspot
    # Configurazioni per i vari componenti
    hotspot_config = list(...),
    module_config = list(...),
    expression_config = list(...),
    microenvironment_config = list(...),
    interaction_config = list(...)
  ),
  random_seed = 42           # Seed per riproducibilità
) {
  # Implementazione:
  # 1. Verifica e preparazione dei dati di input
  # 2. Definizione degli hotspot spaziali
  # 3. Creazione dei moduli genici funzionali
  # 4. Applicazione delle modifiche di espressione negli hotspot
  # 5. Simulazione del microambiente se richiesto
  # 6. Simulazione delle interazioni cellulari se richiesto
  
  # Restituisce un oggetto compatibile con il pipeline esistente:
  # - Matrice di espressione modificata
  # - Metadati degli hotspot per visualizzazione
  # - Statistiche di modifica per validazione
}
```

## Modelli Funzionali e Configurazioni Predefinite

Il modulo includerà configurazioni predefinite per diversi tipi di hotspot, utilizzando indici di geni anziché nomi biologici:

### 1. Configurazione per Microambiente Tumorale

```r
tumor_hotspot_config <- function(n_genes, random_seed = 42) {
  # Genera una configurazione predefinita per simulare un microambiente tumorale
  
  # Usa sottoinsiemi casuali di geni per rappresentare i moduli funzionali
  set.seed(random_seed)
  
  # Configurazione per i moduli genici
  module_config <- list(
    proliferation = list(size = round(n_genes * 0.05), selection = "random"),
    hypoxia = list(size = round(n_genes * 0.03), selection = "random"),
    angiogenesis = list(size = round(n_genes * 0.02), selection = "random"),
    immune = list(size = round(n_genes * 0.04), selection = "random"),
    ecm = list(size = round(n_genes * 0.03), selection = "random"),
    stemness = list(size = round(n_genes * 0.02), selection = "random")
  )
  
  # Configurazione per gli hotspot
  hotspot_config <- list(
    n_hotspots = 2,
    size_range = c(50, 150),
    shape = "irregular",
    intensity_profile = "gaussian",
    overlap_allowed = FALSE
  )
  
  # Configurazione per le modifiche di espressione
  expression_config <- list(
    hotspot_1 = list(
      upregulated = c("proliferation", "hypoxia", "angiogenesis", "stemness"),
      downregulated = c("immune"),
      up_factor = 2.5,
      down_factor = 0.3
    ),
    hotspot_2 = list(
      upregulated = c("proliferation", "stemness", "ecm"),
      downregulated = c("immune"),
      up_factor = 2.0,
      down_factor = 0.4
    )
  )
  
  # Configurazione per il microambiente
  microenvironment_config <- list(
    hypoxia = list(
      enabled = TRUE,
      core_intensity = 0.8,
      gradient_range = 0.6,
      target_modules = c("hypoxia", "angiogenesis")
    ),
    immune_infiltration = list(
      enabled = TRUE,
      intensity = 0.6,
      location = "border",
      target_modules = c("immune")
    ),
    ecm_remodeling = list(
      enabled = TRUE,
      intensity = 0.7,
      target_modules = c("ecm")
    )
  )
  
  # Configurazione per le interazioni
  interaction_config <- list(
    interaction_range = 80,
    interaction_pairs = list(
      pair_1 = list(
        source_module = "hypoxia",
        target_module = "angiogenesis",
        effect = "up",
        strength = 0.8
      ),
      pair_2 = list(
        source_module = "proliferation",
        target_module = "immune",
        effect = "down",
        strength = 0.6
      )
    ),
    decay_with_distance = TRUE,
    decay_factor = 0.05
  )
  
  # Restituisci configurazione completa
  return(list(
    module_config = module_config,
    hotspot_config = hotspot_config,
    expression_config = expression_config,
    microenvironment_config = microenvironment_config,
    interaction_config = interaction_config
  ))
}
```

### 2. Configurazione per Nicchie di Cellule Staminali

```r
stem_cell_niche_config <- function(n_genes, random_seed = 42) {
  # Configurazione simile a quella tumorale, ma con parametri adattati alle
  # nicchie di cellule staminali e appropriati sottoinsiemi di geni
  
  # ...implementazione...
  
  return(config)
}
```

### 3. Configurazione per Siti Infiammatori

```r
inflammatory_site_config <- function(n_genes, random_seed = 42) {
  # Configurazione per siti infiammatori con appropriati sottoinsiemi di geni
  # per simulare risposte infiammatorie
  
  # ...implementazione...
  
  return(config)
}
```

## Considerazioni sulla Validazione

Il modulo includerà funzioni di validazione specifiche per valutare la plausibilità biologica degli hotspot generati:

```r
validate_hotspot_simulation <- function(
  original_expression,      # Matrice originale prima delle modifiche
  modified_expression,      # Matrice dopo l'applicazione degli hotspot
  hotspot_metadata,         # Metadati degli hotspot generati
  gene_modules,             # Moduli genici utilizzati
  coordinates,              # Coordinate spaziali
  validation_plots = TRUE,  # Generare plot di validazione
  output_dir = "plots/hotspot_validation" # Directory output
) {
  # Implementa validazioni multiple:
  
  # 1. Calcolo delle differenze di espressione per hotspot e per modulo
  #    - Verifica che i cambiamenti di espressione siano coerenti con la configurazione
  
  # 2. Analisi dell'autocorrelazione spaziale Moran's I
  #    - Verifica la presenza di clustering spaziale
  
  # 3. Analisi del gradiente
  #    - Verifica la presenza di gradienti biologicamente plausibili
  
  # 4. Verifica della coerenza statistica
  #    - Controlla che le distribuzioni di espressione rimangano realistiche
  
  # 5. Generazione di plot di validazione
  #    - Mappe spaziali degli hotspot
  #    - Heatmap differenziali di espressione
  #    - Gradienti di espressione su distanza
  #    - Correlogrammi spaziali
  
  # Restituisce metriche di validazione e genera plot se richiesto
}
```

## Integrazione con la Pipeline Esistente

Il modulo si integrerà con il framework esistente in due punti principali:

### 1. Integrazione con la Generazione dei Profili di Espressione

Estendere il file `06k_expression_profiles_wrapper.R` per supportare il modulo hotspot:

```r
# Aggiunta al wrapper esistente
generate_expression_profiles <- function(..., hotspot_params = NULL, ...) {
  # Codice esistente per generare i profili di base
  
  # Aggiungi supporto per hotspot se abilitato
  if (!is.null(hotspot_params) && hotspot_params$enabled) {
    # Applica il modulo hotspot
    result <- generate_cellular_hotspots(
      expression_data = result,
      n_genes = n_genes,
      hotspot_params = hotspot_params,
      random_seed = random_seed
    )
  }
  
  return(result)
}
```

### 2. Supporto nei Parametri di Difficoltà

Estendere `07b_difficulty_setup.R` per includere parametri degli hotspot per ciascun livello di difficoltà:

```r
# Estensione dei livelli di difficoltà con parametri hotspot
configure_difficulty_level <- function(level) {
  # Codice esistente
  
  # Aggiungi configurazione hotspot per ogni livello
  if (level == "easy") {
    diff_cfg$hotspot_params <- list(
      enabled = TRUE,
      n_hotspots = 1,
      hotspot_config = list(
        size_range = c(100, 150),
        shape = "circular",
        intensity_profile = "gaussian"
      ),
      # Configurazioni semplificate per moduli e microambiente
      module_config = list(
        proliferation = list(size = 15, selection = "sequential"),
        immune = list(size = 10, selection = "sequential"),
        ecm = list(size = 5, selection = "sequential")
      ),
      # Configurazioni di espressione semplificate
      expression_config = list(
        hotspot_1 = list(
          upregulated = c("proliferation"),
          downregulated = c("immune"),
          up_factor = 1.5,
          down_factor = 0.5
        )
      ),
      # Configurazioni semplificate per altri componenti
      microenvironment_config = list(simple = TRUE),
      interaction_config = list(enabled = FALSE)
    )
  } else if (level == "medium") {
    # Configurazione media complessità
    diff_cfg$hotspot_params <- list(
      enabled = TRUE,
      n_hotspots = 2,
      # ...configurazione più complessa
    )
  } else if (level == "hard") {
    # Configurazione alta complessità
    diff_cfg$hotspot_params <- list(
      enabled = TRUE,
      n_hotspots = 4,
      # ...configurazione molto complessa
    )
  }
  
  return(diff_cfg)
}
```

## Considerazioni sulla Creazione di Moduli Derivati

Il framework di hotspot cellulari supporterà l'estensione per scenari specifici:

1. **Progressione tumorale**: Estensione temporale che simula l'evoluzione degli hotspot nel tempo

2. **Risposta al trattamento**: Modulo che simula la risposta differenziale degli hotspot a trattamenti terapeutici

3. **Pattern multicellulari complessi**: Supporto per pattern di hotspot con strutture gerarchiche annidate

4. **Integrazioni multi-omiche**: Estensione per integrare dati epigenetici, proteomici o metabolomici con i pattern di espressione

## Implementazione Progressiva

L'implementazione del modulo seguirà queste fasi:

1. **Fase 1**: Implementazione delle funzioni core per la definizione spaziale degli hotspot e la creazione dei moduli genici
   - Output: Funzionalità base per creare hotspot e organizzare i geni in moduli funzionali
   - Test: Validazione spaziale e organizzazione modulare

2. **Fase 2**: Implementazione delle modifiche di espressione basate sugli hotspot
   - Output: Capacità di alterare l'espressione genica in specifiche regioni spaziali
   - Test: Verifica dei pattern di espressione differenziale

3. **Fase 3**: Implementazione del microambiente e delle interazioni cellula-cellula
   - Output: Simulazione di interazioni complesse tra celle negli hotspot
   - Test: Validazione dei pattern di co-espressione e correlazione spaziale

4. **Fase 4**: Integrazione completa con il pipeline e creazione di visualizzazioni avanzate
   - Output: Modulo completamente integrato con supporto alla validazione
   - Test: Simulazione end-to-end con hotspot complessi

## Conclusione

Questo modulo hotspot fornirà un potente strumento per simulare nicchie cellulari spazialmente definite all'interno del framework di simulazione di dati spaziali trascrittomici. Lavorando esclusivamente con indici di geni anziché nomi biologici, il modulo manterrà la coerenza con l'approccio esistente, consentendo comunque di simulare pattern realistici come microambienti tumorali, nicchie staminali e siti infiammatori.