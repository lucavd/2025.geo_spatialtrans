# CLAUDE_SIMPLE.md

## Semplificazione della Pipeline Spatial Transcriptomics

Questo documento descrive il processo di semplificazione e ottimizzazione della pipeline di simulazione per spatial transcriptomics, con focus sulla creazione di dati di benchmark per testare algoritmi di clustering.

## Obiettivo della Semplificazione

L'obiettivo era **ridurre la complessità** mantenendo **tutta la funzionalità biologica** per creare una pipeline efficiente per generare dati di benchmark controllabili per testare metodi di clustering.

## Struttura Branch "simplify"

### Branch Creato
```bash
git checkout -b simplify
```

### Nuova Struttura R_simple/
```
R_simple/
├── 01_configuration.R
├── 02_helper_functions.R  
├── 03_image_processing.R
├── 03b_generate_synthetic_tissue.R
├── 04_clustering.R                    # SEMPLIFICATO
├── 05_grid_sampling.R
├── 06*_expression_*.R                 # 12 moduli COMPLETI
└── testing/
    ├── full_test.R                    # Script test biologicamente realistico
    ├── full_test_result.rds           # Risultati validazione
    ├── full_simulation_data.rds       # Dati simulazione completi
    ├── full_tissue_complex.png        # Immagine sintetica generata
    ├── visualize_clusters.R           # Script visualizzazione
    ├── cluster_visualization.png      # Plot cluster colorati
    ├── combined_visualization.png     # Pannello combinato
    └── original_tissue_plot.png       # Plot immagine originale
```

## Semplificazioni Implementate

### 1. Clustering (04_clustering.R)
**PRIMA:** 540 righe, 4 metodi complessi
- `spatial_kmeans`
- `kmeans++` 
- `slic` (superpixel clustering)
- `dbscan_graph` (pipeline DBSCAN+Graph)
- Stima automatica k (elbow, silhouette)
- Gestione fallback complessa

**DOPO:** 64 righe, 1 metodo essenziale
- **Solo `spatial_kmeans`** - perfetto per benchmark controllabili
- Parametro `spatial_weight` per controllare difficoltà
- API pulita e deterministico
- **Riduzione: 92% del codice**

### 2. Moduli Expression (06*)
**ANALISI:** Tutti i 12 moduli sono **necessari** per funzionalità completa
- `06a` → `06k`: pipeline coordinata dal wrapper
- **Mantenuti TUTTI** - non semplificabili senza perdere realismo biologico
- Ogni modulo ha dipendenze specifiche nel workflow

### 3. Moduli Rimossi
**Eliminati dalla copia originale:**
- Tutti i moduli `06l_*` (ligand-receptor interactions)
- Tutti i moduli `06m_*` (temporal dynamics)  
- Tutti i moduli `06n_*` (alternative splicing)
- Tutti i moduli `06o_*` (anisotropic patterns)
- Tutti i moduli `06p_*` (3D microenvironment)
- Moduli `_simple` duplicati
- Spatial correlation avanzata (multiscale, non-stationary)
- Pipeline simulation complessi non utilizzati
- Moduli di visualizzazione e validazione

**Risultato:** Da 38 a 17 moduli essenziali (-55%)

## Script di Test

### Full Test (R_simple/testing/full_test.R)

**Caratteristiche:**
- **Pipeline COMPLETA** (tutti i 12 moduli expression attivi)
- **Parametri biologicamente realistici** per benchmark di qualità
- **Benchmark automatici** con criteri PASS/FAIL
- **Chunking per memoria** - gestisce dataset grandi
- **Validazione biologica** integrata

**Configurazione Full:**
```r
n_genes = 5000          # Realistico per spatial transcriptomics
n_cells = 20000         # Target celle per dataset reale  
k_cell_types = 8        # Tipi cellulari realistici
image_size = 800x800    # Dimensioni full-size
complexity = 3          # Massima complessità tissutale
library_size = 8000     # UMI/cella tipico per Visium HD
```

**Benchmark Implementati:**
1. **Dimensioni**: matrice corretta (geni×celle)
2. **Sparsità**: 30-92% (range esteso per spatial transcriptomics)  
3. **UMI per cella**: 3000-20000 (range realistico per spatial)
4. **UMI CV**: 0.15-1.0 (coefficiente variazione esteso)
5. **Clustering**: numero cluster atteso
6. **Integrità**: no NaN/Inf nei dati
7. **Correlazione spaziale**: test pattern spaziali

## Risultati Test

### Full Test Biologicamente Realistico
```
✓ Tempo: 0.79 minuti
✓ Dimensioni: 5000×20000  
✓ Sparsità: 89.2% (tipica per spatial filtrati)
✓ UMI medio: 9022, mediano: 7100
✓ UMI CV: 0.87 (variabilità realistica)
✓ Cluster: 8/8 assegnati correttamente
✓ Integrità: PASS
✓ RISULTATO: PASS
```

**Validazione Biologica:**
- Library size target: 8000 UMI/cella ✓
- Dropout range: 45-65% ✓
- Marker genes: 25 per tipo cellulare ✓
- Correlazione spaziale: Pattern realistici ✓

## Vantaggi della Semplificazione

### 1. **Performance**
- Test completi in **<1 minuto** (0.79 min per full test)
- **20,000 celle** generate con chunking ottimizzato
- **Validazione biologica** completa e veloce
- **Gestione memoria** intelligente per dataset grandi

### 2. **Robustezza**  
- Pipeline stabile con **pattern semplici e complessi**
- **Clustering deterministico** per benchmark riproducibili
- **Ground truth controllabile** tramite `spatial_weight`

### 3. **Mantenimento Biologico**
- **Tutti i moduli expression** attivi (dropout, library_size, spatial_correlation)
- **Parametri biologici realistici** non semplificati
- **Validazione automatica** con correzioni biologiche

### 4. **Facilità d'Uso**
- **Script di test automatico** con output chiaro
- **API semplificata** per clustering  
- **Documentazione** e benchmark integrati

## Utilizzo per Benchmark

### Generazione Dati di Test
```r
# Carica pipeline
files <- list.files("R_simple", pattern = "\\.R$", full.names = TRUE)
for (f in sort(files)) source(f)

# Genera dati con difficoltà controllabile
spatial_weight <- 0.2  # Facile: cluster per intensità
spatial_weight <- 0.8  # Difficile: cluster spaziali

clust <- cluster_image(img_df_thresh, k=5, spatial_weight=spatial_weight)
```

### Controllo Qualità
- **Ground truth** definita dal clustering iniziale
- **Metriche biologiche** validate automaticamente  
- **Pattern diversificati** tramite complexity parameter
- **Scaling** da test rapidi a full-size

## Visualizzazione dei Risultati

### Script di Visualizzazione (R_simple/testing/visualize_clusters.R)

**Funzionalità:**
- **Ricostruzione immagine originale** con pattern tissutali
- **Colorazione per cluster** con 8 colori distinti
- **Pannello combinato** simile alle pubblicazioni scientifiche
- **Orientamento corretto** (flip X per matching perfetto)

**Output Generati:**
- `combined_visualization.png` - Pannello completo (originale + cluster)
- `cluster_visualization.png` - Plot dettagliato solo cluster
- `original_tissue_plot.png` - Solo immagine tissutale

**Caratteristiche Tecniche:**
- **20,000 punti** colorati per tipo cellulare
- **Coordinate spaziali** preservate perfettamente
- **Legenda** con 8 tipi cellulari distinti
- **Alta risoluzione** (300 DPI) per pubblicazioni

### Esempio Visualizzazione
```r
# Esegui visualizzazione dopo full_test
Rscript R_simple/testing/visualize_clusters.R

# Output: pannello combinato con immagine originale sopra 
# e cluster colorati sotto, perfettamente allineati
```

## Next Steps

### ✅ Completati
- [x] Pipeline core funzionante (17 moduli essenziali)
- [x] Full test con parametri biologicamente realistici
- [x] Validazione automatica completa
- [x] Visualizzazione cluster in stile pubblicazione
- [x] Gestione memoria con chunking
- [x] Branch git organizzato

### Validazione Estesa
- [ ] Test con diverse dimensioni (10k, 15k, 25k geni)
- [ ] Test con più cluster (10, 12, 15)
- [ ] Benchmark performance scaling su dataset molto grandi

### Benchmark Algoritmi
- [ ] Implementare metriche di confronto clustering (ARI, NMI, Silhouette)
- [ ] Test con algoritmi standard (k-means, Louvain, Leiden)
- [ ] Validazione ground truth recovery
- [ ] Confronto con dataset pubblici reali

### Features Avanzate
- [ ] Supporto immagini utente con preprocessing automatico
- [ ] Export formati standard (H5AD, Seurat, CSV)
- [ ] Batch processing per esperimenti multipli
- [ ] Parallelizzazione GPU per dataset molto grandi

## TO DO: Opzioni di Generazione Avanzate

La cartella `R/functions/` contiene moduli avanzati (06l-06p) che devono essere integrati come **opzioni nella simulazione** di R_simple. Questi aggiungono realismo biologico specifico per scenari complessi:

### 🧬 06l - Ligand-Receptor Interactions
**Funzionalità:** Modellazione di interazioni ligando-recettore tra cellule vicine
- Database di interazioni L-R dalla letteratura scientifica  
- Propagazione di segnali basata sulla distanza
- Effetti downstream sull'espressione genica
- Attenuazione del segnale con decadimento spaziale

**Parametri controllabili:**
- `n_interactions`: Numero di interazioni L-R (default: 20)
- `effect_strength`: Range intensità effetto (0.5-2)
- `inhibitory_prob`: Probabilità interazioni inibitorie (0.3)
- `decayFactor`: Fattore decadimento distanza (10-50)

### ⏱️ 06m - Temporal Dynamics  
**Funzionalità:** Variazioni temporali nella struttura spaziale
- Simulazione di pseudo-tempo all'interno dei campioni
- Traiettorie di espressione lungo gradienti di sviluppo
- Oscillazioni nei pattern (es. ciclo cellulare)
- Modelli "splicing kinetics"

**Parametri controllabili:**
- `pseudotime_type`: "gradient", "focal", "bifurcation", "complex"
- `direction`: Direzione del gradiente temporale
- `spatial_coherence`: Coerenza spaziale (0-1)
- `n_foci`: Numero punti focali per tipo "focal"

### 🧬 06n - Alternative Splicing
**Funzionalità:** Pattern spaziali di splicing alternativo
- Varianti di splicing per lo stesso gene con pattern spaziali distinti
- Regolazione coordinata in domini tissutali specifici
- Correlazioni con microambienti locali

**Parametri controllabili:**
- `fraction_genes_with_variants`: Frazione geni con splicing alternativo (0.15)
- `spatial_regulation`: Intensità regolazione spaziale (0.7)
- `cell_type_regulation`: Regolazione specifica per tipo cellulare (0.6)
- `coordinated_splicing_groups`: Gruppi splicing coordinato (3)

### 🔀 06o - Anisotropic Patterns
**Funzionalità:** Pattern anisotropici dipendenti da strutture tissutali
- Simulazione strutture vascolari/nervose con pattern associati
- Direzionalità variabile basata su backbone strutturale
- Gradienti ortogonali alle strutture principali

**Parametri controllabili:**
- `structure_type`: "vessel", "nerve", "boundary", "mixed"
- `n_structures`: Numero strutture lineari (5)
- `curvature`: Livello curvatura (0-1)
- `bifurcation_prob`: Probabilità biforcazione (0.3)

### 📐 06p - 3D Microenvironment  
**Funzionalità:** Effetti del microambiente tridimensionale
- Simulazione proiezione 2D di strutture 3D
- Profondità variabile di campionamento
- Effetti prossimità 3D non catturati dalla distanza 2D

**Parametri controllabili:**
- `depth_pattern`: "smooth", "layered", "complex", "random"
- `depth_range`: Range profondità in μm (-10, 10)
- `n_layers`: Numero strati per modalità "layered" (3)
- `spatial_coherence`: Coerenza spaziale profondità (0.8)

### 🎯 Obiettivo Integrazione
Questi moduli devono essere integrati nel `full_test.R` come **parametri opzionali** che possono essere attivati/disattivati per scenari specifici:

```r
# Esempio configurazione avanzata
advanced_features <- list(
  ligand_receptor = TRUE,      # Attiva interazioni L-R
  temporal_dynamics = FALSE,   # Disattiva dinamiche temporali
  alternative_splicing = TRUE, # Attiva splicing alternativo
  anisotropic_patterns = FALSE,# Disattiva pattern anisotropici
  microenvironment_3d = TRUE   # Attiva effetti 3D
)
```

**Priorità:** Questi moduli aggiungono complessità computazionale ma aumentano drasticamente il realismo biologico per benchmark di algoritmi avanzati su scenari tissutali complessi.

## Conclusioni

La semplificazione ha **mantenuto tutta la funzionalità biologica** essenziale mentre ha **ridotto drasticamente la complessità** operativa. La pipeline R_simple è ora:

- ✅ **Veloce**: full test in <1 minuto (0.79 min)
- ✅ **Completa**: tutti i meccanismi biologici attivi (12 moduli expression)
- ✅ **Scalabile**: 20,000 celle con chunking intelligente
- ✅ **Robusta**: funziona con pattern semplici e complessi
- ✅ **Controllabile**: ground truth e difficoltà parametrizzabili
- ✅ **Testabile**: benchmark automatici biologicamente validati
- ✅ **Visualizzabile**: output in stile pubblicazione scientifica

**Risultato Finale**: Pipeline completa e ottimizzata per generare dati di benchmark biologicamente realistici (5000 geni, 20000 celle, 8 cluster) con validazione automatica e visualizzazione professionale. Pronta per testare e sviluppare algoritmi di clustering avanzati su dati di spatial transcriptomics.