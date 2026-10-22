| frase del riassunto | fonte (tabella o file, riga) | valore letto |
|---|---|---|
| 240 000 nuclei veri delle 30 aree | R6_refs_type_ph1.csv, somma di n | 240 775 |
| cellule di Paneth chiamate «endoteliali» (5 260 nuclei) | labels/A1_exp_..._clusters.csv, cluster 13 (Defa21, Lyz1, Itln1) | 5260 nuclei, type = endoteliale |
| tre gruppi di cardiomiociti senza nome | labels/A6_exp_..._clusters.csv, type = non risolto | 3 cluster, 981 nuclei |
| accordo con Space Ranger 91 % e 79 % (corretto), 88 % e 77 % (pre-registrato) | R6_predictions_ph1.csv e R6_predictions.csv, P-6 A5/A6 | 0.907, 0.794; 0.881, 0.774 |
| impronte spostate di 20 µm: accordo al livello del caso | R6_posthoc_kappa.csv (PH-2), kappa primaria contro spostata | A5 primary: 0.83; A5 shifted: 0.13; A6 primary: 0.29; A6 shifted: 0.10 |
| neuroni circa 3 volte la glia (3,2 corretto; 2,5 pre-registrato, al limite) | R6_predictions.csv / _ph1.csv, P-1 | 2.50 [2.10, 2.95]; 3.22 [2.54, 3.68] |
| neuroni 111 µm², glia 34 (astrociti 50, oligodendrociti 31, microglia 23) | R6_refs_type_ph1.csv A5; R6_posthoc.csv glia pooled | 111; 34; 50, 31, 23 |
| neuroni due nuclei su tre | R6_refs_type_ph1.csv A5 neuroni fraction | 0.65 |
| oligodendrociti dall'8 % al 36 % fra le aree | R6_refs_roi_type_ph1.csv A5 oligodendrociti | 0.08 - 0.36 |
| differenze fra tipi del 24-42 % in intestino, tumore, stroma | R6_predictions_ph1.csv P-2 A1, A2, A3 | A1 1.24, A2 1.39, A3 1.42 |
| linfonodo: B e T differiscono del 5 % (16,8 contro 16,0 µm²) | R6_refs_type_ph1.csv A4; P-2 A4 | 16.8, 16.0; P-2 1.047 |
| cuore: differenza fra tipi 11 %; cardiomiociti da 16 a 28 µm² fra le aree | P-2 A6 (_ph1); R6_refs_roi_type_ph1.csv A6 cardiomiociti | 1.11; 16.1 - 28.3 |
| tumore del colon: dispersione dentro il solo epitelio 0,56 (tutte le cellule 0,55) | R6_posthoc.csv PH-6 A2 epiteliale normale; BL-138 (sigma 0.55) | 0.56 |
| cardiomiociti non più allungati: 0,82 contro 0,83 | R6_refs_type_ph1.csv A6 median_ecc_nuc; P-3 | 0.817, fibroblasti 0.830, immunitari 0.824; P-3 -0.011 [-0.043, 0.013] |
| intestino: epitelio 68 %, stroma e muscolo 21 %, immunitarie 10 % | R6_refs_type_ph1.csv A1 | 0.68, 0.21, 0.10 |
| tumore del colon: epitelio 71 %, fibroblasti 20 % | R6_posthoc.csv PH-4 A2; R6_refs_type_ph1.csv A2 | 0.71 [0.63, 0.79]; 0.20 |
| stroma: fibroblasti 32 %, epitelio 41 %, immunitarie 19 % | R6_refs_type_ph1.csv A3; PH-4 A3 | 0.32, 0.41, 0.19 |
| linfonodo: B 61 % nei follicoli, 42 % nella paracorticale; 46 % senza tipo; 84 trascritti per nucleo | R6_posthoc.csv PH-5; R6_refs_type_ph1.csv A4 non risolto; median_umi | 0.61, 0.42; 0.46; umi mediano A4 84 |
| cuore: cardiomiociti 76 %, fibroblasti 13 %, immunitarie e vascolari 11 % | R6_refs_type_ph1.csv A6 | 0.76, 0.13, 0.11 |
