---
type: Rule
title: "S3c, stimatore appaiato del tetto: T e P appaiati a V per area, bootstrap per coppie, controllo sintetico, regola dell'intervallo"
description: "Regola scritta e committata prima di calcolare D con lo stimatore appaiato: stimatore, incertezza, controlli, controprova, attese ed esito. Nessuna corsa nuova di catena: tabelle di S3b."
status: stable
generated: { by: geo-agent/progetto-master, at: 2026-11-23T23:30:00+01:00 }
sources:
  - id: s3
    path: docs/sonda/REGOLA_S3.md
  - id: s3b
    path: docs/sonda/REGOLA_S3b.md
  - id: diag
    path: results/sonda/S3b/diag_paired.csv
  - id: diagbp
    path: results/sonda/S3b/diag_paired_BP.csv
---
# Perché

In S3b (`REGOLA_S3b.md`[^s3b]) il tetto D(T) è uscito −0,769 h. Lo stimatore di `tools/riesame_D.py` sottrae come rumore s²/n_V + s²/n_T, cioè tratta T e V come aree indipendenti. Le aree di T sono invece le immagini registrate delle aree di V (stesso k): correlazione mediana fra T_k e V_k 0,94 sulle 12 voci (`diag_paired.csv`[^diag]). Con correlazione positiva il rumore vero della differenza è più piccolo, la correzione è troppo grande e D(T) è distorto verso il basso. Questa regola cambia solo come si stima il rumore. La quantità stimata resta la stessa: la media sulle 12 voci dello scarto sistematico al quadrato fra braccio e V, in unità h.

# Già noto prima del commit

Dalla diagnosi di S3b (`diag_paired.csv`[^diag], `diag_paired_BP.csv`[^diagbp]), con lo stesso stimatore puntuale di questa regola: D(T) appaiato 0,336; D(B) 3,15 in entrambi i modi (correlazione mediana B-V 0,11); D(P) appaiato 4,02. I valori puntuali non sono quindi una previsione. Sono nuovi: intervalli, controllo sintetico appaiato, scelta dell'intervallo, C-4 e C-5 con lo stimatore nuovo, controprova sulle permutazioni.

# Dati

Tabella `results/sonda/S3b/area_values.csv` e margini `results/sonda/S3b/h.csv`, invariati (md5 prima e dopo). Bracci: V (H2, 5 aree per tessuto), T (H1 omologo, stesso k), B (media dei 2 semi per area k), P (controllo positivo di REGOLA_S3[^s3]), B0 e B1 (B con il solo seme 0 o 1). M-4 si toglie per un tessuto solo se la media dei singoletti di V è ≥ 0,99 (regola di S3, invariata).

# Stimatore appaiato

Per braccio a, tessuto e misura (una cella): d_k = a_k − V_k sulle 5 aree appaiate per k; q = (media(d)² − var(d)/n) / h², var con ddof = 1, n = 5. D(a) = radice con segno della media di q sulle celle. È non distorto per la media di μ²/h² qualunque sia la correlazione fra a_k e V_k, se le coppie sono indipendenti fra loro. Si applica a tutti i bracci (T, B, P, B0, B1), così il metro è uno solo. Codice nuovo: `tools/sonda/S3c_D.py`; `tools/riesame_D.py` e `tools/sonda/S3_D.py` non cambiano (C-2 li usa).

# Incertezza

Unità indipendente: l'area. Bootstrap per coppie stratificato per tessuto: si ricampionano con reinserimento 5 indici k fra le aree di tumore e, indipendentemente, 5 fra quelle di stroma; gli stessi indici valgono per V e per tutti i bracci, quindi restano le correlazioni fra misure e fra bracci. B = 2000, seme 20261123. Si ricalcola q con media e varianza del campione ricampionato. Intervallo base come in S3: percentili 2,5 e 97,5 spostati della distorsione mediana; w = semi-ampiezza al 95 % della distribuzione di D(a) − D(b); se = deviazione standard bootstrap.

Due candidati per l'intervallo:

- (i) bootstrap per coppie come sopra;
- (ii) lo stesso, con semi-ampiezze ed se moltiplicate per f = √(n/(n−1)) · t(0,975; n−1) / 1,96 = 1,584 (n = 5; calcolato nel codice con scipy): corregge il restringimento del bootstrap con poche unità.

Regola di scelta, decisa dai dati sintetici (sotto): si usa (i) se la sua copertura è ≥ 0,90 in tutti i casi valutati; altrimenti (ii). Se anche (ii) sta sotto 0,90 in un caso, si usa (ii) e lo si dichiara come dubbio aperto.

# Controlli (cancelli; se uno fallisce D non si legge e nulla entra in `progetto.json`)

| id | controllo | esito richiesto |
|---|---|---|
| G-ID | identità algebrica: la media di q appaiato sulle 120 permutazioni di k di T (V fermo) è uguale a q indipendente di `tools/riesame_D.py` (stessa n, deviazione combinata) per ogni cella | scarto relativo ≤ 1e-9 |
| C-1p | sintetico appaiato: per tessuto, normale multivariata sui 24 valori per area (6 misure × V, T, B, P) con media = media osservata di V più lo spostamento e covarianza Σ = 0,95 S + 0,05 diag(S), S la covarianza osservata fra le 5 aree (rango 4, così diventa definita positiva); spostamento 0, 1, 2 h per ogni cella dei bracci; 1000 repliche, seme 20261127 | quadrato medio a 0 entro ±0,15; 1 h entro 1 ± 0,25; 2 h entro 2 ± 0,25 (tolleranze di C-1 di S3) |
| COV | copertura al 95 % dei candidati (i) e (ii), stessi dati sintetici: D(T), D(B), D(P) negli spostamenti 1 h e 2 h e nel caso «osservato» (medie dei bracci = medie osservate, verità = radice della media di (μ_a − μ_V)²/h²); differenze D(P) − D(T) e D(B) − D(T) nel caso osservato; 1000 repliche × 2000 bootstrap | decide l'intervallo (regola sopra); errore Monte Carlo ≈ 0,01 |
| C-2 | identità con il colon di `tools/sonda/S3_D.py`, invariata | come REGOLA_S3 |
| C-4 | controllo positivo P con lo stimatore appaiato e l'intervallo scelto | D(P) − D(T) «diverso», P peggio |
| C-5 | controllo negativo B0 contro B1, appaiati a V | non «diverso» |

C-3 (catena ridotta identica) resta nel check smoke, invariato.

# Controprova

CP-perm (l'appaiamento è vero): per tessuto, rumore R(π) = somma sulle 6 misure di var(T_π(k) − V_k)/h². La permutazione identica deve stare fra le più basse: rango ≤ 6 su 120 (5 %) in tumore e in stroma. Se fallisce, lo stimatore resta non distorto ma l'omologia delle aree non coglie la struttura fra aree: si scrive, e la diagnosi di S3b va rivista.

# Attese (scritte prima)

- G-ID, C-1p, C-2 passano (algebra e codice invariato).
- COV: con 5 aree il bootstrap (i) copre meno del 95 %; mi aspetto 0,80-0,90, quindi probabilmente si sceglie (ii).
- C-4 passa: P − T circa 3,7 h, w fra 1,5 e 2,7. C-5 passa.
- CP-perm: rango 1 su 120 in entrambi i tessuti.
- Pavimento − tetto circa 2,8 h > δ: G1 non chiude per la prima condizione.

# Esito e conseguenze

- Tutti i cancelli passano: pavimento = D(B) appaiato, tetto = D(T) appaiato, soglia = tetto + 1,0 h (sensibilità 0,5 e 2) in `progetto.json`; `results/misura.json` con A1 = D(B) e se dell'intervallo scelto; il check smoke usa `tools/sonda/S3c_D.py` con identità su `results/sonda/S3c/D.json`; bilancio dell'errore rifatto con T appaiato (`results/sonda/S3c/per_measure.csv`). S3c si chiude come misura su A1.
- Un cancello fallisce: numeri provvisori, `progetto.json` invariato; la sessione lo dichiara.

[^s3]: docs/sonda/REGOLA_S3.md
[^s3b]: docs/sonda/REGOLA_S3b.md
[^diag]: results/sonda/S3b/diag_paired.csv
[^diagbp]: results/sonda/S3b/diag_paired_BP.csv
