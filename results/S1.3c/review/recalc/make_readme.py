
import pandas as pd, numpy as np, json, os, glob
V=os.path.expanduser('~/2025.geo_spatialtrans/results/S1.3c/review/recalc'); R=os.path.dirname(os.path.dirname(V))
T=pd.read_csv(V+'/recalc_vs_reported.csv')
def f(v, q=''):
    if v is None or (isinstance(v,float) and np.isnan(v)): return '—'
    if float(v).is_integer() and abs(v)>=1 and not q.startswith(('phi','E','BA','R ','umi','frac')): return f'{int(v)}'
    return f'{v:.6g}'
L=[]
L.append('# S1.3c — revisione avversariale «Ricalcolo» (recalc)\n')
L.append('Revisore: sotto-agente «Ricalcolo», 2026-10-12. Codice indipendente: nessun import da `tools/S1.3c_ecm.py`, nessuna funzione copiata; '
         'il codice della pipeline è stato **letto solo dopo** il primo confronto, per trovare la causa delle discrepanze. Commit di riferimento `5fa0294`.\n')
L.append('## File\n- `recalc_umi_bg.py` — UMI per bin (geni della filtrata letti dalla grezza; tutti i geni), C-3, soglie di sfondo q95 (k = 3, 5) per i 5 dataset → `recalc_bg_<dataset>.json`.\n'
         '- `recalc_roi.py` — per ROI: n bin, bin in tessuto valido, territori (generatore R3 più vicino, cKDTree), d_nuc (EDT nativa sui nuclei `keep`, scala 1), finestre, φ_obs, nullo (permutazione **mia**, seed 9001–9030, 30 permutazioni; 200 per A2 e A6), E, r*, BA, R, CP-2; `ctrl` = controllo positivo tipo A su A3 r1. Variabile `RECALC_MAP`: `rint_noclip` (prima versione del revisore), `floor_clip` (mappatura esatta della maschera R2), `px4` (mappatura della pipeline, per isolare la causa).\n'
         '- `compare_bg.py`, `compare_roi.py`, `summarise.py` — confronti; `recalc_vs_reported.csv` = tabella completa (198 righe).\n'
         '- Intermedi grandi (bin per dataset) in `/tmp/recalc_s13c/` sul server (fuori repo, non versionati).\n')
# 1. background
L.append('## 1. UMI e sfondo (tutti e 5 i dataset; richiesto: Lymph Node)\n')
L.append('| dataset | C-3 bin in tessuto | n filtrata | bin con UMI(geni comuni, grezza) ≠ filtrata | idem tutti i geni | soglia k=3 (rip./ric.) | n finestre k=3 (rip./ric.) | soglia k=5 (rip./ric.) | n finestre k=5 (rip./ric.) |')
L.append('|---|---|---|---|---|---|---|---|---|')
for p in sorted(glob.glob(V+'/recalc_bg_*.json')):
    ds=os.path.basename(p)[10:-5]; me=json.load(open(p)); rep=json.load(open(f'{R}/S1.3c_bg_{ds}_bg.json'))
    b3m,b3r=me['bg']['com']['3'],rep['bg']['com']['3']; b5m,b5r=me['bg']['com']['5'],rep['bg']['com']['5']
    L.append(f"| {ds} | {me['C3']['n_in_tissue']} | {me['C3']['n_filtered']} | {me['C3']['n_diff_com']} (rip. {rep['C3']['common_genes']['n_diff']}) | {me['C3']['n_diff_all']} (rip. {rep['C3']['prereg_all_genes']['n_diff']}) | {b3r['thr']:g} / {b3m['thr']:g} | {b3r['n_windows']} / {b3m['n_windows']} | {b5r['thr']:g} / {b5m['thr']:g} | {b5r['n_windows']} / {b5m['n_windows']} |")
si=json.load(open(V+'/recalc_bg_Visium_HD_Mouse_Small_Intestine.json'))
L.append(f"\nAnche media e frazione di zeri delle finestre di sfondo (com/all, k = 3/5) coincidono entro 1e-12 in tutti i dataset. "
         f"Small Intestine (Space Ranger 3.0.0): {si['n_umi_missing']} barcode di `tissue_positions` assenti dalla matrice grezza, di cui {si['n_umi_missing_in_tissue']} con `in_tissue = 1`; vanno contati come 0 UMI (come fa la pipeline). La mia prima versione li lasciava NaN → soglia NaN: errore **mio**, corretto.\n")
# 2. ROI table k=3
L.append('## 2. ROI (variante primaria, geni comuni)\n')
L.append('Colonne: **rip.** = riportato; **ric. R2** = mio ricalcolo con la mappatura bin→maschera ds4 esatta di R2 (`floor(px/4)`, clip all’ultimo pixel); **ric. pipeline** = mio ricalcolo con la mappatura della pipeline (`px4 = side_um/913`); per φ_null ed E il ricalcolo usa le mie 30 permutazioni, quindi si riporta anche z = (ric. − rip.)/SE, SE = sd·√(1/20+1/30).\n')
qs=['n_bins_roi (C-1)','n_tissue_bins (PH1)','n_gen','R [4,8)/[0,2)','n_lumen_bins CP-2 (PH1)','umi_lumen_mean (PH1)']
qk=['thr','n_tissue_windows','n_empty','phi_obs','phi_null (media)','phi_null_min','phi_null_max','E','r_star','BA','CP-2 n_windows','CP-2 n_empty']
for k in (3,5):
    L.append(f'\n### k = {k} (s = {2*k} µm)\n')
    L.append('| ROI | grandezza | rip. | ric. R2 | ric. pipeline | ric.R2 − rip. | z (nullo) |'); L.append('|---|---|---|---|---|---|---|')
    for A in ['A1_r1','A2_r1','A3_r1','A4_f1','A5_r1','A6_r1']:
        rows=[]
        if k==3: rows+= [T[(T.unit==A)&(T.quantity==q)].iloc[0] for q in qs]
        rows+= [T[(T.unit==f'{A}_k{k}')&(T.quantity==q)].iloc[0] for q in qk]
        for x in rows:
            z='' if pd.isna(x.seed_se) or x.seed_se==0 else f'{x.diff_px4/x.seed_se:.2f}'
            L.append(f'| {A} | {x.quantity} | {f(x.reported,x.quantity)} | {f(x.recalc_floor,x.quantity)} | {f(x.recalc_px4,x.quantity)} | {f(x.diff_floor,"phi")} | {z} |')
L.append('\nIntervallo del nullo: il mio min–max su 30 permutazioni si sovrappone a quello riportato (20 permutazioni) in tutti i ROI; con 200 permutazioni (A2, A6) il riportato cade sempre dentro il mio min–max (vedi rilievo RC-5).\n')
# 3. ctrl
L.append('## 3. Controllo positivo tipo A su A3 r1 (k = 3)\n')
L.append('| modo | grandezza | rip. (C5 / PH-3 f=0.10) | ric. R2 | ric. pipeline | z |'); L.append('|---|---|---|---|---|---|')
for x in T[T.unit.str.startswith('A3_r1_ctrlA')].itertuples():
    z='' if pd.isna(x.seed_se) else f'{x.diff_px4/x.seed_se:.2f}'
    L.append(f"| {x.unit.split('_')[-1]} | {x.quantity} | {f(x.reported,'phi')} | {f(x.recalc_floor,'phi')} | {f(x.recalc_px4,'phi')} | {z} |")
L.append('\n`thr` = azzerati i bin di tessuto con d_nuc > r₉₀ (definizione della pipeline, nota 12); `rank` = esattamente il 10 % dei bin con d_nuc più alto (pareggi rotti a caso, seed 777), lettura letterale del compito. C5 `A_E` e PH-3 (A, f = 0.10) sono lo stesso numero (0.0688495) perché escono dallo stesso codice: non sono due conferme indipendenti.\n')
open(V+'/README_tables.md','w').write('\n'.join(L)+'\n')
print(len(L))
