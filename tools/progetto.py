#!/usr/bin/env python3
"""progetto.py - apertura, controllo e chiusura di una sessione (skill progetto-master).

Solo libreria standard di Python 3. Si copia in tools/ del repository.
  python3 tools/progetto.py init
  python3 tools/progetto.py apri
  python3 tools/progetto.py quadro
  python3 tools/progetto.py check smoke|full
  python3 tools/progetto.py chiudi --sessione ID --ipotesi TESTO --affermazione ID [--atteso A] [--infra | --riesame] [--nota TESTO]
  python3 tools/progetto.py cancello --fase FASE --consenso TESTO [--budget N]
  python3 tools/progetto.py proroga --sessioni N --motivo TESTO --atteso X --altrimenti TESTO --consenso TESTO

Stato: progetto.json (contratto in numeri, cancelli, proroghe), results/misure.tsv (una riga per sessione chiusa),
results/quadro.tsv (a ogni chiusura, il valore di ogni affermazione).
Il check scrive results/misura.json: {"valore": x, "se": s, "voci": {"nome": v, ...}, "affermazioni": {"A2": {"valore": y, "se": t}, ...}}.
'voci' e' facoltativo: se la misura principale e' una media, sono le sue parti, e il quadro le mostra dalla piu' grande.
'valore' e 'se' in cima sono la misura principale, cioe' l'affermazione con "principale": true.
Quadro dell'articolo: per ogni affermazione il valore di oggi contro il numero che serve. 'regge' se l'intervallo
valore +/- 2 se sta tutto dalla parte giusta, 'non regge' se sta tutto dalla parte sbagliata, altrimenti 'incerta'.
Esito di una sessione: in costruzione 'avanza' se l'affermazione servita batte il suo migliore valore tenuto di oltre
2 errori standard e raggiunge l'atteso, altrimenti 'colpo'; nelle altre fasi 'misura'; 'contratto' nella scommessa;
'infra' o 'riesame' se dichiarati.
"""
import argparse
import csv
import datetime as dt
import hashlib
import json
import os
import subprocess
import sys

CFG, MIS, QST, LAST = "progetto.json", "results/misure.tsv", "results/quadro.tsv", ".progetto_check.json"
COLS = ["data", "sessione", "fase", "affermazione", "ipotesi", "atteso", "valore", "se", "esito", "nota"]
QCOLS = ["data", "sessione", "affermazione", "valore", "se", "stato"]
FASI = ["scommessa", "sonda", "costruzione", "conferma", "scrittura"]
DEFAULT = {
    "messaggio": "DA SCRIVERE: il messaggio dell'articolo in una frase",
    "misura": "DA SCRIVERE: misura principale sull'oggetto vero",
    "verso": "alto",
    "pavimento": None, "tetto": None, "soglia": None,
    "affermazioni": [{"id": "A1", "testo": "DA SCRIVERE: affermazione principale", "figura": "Fig. 1", "principale": True}],
    "fase": "scommessa",
    "budget": {"scommessa": 2, "sonda": 3, "costruzione": None, "conferma": 2, "scrittura": None},
    "cancelli": [], "proroghe": [],
    "check": {"smoke": "echo 'DA SCRIVERE: comando smoke' && false", "full": "echo 'DA SCRIVERE: comando full' && false"},
    "file_misura": "results/misura.json",
    "tetti_parole": {"docs/PAPER.md": 500, "docs/STATO.md": 800},
}
ESCLUSI = [":!docs", ":!" + MIS, ":!" + QST, ":!" + LAST, ":!" + CFG]


def sh(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True)


def oggi_iso():
    return dt.date.today().isoformat()


def carica():
    if not os.path.exists(CFG):
        sys.exit(f"Manca {CFG}: esegui 'init' nella radice del repository.")
    with open(CFG, encoding="utf-8") as f:
        return json.load(f)


def salva(c):
    with open(CFG, "w", encoding="utf-8") as f:
        json.dump(c, f, indent=2, ensure_ascii=False)
        f.write("\n")


def leggi_tsv(path):
    if not os.path.exists(path):
        return []
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def scrivi_tsv(path, cols, nuove):
    """Aggiunge righe; se l'intestazione del file e' vecchia, lo riscrive con le colonne di oggi."""
    vecchie = leggi_tsv(path)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for x in vecchie + nuove:
            w.writerow({k: x.get(k, "") for k in cols})


def righe():
    return leggi_tsv(MIS)


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def parole(path):
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8") as f:
        return len(f.read().split())


def impronta():
    """Stato del codice, esclusi documenti e registri: cambia se il codice cambia dopo il check."""
    args = " ".join(f"'{e}'" for e in ESCLUSI)
    testo = sh("git rev-parse HEAD").stdout + sh(f"git status --porcelain -- . {args}").stdout + sh(f"git diff HEAD -- . {args}").stdout
    return hashlib.sha1(testo.encode()).hexdigest()


def segno(c):
    return 1.0 if c.get("verso", "alto") == "alto" else -1.0


def affermazioni(c):
    return c.get("affermazioni") or []


def principale(c):
    affs = affermazioni(c)
    for x in affs:
        if x.get("principale"):
            return x["id"]
    return affs[0]["id"] if affs else ""


def serve(c, aff):
    """Il numero che serve all'affermazione e il suo verso (+1 meglio in alto, -1 meglio in basso)."""
    if aff["id"] == principale(c):
        return num(c.get("soglia")), segno(c)
    return num(aff.get("serve")), (1.0 if aff.get("verso", c.get("verso", "alto")) == "alto" else -1.0)


def misure_oggi(c):
    """I valori scritti dall'ultimo check: {id: (valore, se)}."""
    fm = c.get("file_misura", "results/misura.json")
    if not os.path.exists(fm):
        return {}
    with open(fm, encoding="utf-8") as f:
        m = json.load(f)
    out = {}
    for k, v in (m.get("affermazioni") or {}).items():
        if num(v.get("valore")) is not None:
            out[k] = (float(v["valore"]), float(v.get("se") or 0.0))
    if num(m.get("valore")) is not None and principale(c):
        out[principale(c)] = (float(m["valore"]), float(m.get("se") or 0.0))
    return out


def voci_oggi(c):
    """Le voci di una misura principale che e' una media: {nome: valore}, scritte dal check in 'voci'."""
    fm = c.get("file_misura", "results/misura.json")
    if not os.path.exists(fm):
        return {}
    with open(fm, encoding="utf-8") as f:
        m = json.load(f)
    return {k: float(v) for k, v in (m.get("voci") or {}).items() if num(v) is not None}


def stato(c, aff, oggi):
    n, s = serve(c, aff)
    if aff["id"] not in oggi:
        return "non misurata"
    if n is None:
        return "senza numero"
    v, e = oggi[aff["id"]]
    if (v - s * 2 * e - n) * s >= 0:
        return "regge"
    if (v + s * 2 * e - n) * s < 0:
        return "non regge"
    return "incerta"


def quadro(c):
    oggi, p, conta = misure_oggi(c), principale(c), {}
    out = [f"QUADRO DELL'ARTICOLO: {c.get('messaggio', '(messaggio non scritto)')}"]
    for aff in affermazioni(c):
        st = stato(c, aff, oggi)
        conta[st] = conta.get(st, 0) + 1
        n, s = serve(c, aff)
        val = "-" if aff["id"] not in oggi else f"{oggi[aff['id']][0]:g} +/- {oggi[aff['id']][1]:g}"
        req = "?" if n is None else f"{'>=' if s > 0 else '<='} {n:g}"
        out.append(f"  {aff['id']}{'*' if aff['id'] == p else ' '} {st:<12} oggi {val:<16} serve {req:<8} {aff.get('figura', '-'):<7} {aff.get('testo', '')}")
    if not affermazioni(c):
        out.append("  (nessuna affermazione nel contratto)")
    else:
        out.append("  in tutto: " + ", ".join(f"{v} {k}" for k, v in conta.items()) + " (* = misura principale)")
    voci = voci_oggi(c)
    if voci:
        ordinate = sorted(voci.items(), key=lambda kv: -abs(kv[1]))[:8]
        out.append("  voci della misura principale, dalla piu' grande: " + ", ".join(f"{k} {v:g}" for k, v in ordinate))
    return out


def riferimento(c, r, aff):
    """Il valore da battere per un'affermazione: il migliore fra quelli tenuti (sessioni 'misura' e 'avanza'), altrimenti il pavimento per la principale."""
    per_id = {x["id"]: x for x in affermazioni(c)}
    s = serve(c, per_id[aff])[1] if aff in per_id else segno(c)
    tenute_sess = {x["sessione"] for x in r if x["esito"] in ("misura", "avanza")}
    tenute = [(num(q["valore"]), num(q["se"]) or 0.0) for q in leggi_tsv(QST)
              if q["affermazione"] == aff and q["sessione"] in tenute_sess and num(q["valore"]) is not None]
    if tenute:
        return max(tenute, key=lambda t: t[0] * s), s
    pav = num(c.get("pavimento"))
    if aff == principale(c) and pav is not None:
        return (pav, 0.0), s
    return None, s


def regressioni():
    """Affermazioni che reggevano alla penultima chiusura misurata e non reggono piu' all'ultima."""
    q = leggi_tsv(QST)
    sess = []
    for x in q:
        if x["sessione"] not in sess:
            sess.append(x["sessione"])
    if len(sess) < 2:
        return []
    prima = {x["affermazione"]: x["stato"] for x in q if x["sessione"] == sess[-2]}
    dopo = {x["affermazione"]: x["stato"] for x in q if x["sessione"] == sess[-1]}
    return [f"{k} reggeva alla sessione {sess[-2]} e alla sessione {sess[-1]} e' '{v}'" for k, v in dopo.items() if prima.get(k) == "regge" and v != "regge"]


def sessioni_spese(r, fase):
    """Le sessioni che contano nel budget della fase: tutte tranne i riesami."""
    return sum(1 for x in r if x["fase"] == fase and x["esito"] != "riesame")


def dopo_riesame(r):
    """Le righe dopo l'ultimo riesame: i colpi si contano da li'."""
    idx = max((i for i, x in enumerate(r) if x["esito"] == "riesame"), default=-1)
    return r[idx + 1:]


def allarmi(c, r):
    out, fase, s = [], c["fase"], segno(c)
    nella_fase = [x for x in r if x["fase"] == fase]
    bud, spese = c["budget"].get(fase), sessioni_spese(r, fase)
    if bud is None:
        out.append(f"BUDGET non firmato per la fase '{fase}'. Proponilo al PI con il motivo; lo registra 'cancello --budget'.")
    elif spese >= bud:
        n_pro = sum(1 for x in c.get("proroghe", []) if x["fase"] == fase)
        coda = "Serve il cancello del PI, non un'altra sessione." if n_pro == 0 else f"Fase gia' prorogata {n_pro} volta/e: un'altra proroga chiede prima un riesame."
        out.append(f"BUDGET: {spese} sessioni su {bud} nella fase '{fase}'. " + coda)
    pav, tet, sog = num(c.get("pavimento")), num(c.get("tetto")), num(c.get("soglia"))
    fine_sonda = fase not in ("scommessa", "sonda") or (fase == "sonda" and len(nella_fase) >= (bud or 3))
    if fine_sonda and None in (pav, tet, sog):
        out.append("CRITERIO: pavimento, tetto o soglia non misurati. Senza i tre numeri non si costruisce.")
    if None not in (pav, tet, sog):
        if (sog - tet) * s > 0:
            out.append(f"CRITERIO: la soglia ({sog}) sta oltre il tetto ({tet}). Il criterio non e' raggiungibile: riposiziona.")
        if (sog - pav) * s <= 0:
            out.append(f"CRITERIO: la soglia ({sog}) non supera il pavimento ({pav}). Non c'e' niente da dimostrare.")
    affs, oggi, p = affermazioni(c), misure_oggi(c), principale(c)
    if fine_sonda:
        if not affs:
            out.append(f"QUADRO: nessuna affermazione in {CFG}. Senza affermazioni non c'e' un articolo da costruire.")
        mancanti = [a["id"] for a in affs if a["id"] not in oggi]
        senza = [a["id"] for a in affs if serve(c, a)[0] is None]
        if mancanti:
            out.append(f"QUADRO: affermazioni non misurate: {', '.join(mancanti)}. La catena deve produrre in forma grezza tutte le figure del contratto.")
        if senza:
            out.append(f"QUADRO: affermazioni senza il numero che serve: {', '.join(senza)}.")
    for x in regressioni():
        out.append(f"REGRESSIONE: {x}. Ripara o torna indietro prima di ogni altro lavoro.")
    lavoro = [x for x in dopo_riesame(nella_fase) if x["esito"] in ("colpo", "avanza")]
    if len(lavoro) >= 2 and lavoro[-1]["esito"] == lavoro[-2]["esito"] == "colpo":
        if lavoro[-1]["ipotesi"] == lavoro[-2]["ipotesi"]:
            out.append(f"DUE COLPI sulla stessa ipotesi ('{lavoro[-1]['ipotesi']}'). Strada chiusa: niente terza prova. Altra voce del bilancio o riesame.")
        else:
            out.append("DUE COLPI di fila su ipotesi diverse. Rifai il bilancio dell'errore prima di lavorare.")
    if len(lavoro) >= 3 and all(x["esito"] == "colpo" for x in lavoro[-3:]):
        out.append("TRE COLPI di fila. La prossima sessione e' un riesame.")
    ult3 = [x for x in r if x["esito"] not in ("riesame", "contratto")][-3:]
    if sum(1 for x in ult3 if x["esito"] == "infra") > 1:
        out.append("INFRASTRUTTURA: piu' di una sessione su tre senza misura.")
    if p in oggi and None not in (tet, sog):
        v, e = oggi[p]
        if (v - tet) * s >= 0:
            out.append(f"TETTO: la misura principale ({v:g}) e' al livello del tetto ({tet:g}). Altro lavoro su questa misura non si distingue: vai al cancello.")
        elif (v - s * 2 * e - sog) * s >= 0 and fase == "costruzione":
            altre = [a["id"] for a in affs if a["id"] != p and stato(c, a, oggi) != "regge"]
            if altre:
                out.append(f"SOGLIA raggiunta sulla misura principale ({v:g} +/- {e:g} contro {sog:g}), ma non reggono ancora: {', '.join(altre)}. Lavora su quelle, o proponi al PI di restringere il messaggio.")
            else:
                out.append(f"SOGLIA raggiunta ({v:g} +/- {e:g} contro {sog:g}) e tutte le affermazioni reggono. Proponi il cancello verso la conferma.")
    for path, limite in c.get("tetti_parole", {}).items():
        n = parole(path)
        if n is None:
            out.append(f"DOCUMENTO mancante: {path}")
        elif n > limite:
            out.append(f"DOCUMENTO troppo lungo: {path} ha {n} parole, limite {limite}.")
    return out


def cmd_init(_):
    if not os.path.exists(CFG):
        salva(DEFAULT)
    if not os.path.exists(MIS):
        scrivi_tsv(MIS, COLS, [])
    os.makedirs("docs", exist_ok=True)
    ignora = [LAST, "results/check_*.log"]
    presenti = open(".gitignore", encoding="utf-8").read().split() if os.path.exists(".gitignore") else []
    with open(".gitignore", "a", encoding="utf-8") as f:
        f.writelines(x + "\n" for x in ignora if x not in presenti)
    print(f"Creati {CFG} e {MIS}. Compila {CFG} (messaggio, misura, affermazioni), poi scrivi docs/PAPER.md e docs/STATO.md.")


def cmd_quadro(_):
    print("\n".join(quadro(carica())))


def cmd_apri(_):
    c, r = carica(), righe()
    fase = c["fase"]
    n = sessioni_spese(r, fase)
    pro = [x for x in c.get("proroghe", []) if x["fase"] == fase]
    bud = c["budget"].get(fase)
    coda = "" if not pro else f" (prorogata {len(pro)} volta/e, da {pro[0]['da']} a {pro[-1]['a']})"
    print(f"progetto-master | fase {fase} | sessione {n + 1} di {'budget non firmato' if bud is None else bud}")
    print(f"fase: {fase} | sessioni nella fase: {n}/{'non firmato' if bud is None else bud}{coda} | sessioni totali: {len(r)}")
    verso_txt = "il basso" if c.get("verso") == "basso" else "l'alto"
    print(f"misura principale: {c['misura']} (meglio verso {verso_txt}) | pavimento {c.get('pavimento')} | soglia {c.get('soglia')} | tetto {c.get('tetto')}")
    print("\n".join(quadro(c)))
    for x in r[-5:]:
        print(f"  {x['data']} {x['sessione']:<8} {x['fase']:<11} {x.get('affermazione', ''):<3} {x['esito']:<9} valore {x['valore'] or '-'} +/- {x['se'] or '-'} (atteso {x['atteso'] or '-'}) | {x['ipotesi']}")
    sporco = sh("git status --porcelain").stdout.strip()
    print("git: " + ("pulito" if not sporco else f"LAVORO NON COMMITTATO ({len(sporco.splitlines())} file): un'altra chat e' attiva? Fermati e chiedi."))
    if os.path.exists(LAST):
        with open(LAST, encoding="utf-8") as f:
            k = json.load(f)
        print(f"ultimo check: {k['livello']} rc={k['rc']} il {k['quando']}" + ("" if k["impronta"] == impronta() else " (il codice e' cambiato dopo)"))
    al = allarmi(c, r)
    print("ALLARMI: " + ("nessuno" if not al else ""))
    for a in al:
        print("  - " + a)


def cmd_check(a):
    c = carica()
    comando = c["check"][a.livello]
    os.makedirs("results", exist_ok=True)
    log = f"results/check_{a.livello}.log"
    p = sh(comando)
    with open(log, "w", encoding="utf-8") as f:
        f.write(p.stdout + p.stderr)
    with open(LAST, "w", encoding="utf-8") as f:
        json.dump({"livello": a.livello, "rc": p.returncode, "impronta": impronta(), "quando": dt.datetime.now().isoformat(timespec="seconds")}, f)
    coda = (p.stdout + p.stderr).strip().splitlines()[-15:]
    print("\n".join(coda))
    print(f"check {a.livello}: {'VERDE' if p.returncode == 0 else 'ROSSO'} (rc={p.returncode}), log completo in {log}")
    print("\n".join(quadro(c)))
    sys.exit(p.returncode)


def cmd_chiudi(a):
    c, r = carica(), righe()
    fase, per_id = c["fase"], {x["id"]: x for x in affermazioni(c)}
    if any(x["sessione"] == a.sessione for x in r):
        sys.exit(f"CHIUSURA RIFIUTATA: la sessione '{a.sessione}' e' gia' nel registro")
    senza_misura = a.infra or a.riesame or fase == "scommessa"
    rifiuti = []
    if not (a.riesame or fase == "scommessa"):
        if not os.path.exists(LAST):
            rifiuti.append("nessun check eseguito")
        else:
            with open(LAST, encoding="utf-8") as f:
                k = json.load(f)
            if k["rc"] != 0:
                rifiuti.append("l'ultimo check e' rosso")
            elif k["impronta"] != impronta():
                rifiuti.append("il codice e' cambiato dopo l'ultimo check")
    for path, limite in c.get("tetti_parole", {}).items():
        n = parole(path)
        if n is None or n > limite:
            rifiuti.append(f"{path}: {'manca' if n is None else str(n) + ' parole, limite ' + str(limite)}")
    aff, oggi, valore, se = a.affermazione or "", misure_oggi(c), None, None
    if not senza_misura:
        if not aff:
            rifiuti.append("manca --affermazione: ogni sessione serve un'affermazione del contratto")
        elif aff not in per_id:
            rifiuti.append(f"l'affermazione '{aff}' non e' nel contratto ({', '.join(per_id) or 'nessuna'})")
        elif aff not in oggi:
            rifiuti.append(f"manca la misura di {aff}: il check deve scriverla in {c.get('file_misura', 'results/misura.json')} (oppure dichiara --infra)")
        else:
            valore, se = oggi[aff]
    if rifiuti:
        sys.exit("CHIUSURA RIFIUTATA: " + "; ".join(rifiuti))
    if a.riesame:
        esito = "riesame"
    elif fase == "scommessa":
        esito = "contratto"
    elif a.infra:
        esito = "infra"
    else:
        rif, s = riferimento(c, r, aff)
        if fase != "costruzione" or rif is None:
            esito = "misura"
        else:
            guadagno = (valore - rif[0]) * s
            margine = 2 * max(se or 0.0, rif[1])
            raggiunto = a.atteso is None or (valore - a.atteso) * s >= 0
            esito = "avanza" if guadagno > margine and raggiunto else "colpo"
    riga = {"data": oggi_iso(), "sessione": a.sessione, "fase": fase, "affermazione": aff, "ipotesi": a.ipotesi,
            "atteso": "" if a.atteso is None else a.atteso, "valore": "" if valore is None else valore,
            "se": "" if se is None else se, "esito": esito, "nota": a.nota}
    scrivi_tsv(MIS, COLS, [riga])
    if not senza_misura:
        scrivi_tsv(QST, QCOLS, [{"data": oggi_iso(), "sessione": a.sessione, "affermazione": x["id"],
                                 "valore": oggi[x["id"]][0], "se": oggi[x["id"]][1], "stato": stato(c, x, oggi)}
                                for x in affermazioni(c) if x["id"] in oggi])
    print(f"sessione {a.sessione} chiusa: esito '{esito}'. Ora commit e push.")
    print("\n".join(quadro(c)))
    for x in allarmi(c, righe()):
        print("  ALLARME: " + x)


def cmd_cancello(a):
    c, r = carica(), righe()
    rifiuti = []
    if a.fase == "costruzione":
        oggi, affs = misure_oggi(c), affermazioni(c)
        if None in (num(c.get("pavimento")), num(c.get("tetto")), num(c.get("soglia"))):
            rifiuti.append("pavimento, tetto o soglia non misurati")
        if not affs:
            rifiuti.append("nessuna affermazione nel contratto")
        mancanti = [x["id"] for x in affs if x["id"] not in oggi]
        senza = [x["id"] for x in affs if serve(c, x)[0] is None]
        if mancanti:
            rifiuti.append(f"affermazioni non misurate: {', '.join(mancanti)}")
        if senza:
            rifiuti.append(f"affermazioni senza il numero che serve: {', '.join(senza)}")
    bud = a.budget if a.budget is not None else c["budget"].get(a.fase)
    if bud is None:
        rifiuti.append(f"il budget della fase '{a.fase}' non e' firmato: usa --budget N")
    if rifiuti:
        sys.exit("CANCELLO RIFIUTATO: " + "; ".join(rifiuti))
    c.setdefault("cancelli", []).append({"data": oggi_iso(), "da": c["fase"], "a": a.fase, "budget": bud, "consenso": a.consenso, "righe": len(r)})
    c["budget"][a.fase], c["fase"] = bud, a.fase
    salva(c)
    print(f"cancello registrato: fase '{a.fase}', budget {bud} sessioni. Consenso del PI: {a.consenso}")


def cmd_proroga(a):
    c, r = carica(), righe()
    fase = c["fase"]
    prec = [x for x in c.get("proroghe", []) if x["fase"] == fase]
    if prec and not any(x["esito"] == "riesame" for x in r[prec[-1]["righe"]:]):
        sys.exit(f"PROROGA RIFIUTATA: la fase '{fase}' e' gia' stata prorogata il {prec[-1]['data']}. Prima serve un riesame (chiudi --riesame), poi decide il PI.")
    da = c["budget"].get(fase) or 0
    p = principale(c)
    c.setdefault("proroghe", []).append({
        "data": oggi_iso(), "fase": fase, "da": da, "a": da + a.sessioni, "misura_oggi": misure_oggi(c).get(p, [None])[0],
        "motivo": a.motivo, "atteso": a.atteso, "altrimenti": a.altrimenti, "consenso": a.consenso, "righe": len(r)})
    c["budget"][fase] = da + a.sessioni
    salva(c)
    print(f"proroga registrata: fase '{fase}' da {da} a {da + a.sessioni} sessioni. Atteso {a.atteso:g}; altrimenti: {a.altrimenti}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("init").set_defaults(f=cmd_init)
    sub.add_parser("apri").set_defaults(f=cmd_apri)
    sub.add_parser("quadro").set_defaults(f=cmd_quadro)
    k = sub.add_parser("check")
    k.add_argument("livello", choices=["smoke", "full"])
    k.set_defaults(f=cmd_check)
    z = sub.add_parser("chiudi")
    z.add_argument("--sessione", required=True)
    z.add_argument("--ipotesi", required=True)
    z.add_argument("--affermazione")
    z.add_argument("--atteso", type=float)
    z.add_argument("--infra", action="store_true")
    z.add_argument("--riesame", action="store_true")
    z.add_argument("--nota", default="")
    z.set_defaults(f=cmd_chiudi)
    g = sub.add_parser("cancello")
    g.add_argument("--fase", required=True, choices=FASI)
    g.add_argument("--consenso", required=True)
    g.add_argument("--budget", type=int)
    g.set_defaults(f=cmd_cancello)
    e = sub.add_parser("proroga")
    e.add_argument("--sessioni", required=True, type=int)
    e.add_argument("--motivo", required=True)
    e.add_argument("--atteso", required=True, type=float)
    e.add_argument("--altrimenti", required=True)
    e.add_argument("--consenso", required=True)
    e.set_defaults(f=cmd_proroga)
    a = p.parse_args()
    a.f(a)


if __name__ == "__main__":
    main()
