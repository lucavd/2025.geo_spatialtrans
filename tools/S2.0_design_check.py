"""tools/S2.0_design_check.py — S2.0, cancello G-3 (docs/steps/S2.0.md, afef06b): il documento dello Step 2 copre tutto ed e' tracciabile.
1. ogni voce dei dubbi ricevuti, ogni elemento dello Step 2 della revisione di Daniele (3.2.A-E, 4.1-4.5) e ogni decisione che tocca lo Step 2
   (DEC-034, DEC-035, DEC-037 e le decisioni con step: "S2.0") e' nominata in una sezione del documento;
2. ogni riga della tabella dei parametri (sezione '## Tabella dei parametri') ha una fonte che esiste: riga B-nnn di BIO_REFERENCES, file DEC,
   file del repository, oppure 'da misurare in <step>' con un id di step;
3. due mutanti (una sezione tolta, un B-nnn inventato) sono rilevati.
Uso: python3 -B tools/S2.0_design_check.py [--doc <file>]  (PASS/FAIL per asserzione)"""
import os, re, sys, glob
DOC = "docs/superpowers/specs/2026-11-03-step2-visium-grid-design.md"
BL = ["BL-197", "BL-192", "BL-194", "BL-195", "BL-001", "BL-005", "BL-010", "BL-146", "BL-201", "BL-077", "BL-086", "BL-087",
      "BL-045", "BL-057", "BL-072", "BL-144"]
DANIELE = ["3.2.A", "3.2.B", "3.2.C", "3.2.D", "3.2.E", "4.1", "4.2", "4.3", "4.4", "4.5"]
RES = []

def check(name, ok, val=""):
    RES.append(ok); print(f"{'PASS' if ok else 'FAIL'} {name} {val}")

def decs_s20():
    out = ["DEC-034", "DEC-035", "DEC-037"]
    for f in sorted(glob.glob("docs/decisions/DEC-*.md")):
        if re.search(r'^step:\s*"?S2\.0"?\s*$', open(f).read(), re.M):
            out.append(os.path.basename(f)[:-3])
    return out

def sections(text):
    parts = re.split(r"(?m)^(#{1,4} .*)$", text)
    return [(parts[i], parts[i + 1]) for i in range(1, len(parts) - 1, 2)]

def bio_ids():
    return set(re.findall(r"(?m)^\| (B-\d{3}) \|", open("docs/BIO_REFERENCES.md").read()))

STEP_RE = re.compile(r"\b(S\d(?:\.\d)?[a-z]?|R\d[a-z]?|F3\.\d|DOC1)\b")
def source_ok(cell, bio):
    ok_any, bad = False, []
    for b in re.findall(r"\bB-\d{3}\b", cell):
        if b in bio: ok_any = True
        else: bad.append(b)
    for d in re.findall(r"\bDEC-\d{3}\b", cell):
        if os.path.exists(f"docs/decisions/{d}.md"): ok_any = True
        else: bad.append(d)
    for p in re.findall(r"`([^`]+)`", cell):
        p0 = p.split(":")[0].split(" ")[0]
        if "/" in p0:
            if os.path.exists(p0): ok_any = True
            else: bad.append(p0)
    if re.search(r"da misurare in", cell) and STEP_RE.search(cell.split("da misurare in", 1)[1]):
        ok_any = True
    return ok_any and not bad, bad

def run(text, label=""):
    n0 = len(RES)
    secs = sections(text); body = "\n".join(h + b for h, b in secs)
    for t in BL + DANIELE + decs_s20():
        hit = [h.strip() for h, b in secs if re.search(re.escape(t) + r"(?![\d])", h + b)]
        check(f"{label}cov_{t}", len(hit) > 0, f"sezioni={len(hit)}")
    # tabella dei parametri
    m = re.search(r"(?ms)^## Tabella dei parametri.*?(?=^## )", text + "\n## FINE\n")
    check(f"{label}param_table_present", m is not None)
    if m:
        rows = [l for l in m.group(0).splitlines() if l.startswith("|") and not re.match(r"^\|\s*-", l)]
        hdr = [c.strip().lower() for c in rows[0].strip("|").split("|")]
        j = [i for i, h in enumerate(hdr) if h.startswith("fonte")]
        check(f"{label}param_table_has_source_column", len(j) == 1, str(hdr))
        bio = bio_ids(); nbad = 0
        for l in rows[1:]:
            cells = [c.strip() for c in l.strip("|").split("|")]
            ok, bad = source_ok(cells[j[0]], bio) if j else (False, [])
            if not ok:
                nbad += 1; print(f"   {label}riga senza fonte valida: {cells[0][:60]} | {cells[j[0]][:80] if j else ''} | {bad}")
        check(f"{label}param_rows_sourced", nbad == 0, f"righe={len(rows) - 1} senza_fonte={nbad}")
    return all(RES[n0:])

if __name__ == "__main__":
    doc = sys.argv[sys.argv.index("--doc") + 1] if "--doc" in sys.argv else DOC
    text = open(doc).read()
    ok = run(text)
    # mutanti
    secs = sections(text)
    tgt = [i for i, (h, b) in enumerate(secs) if "BL-197" in h + b]
    mut1 = "".join(h + b for i, (h, b) in enumerate(secs) if i not in tgt) if tgt else text
    n0 = len(RES); r1 = run(mut1, "M1_"); del RES[n0:]
    mut2 = re.sub(r"\bB-(\d{3})\b", "B-999", text, count=0) if False else None
    m = re.search(r"(?ms)^## Tabella dei parametri.*?(?=^## )", text)
    mut2 = text[:m.start()] + m.group(0).replace("B-3", "B-9", 1) + text[m.end():] if m else text
    n0 = len(RES); r2 = run(mut2, "M2_"); del RES[n0:]
    check("mutant_M1_section_removed_detected", not r1)
    check("mutant_M2_invented_B_detected", not r2)
    print(f"S2.0_design_check: {sum(RES)}/{len(RES)} PASS")
