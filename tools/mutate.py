"""tools/mutate.py - S1.6 (BL-155): mutanti come copie modificate del codice, senza interruttori nel codice di produzione.
Una specifica (file Python) definisce MUTANTS = {nome: [(file relativo, testo esatto, testo sostitutivo), ...]}; ogni testo deve comparire
esattamente una volta nel file originale, altrimenti il mutante e' rifiutato (exit 3: sostituzione non applicata).
Uso: python3 -B tools/mutate.py <specifica.py> <mutante> <cartella di destinazione>   (copia R/*.R e i file toccati, poi applica)
     python3 -B tools/mutate.py <specifica.py> --list"""
import sys, os, shutil, glob, runpy
spec = runpy.run_path(sys.argv[1])["MUTANTS"]
if sys.argv[2] == "--list": print(" ".join(spec)); sys.exit(0)
name, dest = sys.argv[2], sys.argv[3]
os.makedirs(os.path.join(dest, "R"), exist_ok=True)
for f in glob.glob("R/*.R"): shutil.copy(f, os.path.join(dest, "R"))
for f, old, new in spec[name]:
    p = os.path.join(dest, f); os.makedirs(os.path.dirname(p), exist_ok=True)
    if not os.path.exists(p): shutil.copy(f, p)
    s = open(p).read(); n = s.count(old)
    if n != 1: print(f"mutante {name}: testo trovato {n} volte in {f}"); sys.exit(3)
    open(p, "w").write(s.replace(old, new, 1))
print(f"mutante {name}: {len(spec[name])} sostituzioni in {dest}")
