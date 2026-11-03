#!/usr/bin/env bash
# tools/R5_repro2.sh <A> <B> — G-4 (dopo la revisione): confronta per md5 l'unione delle tabelle (calcolo e post hoc), kernels.json, params.json
# di due esecuzioni complete; un file presente in una sola cartella e' MISSING. Scrive results/R5/R5_G4.txt.
set -uo pipefail
cd "$(dirname "$0")/.."
A=$1; B=$2; OUTF=results/R5/R5_G4.txt
{ for f in $( (cd "$A/tables" && ls *.csv; cd "$B/tables" && ls *.csv) | sort -u) kernels.json params.json; do
    pa=$A/tables/$f; pb=$B/tables/$f; [ -f "$pa" ] || [ -f "$pb" ] || { pa=$A/$f; pb=$B/$f; }
    if [ -f "$pa" ] && [ -f "$pb" ]; then ma=$(md5sum < "$pa" | cut -c1-32); mb=$(md5sum < "$pb" | cut -c1-32)
        echo "$f $ma $mb $([ "$ma" = "$mb" ] && echo identical || echo DIFFERENT)"
    else echo "$f MISSING"; fi; done; } > "$OUTF"
d=$(grep -c 'DIFFERENT\|MISSING' "$OUTF" || true); t=$(wc -l < "$OUTF")
echo "G-4: $((t - d))/$t files identical ($A vs $B)" >> "$OUTF"; tail -1 "$OUTF"
