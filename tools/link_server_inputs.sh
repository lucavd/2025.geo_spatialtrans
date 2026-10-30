#!/usr/bin/env bash
# tools/link_server_inputs.sh - S1.6b: in un clone nuovo, collega gli input dei test che non sono in git (file derivati ignorati da .gitignore:
# maschere e nuclei di R2, finestre di R2b, ambienti Python, link ai dati reali) alla copia di lavoro del server del laboratorio.
# Uso (radice del clone): bash tools/link_server_inputs.sh [<copia di lavoro sorgente>]   (predefinita: ~/2025.geo_spatialtrans)
# Esclusi: libreria R (si ricostruisce con tools/restore_library.R), log, cache, file dell'editor. Un percorso gia' presente non si tocca.
# I collegamenti puntano alla copia di lavoro sorgente: i test li leggono soltanto (controllo «nessun file versionato modificato»).
set -u
SRC=${1:-$HOME/2025.geo_spatialtrans}
[ -d "$SRC/.git" ] || { echo "FAIL $SRC non e' una copia di lavoro git"; exit 1; }
[ "$(cd "$SRC" && pwd -P)" = "$(pwd -P)" ] && { echo "INFO questa e' gia' la copia di lavoro sorgente: nulla da collegare"; exit 0; }
n=0
while IFS= read -r p; do
  p=${p%/}
  [ -e "$p" ] || [ -L "$p" ] && continue
  mkdir -p "$(dirname "$p")" && ln -s "$SRC/$p" "$p" && n=$((n + 1)) && echo "link $p"
done < <(git -C "$SRC" status --ignored --porcelain | sed -n 's/^!! //p' |
         grep -vE '^(renv/library/|renv/staging/|logs/|\.Rproj\.user/|\.Rhistory$|\.aider|\.claude/|reports/.*_files/|reports/\.quarto/|\.quarto/)' |
         grep -v '__pycache__')
echo "INFO collegati $n percorsi da $SRC"
