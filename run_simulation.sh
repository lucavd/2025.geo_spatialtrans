#!/bin/bash

# Nome della sessione screen
SCREEN_NAME="r_simulation"

# Percorso del progetto
PROJECT_DIR="$HOME/2025.geo_spatialtrans"

# Script R da eseguire
SCRIPT_PATH="R/daniele/image_simulations_2_claude.R"

# Cartella per i log
LOG_DIR="$PROJECT_DIR/logs"
mkdir -p "$LOG_DIR"

# Nome del file di log con timestamp
LOG_FILE="$LOG_DIR/simulation_output_$(date +%Y%m%d_%H%M%S).log"

# Avvia screen e esegui lo script R catturando l'output
screen -dmS "$SCREEN_NAME" bash -c "cd $PROJECT_DIR && Rscript $SCRIPT_PATH | tee $LOG_FILE"

echo "Script avviato in screen con nome '$SCREEN_NAME'."
echo "Per rientrare nella sessione: screen -r $SCREEN_NAME"
echo "Per vedere il log: tail -f $LOG_FILE"
