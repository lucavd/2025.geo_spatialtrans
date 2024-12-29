#!/bin/bash
TARGET_SESSION="3231775.model_MMM"

echo "Waiting for session $TARGET_SESSION to finish..."
while screen -ls | grep -q "$TARGET_SESSION"; do
    sleep 60
done

echo "Session model_MMM finished, starting new targets pipeline..."
screen -dm -S "next_model" bash -c "cd ~/Projects/2025.geo_spatialtrans && R -e 'targets::tar_make()'"
