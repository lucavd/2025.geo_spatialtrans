#!/bin/bash
# R/testing/test_R2b.sh — un comando, PASS/FAIL per asserzione (sessione R2b)
cd "$(dirname "$0")/../.." && .venv/bin/python R/testing/test_R2b.py
