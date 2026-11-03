#!/usr/bin/env bash
# tools/R5_G2.sh — G-2: confronta results/R5/R5_md5_before.txt e R5_md5_after.txt (stessi file, stessi md5) -> results/R5/R5_G2.txt
cd "$(dirname "$0")/.."
B=results/R5/R5_md5_before.txt; A=results/R5/R5_md5_after.txt; O=results/R5/R5_G2.txt
nb=$(wc -l < $B); na=$(wc -l < $A); d=$(diff <(sort -k2 $B) <(sort -k2 $A) | grep -c '^[<>]')
{ echo "before: $nb files; after: $na files; differing lines: $d"; git diff --quiet HEAD -- results/R4 results/R2 results/R6 && echo "results/R2, R4, R6 unchanged vs git HEAD" || echo "results/R2, R4, R6 CHANGED vs git HEAD";
  [ "$d" = 0 ] && [ "$nb" = "$na" ] && echo "G-2: PASS ($nb files identical before and after)" || echo "G-2: FAIL"; } > $O
cat $O
