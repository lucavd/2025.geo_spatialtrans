Lo step risponde a due domande, scelte da Luca all'apertura.

1. **Zone comuni** ([step sulle aree sovrapposte prima di S2.5 (DEC-048)](../docs/decisions/DEC-048.md), voce BL-215): le 5 aree di 1 mm² di R2 dello stesso tessuto si sovrappongono, quindi non sono unità indipendenti.
2. **Poche unità** (voci BL-186 e BL-128): con 5 aree (o 4 finestre contate a mano) l'intervallo bootstrap percentile della media è più stretto del dovuto anche senza zone comuni.

Parole di Luca (chat, 2026-11-08), alla domanda se includere le poche unità in R7: «Sì, R7 copre sovrapposizioni e poche aree (Recommended)».

Fuori dal perimetro, dichiarato: la struttura della sezione fuori dalle aree (esplorazione E-1: le aree di R2 sono dentro un compartimento e il tessuto intorno è più variabile); le validazioni lasciando fuori un'area si misurano (quota di nuclei condivisi), ma non si rifanno.
