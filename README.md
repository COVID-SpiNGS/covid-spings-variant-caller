# Live Variant Caller

## Installation
1. Clone repo and enter directory
2. Install required dependencies with ```pipenv install -r requirements.txt```
3. A) Run only variant caller with ```python -m src/main.py```
3. B) Run variant caller with whole infrastructure (client/server) with ``` scripts/run_full.sh```




## Notizen

POS = Position in der Referenz Sequenz  \\
REF = Base in der Referenz Sequenz (bisher noch nicht rausgefunden wie ich die bekomme) \\
DEPTH = Anzahl beruecksichtiger Reads \\
A = Anzahl As \\
A% = Prozentanteil von As im Verhaeltnis zur Anzahl der Reads. Oder Ergebnis der Formel zur Berechnung der Genetype lilelihood aus dem Paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3198575/ \\


...Rest so wie bei A un A%


