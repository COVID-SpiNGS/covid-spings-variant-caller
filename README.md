# Live Variant Caller

This repository is a live variant caller that let's you run variant calling in real-time.

## Installation
1. Clone repo and enter directory
2. Install required dependencies with ```pipenv install -r requirements.txt```
3. A) Run only variant caller with ```python -m src/main.py```
3. B) Run variant caller with whole infrastructure (client/server) with ``` scripts/run_full.sh```

#### Make sure, you have an in- and output directory in the root of the project. The input dir should contain a reference FASTA or FASTQ as well as an binary alignment file (BAM).

#### Alternatively, you can use ```covid-spings-pipeline``` simulation as the variant callers output. Just make sure to provide the correct paths.


## Notizen

POS = Position in der Referenz Sequenz  \\
REF = Base in der Referenz Sequenz (bisher noch nicht rausgefunden wie ich die bekomme) \\
DEPTH = Anzahl beruecksichtiger Reads \\
A = Anzahl As \\
A% = Prozentanteil von As im Verhaeltnis zur Anzahl der Reads. Oder Ergebnis der Formel zur Berechnung der Genetype lilelihood aus dem Paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3198575/ \\


...Rest so wie bei A un A%


