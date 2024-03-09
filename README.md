# Live Variant Caller

This repository is a live variant caller that let's you run variant calling in real-time. 

## Installation
0. If needed - install Pipenv.
1. Clone repo and enter directory
2. Install required dependencies with ```pipenv install -r requirements.txt```
3. A) Run variant caller only once with ```python -m src/main.py```
3. B) Run live variant caller with whole infrastructure (client/server) with ``` scripts/run_full.sh [YOUR PATH HERE]``` - where you add your preferred path - it should be your sequencing output path.

For the input, make sure you have a reference FASTA or FASTQ for your genome. you _can_ also provide a BAM-file, but this is optional.

#### Alternatively, you can use ```covid-spings-pipeline``` simulation as the variant callers output. Just provide the path to your simulation output directory.
