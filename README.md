# Live Variant Caller

This repository is a live variant caller that let's you run variant calling in real-time. 

## Installation
1. Clone repo and enter directory
2. Install required dependencies with ```pipenv install -r requirements.txt```
2a. If needed - install Pipenv.
3. A) Run variant caller only once with ```python -m src/main.py```
3. B) Run variant caller with whole infrastructure (client/server) with ``` scripts/run_full.sh [YOUR PATH HERE]``` - where you add your preferred path - it should be your sequencing output path.

#### Make sure you have a reference FASTA or FASTQ for your genome. you can also provide a BAM-file, but this is optional.

#### Alternatively, you can use ```covid-spings-pipeline``` simulation as the variant callers output. Refer to the Installation section.
