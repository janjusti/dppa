DPPA (Deep Protein Polarity Analyser)

Installation
=================
```bash
pip install dppa
```

Usage
=================
```bash
usage: dppa [-h] -t TARGET -r {csv,xls,all} [--debug]

Analyse all protein alignment .fasta files from a target.

optional arguments:
  -h, --help            show this help message and exit
  -t TARGET, --target TARGET
                        Target .fasta file to be analysed.
  -r {csv,xls,all}, --report {csv,xls,all}
                        Output report file type.
  --debug               Turn debug messages on.
```