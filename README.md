DPPA (Deep Protein Polarity Analyser)
=================

Installation
-----------------
```bash
pip install dppa
```

Usage
----------------
Terminal:
```bash
usage: rundppa [-h] -t TARGET -r {csv,xls,all} [--debug]

Analyse all protein alignment .fasta files from a target.

optional arguments:
  -h, --help            show this help message and exit
  -t TARGET, --target TARGET
                        Target .fasta file to be analysed.
  -r {csv,xls,all}, --report {csv,xls,all}
                        Output report file type.
  --debug               Turn debug messages on.
```
Using Python:

```python
import dppa

target_name = 'example.fasta'
report_name = 'dppa-report'
report_type = 'xls'

dppa.set_debug_mode(True) # optional
[results_dataframe, alerts_dataframe] = dppa.run(target_name)
dppa.export(report_name, report_type, results_dataframe, alerts_dataframe)
```