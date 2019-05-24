DPPA (Deep Protein Polarity Analyser)
=================

Installation
-----------------
```bash
pip install dppa
```

Usage
----------------
Bash:
```bash
usage: run-dppa [-h] [--reportName REPORTNAME] [--debug] TARGET REPORTTYPE

Analyse all protein alignment .fasta files from a target.

positional arguments:
  TARGET                Target .fasta file to be analysed.
  REPORTTYPE            Output report file type.

optional arguments:
  -h, --help            show this help message and exit
  --reportName REPORTNAME
                        Output report file name.
  --debug               Turn debug messages on.
```
Python:

```python
import dppa

target_name = 'example.fasta'
report_name = 'dppa-report'
report_type = 'xls'

dppa.set_debug_mode(True) # optional

# results[0] -> polarity results dataframe
# results[1] -> alerts dataframe
results = dppa.run(target_name)

dppa.export(report_name, report_type, results)
```