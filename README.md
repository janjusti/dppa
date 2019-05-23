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
usage: rundppa [-h] --target TARGET --reportType {csv,xls,all}
               [--reportName REPORTNAME] [--debug]

Analyse all protein alignment .fasta files from a target.

optional arguments:
  -h, --help            show this help message and exit
  --target TARGET       Target .fasta file to be analysed.
  --reportType {csv,xls,all}
                        Output report file type.
  --reportName REPORTNAME
                        Output report file name.
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