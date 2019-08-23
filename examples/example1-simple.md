# Example 1 - Simple Execution

Open terminal and create a new folder.

```bash
mkdir dppa-examples
```

Get inside `dppa-examples` and download `ExSimple` module.

```bash
cd dppa-examples
wget https://github.com/janjusti/dppa/raw/master/examples/ExSimple.py
```

Open Python3 from terminal.

```bash
python3
Python 3.6.7 (default, Oct 22 2018, 11:32:17) 
[GCC 8.2.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> 
```

Import all functions from `ExSimple` module.

```python
>>> from ExSimple import *
```

This example will use four different genes of HPV 16 (E6), extracted from [GenBank](https://www.ncbi.nlm.nih.gov/genbank/) and identified as [AGM34409.1](https://www.ncbi.nlm.nih.gov/protein/AGM34409.1), [AGM34408.1](https://www.ncbi.nlm.nih.gov/protein/AGM34408.1), [AIQ82784.1](https://www.ncbi.nlm.nih.gov/protein/AIQ82784.1) and [AIQ82783.1](https://www.ncbi.nlm.nih.gov/protein/AIQ82783.1).

Generate .fasta file from ID list.

```python
>>> id_list = ['AGM34409.1', 'AGM34408.1', 'AIQ82784.1', 'AIQ82783.1']
>>> generate_fasta_from_ids(id_list, 'simple-unaligned.fasta')
KC662553.1 successfully fetched.
KC662552.1 successfully fetched.
KM058604.1 successfully fetched.
KM058603.1 successfully fetched.
simple-unaligned.fasta sucessfully created.
```

Align all sequences with any Multiple Sequence Alignment (MSA) software. In this case, [MUSCLE](https://www.drive5.com/muscle/) is used.

```python
>>> align_via_muscle('simple-unaligned.fasta', 'simple.fasta')
Alignment simple-unaligned.fasta > simple.fasta done.
```

Finally, analysis results are obtained from DPPA's solver.

```python
>>> results = solver.run('simple.fasta')
>>> solver.export(results, 'csv', 'simple')
```

## Output ([More details](../docs/report-exp.md))

### `simple-pols.csv`

| ColNum | PossibleAminos               | PossiblePols                   | PolScore |
|--------|------------------------------|--------------------------------|----------|
| 32     | "\{'H': 0\.25, 'D': 0\.75\}" | "\{'Pp': 0\.25, 'Pn': 0\.75\}" | 7\.45    |
| 90     | "\{'L': 0\.75, 'V': 0\.25\}" | \{'Np': 1\.0\}                 | 0\.0     |


### `simple-alerts.csv`

All gaps and unidentified amino acids are listed in this file.

| SeqName                                                    | ColNum | AlertType |
|------------------------------------------------------------|--------|-----------|
| "AGM34408\.1 E6, partial \[Human papillomavirus type 16\]" | 1      | Pure \-   |
| "AGM34409\.1 E6, partial \[Human papillomavirus type 16\]" | 1      | Pure \-   |
| "AGM34408\.1 E6, partial \[Human papillomavirus type 16\]" | 2      | Pure \-   |
| "AGM34409\.1 E6, partial \[Human papillomavirus type 16\]" | 2      | Pure \-   |
| \.\.\.                                                     | \.\.\. | \.\.\.    |
| "AGM34408\.1 E6, partial \[Human papillomavirus type 16\]" | 157    | Pure \-   |
| "AGM34409\.1 E6, partial \[Human papillomavirus type 16\]" | 157    | Pure \-   |
| "AGM34408\.1 E6, partial \[Human papillomavirus type 16\]" | 158    | Pure \-   |
| "AGM34409\.1 E6, partial \[Human papillomavirus type 16\]" | 158    | Pure \-   |