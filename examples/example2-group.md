# Example 2 - Group Execution

Open terminal and create a new folder.

```bash
mkdir dppa-examples
```

Get inside `dppa-examples` and download `ExGroup` module.

```bash
cd dppa-examples
wget https://github.com/janjusti/dppa/raw/master/example/ExGroup.py
```

Open Python3 from terminal.

```bash
python3
Python 3.6.7 (default, Oct 22 2018, 11:32:17) 
[GCC 8.2.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> 
```

Import all functions from `ExGroup` module.

```python
>>> from ExGroup import *
```

This example will use four different genes of HPV 16 (E6), extracted from [GenBank](https://www.ncbi.nlm.nih.gov/genbank/) and identified arbitrarily as Group A ([AGM34409.1](https://www.ncbi.nlm.nih.gov/protein/AGM34409.1) and [AGM34408.1](https://www.ncbi.nlm.nih.gov/protein/AGM34408.1)) and Group B ([AIQ82784.1](https://www.ncbi.nlm.nih.gov/protein/AIQ82784.1) and [AIQ82783.1](https://www.ncbi.nlm.nih.gov/protein/AIQ82783.1)).

Generate .fasta files from ID list for each group.

```python
>>> groupA_ids = ['AGM34409.1', 'AGM34408.1']
>>> groupB_ids = ['AIQ82784.1', 'AIQ82783.1']
>>> generate_fasta_from_ids(groupA_ids, 'groupA-unaligned.fasta')
AGM34409.1 successfully fetched.
AGM34408.1 successfully fetched.
groupA-unaligned.fasta sucessfully created.
>>> generate_fasta_from_ids(groupB_ids, 'groupB-unaligned.fasta')
AIQ82784.1 successfully fetched.
AIQ82783.1 successfully fetched.
groupB-unaligned.fasta sucessfully created.
```

Align all sequences with any Multiple Sequence Alignment (MSA) software. In this case, [MUSCLE](https://www.drive5.com/muscle/) is used.

```python
>>> align_via_muscle('groupA-unaligned.fasta', 'groupA.fasta')
Alignment groupA-unaligned.fasta > groupA.fasta done.
>>> align_via_muscle('groupB-unaligned.fasta', 'groupB.fasta')
Alignment groupB-unaligned.fasta > groupB.fasta done.
```

Get consensus sequence for each group and generate a new .fasta file with both of them.

```python
>>> seqA = get_consensus_seq('groupA.fasta', 'groupA_any_sentence')
'groupA_any_sentence' consensus sequence generated (original: groupA.fasta)
>>> seqB = get_consensus_seq('groupB.fasta', 'groupB_any_sentence')
'groupB_any_sentence' consensus sequence generated (original: groupB.fasta)
>>> write_fasta('groupAB-cons.fasta', [seqA, seqB])
groupAB-cons.fasta successfully written.
```

Align all consensus sequences from previous step.

```python
>>> align_via_muscle('groupAB-cons.fasta', 'groups-target.fasta')
Alignment groupAB-cons.fasta > groups-target.fasta done.
```

In this example, "any sentence" is being used as *searchable keyphrase*. If there is any unrecognised amino acid (any symbol [not recognised by IUPAC](https://www.bioinformatics.org/sms2/iupac.html)) in any consensus sequence, DPPA will search for matching amino acids in original files. 

*e.g.*: an unrecognised amino acid (X, position 32) is found in `groupB_any_sentence` sequence. DPPA will look into `groupB.fasta` to check which amino acids X could be.

This feature works for any amount of levels, as long as the sequence's name containing an unrecognised amino acid has a searchable keyphrase settled. If there is not any other searchable sequence to check, this unrecognised amino acid will be considered as "pure" in `alerts`.

Finally, analysis results are obtained from dppa's solver.

```python
>>> results = solver.run('groups-target.fasta', searchable_keyphrase='any sentence')
>>> solver.export(results, 'csv', 'groups')
```

## Output ([More details](../docs/report-exp.md))

### `groups-pols.csv`

| ColNum | PossibleAminos               | PossiblePols                   | PolScore |
|--------|------------------------------|--------------------------------|----------|
| 32     | "\{'D': 0\.75, 'H': 0\.25\}" | "\{'Pp': 0\.25, 'Pn': 0\.75\}" | 7\.45    |
| 90     | "\{'V': 0\.25, 'L': 0\.75\}" | \{'Np': 1\.0\}                 | 0\.0     |

### `groups-alerts.csv`

All gaps and unidentified amino acids are listed in this file.

| SeqName             | ColNum | AlertType |
|---------------------|--------|-----------|
| groupA any sentence | 1      | Pure \-   |
| groupA any sentence | 2      | Pure \-   |
| groupA any sentence | 3      | Pure \-   |
| groupA any sentence | 4      | Pure \-   |
| \.\.\.              | \.\.\. | \.\.\.    |
| groupA any sentence | 155    | Pure \-   |
| groupA any sentence | 156    | Pure \-   |
| groupA any sentence | 157    | Pure \-   |
| groupA any sentence | 158    | Pure \-   |