#### `ColNum`

Position of mutated amino acid, starting from 1. For example, `ColNum` 16 would be related to every amino acid in position 16 in .fasta target file.

#### `PossibleAminos`

Percentage values dictionary showing the proportion of each amino acid on `ColNum`. Only recognized amino acids are listed in `PossibleAminos`; gaps and unidentified amino acids are also considered, but not shown (see `alerts` instead).

#### `PossiblePols`

Percentage values dictionary showing the proportion of polarity from each amino acid listed on `PossibleAminos`. Score values are used to calculate `PolScore`.

| Type | Name               | Score |
|------|--------------------|-------|
| Pp   | Positive Polar     | 1     |
| Pn   | Negative Polar     | 1     |
| Nc   | Neutral Polar      | 0     |
| Np   | Neutral Non\-polar | 0     |

#### `PolScore`

Polarity score is calculated using score's sum of all existent polarities on `PossiblePols` (`SumScores`), length of `PossiblePols` size (`PolListSize`) and smallest percentage value between all existent polarities (`MinPolsPerc`).

```math
PolScore = SumScores*3 + PolListSize*0.6 + MinPolsPerc
```