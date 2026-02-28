# Needleman-Wunsch Sequence Alignment

A Python implementation of the Needleman-Wunsch algorithm for global sequence alignment, provided in two versions: a simple version and an advanced version with step-by-step visualization.

---

## Table of Contents

- [Algorithm Description](#algorithm-description)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Scoring Parameters](#scoring-parameters)
- [Example Output](#example-output)
- [Applications](#applications)
- [Contributing](#contributing)
- [License](#license)

---

## Algorithm Description

The Needleman-Wunsch algorithm is a dynamic programming method for global sequence alignment. It works by:

1. Building a score matrix using dynamic programming
2. Tracking back through the matrix to find the optimal path
3. Constructing the final alignment from the traceback path

### Default Scoring Scheme

| Event | Score |
|-------|-------|
| Match | +1 |
| Mismatch | -1 |
| Gap | -2 |

---

## Features

### Simple Version
- User input for two sequences
- Score matrix display with sequence headers
- Optimal alignment with match/mismatch indicators
- Alignment score calculation

### Advanced Version
- All features from the simple version, plus:
- Step-by-step matrix filling visualization
- Trace matrix showing path directions
- Detailed backtrace path listing
- Alignment statistics (matches, mismatches, gaps, identity %)
- Customizable scoring parameters
- Input validation

---

## Installation

### Prerequisites
- Python 3.x
- NumPy

```bash
# Clone the repository
git clone https://github.com/yourusername/needleman-wunsch-alignment.git

# Navigate to the project directory
cd needleman-wunsch-alignment

# Install dependencies
pip install numpy
```

---

## Usage

### Simple Version
```bash
python needleman_wunsch_simple.py
```

### Advanced Version
```bash
python needleman_wunsch_advanced.py
```

---

## Scoring Parameters

Both versions support customizable scoring. The advanced version prompts you to enter custom values at runtime.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `match` | `+1` | Score awarded for matching characters |
| `mismatch` | `-1` | Penalty for mismatching characters |
| `gap` | `-2` | Penalty for inserting a gap |

---

## Example Output

### Simple Version

```
NEEDLEMAN-WUNSCH ALIGNMENT
------------------------------
Enter first sequence: ACCG
Enter second sequence: ACTG

SCORE MATRIX:
      A  C  T  G
    0 -2 -4 -6 -8
  A -2  1 -1 -3 -5
  C -4 -1  2  0 -2
  C -6 -3  0  1 -1
  G -8 -5 -2 -1  2

========================================
RESULTS:
========================================
Sequence 1: ACCG
Sequence 2: ACTG

Alignment score: 2

Optimal alignment:
Seq1: ACCG
Seq2: ACTG
      || *|
```

### Advanced Version

The advanced version additionally outputs a step-by-step matrix fill log, a trace matrix, and a detailed backtrace path, followed by alignment statistics:

```
Statistics:
  Length: 4
  Matches: 3
  Mismatches: 1
  Gaps: 0
  Identity: 75.0%
```

---

## Applications

- **Bioinformatics** — DNA, RNA, and protein sequence alignment
- **Natural Language Processing** — String similarity and comparison
- **Plagiarism Detection** — Document similarity analysis
- **Version Control** — File difference analysis
- **Phylogenetics** — Evolutionary relationship analysis

---

## Contributing

Contributions are welcome! Please feel free to open an issue or submit a pull request.


---

## Acknowledgments

- **Saul B. Needleman and Christian D. Wunsch** (1970) for the original algorithm
- The bioinformatics community for continued development and applications