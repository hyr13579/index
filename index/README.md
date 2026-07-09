# ZA Primary Enzyme Index & MycoDeg Phylogeny Index Calculator

GAI-based indices for fungal enzyme profile characterization, based on the method by Lai et al. (2023).

## Indices

### Formula (2-1) ZA Primary Enzyme Index
Measures whether a strain's enzyme profile aligns more with dikaryon (MZ) or non-dikaryon (MM) subkingdoms.

```
1. Convert to relative abundance: rel = data / rowSums(data)
2. p_num = count of MZ families present (>0) in sample
   n_num = count of MM families present (>0) in sample
3. numerator   = Σ(rel[j] × p_num/|MZ|) for j∈MZ + β
   denominator = Σ(rel[j] × n_num/|MM|) for j∈MM + β
4. ZA_score = log10(numerator / denominator)
```

### Formula (2-2) MycoDeg Phylogeny Index
Measures alignment with Ascomycota (MA) vs Basidiomycota (MS).

```
Same logic with MA (positive) and MS (negative) groups
```

- β = 1e-3 (smoothing factor to avoid log(0))

## Installation

### 1. Install Python 3.8+

Download from https://www.python.org/downloads/

### 2. Install dependencies

```bash
pip install -r requirements.txt
```

## Usage

### Basic (single file)

```bash
python za_index_calculator.py data/test1.txt
```

### Multiple files

```bash
python za_index_calculator.py data/test1.txt data/test2.txt
```

### Custom directories

```bash
python za_index_calculator.py your_input.txt --data-dir path/to/reference --output-dir path/to/output
```

### Default directories

| Item | Default path |
|------|-------------|
| Reference data | `data/` (next to script) |
| Output | `result/` (next to script) |

## Input Format

### Your input file (tab-separated)
| Column 1 | Column 2+ |
|-----------|-----------|
| Strain name | Gene family counts (e.g. GH3, GH5, AA3, ...) |

Example:
```
strain_name	GH3	GH5	AA3	...
MyStrain	12	5	8	...
```

### Required reference files in `data/`

| File | Description |
|------|-------------|
| `主酶赋值.txt` | Gene family weights for Formula 2-1 (56 families, Estimate = ±0.1) |
| `系统发育赋值.txt` | Gene family weights for Formula 2-2 (17 families, Estimate = ±0.1) |
| `双核.txt` | Dikaryon reference strains (105 strains, with class column A/B/C/D) |
| `非双核.txt` | Non-dikaryon reference strains (50 strains) |

## Output Structure

```
result/
├── strain_name/
│   ├── figures/
│   │   ├── formula_2_1_ranking.png    # ZA index ranking plot
│   │   └── formula_2_2_ranking.png    # MycoDeg phylogeny ranking plot
│   └── strain_name_summary.txt        # Per-strain summary
├── overall_summary.txt                 # Combined text summary
└── overall_summary.csv                 # Combined CSV summary
```

## Output Interpretation

### Formula (2-1) ZA Primary Enzyme Index
- **Score > 0**: Enzyme profile aligns more with dikaryon (MZ) subkingdom
- **Score < 0**: Enzyme profile aligns more with non-dikaryon (MM) subkingdom
- **"Possibly lacking primary enzymes"**: Score at or below the minimum of all dikaryon reference strains

### Formula (2-2) MycoDeg Phylogeny Index
- **Score > 0**: Profile closer to Ascomycota (MA)
- **Score < 0**: Profile closer to Basidiomycota (MS)
- **Closest class**: A/B/C/D phylogenetic class with nearest mean score

## Special Handling: `*_t` Families

The following 7 gene families ending in `_t` are NOT single enzymes but totals of all sub-families:

- `GH55_t`, `PL1_t`, `AA3_t`, `GH51_t`, `GH43_t`, `GH13_t`, `GH16_t`

When a `*_t` column is absent or zero, the calculator sums all columns matching the prefix (e.g., `GH55_t = GH55 + GH55_1 + GH55_2 + ...`).

## Dependencies

- Python >= 3.8
- pandas
- numpy
- matplotlib

## Reference

Lai et al., 2023. GAI (Gene Abundance Index) method for microbial community analysis.
