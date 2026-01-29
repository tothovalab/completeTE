# CompleteTE Copilot Instructions

## Project Overview

**completeTE** is a bioinformatics tool that curates complete proviral retrotransposon elements by identifying internal regions (INT) flanked by matching Long Terminal Repeats (LTRs) from UCSC RepeatMasker annotations.

### Core Concept
The tool searches for proviral structures in genomic data following the **LTR–INT–LTR** pattern:
- **LTR (upstream)**: Long Terminal Repeat matching a target subfamily
- **INT (internal)**: Internal region between LTRs  
- **LTR (downstream)**: Long Terminal Repeat matching a target subfamily on the same strand

## Architecture

### Single Entry Point: `scripts/curate_complete_elements.py`

The codebase is intentionally minimal with one primary script handling the complete workflow:

1. **Data Loading** (`load_ucsc_rmsk`): Loads UCSC RepeatMasker tab-delimited annotation files, filters repetitive noise (Simple_repeat, Low_complexity, Satellite), and normalizes to BED-like DataFrame
2. **Element Matching**: For each INT element, searches for flanking LTRs within a tolerance window (default 10,000 bp) on the same strand
3. **Proximity Scoring**: Selects nearest LTR on each side when multiple candidates exist
4. **Aggregation**: Groups overlapping INT regions by matching LTR boundaries
5. **Strand Normalization**: Converts 5' and 3' coordinates based on strand (+/-)
6. **Output**: Exports complete elements as CSV with columns: chr, strand, ltr_5_start, ltr_5_end, int_start, int_end, ltr_3_start, ltr_3_end

## Key Patterns & Conventions

### Pandas-First Data Processing
- All genomic regions handled as pandas DataFrames with columns: chr, start, end, subfamily, class, strand
- **Strand symmetry**: Always filter for valid strands (`+ or -`) early; INT/LTR must share same strand
- **Chromosome sorting**: Use `sort_chromosomes()` function to standardize genomic order (autosomes 1-22, X, Y, MT, then unplaced contigs)

### Genomic Coordinate Conventions
- **0-based half-open intervals**: Coordinates follow BED format (start inclusive, end exclusive)
- **Strand-aware operations**: Upstream/downstream LTRs depend on strand direction:
  - **Forward strand (+)**: ltr_up is upstream, ltr_down is downstream
  - **Reverse strand (−)**: positions reverse, ltr_down becomes 5', ltr_up becomes 3'
- **Tolerance parameter**: Controls acceptable gap/overlap between elements (default 10,000 bp); controls genomic "fuzziness"

### Column Naming Scheme
- Input from RepeatMasker: bin, swScore, milliDiv, milliDel, milliIns, chr, start, end, genoLeft, strand, subfamily, family, class, repStart, repEnd, repLeft, id
- Intermediate: ltr_up/ltr_down (relative to sequence direction)
- Output: ltr_5/ltr_3 (biological polarity) + int_start, int_end

## Critical Command-Line Interface

```bash
python curate_complete_elements.py \
  -r <UCSC_RMSK_FILE> \
  -i <INT_SUBFAMILY> \
  -l <LTR_SUBFAMILY> \
  -t <TOLERANCE_BP>  # optional, default=10000
```

**Arguments are required** (except tolerance); script will fail silently or with cryptic pandas errors if missing output file parameter.

## Integration Notes

### Input Data Requirements
- RepeatMasker output must be UCSC-format (tab-delimited with standard column order)
- Subfamily names must exactly match RepeatMasker subfamily annotations (e.g., "HERV-K", "ERVK")

### Missing Pieces (Implementation Opportunities)
- ⚠️ **Output file parameter not implemented**: Script references undefined `out_file` variable at end—must add `-o/--output` argument
- ⚠️ **No logging mechanism**: Print statements should route to log file in output directory (see commented TODO)
- ⚠️ **Missing import**: `re` module used in `sort_chromosomes()` but not imported
- ⚠️ **No numpy import**: Script uses `np.where()` without importing numpy

## Common Development Tasks

### Running the Tool
```bash
# Minimal example (will error due to missing output param):
python scripts/curate_complete_elements.py -r rmsk.txt -i HERVK-int -l HERVK-LTR

# After fixes:
python scripts/curate_complete_elements.py -r rmsk.txt -i HERVK-int -l HERVK-LTR -o output.csv -t 5000
```

### Testing Approach
- RepeatMasker output can be validated by checking: (1) correct column count, (2) strand values are only +/-, (3) chromosome format (chr1, chrX, etc.)
- Element matching logic depends on sorted dataframes—verify `sort_chromosomes()` produces correct genomic order
- Check edge cases: overlapping LTRs, zero-width intervals, elements at chromosome boundaries

## Project Structure
```
completeTE/
├── .github/
│   └── copilot-instructions.md
├── scripts/
│   └── curate_complete_elements.py    # Single entry point
├── examples/                          # Currently empty
├── README.md                          # Brief project summary
└── .git/
```

## Immediate Next Steps for Contributors
1. Fix missing imports (re, numpy) and undefined variables (out_file)
2. Add `-o/--output` argument to argument parser
3. Implement log file output (currently just print statements)
4. Add validation for subfamily name matching against input data
5. Consider test cases with small example RepeatMasker subsets
