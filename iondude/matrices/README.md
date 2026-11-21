# Channel matrices

This folder has no matrices now. The matrices will be estimated again with the method of
Song et al. 2017, *Comparison of error correction algorithms for Ion Torrent PGM data:
application to hepatitis B virus*, Sci Rep 7:8106, doi:[10.1038/s41598-017-08139-y](https://doi.org/10.1038/s41598-017-08139-y).

The estimator writes the matrices and iondude reads them. The file format is JSON, described below.
`iondude/unit_test/matrices.json` is a small valid example.

## File format, version 1

```json
{
  "format": "iondude-matrices",
  "format_version": 1,
  "source": {
    "tool": "iondude-estimator",
    "version": "0.4",
    "date": "2026-11-20",
    "input": "B22.bam"
  },
  "estimation": {
    "max_read_length": 200,
    "range": 20,
    "min_mapping_quality": 20
  },
  "slot_type": "position",
  "range": 20,
  "max_length": 200,
  "substitutions": {
    "symbols": ["A", "C", "G", "T"],
    "slots": [
      {"min": 0, "max": 20,
       "counts":  [[9970, 12, 9, 9], ...],
       "channel": [[0.9970, 0.0012, 0.0009, 0.0009], ...],
       "inverse": [[1.0030, -0.0012, -0.0009, -0.0009], ...]}
    ]
  },
  "homopolymers": {
    "symbols": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    "slots": [ ... the same shape, 10 x 10 ... ]
  }
}
```

## Fields

| Field | Required | Rule |
|---|---|---|
| `format` | yes | must be `"iondude-matrices"` |
| `format_version` | yes | must be `1`. A reader stops with an error for other values |
| `source` | yes | `tool`, `version`, `date` (ISO 8601, `YYYY-MM-DD`), `input` (the SAM/BAM path as given). Readers do not use it |
| `estimation` | yes | all the estimator CLI settings, by name. Readers do not use it |
| `slot_type` | yes | `"position"` (a slot is a read position range) or `"quality"` (a slot is one Phred score) |
| `range` | yes | width of each slot. `1` for `quality` |
| `max_length` | yes | the max read length (`position`) or the max quality score + 1 (`quality`) |
| `substitutions`, `homopolymers` | yes | both sections |
| `symbols` | yes | the row and column names: `["A", "C", "G", "T"]` for substitutions, `[1, ..., 10]` (homopolymer length) for homopolymers |
| `slots` | yes | one slot for each range from `0` to `max_length`, ordered by `min`. No gaps, no overlaps |
| `slots[].min`, `slots[].max` | yes | `max - min == range`. The first slot starts at `0` |
| `slots[].counts` | yes | integers, n x n. The errors, and the correct calls on the diagonal. iondude does not use them |
| `slots[].channel` | yes | numbers, n x n. Rows are the true symbol, columns the observed symbol. Each row sums to 1 |
| `slots[].inverse` | yes | numbers, n x n. `channel x inverse` must be the identity matrix |

Notes:

- The last slot can go past `max_length`. For example, `max_length` 150 and `range` 20 give 8 slots, the last one is 140-160.
- iondude uses the last slot for the read positions past `max_length`.
- With `slot_type` `quality`, iondude does not denoise the bases with a quality score of `max_length` or more.
- Each matrix row is on one line, so that a diff of two files is readable.

## Read the file in Python

```python
from estimator.matrices import get_slot, read_matrices

matrices = read_matrices("matrices.json")
channel = get_slot(matrices, "substitutions", 25)["channel"]   # the slot 20-40
```
