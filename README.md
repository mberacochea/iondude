# iondude

Implementation of DUDE ([Discrete Universal DEnoiser](https://en.wikipedia.org/wiki/Discrete_Universal_Denoiser)) to correct Ion Torrent reads.

iondude corrects substitutions and homopolymer length errors. It uses channel matrices: the probability that the sequencer reads each base (or homopolymer length) as another one. The estimator calculates these matrices from reads mapped against a reference.

## Repository layout

| Folder | Content |
|---|---|
| `iondude/` | The denoiser, writen in C |
| `iondude/matrices/` | The description of the channel matrices file (JSON) |
| `estimator/` | An error estimator, written Python |

## IonDUDE

### Build

Requirements: a C99 compiler, `make` and [zlib](https://www.zlib.net/).

```bash
cd iondude
make        # builds build/iondude
make test   # runs the iondude unit tests
```

For a debug build: `make clean && make OPT="-O0 -g"`.

### Run

```bash
$ iondude -h

iondude 0.4 (branch main, commit none yet)

DUDE denoiser for Ion Torrent reads.

Usage: iondude -q reads.fastq -m matrices.json -o denoised [options]

Required:
  -q <file>    reads, FASTQ (plain or gzip)
  -m <file>    channel matrices made by the estimator (JSON)
  -o <prefix>  output prefix: <prefix>_subs.fastq and <prefix>_hom.fastq

Options:
  -k <int>     context length for substitutions                  [5]
  -i <int>     context length for homopolymers                   [2]
  -c <int>     correction type                                   [0]
                 0  substitutions only
                 1  homopolymers only
                 2  substitutions, then homopolymers
                 3  homopolymers, then substitutions
  -t           run the self-tests
  -v           print the version
  -h           print this help
```

```bash
./build/iondude -q reads.fastq -m matrices.json -k 5 -i 2 -c 2 -o denoised
```

This writes `denoised_subs.fastq` and `denoised_hom.fastq`.

iondude reads the slot type (read position or quality score), the width of the position ranges and
the max read length from the matrices file. The format is in [iondude/matrices/README.md](iondude/matrices/README.md).

## Estimator

Requirements: [uv](https://docs.astral.sh/uv/).

```bash
cd estimator
make init   # install the dependencies
make test   # runs the tests
```

Estimate the matrices from a SAM/BAM file:

```bash
uv run iondude-estimator mapped.sam -l 200 -r 20 -o matrices.json
```

See [estimator/README.md](estimator/README.md) for the details of the estimation.

## License

Apache License 2.0, see [LICENSE](LICENSE).
The bundled `klib` (MIT), `cJSON` (MIT) and `zlib` (zlib license) keep their own licenses.
