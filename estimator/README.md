## Estimator

[SAM/BAM spaces](https://samtools.github.io/hts-specs/SAMv1.pdf)
[pysam](https://pysam.readthedocs.org/en/latest/index.html)

# Estimator

It is a program that is part of Iondude, written in Python.

The only requirement is pysam.

Its function is to estimate sequencing errors from reads mapped against a reference.

It uses a .sam file as input.

Short description of the SAM format. For more information see the [specification](https://samtools.github.io/hts-specs/SAMv1.pdf).

## SAM

It is used to store aligned/mapped sequences. SAM means: Sequence Alignment/Map format.

It is a TAB-delimited text file made of:
- header (optional)
- alignment

The header lines start with @, the alignment lines do not.
Each alignment has 11 mandatory fields.

These fields are

 Name   | Description
--------|------------------------------------------------------
 QNAME  | Query name
 FLAG   | Bitwise FLAG
 RNAME  | Reference name
 POS    | Mapping position, from the left, starting at 1
 MAPQ   | Mapping quality (a value of 255 means not available)
 CIGAR  | CIGAR string
 RNEXT  | Name of the next mapped read (of the mate pair)
 PNEXT  | Position of the next mapped read (of the mate pair)
 TLEN   | Total length of the reads, for paired reads
 SEQ    | Read sequence
 QUAL   | Read quality

*CIGAR*:
It is a string encoded with the format \*|([0-9]+[MIDNSHPX=])+

 Value | Description
-------|-------------------------------------------------
 M     | Match (match or mismatch)
 I     | Insertion to the reference
 D     | Deletion from the reference
 N     | Skipped region
 S     | Soft-clip
 H     | Hard-clip
 P     | Padding (silent deletion from padded reference)
 =     | Exact match
 X     | Mismatch

_Example_

10M2I11M1D5M

- 10 matches
- 2 insertions
- 11 matches
- 1 deletion
- 5 matches

SAM also defines a set of optional TAGs (TAG:TYPE:VALUE). The ones of interest for the
Estimator are:

*MD*: String of the mismatches. It is a string that matches the following regex:
[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*

This tag only stores information about the reference. This makes it possible to find SNPs/indels without
the reference.

_Example_

10A0A5^AC6

From the left of the mapping we have:
- 10 matches followed by a substitution of an A by the base that is at that
position in the read.
- At 11 there is another A in the reference, the 0 between the A's indicates that they are consecutive
- 5 matches, from position 12 (where the A is)
- 2 bases that are in the reference and not in the read, the dinucleotide AC
- ends with 6 matches

_Example_

READ: TGCTCGGCCATCAGCATCGGGGCGAATTCGATCACATGCGTTTCGAC-ACGCGTTG ... AGC<br>
                                                      -  - - - 319 -<br>
REF : TGCTCGGCCATCAGCATCGGGGCGAATTCGATCACATGCGTTTCGACGCCGAGGTT ... TCG <br>

MD  : 47^G0C2A1G1T319T2

C=>A , A=>C. G=>T, T=>G and T=>A

*NM*: Edit distance.

*Note* The insertions to the reference are not stored in MD, because they represent a loss of information about the reference.

When the MD TAG is used, the position that it stores does not take insertions into account, so if e (this note is incomplete in the original text)

## Model

The estimator starts from a mapping with BWA (indel tolerant).

The .sam file that BWA produces is processed with pysam.

The output is one JSON file with the channel matrices, for each range of read positions:

- substitutions: 4 x 4, `[true base][observed base]`
- homopolymers: 10 x 10, `[true length][observed length]`. The max length is 10, because the
  probability of longer homopolymers is very low

Each matrix has the counts, the frequencies (`channel`) and the inverted frequencies (`inverse`).
The file also records the estimator version, the input file and the settings of the estimation.
The format is in [iondude/matrices/README.md](../iondude/matrices/README.md).
Use `estimator.matrices.read_matrices` to read the file in Python.
