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

The output of the program will be these files:
- substitutions_frecuency.csv
- homopolymers_frecuency.csv
- log.txt

_Substitutions_frecuency.csv_

 - | A       | C       | G       | T
---|---------|---------|---------|---------
 A | -       | f(A->C) | f(A->G) | f(A->T)
 C | f(C->A) | -       | f(C->G) | f(C->T)
 G | f(G->A) | f(G->C) | -       | f(G->T)
 T | f(T->A) | f(T->C) | f(T->G) | -

Where f(A->C) is the frequency with which a base A changes to C.

_homopolymers_frecuency.csv_

 - | 1       | 2       | 3
---|---------|---------|---------
 1 | -       | f(1->2) | f(1->3)
 2 | f(2->1  | -       | f(2->3)
 3 | f(3->1) | f(3->2) | -

Where f(1->2) is the frequency with which a homopolymer of length 1 changes to a homopolymer of length 2.
*It is not limited to 3, the max considered is 10 because the probability of longer homopolymers is very low.*

_log.txt_

Extra information. To be defined.
