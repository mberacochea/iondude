===========================
Validation of the Estimator
===========================

**26/02/2016**

*Original reads and reference files: https://github.com/nickloman/benchtop-sequencing-comparison*

Data used for this purpose:

- Reference produced by Loman et al. 2012 (/assemblies/newbler)
    - Scaffolds from isolate E. coli O104:H4
        - 454Scafolds.fna
        - 454Scafolds.qual
        - 454Scafolds.txt
- Reads produced by Loman et al. 2012 (/reads)
    - IonTorrent
        - IonTorrent-1.fastq

Reads where mapped with bwa 0.7.12-r1039::

 bwa index 454Scaffolds.fna
 bwa mem 454Scaffolds.fna iontorrent_1.fastq > loman_ref_iontorrent1.sam

Estimation of the errors with Estimator::

 python run.py validation/loman_ref_iontorrent1.sam 2> error.log

There are many reads with substitutions for the same base, error
In order to correct this, 2 reads that presented this error were used to build 2 more tests for the estimator.
Reads:
- 9TBRU:4:92
- 9TBRU:4:774

Those reads were extracted with seqtk::

 seqtk subseq IonTorrent-1.fastq reads_error_to_tests.txt > reads_error.fq

Then were mapped with bwa and the changes counter manually with Tablet::

 bwa mem 454Scaffolds.fna estimator_correction/reads_error.fq  > estimator_correction/mapped.sam

**27/02/2016**

I have replaced the code to count the substitutions, the code uses pysam built-code functionality.

Following the steps from 26/2 to create some more tests and validate the output.

Final validation against the IonTorrent reads set 1::

 python run.py validation/loman_ref_iontorrent1.sam