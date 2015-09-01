# HIV WHOLE-GENOME LONGITUDINAL DEEP SEQUENCING
- Authors: F. Zanini, R. Neher (+ other authors for the non-coding part), 2013-2015
- License: MIT

## OVERVIEW
The analysis of HIV whole-genome longitudinal sequences consists of a set of
scripts and modules. There are two levels of analysis:

1. SAMPLE: Reads from each sequenced sample must be cleaned and assigned to
           regions. This step is performed on non-patient samples too. The
           scripts and modules are located in the "sequencing" folder.

2. PATIENT: Cleaned reads are mapped against a patient-specific reference (the
            majotiry sequence of the first time point), allele counts and linkage
            data structures are precomputed. The scripts and modules are located
            in the subfolder "store".

Scripts exist for checking the status of the analysis, both for the SAMPLE and 
for the PATIENT parts; see README in the respective folders.

A standard control for sequencing quality is co-sequencing of PhiX, a plasmid of known
sequence. The analysis pipeline on those reads is located in the README, "phix" folder.


## 1. MAPPING/FILTERING, SAMPLE BY SAMPLE
See [sequencing/README](hivwholeseq/sequencing/README).

![scheme of sample by sample pipeline](scheme.png)

## 2. MAPPING/FILTERING + INTERMEDIATE DATA STRUCTURES, PATIENT BY PATIENT
To be done after 1. See [store/README](hivwholeseq/store/README).

![scheme of patient by patient pipeline](scheme2.png)

## 3. ANALYSIS OF CLEAN READS
See `patients/patients.py` and `patients/samples.py` for general classes, and the scripts in `patients` for typical analyses.

At this stage, the [HIVEVO_access](https://github.com/neherlab/HIVEVO_access) repository is recommended. the codebase is much smaller and personal data on patients/samples should be hidden.

##DIAGRAMS
```
1. START -> TABLE -> SYMLINK -> PREMAP -> DIVIDE -> CONSENSUS -> MAP -> FILTER

2. REFERENCE -> MAP -> FILTER -> ALLELES -> TRAJECTORIES
                                         -> PAIRS
                                         -> LOCAL HAPLOTYPES
```
