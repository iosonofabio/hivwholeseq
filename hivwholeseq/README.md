# HIV WHOLE-GENOME LONGITUDINAL DEEP SEQUENCING
- Authors: F. Zanini, R. Neher (+ other authors for the non-coding part), 2013-2015
- License: MIT

## OVERVIEW
The analysis of HIV whole-genome longitudinal sequences consists of a set of
scripts and modules. There are two levels of analysis:

1. MAPPING: Reads from each sequenced sample must be cleaned and assigned to
            regions. This step is performed on non-patient samples too. The
            scripts and modules are located in the "sequencing" folder.

2. PATIENT: Cleaned reads are mapped against a patient-specific reference (the
            majotiry sequence of the first time point), allele counts and linkage
            data structures are precomputed, and population genetics is
            studied. The scripts and modules are located in the subfolder "patients"

Scripts exist for checking the status of the analysis, both for the SEQUENCING and 
for the PATIENT parts; see README in the respective folders.

A standard control for sequencing quality is co-sequencing of PhiX, a plasmid of known
sequence. The analysis pipeline on those reads is located in the README, "phix" folder.


## MAPPING PIPELINE FOR HIV SAMPLES
See "sequencing/README".


## POST-MAPPING PIPELINE FOR HIV SAMPLES
At this point, reads are mapped against the patient-specific reference. See README in
the "patients" folder for details on that part of the analysis.

Additional sample-by-sample analyses are possible:

- Calculate useful observables: coverage, allele frequencies, mutations.

- Various popgen measures.


##DIAGRAMS
```
1. MAPPING: START -> TABLE -> SYMLINK -> PREMAP -> DIVIDE -> CONSENSUS -> MAP -> FILTER

2.1 PATIENT: REFERENCE -> MAP -> FILTER -> ALLELES -> TRAJECTORIES
                                        -> PAIRS
                                        -> LOCAL HAPLOTYPES

2.2 CROSS SECTIONAL: CONSENSUS TREE
```
