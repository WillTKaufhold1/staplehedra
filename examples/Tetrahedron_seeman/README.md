# Replicating a Seeman tetrahedron using Staplehedra

## Design

Our objective here was to build a bog-standard DNA tetrahedron with no overhangs, and no single stranded DNA between edges. This is a test case to investigate folding efficiency, melting temperature and hydrodynamic radius of a simple polyhedron. Any degrees of freedom in the choice of DNA have been chosen by using the M13mp18 sequence.

## Experimental Characterization

### Oligo synthesis

DNA oligos are synthesized predominantly via phosphoramidite chemistry. The direction of synthesis is from 3' to 5'. At a coupling efficiency of 99.25%, an oligo of length 90 nts would be expected to have a full length product fraction of 0.9925^90 = 51%. So about half of all DNA oligos will be missing some fraction of their 5' ends. Incomplete formation of the desired polyhedron is therefore expected to some extent...

<img src = https://sfvideo.blob.core.windows.net/sitefinity/images/default-source/product-page-images/custom-dna-rna/coupling-efficiency-curves-dna-oligos-only.png?sfvrsn=f2fd3407_8 width = 600/>
*Figure 1*: IDT's stated coupling efficiency indicate that at 90 nts, an oligo will have a full length fraction of about 50%.

### Molecular simulation

Sequence generation and molecular simulation can be generated with the following command, run in the simulation directory.

```bash
 ~/anaconda3/bin/python ../../../bin6/polymebabyonemoretime.py ../plyfiles/regular_tetrahedron.ply --coarse-steps 1e5 --fine-steps 1e5 --oxdna-steps 1e7 -d out --spacers 0 --bp 20 --nicks 1
```

### DLS

<img src = https://imgur.com/VKQkLFW.png width = 300>
*Figure 2* : DLS (Zetasizer Ultra) data at 588 nM structure concentration. Indication of a  


### A260 curves

### qPCR 
