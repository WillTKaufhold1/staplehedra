# Designing amphiphilic DNA polyhedra which crystallize under annealing

## Introduction

This example concerns designing crystals which are similar to those developed by Ryan Brady in 2017, in [Crystallization of Amphiphilic DNA C-Stars](
https://doi.org/10.1021/acs.nanolett.7b00980). The principle idea of this paper is to use cholesterol conjugated nanostars to form crystals. However, the geometric structure of the placement of cholesterol groups can be achieved using a DNA polyhedra with carefully placed overhangs. These overhangs can then bind to a cholesterol bearing strand when annealing the crystal.

Our expectation is that the melting temperature of the polyhedron will be substantially higher than the melting temperature of the cholesterol bearing strand binding to the polyhedron. Therefore when annealing slowly, the polyhedron will form first.   

## DNA sequences

The DNA sequence of the cholesterol bearing strand is ```5'-GCGTCACGGCGTCG AA 3' - TEG Cholesterol```. So a reasonable overhang sequence is given by ```5'-CGACGCCGTGACGC-3' - tetrahedron```, with the ```5'``` end facing outwards. As for the geometry of the tetrahedron, we have chosen three kinds:

1. A regular tetrahedron, side length 21 bp.
2. A tetrahedron whose geometry is consistent with the proposed nanostar structure in the paper. This tetrahedrdon has side lengths which are in the ratio 2:SQRT(2). Since the side length of one is 21, the other should be of length 14.8. But we'd like both of these if possible to be close to integer multiples of 10.5, so as to match the helicity of the DNA. One approach is to lengthen the longer side to 31 bp, so that the other should be 22, and we come back fairly close to obeying our helicity condition.
3. A tetrahedron which locates the cholesterol bearing strands at locations of where they should be in the crystal. This is very slightly different to our previous criterion, in that here, we also need to recognize that the double stranded DNA in the overhang contributes to the location of the cholesterol group...

## Running the simulations 

To run the regular tetrahedron simulation, which will also generate DNA strands, I can navigate to the tetrahedron directory, and run the catchy command:

```bash
 ~/anaconda3/bin/python ~/Dropbox/PROJECTS/Imperial_projects/poly/dev/bin6/polymebabyonemoretime.py ../plyfiles/regular_tetrahedron.ply -d out --coarse-steps 1e5 --fine-steps 1e5 --oxdna-steps 1e5 --spacers 0 --overhangfile tetrahedron_sequence.overhang --nicks 0 --bp 21
```
But I want to do some in silico optimization, just to minimize unusual strains on the molecule as a consequence of side lengths not quite matching the DNA helicity.









