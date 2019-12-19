# Seeman vs Turberfield nick patterns

## Intro

Since synthetic oligonucleotides are linear (as opposed to circular) molecules, a polyhedron which is made of these oligos must have nicks. 

In the latter part of the 20th century, two design patterns were developed for deciding nick locations in single helix staple-only origami. One option, as developed by Ned Seeman is to have strand nicks in the center of edges, and this heuristic was later developed for the vHelix design software. 

Another option, as developed by Andrew Turberfield, is to have nicks at the vertex. Every vertex in a polyhedra is an intersection of at least three edges. In the Turberfield design pattern, each of these edges will necessarily consist of a different strand of DNA, and a nick should occur at only one strand (or equivalently edge) per vertex.

In Staplehedra, there are options for either Seeman style nicks (```--nicks 1```), or Turberfield nicks (```--ldmoverhangfile test.overhangs```). If the Seeman option is chosen (which is the default), then nick locations are decided by the software, by assigning each face (which corresponds to a strand) to a unique edge which borders it (which is where the strand will be broken), using a graph traversal algorithm. Alternatively, the Turberfield nick system requires an additional input file, which describes on which vertex each face is nicked. 

The objective of this exercise is to investigate whether the choice between Seeman nicks (at edges), and Turberfield nicks (at vertices) influences the structure and stability of the DNA tetrahedron.

## Running simulations

I suspect that single stranded DNA between edges at vertices affects structure and fluctuations. So I will run simulations with 0, 2 and 4 ssDNA nts between edges for both Turberfield and Seeman structures.

Simulations for the Seeman tetrahedron are run with the command

```bash
for i in 0 2 4; do (mkdir Seeman_$i && cd Seeman_$i && staplehedra ../plyfiles/tetrahedron.ply --spacers $i --coarse-steps 1e4 --fine-steps 1e4 --oxdna-steps 1e6 -d out --nicks 1 ) & done
```
Simulations for the Turberfield tetrahedron are run with the command

```bash
for i in 0 2 4; do (mkdir Turberfield_$i && cd Turberfield_$i && staplehedra ../plyfiles/tetrahedron.ply --spacers $i --coarse-steps 1e4 --fine-steps 1e4 --oxdna-steps 1e6 -d out --nicks 0 --ldmoverhangfile ../../../overhang_files/ldm_tet_seq.overhang) & done
```

Since these are quite small simulations, they can be run simultaneously on the same GPU. I imagine this substantially impacts performance, but I haven't ever done an explicit comparison.

## Post simulation analysis

Let's look first at the Seeman style edge nicks.

<img src = "https://imgur.com/PRDpFfY.gif"/>
Figure 1: Seeman nicks, 0 spacers. Reasonable tetrahedral shape.


<img src = "https://imgur.com/kRWdZWF.gif"/>
Figure 2: Seeman nicks, 2 spacers. Additional spacers allow deformation of molecule into structure with maximum stacking. This requires the nick to also correspond to a vertex, so the end resulting structures is no longer tetrahedral.

<img src = "https://imgur.com/gwzKznO.gif"/>
Figure 3: Seeman nicks, 4 spacers. Now with 4 spacers, stacking of the duplexes at edges is no longer particularly favorable, so well structured tetrahedron is recovered.

<img src = "https://imgur.com/5Pa2BaA.gif"/>
Figure 4: Turberfield nicks, 0 spacers.  While topologically a tetrahedron, this doesn't particularly resemble the original design, as at vertices, the edges are free to rotate due to the nick. 

<img src = "https://imgur.com/0FPDFwz.gif"/>
Figure 5: Turberfield nicks, 2 spacers. Similar to Figure 4, except that now with additional spacers, the vertex is more symmetric: while the nick does lead to rotation and opening of the vertex, the edge structure on the vertex is essentially still symmetric. This is because additional ssDNA spacers between adjacent edges add some bulk in that part of the vertex where the strand is not nicked.

<img src = "https://imgur.com/S0iiEnk.gif"/>
Figure 6: Turberfield nicks, 4 spacers. Essentially the same as Figure 5, but with enhanced vertex size.























