# Staplehedra

## Release 0.1

N.B: this script is barely functional at the moment, and it shouldn't really be used by anyone for anything... I also haven't got round to cryo-EM-ing structures, so I don't even know if the designed staplehedra fold correctly.  

## DNA polyhedra generation from staples

A Python3 package for the design and simulation of single-stranded DNA polyhedra constructed from staples only. The package takes as an input a ply file, denoting a faceted structure which is expandable to a sphere. It will return sequence choice for oligos, and files for simulation in MRDNA (Multi-Resolution DNA Simulations). Through the MRDNA framework it can perform hierarchical MD simulations, and output simulation files for the coarse grained force field oxDNA, and for atomistic simulation in NAMD.

## The design

I wrote this program so I could design and simulate small scale polyhedra (<20 nm in diameter) which would fold at high efficiency at vast concentrations (>5 uM) with 1:1 stoichiometry of strands. 

The idea is that you might be planning material synthesis from DNA (for example, for 3D crystals or hydrogels). But using scaffolded origami (e.g. m13mp18) starts to become super expensive if you are, say, creating visible quantities of product, and don't have an easy way of harvesting phage genome. On the other hand, staple-only DNA structures can be fairly affordable, even for the synthesis of macroscopic amounds of product.`

The current (tentative) design strategy is to use one oligonucleotide per face, as illustrated by the lighter colored strand in the tetrahedron below, and place the strand break in that oligonucleotide at a location midway between two vertices. Then repeat for all others faces.

I have also included an option to have a variable number of single stranded DNAs between adjacent faces, which I imagine would control the flexibility of the structure.

<img src="https://i.imgur.com/FOQmBUB.png" width="1000"/>

Figure 1: various DNA Platonic Solids generated with Staplehedra. Each face is a separate DNA strand, and for each polyhedron, one DNA strand has been highlighted. Note that some nucleotides are displaced -- these can be relaxed into their proper positions using gradient descent with atomistic simulation.

## Relationship to other software

There are several projects that are super close to doing exactly what I need. vHelix creates scaffolded polyhedra from single helices, so in principle strategic rebreaking/merging of strands could be used to create DNA polyhedra. However, this hacky solution would perform suboptimally when designing polyhedra whose connectivity of vertices is not an Eulerian Graph. Also the vHelix router has been written with Windows in mind, so it doesn't fit well my existing simulation/design workflows in Linux.

Alternatively, TALOS and DAEDALUS are great scaffold design softwares for DX tile / 4 helix bundle generation of polhedra. But I suspect that DX tiles will not form at high efficiency from 1:1 stoichiometry where the scaffold has been broken up into smaller sections. Also, I have found (via oxDNA simulation) that the vertices of DX tile polyhedra are too deformable for design of DNA materials with controllable Young's Modulus / Poisson ratio.

Staplehedra is essentially a very specific front end for the MRDNA package/project, which automates lots of the awkward / difficult geometrical descriptions of DNA. MRDNA also has a hierarchical coarse-grained MD relaxation proceedure which allows Staplehedra to be vague when describing the overarching geometry of DNA nanostructures, and to efficiency relax to an atomistic representation.

## Dependencies

* Linux operating system 
* MRDNA
* Python >= 3.5
** plyfile 

## Stuff that works now

* Generation of PBBs of polyhedra from .ply files, made from circular DNA
* Relaxation in MRDNA.
* Controllable numbers of unpaired DNA nucleotides between vertices (1 or more...)

## Stuff that will (hopefully) work tomorrow (evening of 26/11/19 GMT)

* Strategic placing of strand breaks to create polyhedra from linear DNA
* Export of sequences to a .csv file
* Automatic atomistic / oxDNA simulation 
* Eliminating lower limit of number of unpaired DNA nucleotides between vertices

## Plans for the immediate future (like early December)

* Simulation informed optimization of structures like in vHelix
* Optional addition of customizable overhangs at vertices/edges for aiding material synthesis, fluorophore / nanoparticle positioning for optoelectronics etc.

## Plans for early next year

* Experimental synthesis of structures
* Protocol optimization for folding
* Gels / AFM / CryoEM 
* Very cheap webserver for job submission if these experimental structures look okay

## Usage

You can call the script with:

```bash
bin/polymebabyonemoretime.py 10 ../ply_files/dodec.ply sim1 1
```

Here, the length of the smallest edge is 10 nucleotides (other distances will be normalized to take account of this). At the moment, it makes a lot of sense to have this number as an integer multiple of 10.5, so as to match the helicity of DNA. The ply file is ../ply_files/dodec.ply -- I have included a tiny library of ply files in that directory. sim1 is the output directory name. Here, I have 1 unparied strand of DNA between edges. For reasons releated to my laziness, this can't be reduced to 0 yet...
