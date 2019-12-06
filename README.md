# Staplehedra

## Release 0.1

N.B: this script is barely functional at the moment, and it shouldn't really be used by anyone for anything... I also haven't got round to cryo-EM-ing structures, so I don't even know if the designed staplehedra fold correctly.  

## DNA polyhedra generation from staples

A Python3 package for the design and simulation of single-stranded DNA polyhedra constructed from staples only. The package takes as an input a ply file, denoting a faceted structure which is expandable to a sphere. It will return sequence choice for oligos, and files for simulation in MRDNA (Multi-Resolution DNA Simulations). Through the MRDNA framework it can perform hierarchical MD simulations, and output simulation files for the coarse grained force field oxDNA, and for atomistic simulation in NAMD.

## The design

I wrote this program so I could design and simulate small scale polyhedra (<20 nm in diameter) which would fold at high efficiency at vast concentrations (>5 uM) with 1:1 stoichiometry of strands. 

The idea is that you might be planning material synthesis from DNA (for example, for 3D crystals or hydrogels). But using scaffolded origami (e.g. m13mp18) starts to become super expensive if you are, say, creating visible quantities of product, and don't have an easy way of harvesting phage genome. On the other hand, staple-only DNA structures can be fairly affordable, even for the synthesis of macroscopic amounds of product.`

The current (tentative) design strategy is to use one oligonucleotide per face, as illustrated by the lighter colored strand in the tetrahedron below, and place the strand break in that oligonucleotide at a location midway between two vertices. Then repeat for all others faces.

<img src="https://i.imgur.com/R2x5KzT.png" width="1000"/>

**Figure 1**: Arbitrary polyhedra can be designed with staplehedra. Here, various DNA Platonic Solids are illustrated. Each face is a separate DNA strand, and for each polyhedron, one DNA strand has been highlighted.


### Controllable Flexibility

I have also included an option to have a variable number of single stranded DNAs between adjacent faces, which I imagine would control the flexibility of the structure. This may be useful if you are creating macromolecular structures where some feature has an unknown size -- for example, the size of cholesterol micelles in an amphiphilic structure. Alteratively, such a feature may be useful if studying patchy colloids with varying size of patch -- a more flexible DNA polyhedron would correspond to a patchy colloid with a larger patch.

<img src = "https://imgur.com/2TJtapN.gif" width=300/>
<img src = "https://imgur.com/FWDp1TM.gif" width=300/>
<img src = "https://imgur.com/zKYRxJr.gif" width=300/>

**Figure 2**: The flexibility of staplehedra can be tuned by controlling the number of unpaired bases between adjacent helices. Here are two cubes with 0 unpaired nucleotides (top), and 6 unparied nucleotides (middle), and 9 unpaired nucleotides (bottom), simulated with oxDNA using the MRDNA interface.

￼   ￼
<img src = https://imgur.com/6Kvx26m.gif width=300/>
<img src = https://imgur.com/SNVV66E.gif width=300/>
<img src = https://imgur.com/aVSEpOm.gif width=300/>

**Figure 3**: Control of flexibility is even more apparent in larger staplehedra, here, the icosahedron. Shown decending are short MD simulations for 0, 3, and 9, free nucleotides between adjacent edges. Note that with 3 free nucleotides, the formation of a re-entrant polyhedron is possible.

### Overhangs

Since the purpose of staplehedra is for the development of DNA polyhedron based materials, individual polyhedra need to be able to interact with each other. Also, most modified DNA (from modified phosphoramidite building blocks), like fluorophores / quenchers / chemical groups for conjugation to nanoparticles / cholesterol / biotin / digoxigenin etc. are cheapest if they occur at the terminus of a DNA strand. So it makes sense to write the abilility to create overhangs into Staplehedra. 

ssDNA overhangs are implemented as illustrated below in the oxDNA graphic. 

<img src = "https://i.imgur.com/dtxN7iH.png" width = 300>

**Figure 4**: Single stranded overhangs are implemented as illustrated: the purple strand begins at one vertex, wraps around the entire face and subsequenctly binds to itself. Strand breaks can also be implemented if this results in an excessively long strand of DNA for low error synthesis. Notice that here, the double stranded region of the overhang stacks, so the overhang position is not particularly flexible.

The stacking of the overhanging double stranded region of the sdDNA with the DNA in the frame of the structure may be problematic. It can be very slightly reduced with the use of ssDNA overhangs between the frame and the overhangs, although I suspect excessive sized overhangs will start to disrupt the structure. Ideally one would just have the connectivity without the stacking, using an abasic site -- either an IDT [dSpacer](https://eu.idtdna.com/site/Catalog/Modifications/Product/1202), or an IDT [C3 Spacer](https://eu.idtdna.com/site/Catalog/Modifications/Product/1056). 


<img src = ...>
**Figure 6**: 


For addition of overhangs in the design, an overhang file can be provied as demonstrated below.

```bash 
--overhangfile overhangfile.overhangs 
```

```
face,overhang,ds_length,ss_length,out_side,ss_extra
0,0,10,10,5,5
1,1,10,5,3,5
```

* face: the index of the face, whose DNA strand overhangs. The index is consistent with the ply file
* overhang: the index of the vertex where the overhang overhangs
* ds_length: the length of the double stranded region for the overhang
* ss_length: the length of the single stranded region in the overhang 
* out_side: whether the 5' or 3' side of the ssDNA faces outward
* ss_extra: any additional single stranded DNA between the overhanging region, and the DNA which is part of the face.


## Relationship to other software

There are several projects that are super close to doing exactly what I need. vHelix creates scaffolded polyhedra from single helices, so in principle strategic rebreaking/merging of strands could be used to create DNA polyhedra. However, this hacky solution would perform suboptimally when designing polyhedra whose connectivity of vertices is not an Eulerian Graph. Also the vHelix router has been written with Windows in mind, so it doesn't fit well my existing simulation/design workflows in Linux.

Alternatively, TALOS and DAEDALUS are great scaffold design softwares for DX tile / 4 helix bundle generation of polhedra. But I suspect that DX tiles will not form at high efficiency from 1:1 stoichiometry where the scaffold has been broken up into smaller sections. Also, I have found (via oxDNA simulation) that the vertices of DX tile polyhedra are too deformable for design of DNA materials with controllable Young's Modulus / Poisson ratio.

Staplehedra is essentially a very specific front end for the MRDNA package/project, which automates lots of the awkward / difficult geometrical descriptions of DNA. MRDNA also has a hierarchical coarse-grained MD relaxation proceedure which allows Staplehedra to be vague when describing the overarching geometry of DNA nanostructures, and to efficiency relax to an atomistic representation.

## Dependencies

* **Linux operating system** (may also work with OS X, or possibly but not probably the Ubuntu app in the Windows Store)
* MRDNA
* Python >= 3.5
* plyfile module for Python

## Stuff that works now

* Generation of PBBs of polyhedra from .ply files, made from circular DNA
* Relaxation in MRDNA.
* Controllable numbers of unpaired DNA nucleotides between vertices 
* Strategic placing of strand breaks to create polyhedra from linear DNA
* Export of sequences to a .csv file

## Plans for the immediate future (like early December)

* Simulation informed optimization of structures like in vHelix
* Optional addition of customizable overhangs at vertices/edges for aiding material synthesis, fluorophore / nanoparticle positioning for optoelectronics etc.

## Plans for early next year

* Experimental synthesis of structures
* Protocol optimization for folding
* Gels / AFM / CryoEM 
* Very cheap webserver for job submission if these experimental structures look okay

## Usage

You can call the script using all of the MRDNA arguments (including --coarse-steps --fine-steps --oxdna-steps -d etc.), with the additional options of --spacers, --bp, which are specific for the creation of the polyhedra.

```bash
bin/polymebabyonemoretime.py --bp 21 -d out --spacers 1 --coarse-steps 1e5 --fine-steps 1e5 ../ply_files/dodec.ply
```

Here, the length of the shortest edge is 21 nucleotides (other distances will be normalized to take account of this). At the moment, it makes a lot of sense to have this number as an integer multiple of 10.5, so as to match the helicity of DNA. The ply file is ../ply_files/dodec.ply -- I have included a tiny library of ply files in that directory. out is the output directory name. Here, I have 1 unparied strand of DNA between edges.
