# Staplehedra


## Release 0.0.

N.B: this script is barely function, and I'm shouldn't really be used by anyone for anything... I also haven't got round to cryo-EM-ing structures, so I don't even know if the designed structures fold correctly...

## DNA polyhedra generation from staples

A Python3 package for the design and simulation of single-stranded DNA polyhedra constructed from staples only. The package takes as an input a ply file, denoting a faceted structure which is expandable to a sphere. It can return sequence choice for oligos, and files for simulation in MRDNA (Multi-Resolution DNA Simulations). Through the MRDNA framework it can return simulation files for the coarse grained force field oxDNA, and for atomistic simulation. 

## The design

I wrote this program so I could design small scale polyhedra (<10 nm in diamter) which would fold at high efficiency at vast concentrations (>1 uM) with 1:1 stoichiometry of strands.

The idea is that if you're planning material synthesis from DNA (for example, for 3D crystals or gels), using scaffolded origami (e.g. m13mp18) starts to become super expensive for creating visible quantities of product. But staple-only DNA structures can be fairly affordable, even for the synthesis of macroscopic amounds of product.`

The current tentative design strategy is to take use one oligonucleotide per face, 

<img src="https://i.imgur.com/WlF92LX.png" width="200"/>
Fig1: a DNA 

## Dependencies

* Linux operating system 
* MRDNA
* Python >= 3.5
** plyfile 

## Usage

```bash
bin/polymebabyonemoretime.py 10 ../ply_files/dodec.ply  
```
