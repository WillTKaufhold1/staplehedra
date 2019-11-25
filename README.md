# Staplehedra

![bicolored tetrahedron](https://imgur.com/WlF92LX)
Format: ![Alt Text](url)

## Release 0.0.

N.B: this script is barely function, and I'm shouldn't really be used by anyone for anything... I also haven't got round to cryo-EM-ing structures, so I don't even know if the designed structures fold correctly...

## DNA polyhedra generation from staples

A Python3 package for the design and simulation of single-stranded DNA polyhedra constructed from staples only. The package takes as an input a ply file, denoting a faceted structure which is expandable to a sphere. It can return sequence choice for oligos, and files for simulation in MRDNA (Multi-Resolution DNA Simulations). Through the MRDNA framework it can return simulation files for the coarse grained force field oxDNA, and for atomistic simulation. 

## Dependencies

* Linux operating system 
* MRDNA
* Python >= 3.5
** plyfile 

## Usage

```bash
bin/polymebabyonemoretime.py 10 ../ply_files/dodec.ply  
```
