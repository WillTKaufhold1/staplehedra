# Staplehedra

## Release -1.0...

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
