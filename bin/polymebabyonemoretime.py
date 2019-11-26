#! /usr/bin/env python3
import numpy as np
import mrdna
import sys
import argparse
from copy import deepcopy as copy

from mrdna.coords import readArbdCoords
import segment_maker

def make_parser():

    parser = argparse.ArgumentParser(description = "Make poly")
    parser.add_argument('plyfile',type=str,help="take in a .ply format. Remember to have consistent facing normals!")

    parser.add_argument('--bp',type=int,default = 21,help="Length of shortest side in nucleotides")
    parser.add_argument('--dirname',type=str, default="out",help="Simulation directory, created")
    parser.add_argument('--spacers',type=int,default=1,help="ssDNA between adjacent edges")

    parser.add_argument('--atoms',action="store_true",help="generate files for NAMD simulation")
    parser.add_argument('--oxdna',action="store_true",help="generate oxDNA files")

    return parser

def __main__():

    #TODO: add sane method of argparser

    SIM_NAME = "DNA"

    parser = make_parser()
    args = parser.parse_args()

    generate_atomic = bool(args.atoms)
    generate_oxdna = bool(args.oxdna)
    LENGTH_OF_SMALLEST = int(args.bp)
    FNAME = str(args.plyfile)
    DIRNAME = str(args.dirname)
    SPACERS = int(args.spacers)
    
    segs_list = segment_maker.get_segments(FNAME, LENGTH_OF_SMALLEST, SPACERS)

    model = mrdna.SegmentModel(
                        segs_list,
                        local_twist = True,
                        max_basepairs_per_bead = 5,
                        max_nucleotides_per_bead = 5,
                        dimensions=(5000,5000,5000),
                    )

    model.simulate(output_name = SIM_NAME, directory=DIRNAME,
               num_steps=1e4, output_period=1e3)

    coords = readArbdCoords('%s/output/%s.restart'%(DIRNAME,SIM_NAME))
    model.update_splines(coords) 

    from mrdna.model.dna_sequence import read_sequence_file
    seq = read_sequence_file() 

    model.set_sequence(seq)   

    if generate_atomic:
        model.generate_atomic_model()
        model.write_atomic_ENM( output_name = "%s/%s"%(DIRNAME,SIM_NAME+'-atomic',))  
        model.atomic_simulate( output_name = "%s/%s"%(DIRNAME,SIM_NAME+'-atomic') )
    if generate_oxdna:
        model.generate_oxdna_model()
        model._write_oxdna_configuration('%s/prova.conf'%(DIRNAME,))
        model._write_oxdna_topology('%s/top.top'%(DIRNAME,))
        model._write_oxdna_input('%s/in'%(DIRNAME,),'%s/top.top'%(DIRNAME,),'%s/prova.conf'%(DIRNAME,),'%s/trajectory.dat'%(DIRNAME,),'%s/last_conf.dat'%(DIRNAME),
                            '%s/log'%(DIRNAME,))

__main__()
