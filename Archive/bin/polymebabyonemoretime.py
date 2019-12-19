#! /usr/bin/env python3
import numpy as np
import mrdna
import sys
import argparse
from copy import deepcopy as copy

from mrdna.coords import readArbdCoords
import segment_maker

from mrdna.simulate import multiresolution_simulation as simulate

def make_parser():

    #specific options for staplehedra

    parser = argparse.ArgumentParser(description = "Make poly")
    parser.add_argument('plyfile',type=str,help="take in a .ply format. Remember to have consistent facing normals!")

    parser.add_argument('--bp',type=int,default = 21,help="Length of shortest side in nucleotides")
    parser.add_argument('--spacers',type=int,default=1,help="ssDNA between adjacent edges")

    #code lifted from MRDNA for consistency in inputs ... I should move this to a module...
    #or reformat the entire thing to be a part of MRDNA? 

    parser.add_argument('-o','--output-prefix', type=str, default=None,
                        help="Name for your job's output")
    parser.add_argument('-d','--directory', type=str,  default=None,
                        help='Directory for simulation; does not need to exist yet')
    parser.add_argument('-g','--gpu', type=int, default=0,
                        help='GPU for simulation; check nvidia-smi for availability')
    parser.add_argument('--debye-length', type=float, default=None,
                        help='Adjust the electrostatic repulsion matching osmotic pressure data from Rau and Parsegian under 25 mM MgCl2 condition with a Debye-Hueckel correction from the default 11.1 Angstroms.')

    parser.add_argument('--temperature', type=float, default=295,
                        help='Temperature in Kelvin.')

    parser.add_argument('--sequence-file', type=str, default=None,
                        help='Sequence of longest strand.')

    parser.add_argument('--output-period', type=float, default=1e4,
                        help='Simulation steps between DCD frames')
    # parser.add_argument('--minimization-steps', type=float, default=0,
    #                     help='Simulation steps between DCD frames')
    parser.add_argument('--coarse-local-twist', action='store_true',
                        help='Use local twist for coarsest representation?')
    parser.add_argument('--fix-linking-number', action='store_true',
                        help='Fix the linking number so that it cannot change once local twist is introduced?')
    parser.add_argument('--coarse-steps', type=float, default=1e7,
                        help='Simulation steps for coarse model (200 fs/step)')
    parser.add_argument('--fine-steps', type=float, default=1e7,
                        help='Simulation steps for fine model (50 fs/step)')

    parser.add_argument('--oxdna-steps', type=float, default=None,
                        help='Perform an oxDNA simulation instead of creating atomic model')
    parser.add_argument('--oxdna-output-period', type=float, default=1e4,
                        help='Simulation steps between oxDNA configuration and energy output')

    parser.add_argument('--backbone-scale', type=float, default=1.0,
                        help='Factor to scale DNA backbone in atomic model; try 0.25 to avoid clashes for atomistic simulations')
    parser.add_argument('--debug', action='store_true',
                        help='Run through the python debugger?')

    parser.add_argument('--draw-cylinders', action='store_true',
                        help='Whether or not to draw the cylinders')
    parser.add_argument('--draw-tubes', action='store_true',
                        help='Whether or not to draw the tubes')

    return parser


def m13seq():
    from functools import reduce
    with open('../m13mp18.dat','r') as f:
        data = f.readlines()
    filtered = list(filter(lambda x : x[0] != ';',data))
    stripped = list(map(lambda x : x.replace(' ','').replace('\n',''), filtered))
    accumulalted = str(reduce( lambda x,y : x + y, stripped, ''))
    return accumulalted


def __main__():

    #TODO: add sane method of argparser

    parser = make_parser()
    args = parser.parse_args()

    LENGTH_OF_SMALLEST = int(args.bp)
    FNAME = str(args.plyfile)
    SPACERS = int(args.spacers)
    
    segs_list = segment_maker.get_segments(FNAME, LENGTH_OF_SMALLEST, SPACERS)

    model = mrdna.SegmentModel(
                        segs_list,
                        local_twist = True,
                        dimensions=(5000,5000,5000),
                    )
    #apply sequence here

    model.set_sequence(m13seq(),force=False)


    prefix = "DNA" 


    run_args = dict(
                model = model,
                output_name = prefix,
                job_id = "job-" + prefix,
                directory = args.directory,
                gpu = args.gpu,
                minimization_output_period = int(args.output_period),
                coarse_local_twist = args.coarse_local_twist,
                fix_linking_number = args.fix_linking_number,
                coarse_output_period = int(args.output_period),
                fine_output_period = int(args.output_period),
                minimization_steps = 0, # int(args.minimization_steps),
                coarse_steps = int(args.coarse_steps),
                fine_steps = int(args.fine_steps),
                backbone_scale = args.backbone_scale,
                oxdna_steps = args.oxdna_steps,
                oxdna_output_period = args.oxdna_output_period
            )

    simulate( **run_args )

__main__()
