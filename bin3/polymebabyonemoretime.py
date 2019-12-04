#! /usr/bin/env python3
import numpy as np
import pandas as pd
import mrdna
import sys
import argparse
from copy import deepcopy as copy
from functools import reduce
from mrdna.coords import readArbdCoords
from mrdna.simulate import multiresolution_simulation as simulate

from mrdna_parser import make_parser
import segment_maker

from m13_grabber import m13seq


def read_overhang_file(fname):
    return pd.read_csv(fname)

def get_strand_sequence(strand):
    seg_sequences = []
    for i in strand.strand_segments:
        seg_sequences.append(i.get_sequence())
    return reduce(lambda x,y : x+y, 
               reduce(lambda x,y : x+y, seg_sequences,[]),
            '')

def get_sequences(model):
    seqs = []
    for index,i in enumerate(model.strands):
        print (i)
        print (index)
        seqs.append(get_strand_sequence(i))
        print (get_strand_sequence(i))
    return seqs

def export_sequences(model,fname):
    seqs = get_sequences(model)
    with open(fname,'w+') as f: 
        for i in seqs:
            f.writelines(i+'\n')

def __main__():

    parser = make_parser()
    args = parser.parse_args()

    LENGTH_OF_SMALLEST = int(args.bp)
    FNAME = str(args.plyfile)
    SPACERS = int(args.spacers)

    overhangfilename = args.overhangfile 
    if overhangfilename is not None:
        overhangs = read_overhang_file(overhangfilename) 

    #breakpoint()

    segs_list = segment_maker.get_segments(FNAME, LENGTH_OF_SMALLEST, SPACERS,overhangs)
    

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

    export_sequences(model,args.seqfile)
    
    simulate( **run_args )


__main__()
