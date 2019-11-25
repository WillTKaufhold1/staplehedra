#! /usr/bin/env python3
import numpy as np
import pandas as pd
import mrdna
import sys
import argparse
from plyfile import PlyData
from copy import deepcopy as copy

def read_ply(fname):
    polydata = PlyData.read(fname)
    locs = list(zip(polydata['vertex']['x'],polydata['vertex']['y'],polydata['vertex']['z']))
    faces = polydata['face']['vertex_index']
    return locs,faces

class face():
    def __init__(self,vertices):
        self.v = vertices
        self.e = list(zip(self.v,np.roll(self.v,-1)))
        self.ori = []

        self.connections = list(zip(self.e,
                                    map(tuple,np.roll(self.e,-1,axis = 0))))    

class edge():
    def __init__(self,index,start,stop):
        self.index = index
        self.start = start
        self.stop = stop
        self.rawlen = np.sqrt(np.sum((np.array(stop)-np.array(start))**2))
        
    def normalize(self,min_length,bp):
        
        compress_factor = 0.1
        
        self.nStart = self.start / min_length * bp * 3.4
        self.nStop = self.stop / min_length * bp * 3.4
        
        #a bit of compression to avoid collisions.
        vector = self.nStop-self.nStart
        self.nStart_c = self.nStart + vector * compress_factor 
        self.nStop_c  = self.nStop - vector * compress_factor
        
        self.nLen = np.sqrt(np.sum((np.array(self.nStop)-np.array(self.nStart))**2))
        self.nBp = bp


def __main__():
    
    generate_atomic = True
    generate_oxdna = True

    parser = argparse.ArgumentParser(description = "Make poly")
    parser.add_argument('bp')
    parser.add_argument('plyfile')
    parser.add_argument('dirname')
    parser.add_argument('spacers')

    args = parser.parse_args()

    LENGTH_OF_SMALLEST = int(args.bp)
    FNAME = str(args.plyfile)
    DIRNAME = str(args.dirname)
    SPACERS = int(args.spacers)

    locs,faces = read_ply(FNAME)

    mentioned_edges = set()
    face_data = [face(i) for i in faces]

    for i in face_data:
        for e in i.e:
            if e not in mentioned_edges and e[::-1] not in mentioned_edges:
                i.ori.append(1)
                mentioned_edges.add(e)
            else:
                i.ori.append(-1)

    edges = []

    for i in mentioned_edges:
        edges.append(
            edge(i,locs[i[0]],locs[i[1]])
        )

    min_length = np.min([x.rawlen for x in edges])
    max_length = np.max([x.rawlen for x in edges])

    for i in edges:
        i.normalize(min_length,LENGTH_OF_SMALLEST)

    single_stranded_dna = []

    segs = {}

    for index, e in enumerate(edges):
        segs[tuple(e.index)] = (
            mrdna.DoubleStrandedSegment(name = 'helix%s'%(index,),
                                num_bp = e.nBp,
                                start_position = e.nStart_c,           
                                end_position = e.nStop_c
                                       ))

    for f in face_data:

        for con in f.connections:

            c1,c2 = con
            c1_positive = c1 in segs
            c2_positive = c2 in segs

            if c1_positive and c2_positive:

                ss = copy(mrdna.SingleStrandedSegment("strand2",
                      start_position = segs[c1].end_position ,
                      end_position =   segs[c2].start_position,
                      num_nt = SPACERS))

                segs[c1].connect_end3(ss)
                segs[c2].connect_start5(ss)

            elif c1_positive and not c2_positive:
                ss = copy(mrdna.SingleStrandedSegment("strand2",
                      start_position = segs[c1].end_position ,
                      end_position =   segs[c2[::-1]].end_position,
                      num_nt = SPACERS))

                segs[c1].connect_end3(ss)
                segs[c2[::-1]].connect_end5(ss)

            elif not c1_positive and c2_positive:

                ss = copy(mrdna.SingleStrandedSegment("strand2",
                      start_position = segs[c1[::-1]].start_position ,
                      end_position =   segs[c2].start_position,
                      num_nt = SPACERS))

                segs[c1[::-1]].connect_start3(ss)
                segs[c2].connect_start5(ss)


            elif not c1_positive and not c2_positive:

                ss = copy(mrdna.SingleStrandedSegment("strand2",
                      start_position = segs[c1[::-1]].start_position ,
                      end_position =   segs[c2[::-1]].end_position,
                      num_nt = SPACERS))      

                #segs[c1[::-1]].connect_start3(segs[c2[::-1]].end5)

                segs[c1[::-1]].connect_start3(ss)
                segs[c2[::-1]].connect_end5(ss)

            single_stranded_dna.append(ss)

            '''
            if c1_positive and c2_positive:
                segs[c1].connect_end3(segs[c2].start5)
            elif c1_positive and not c2_positive:
                segs[c1].connect_end3(segs[c2[::-1]].end5)
            elif not c1_positive and c2_positive:
                segs[c1[::-1]].connect_start3(segs[c2].start5)
            elif not c1_positive and not c2_positive:
                segs[c1[::-1]].connect_start3(segs[c2[::-1]].end5)
            else:
                print ('ohfuck') 
            '''

    segs_list = [segs[i] for i in segs] + single_stranded_dna

    model = mrdna.SegmentModel(
                        segs_list,
                        local_twist = True,
                        max_basepairs_per_bead = 5,
                        max_nucleotides_per_bead = 5,
                        dimensions=(5000,5000,5000),
                    )

    model.simulate(output_name = "loop", directory=DIRNAME,
               num_steps=1e6, output_period=1e3)

    from mrdna.coords import readArbdCoords
    coords = readArbdCoords('%s/output/loop.restart'%(DIRNAME,))
    model.update_splines(coords) 

    from mrdna.model.dna_sequence import read_sequence_file
    seq = read_sequence_file() 

    model.set_sequence(seq)   


    if generate_atomic:
        model.generate_atomic_model()
        model.write_atomic_ENM( output_name = "%s/loop-atomic"%(DIRNAME,))  

    if generate_oxdna:
        model.generate_oxdna_model()
        model._write_oxdna_configuration('%s/prova.conf'%(DIRNAME,))
        model._write_oxdna_topology('%s/top.top'%(DIRNAME,))
        model._write_oxdna_input('%s/in'%(DIRNAM,),'%s/top.top'%(DIRNAM,),'%s/prova.conf'%(DIRNAME,),'%s/trajectory.dat'%(DIRNAME,),'%s/last_conf.dat'%(DIRNAME),
                            '%s/log'%(DIRNAME,))

__main__()