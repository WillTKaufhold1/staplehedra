import mrdna
from plyfile import PlyData
import numpy as np
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
    def apply_nick(self,nick):
        self.nick = nick

class edge():
    def __init__(self,index,start,stop):
        self.index = index
        self.start = start
        self.stop = stop
        self.rawlen = np.sqrt(np.sum((np.array(stop)-np.array(start))**2))
        
    def normalize(self,min_length,bp,overhang=False,bp_overhang=None):

        #bp is the number of base pairs in the shortest length.

        compress_factor = 0.1
        
        self.nStart = self.start / min_length * bp * 3.4
        self.nStop = self.stop / min_length * bp * 3.4
        
        #a bit of compression to avoid collisions... TODO: be more strategic...
        vector = self.nStop-self.nStart
        self.nStart_c = self.nStart + vector * compress_factor 
        self.nStop_c  = self.nStop - vector * compress_factor
        
        self.nLen = np.sqrt(np.sum((np.array(self.nStop)-np.array(self.nStart))**2))
        if not overhang:
            self.nBp = bp
        elif overhang:
            self.nBp = bp_overhang
        else:
            print('Bug')

def get_graph(face_data):
    E_vals = [i.e for i in face_data]
    graph = {}
    attached = []

    for face_index,face in enumerate(E_vals):
        graph[face_index] = []
        for edge in face:
            for other_edge_index,other_edge_list in enumerate(E_vals):
                if (edge[::-1] in other_edge_list):
                    graph[face_index].append(other_edge_index)
    return graph

def get_break_index(face_data):
    
    graph = get_graph(face_data)
    breaks = {}
    face = 0
    breaks[face] = graph[face][0]
    graph[face].remove(breaks[face])
    graph[breaks[face]].remove(face)
    leading_edge = copy(graph[face])
    for l in leading_edge:
        breaks[l] = face
        
    while leading_edge != []:        
        for face in copy(leading_edge):
            for otherface in copy(graph[face]):
                if otherface not in breaks.keys():

                    leading_edge.append(copy(otherface))
                    breaks[otherface] = copy(face)

                    graph[face].remove(otherface)
                    graph[otherface].remove(face)

            leading_edge.remove(face)

    return breaks

def assign_nicks(face_data):

    breaks = get_break_index(face_data)

    E_vals = [i.e for i in face_data]

    edge_dict = {}

    for face1,edges_1 in enumerate(E_vals):
        for edge in edges_1:
            for face2,edges_2 in enumerate(E_vals):
                if (edge[::-1] in edges_2 and edge not in edges_2):
                    edge_dict[(face1,face2)] = edge

    for b in breaks:
        face_data[b].apply_nick(edge_dict[(b,breaks[b])])

def break_me(segs,breaks):
    #segs_list[0].add_nick(34,fwd_strand=True)
    pass

def get_segments(FNAME, LENGTH_OF_SMALLEST, SPACERS,overhangs):

    locs,faces = read_ply(FNAME)

    mentioned_edges = set()
    face_data = [face(i) for i in faces]

    #NICKS
    #now we work out where the nicks are so we can choose the correct orientation of our strands.
    assign_nicks(face_data)
    nick_locs = []
    for f in face_data:
        nick_locs.append(f.nick)
    mentioned_edges = set(nick_locs)
    #NICKS

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

    #add the dsDNA overhangs here. the current problem is related to the normalization 
    overhang_ds_edges = []
    for i in range(len(overhangs)):
        vertex = overhangs.iloc[i]['overhang'] 
        vertex_location = locs[vertex] 
        bp = overhangs.iloc[i]['ds_length'] 
        length = bp * min_length / float(LENGTH_OF_SMALLEST) / 3.4 

        end = length * np.array(vertex_location) / np.sqrt(np.sum(np.array(vertex_location)**2)) + np.array(vertex_location)

        print (bp,length,vertex_location,end)
        overhang_ds_edges.append(edge(vertex,vertex_location,end))

    #serious problem with normalization of side lengths...


    for i in edges:
        i.normalize(min_length,LENGTH_OF_SMALLEST)
    
    for i in overhang_ds_edges:
        i.normalize(min_length,LENGTH_OF_SMALLEST,overhang=True,bp_overhang=5) # the 5 is a placeholder...

    single_stranded_dna = []

    segs = {}

    for index, e in enumerate(edges):
        segs[tuple(e.index)] = (
            mrdna.DoubleStrandedSegment(name = 'helix%s'%(index,),
                                num_bp = e.nBp,
                                start_position = e.nStart_c,  
                                end_position = e.nStop_c
                                       ))

    #overhangs
    overhang_segs = {}#you should rethink this datastructure.
    for index, e in enumerate(overhang_ds_edges):
        print (e.nBp, e.nStart_c,e.nStop_c)
        tmp = copy((mrdna.DoubleStrandedSegment(name = 'helix%s'%(index,),
                                num_bp = int(e.nBp),
                                start_position = list(e.nStart_c),  
                                end_position = list(e.nStop_c)
                                       )))
        overhang_segs[e.index] = tmp

    ssDNA_index = 0

    for f in face_data:

        for con in f.connections:

            ssDNA_index += 1

            c1,c2 = con
            c1_positive = c1 in segs
            c2_positive = c2 in segs

            if SPACERS != 0:
                if c1_positive and c2_positive:

                    ss = copy(mrdna.SingleStrandedSegment("strand%s"%(ssDNA_index,),
                        start_position = segs[c1].end_position ,
                        end_position =   segs[c2].start_position,
                        num_nt = SPACERS))

                    segs[c1].connect_end3(ss)
                    segs[c2].connect_start5(ss)

                elif c1_positive and not c2_positive:
                    ss = copy(mrdna.SingleStrandedSegment("strand%s"%(ssDNA_index,),
                        start_position = segs[c1].end_position ,
                        end_position =   segs[c2[::-1]].end_position,
                        num_nt = SPACERS))

                    segs[c1].connect_end3(ss)
                    segs[c2[::-1]].connect_end5(ss)

                elif not c1_positive and c2_positive:

                    ss = copy(mrdna.SingleStrandedSegment("strand%s"%(ssDNA_index,),
                        start_position = segs[c1[::-1]].start_position ,
                        end_position =   segs[c2].start_position,
                        num_nt = SPACERS))

                    segs[c1[::-1]].connect_start3(ss)
                    segs[c2].connect_start5(ss)

                elif not c1_positive and not c2_positive:

                    ss = copy(mrdna.SingleStrandedSegment("strand%s"%(ssDNA_index,),
                        start_position = segs[c1[::-1]].start_position ,
                        end_position =   segs[c2[::-1]].end_position,
                        num_nt = SPACERS))      

                    #segs[c1[::-1]].connect_start3(segs[c2[::-1]].end5)

                    segs[c1[::-1]].connect_start3(ss)
                    segs[c2[::-1]].connect_end5(ss)

                single_stranded_dna.append(ss)
            else:
                if c1_positive and c2_positive:
                    segs[c1].connect_end3(segs[c2].start5,type_="terminal_crossover")
                elif c1_positive and not c2_positive:
                    segs[c1].connect_end3(segs[c2[::-1]].end5,type_="terminal_crossover")
                elif not c1_positive and c2_positive:
                    segs[c1[::-1]].connect_start3(segs[c2].start5,type_="terminal_crossover")
                elif not c1_positive and not c2_positive:
                    segs[c1[::-1]].connect_start3(segs[c2[::-1]].end5,type_="terminal_crossover")

    
    #NICKS
    for f in face_data:
        nick = f.nick
        if nick in segs:
            segs[nick].add_nick(10,on_fwd_strand=True)
        else:
            print ('failure')


    segs_list = [segs[i] for i in segs] + single_stranded_dna + [overhang_segs[i] for i in overhang_segs]

    return segs_list