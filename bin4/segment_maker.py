import mrdna
from plyfile import PlyData
import numpy as np
from copy import deepcopy as copy

#you absolutely have to restructure this program, it's getting totally unmanageable...
#also you need to think systematically about the datastructures

def read_ply(fname):
    polydata = PlyData.read(fname)
    locs = list(zip(polydata['vertex']['x'],polydata['vertex']['y'],polydata['vertex']['z']))
    faces = polydata['face']['vertex_index']
    return locs,faces

class FACE():
    def __init__(self,vertices,index):
        self.v = vertices
        self.e = list(zip(self.v,np.roll(self.v,-1)))
        self.ori = []
        self.index = index

        self.connections = list(zip(self.e,
                                    map(tuple,np.roll(self.e,-1,axis = 0))))    
    def apply_nick(self,nick):
        self.nick = nick

class edge():
    def __init__(self,index,start,stop,overhang=False,bp_overhang=None,out_orientation=None):
        self.index = index
        self.start = start
        self.stop = stop
        self.rawlen = np.sqrt(np.sum((np.array(stop)-np.array(start))**2))
        self.overhang = overhang
        
        if overhang:
            self.bp_overhang = bp_overhang
        
        if out_orientation: #lol, this is literally what inheritance is designed to stop...
            if out_orientation not in [5,3]:
                print ('orientation must be one of 5 or 3')
            self.out_orientation = out_orientation #the ssDNA polarity facing outwards.
        
    def normalize(self,min_length,bp):

        #bp is the number of base pairs in the shortest length.

        compress_factor = 0.1
        
        self.nStart = self.start / min_length * bp * 3.4
        self.nStop = self.stop / min_length * bp * 3.4
        
        #a bit of compression to avoid collisions... TODO: be more strategic...
        vector = self.nStop-self.nStart
        self.nStart_c = self.nStart + vector * compress_factor 
        self.nStop_c  = self.nStop - vector * compress_factor
        
        self.nLen = np.sqrt(np.sum((np.array(self.nStop)-np.array(self.nStart))**2))
        if not self.overhang:
            self.nBp = bp
        elif self.overhang:
            self.nBp = self.bp_overhang
        else:
            print('Bug')

def get_graph(face_data):
    E_vals = [i.e for i in face_data]
    graph = {}

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


def get_segments(FNAME, LENGTH_OF_SMALLEST, SPACERS,overhangs = None):

    locs,faces = read_ply(FNAME)

    mentioned_edges = set()
    face_data = []
    for index,i in enumerate(faces):
        print (FACE)
        face_data.append(FACE(i,index))
    
    #NICKS
    #now we work out where the nicks are so we can choose the correct orientation of our strands.
    assign_nicks(face_data)
    nick_locs = []
    for f in face_data:
        nick_locs.append(f.nick)
    mentioned_edges = set(nick_locs)
    #NICKS

    #choose orientation
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

    #ds overhangs
    #rework this data structure...

    if overhangs is not None:
            
        overhang_ds_edges = []
        overhang_ss_edges = []

        for i in range(len(overhangs)):
            vertex = overhangs.iloc[i]['overhang'] 
            vertex_location = locs[vertex] 
            ##ds overhangs
            bp = overhangs.iloc[i]['ds_length'] 
            length = bp * min_length / float(LENGTH_OF_SMALLEST) / 3.4 
            end = length * np.array(vertex_location) / np.sqrt(np.sum(np.array(vertex_location)**2)) + np.array(vertex_location)
            overhang_ds_edges.append(edge(vertex,vertex_location,end,overhang=True,bp_overhang=bp))

            #ss overhangs        

            bp_ss = overhangs.iloc[i]['ss_length'] 
            length_ss = bp_ss * min_length / float(LENGTH_OF_SMALLEST) / 3.4 
            start_ss = end 
            end_ss = end + length_ss * np.array(vertex_location) / np.sqrt(np.sum(np.array(vertex_location)**2))
            overhang_ss_edges.append(edge(vertex,start_ss,end_ss,overhang=True,bp_overhang=bp_ss,out_orientation=overhangs.iloc[i]['out_side'])) 

        for i in edges + overhang_ds_edges + overhang_ss_edges:
            i.normalize(min_length,LENGTH_OF_SMALLEST)
    else:
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

    if overhangs is not None:
        #overhangs
        overhang_segs = {}
        for index, e in enumerate(overhang_ds_edges):
            tmp = copy((mrdna.DoubleStrandedSegment(name = 'helix%s'%(index+100,),
                                    num_bp = int(e.nBp),
                                    start_position = list(e.nStart_c),  
                                    end_position = list(e.nStop_c)
                                        )))
            overhang_segs[e.index] = tmp

        #overhangs
        ss_overhang_segs = {}
        for index, e in enumerate(overhang_ss_edges):
            if e.out_orientation == 5:
                tmp = copy((mrdna.SingleStrandedSegment(name = 'helix%s'%(index+200,),
                                        num_nt = int(e.nBp),
                                        start_position = list(e.nStart_c),  
                                        end_position = list(e.nStop_c)
                                            )))
            
            elif e.out_orientation == 3:
                tmp = copy((mrdna.SingleStrandedSegment(name = 'helix%s'%(index,),
                                        num_nt = int(e.nBp),
                                        start_position = list(e.nStop_c),
                                        end_position = list(e.nStart_c),  
                                            )))
            else:
                print ('BUG')
            ss_overhang_segs[e.index]= (tmp)

    ssDNA_index = 0

    if overhangs is not None:
    #CONNECT SS_DNA TO OVERHANGS
        for ss in overhang_ss_edges:
            #connected if they have the same index!
            index = ss.index
            if overhangs is not None:
                dsseg = overhang_segs[index]
                ssseg = ss_overhang_segs[index]
                #the dsDNA faces outwards
            
                if overhang_ss_edges[index].out_orientation == 3:
                    overhang_segs[index].connect_end3(ss_overhang_segs[index])
                elif overhang_ss_edges[index].out_orientation == 5:
                    overhang_segs[index].connect_end5(ss_overhang_segs[index])

    def intersection(lst1, lst2): 
        lst3 = [value for value in lst1 if value in lst2] 
        return lst3 
  
    for f in face_data:

        for con in f.connections:
            ### See if we have an overhang 
            overhang_here = False
            if overhangs is not None:
                vertex_index = intersection(con[0],con[1])[0] 
                face_index = f.index    

                for i in range(len(overhangs)):
                    if (overhangs.iloc[i]['face'] == face_index and overhangs.iloc[i]['overhang'] == vertex_index):
                        overhang_here = True

                ###
                if overhang_here:
                    print (f.index)
                    #but obviously there isn't only one overhang?

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
                #let's also add connectivity to the appropriate overhang...
                if overhang_here:
                    print ('adding overhang!')

                    if c1_positive:
                        segs[c1].connect_end3(overhang_segs[vertex_index].start5, type_="terminal_crossover")
                        overhang_segs[vertex_index].connect_start5(segs[c1].end3, type_="terminal_crossover")
                    else:
                        segs[c1[::-1]].connect_start3(overhang_segs[vertex_index].start5, type_="terminal_crossover")
                        overhang_segs[vertex_index].connect_start5(segs[c1[::-1]].start3, type_="terminal_crossover")
                    if c2_positive:
                        segs[c2].connect_start5(overhang_segs[vertex_index].start3, type_="terminal_crossover")
                        overhang_segs[vertex_index].connect_start3(segs[c2].start5, type_="terminal_crossover")
                    else:
                        segs[c2[::-1]].connect_end5(overhang_segs[vertex_index].start3, type_="terminal_crossover")
                        overhang_segs[vertex_index].connect_start3(segs[c2[::-1]].end5, type_="terminal_crossover")
                if not overhang_here:
                    if c1_positive and c2_positive:
                        segs[c1].connect_end3(segs[c2].start5, type_="terminal_crossover")
                    elif c1_positive and not c2_positive:
                        segs[c1].connect_end3(segs[c2[::-1]].end5, type_="terminal_crossover")
                    elif not c1_positive and c2_positive:
                        segs[c1[::-1]].connect_start3(segs[c2].start5, type_="terminal_crossover")
                    elif not c1_positive and not c2_positive:
                        segs[c1[::-1]].connect_start3(segs[c2[::-1]].end5, type_="terminal_crossover")

        
    #NICKS
    nicks = False
    if nicks == True:
        for f in face_data:
            nick = f.nick
            if nick in segs:
                segs[nick].add_nick(10,on_fwd_strand=True)
            else:
                print ('failure')

    if overhangs is not None:
        segs_list = [segs[i] for i in segs] + single_stranded_dna + [overhang_segs[i] for i in overhang_segs] + [ss_overhang_segs[i] for i in ss_overhang_segs]
    else:
        segs_list = [segs[i] for i in segs] + single_stranded_dna

    breakpoint()
    return segs_list
    




























