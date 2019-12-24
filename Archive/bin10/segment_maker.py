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
            self.nBp = int(bp*self.rawlen / min_length)
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

def tn_get_graph(face_data):
    # WORK IN PROGRESS
    # get vertices
    # I don't think this is quite the right thing....
    V_vals = [i.v for i in face_data]
    graph = {}
    for face_index,face in enumerate(V_vals):
        graph[face_index] = []
        for vertex in face:
            for other_vertex_index, other_vertex_list in enumerate(V_vals):
                if (vertex in other_vertex_list):
                    if other_vertex_index not in graph[face_index]:
                        graph[face_index].append(other_vertex_index)
    return graph

def tn_get_break_index(face_data):
    #essentially the same algorithm... as for assigning edge breaks
    graph = tn_get_graph(face_data)
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


def get_face_data(FNAME):

    locs,faces = read_ply(FNAME)

    face_data = []
    for index,i in enumerate(faces):
        face_data.append(FACE(i,index))
    
    return face_data,locs

def get_mentioned_edges(face_data):

    nick_locs = []
    for f in face_data:
        nick_locs.append(f.nick)
    mentioned_edges = set(nick_locs)
    
    return mentioned_edges

def do_orientation_assignment(face_data,mentioned_edges):
    for i in face_data:
        for e in i.e:
            if e not in mentioned_edges and e[::-1] not in mentioned_edges:
                i.ori.append(1)
                mentioned_edges.add(e)
            else:
                i.ori.append(-1)
    return mentioned_edges


def get_overhang_edges(overhangs,locs,min_length,LENGTH_OF_SMALLEST):
    ss_overhang_sequences = {}
        
    overhang_ds_edges = {}
    overhang_ss_edges = {}

    for i in range(len(overhangs)):
        vertex = overhangs.iloc[i]['overhang'] 
        vertex_location = locs[vertex] 
        bp = overhangs.iloc[i]['ds_length'] 
        length = bp * min_length / float(LENGTH_OF_SMALLEST) / 3.4 
        end = length * np.array(vertex_location) / np.sqrt(np.sum(np.array(vertex_location)**2)) + np.array(vertex_location)
        overhang_ds_edges[vertex] = edge(vertex,vertex_location,end,overhang=True,bp_overhang=bp)

        #ss overhangs        

        bp_ss = overhangs.iloc[i]['ss_length'] 
        length_ss = bp_ss * min_length / float(LENGTH_OF_SMALLEST) / 3.4 
        start_ss = end 
        end_ss = end + length_ss * np.array(vertex_location) / np.sqrt(np.sum(np.array(vertex_location)**2))
        overhang_ss_edges[vertex] = edge(vertex,start_ss,end_ss,overhang=True,bp_overhang=bp_ss,out_orientation=overhangs.iloc[i]['out_side'])

        ss_overhang_sequences[vertex] = overhangs.iloc[i]['ss_seq']

    return (overhang_ds_edges, overhang_ss_edges, ss_overhang_sequences)


def get_segs(edges):

    segs = {}

    for index, e in enumerate(edges):
        segs[tuple(e.index)] = (
            mrdna.DoubleStrandedSegment(name = 'helix%s'%(index,),
                                num_bp = e.nBp,
                                start_position = e.nStart_c,  
                                end_position = e.nStop_c
                                       ))

    return segs

def get_overhang_segs(overhang_ds_edges,overhang_ss_edges,ss_overhang_sequences):
    overhang_segs = {}
    for index in overhang_ds_edges:
        e = overhang_ds_edges[index] 
        tmp = copy((mrdna.DoubleStrandedSegment(name = 'helix%s'%(index+100,),
                                num_bp = int(e.nBp),
                                start_position = list(e.nStart_c),  
                                end_position = list(e.nStop_c)
                                    )))
        overhang_segs[e.index] = tmp

    ss_overhang_segs = {}
    for index in overhang_ss_edges:
        e = overhang_ss_edges[index]
        if e.out_orientation == 5:
            tmp = copy((mrdna.SingleStrandedSegment(name = 'helix%s'%(index+200,),
                                    num_nt = int(e.nBp),
                                    start_position = list(e.nStart_c),  
                                    end_position = list(e.nStop_c)
                                        )))
        
        elif e.out_orientation == 3:
            tmp = copy((mrdna.SingleStrandedSegment(name = 'helix%s'%(index+200,),
                                    num_nt = int(e.nBp),
                                    start_position = list(e.nStop_c),
                                    end_position = list(e.nStart_c),  
                                        )))
        tmp.sequence = ss_overhang_sequences[index]
        ss_overhang_segs[e.index]= (tmp)

    return (overhang_segs,ss_overhang_segs)

def connect_ssdna_to_overhangs(overhang_ss_edges,ss_overhang_segs,overhang_ds_edges,overhang_segs):
    for index in overhang_ss_edges:
        #connected if they have the same index!
        ss = overhang_ss_edges[index]
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

def is_there_an_overhang_here(f,con,overhangs):
    
    overhang_here = False

    ss_here = 0
    vertex_index = intersection(con[0],con[1])[0] 
    face_index = f.index    

    if overhangs is not None:

        for i in range(len(overhangs)):
            if (overhangs.iloc[i]['face'] == face_index and overhangs.iloc[i]['overhang'] == vertex_index):
                overhang_here = True
                ss_here = overhangs.iloc[i]['ss_extra']

    return overhang_here,ss_here

def join_overhang(ss_here, c1_positive, c2_positive, overhang_segs, segs, c1, c2, vertex_index):

    single_stranded_dna = []
    
    if ss_here == 0: 
        
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
        
    elif ss_here != 0:
        # add the connectivity here!

        if c1_positive:

            ss = copy(mrdna.SingleStrandedSegment("strand%s"%(np.random.randint(0,high = int(1e10)),),
                start_position = segs[c1].end_position ,
                end_position =   overhang_segs[vertex_index].start_position,
                num_nt = ss_here))

            segs[c1].connect_end3(ss)
            overhang_segs[vertex_index].connect_start5(ss)
            single_stranded_dna.append(ss)

        else:

            ss = copy(mrdna.SingleStrandedSegment("strand%s"%(np.random.randint(0,high = int(1e10)),),
                start_position = segs[c1[::-1]].start_position ,
                end_position =   overhang_segs[vertex_index].start_position,
                num_nt = ss_here))

            segs[c1[::-1]].connect_start3(ss)
            overhang_segs[vertex_index].connect_start5(ss)
            single_stranded_dna.append(ss)

        if c2_positive:

            ss = copy(mrdna.SingleStrandedSegment("strand%s"%(np.random.randint(0,high = int(1e10)),),
                start_position = overhang_segs[vertex_index].end_position ,
                end_position =   segs[c2].start_position,
                num_nt = ss_here))

            overhang_segs[vertex_index].connect_start3(ss)
            segs[c2].connect_start5(ss)
            single_stranded_dna.append(ss)

        else:

            ss = copy(mrdna.SingleStrandedSegment("strand%s"%(np.random.randint(0,high = int(1e10)),),
                start_position = segs[c2[::-1]].end_position ,
                end_position =   overhang_segs[vertex_index].start_position,
                num_nt = ss_here))

            overhang_segs[vertex_index].connect_start3(ss)
            segs[c2[::-1]].connect_end5(ss)
            single_stranded_dna.append(ss)

    return single_stranded_dna

def join_segments(c1_positive, c2_positive, segs, c1, c2, ssDNA_index,SPACERS):

    single_stranded_dna = []


    if SPACERS !=0 :

        r1 = np.random.rand(1)[0] - 0.5
        r2 = np.random.rand(1)[0] - 0.5

        if c1_positive and c2_positive:
            ss = copy(mrdna.SingleStrandedSegment("strand%s"%(ssDNA_index,),
                start_position = segs[c1].end_position + r1,
                end_position =   segs[c2].start_position + r2,
                num_nt = SPACERS))

            segs[c1].connect_end3(ss)
            segs[c2].connect_start5(ss)

        elif c1_positive and not c2_positive:
            ss = copy(mrdna.SingleStrandedSegment("strand%s"%(ssDNA_index,),
                start_position = segs[c1].end_position + r1,
                end_position =   segs[c2[::-1]].end_position + r2,
                num_nt = SPACERS))

            segs[c1].connect_end3(ss)
            segs[c2[::-1]].connect_end5(ss)

        elif not c1_positive and c2_positive:

            ss = copy(mrdna.SingleStrandedSegment("strand%s"%(ssDNA_index,),
                start_position = segs[c1[::-1]].start_position +r1,
                end_position =   segs[c2].start_position + r2,
                num_nt = SPACERS))

            segs[c1[::-1]].connect_start3(ss)
            segs[c2].connect_start5(ss)

        elif not c1_positive and not c2_positive:

            ss = copy(mrdna.SingleStrandedSegment("strand%s"%(ssDNA_index,),
                start_position = segs[c1[::-1]].start_position +r1 ,
                end_position =   segs[c2[::-1]].end_position + r2,
                num_nt = SPACERS))      

            segs[c1[::-1]].connect_start3(ss)
            segs[c2[::-1]].connect_end5(ss)

        single_stranded_dna.append(ss)

    else:
        if c1_positive and c2_positive:
            segs[c1].connect_end3(segs[c2].start5, type_="terminal_crossover")
        elif c1_positive and not c2_positive:
            segs[c1].connect_end3(segs[c2[::-1]].end5, type_="terminal_crossover")
        elif not c1_positive and c2_positive:
            segs[c1[::-1]].connect_start3(segs[c2].start5, type_="terminal_crossover")
        elif not c1_positive and not c2_positive:
            segs[c1[::-1]].connect_start3(segs[c2[::-1]].end5, type_="terminal_crossover")
    
    return single_stranded_dna


def join_spacerless_overhangs(ss_here, c1_positive, c2_positive, overhang_segs, segs, c1, c2, vertex_index):

    single_stranded_dna = []

    if ss_here == 0:
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
    else:

        if c1_positive:

            ss = copy(mrdna.SingleStrandedSegment("strand%s"%(np.random.randint(0,high = int(1e10)),),
                start_position = segs[c1].end_position ,
                end_position =   overhang_segs[vertex_index].start_position,
                num_nt = ss_here))

            segs[c1].connect_end3(ss)
            overhang_segs[vertex_index].connect_start5(ss)
            single_stranded_dna.append(ss)

        else:

            ss = copy(mrdna.SingleStrandedSegment("strand%s"%(np.random.randint(0,high = int(1e10)),),
                start_position = segs[c1[::-1]].start_position ,
                end_position =   overhang_segs[vertex_index].start_position,
                num_nt = ss_here))

            segs[c1[::-1]].connect_start3(ss)
            overhang_segs[vertex_index].connect_start5(ss)
            single_stranded_dna.append(ss)

        if c2_positive:

            ss = copy(mrdna.SingleStrandedSegment("strand%s"%(np.random.randint(0,high = int(1e10)),),
                start_position = overhang_segs[vertex_index].end_position ,
                end_position =   segs[c2].start_position,
                num_nt = ss_here))

            overhang_segs[vertex_index].connect_start3(ss)
            segs[c2].connect_start5(ss)
            single_stranded_dna.append(ss)

        else:

            ss = copy(mrdna.SingleStrandedSegment("strand%s"%(np.random.randint(0,high = int(1e10)),),
                start_position = segs[c2[::-1]].end_position ,
                end_position =   overhang_segs[vertex_index].start_position,
                num_nt = ss_here))

            overhang_segs[vertex_index].connect_start3(ss)
            segs[c2[::-1]].connect_end5(ss)
        
            single_stranded_dna.append(ss)
            
    return single_stranded_dna



def get_segments(FNAME, LENGTH_OF_SMALLEST, SPACERS,nicks=1,overhangs = None,turberfieldnicks = 0,ldmoverhangs = None):

    tn = turberfieldnicks
    lo = ldmoverhangs

    face_data,locs = get_face_data(FNAME)

    
    assign_nicks(face_data)
    mentioned_edges = get_mentioned_edges(face_data)

    mentioned_edges = do_orientation_assignment(face_data,mentioned_edges)

    edges = [edge(i,locs[i[0]],locs[i[1]]) for i in mentioned_edges]

    min_length = np.min([x.rawlen for x in edges])

    if overhangs is not None:
        overhang_ds_edges,overhang_ss_edges,ss_overhang_sequences = get_overhang_edges(overhangs,locs,min_length,LENGTH_OF_SMALLEST)
        for i in edges + list(overhang_ds_edges.values()) + list(overhang_ss_edges.values()):
            i.normalize(min_length,LENGTH_OF_SMALLEST)
   
    else:
        for i in edges:
            i.normalize(min_length,LENGTH_OF_SMALLEST)


    segs = get_segs(edges)

    if overhangs is not None:
        overhang_segs,ss_overhang_segs = get_overhang_segs(overhang_ds_edges,overhang_ss_edges,ss_overhang_sequences)

    if overhangs is not None:
        connect_ssdna_to_overhangs(overhang_ss_edges,ss_overhang_segs,overhang_ds_edges,overhang_segs)

    single_stranded_dna = []
    ssDNA_index = 0


    ldm_breakable_segs = []

    for f in face_data:
        for con in f.connections:

            vertex_index = intersection(con[0],con[1])[0] 
            face_index = f.index    

            overhang_here,ss_here = is_there_an_overhang_here(f,con,overhangs)

            ssDNA_index += 1

            c1,c2 = con
            c1_positive = c1 in segs
            c2_positive = c2 in segs

            ldm_break = False
            if type(ldmoverhangs) != type(None):
                for i in range(len(ldmoverhangs)):
                    if (ldmoverhangs.iloc[i]['face'] == face_index and ldmoverhangs.iloc[i]['overhang'] == vertex_index):
                        ldm_break = True

            if overhang_here:
                single_stranded_dna.extend(join_overhang(ss_here, c1_positive, c2_positive, overhang_segs, segs, c1, c2, vertex_index))

            else:
                if not ldm_break:
                    single_stranded_dna.extend(join_segments(c1_positive, c2_positive, segs, c1, c2,ssDNA_index,SPACERS))
                elif ldm_break:
                    #here we add our specific sequence... maybe just add a sequence, then some nones?
                    #work out which segs need to be broken up
                    segs_broken = []
                    for i in f.e:
                        if vertex_index in i:
                            if i[-1] == vertex_index:
                                segs_broken.append(i)
                            else:
                                segs_broken.append(i[::-1])
                    
                    ldm_breakable_segs.append(segs_broken)
                    #reorient the segs broken appropriately, so the vertex is at the end.


    #break those segs!
    ldm_replacements = []
    
    for pair in ldm_breakable_segs:
        for seg_name in pair:
            if seg_name in segs:
                print ('Woo')
                to_be_broken = segs[seg_name]

                segs.pop(seg_name)

                s = to_be_broken.start_position
                e = to_be_broken.end_position
                nts = to_be_broken.num_nt

                overhang_length = 14 #obviously actually find this from the file.

                # START SIDE.

                start_1 = s
                end_1 = s + (e-s) * (nts - overhang_length) / float(nts)
                nbp_1 = nts - overhang_length

                # OVERHANG SIDE.

                start_2 = end_1
                end_2 = e
                nbp_2 = overhang_length

                start_seg = mrdna.DoubleStrandedSegment(name = 'helix%s'%(np.random.rand()),
                                num_bp = nbp_1,
                                start_position = start_1,  
                                end_position = end_1
                                       )

                end_seg = mrdna.DoubleStrandedSegment(name = 'helix%s'%(np.random.rand()),
                                num_bp = nbp_2,
                                start_position = start_2,
                                end_position = end_2
                                        )

                start_seg.connect_end3(end_seg.start5) #only works if direction of strand is anticlockwise around face
                #start_seg.connect_end5(end_seg.start3) #only works if direction of strand is anticlockwise around face

                #now we reconnect this guy


                breakpoint()
                #TODOTODOTODOTODOTODO

                print (to_be_broken)
                for connection in to_be_broken.connections:
                    A = connection.A 
                    B = connection.B
                    breakpoint()
                    if B.container == to_be_broken:
                        if B.on_fwd_strand:
                            A.container.connections.append(
                                mrdna.segmentmodel.Connection(A,end_seg.start5, type_=connection.type_))
                        else:
                            A.container.connections.append(
                                mrdna.segmentmodel.Connection(A.end3,start_seg, type_=connection.type_))
                        
                    else:
                        pass
                        '''
                        if A.on_fwd_strand:
                            #connect out end
                            B.container.connections.append(
                                mrdna.segmentmodel.Connection(B,end_seg, type_=connection.type_))
                        else:
                            B.container.connections.append(
                                mrdna.segmentmodel.Connection(B,start_seg, type_=connection.type_))
                        '''

                    breakpoint()
               

                #reconnect prexisting connections!



                ldm_replacements.append(start_seg)
                ldm_replacements.append(end_seg)
                #breakpoint()

            else:
                print ('Boo!')


    no_nick_faces = []
    if type(overhangs) != type(None):
    #if overhangs != None:
        no_nick_faces = list(overhangs['face'])

    if nicks:
        print ("adding nicks!")
        for f in face_data:

            nick = f.nick
            if ((nick in segs) and (f.index not in no_nick_faces)):
                segs[nick].add_nick(10,on_fwd_strand=True)

            else:
                print ('failure')

    #make sure that all the ssDNA has the sequence 'TTT...'

    for s in single_stranded_dna:
        s.sequence = s.num_nt * 'T'

    if overhangs is not None:
        segs_list = [segs[i] for i in segs] + single_stranded_dna + [overhang_segs[i] for i in overhang_segs] + [ss_overhang_segs[i] for i in ss_overhang_segs]

    else:
        segs_list = [segs[i] for i in segs] + single_stranded_dna


    return segs_list
    




