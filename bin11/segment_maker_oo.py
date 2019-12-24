import mrdna
from plyfile import PlyData
import numpy as np
from copy import deepcopy as copy

def read_ply(fname):
    polydata = PlyData.read(fname)
    locs = list(zip(polydata['vertex']['x'],polydata['vertex']['y'],polydata['vertex']['z']))
    faces = polydata['face']['vertex_index']
    return locs,faces

class Polyhedron():
    def __init__(self,locs_xyz,faces_index):
        self.edges = {}
        self.faces = {}
        self.connections = {}

        for index,vertices in enumerate(faces_index):
            self.faces[index] = Face(vertices,index)

        for face in self.faces:
            for vertex_0,vertex_1 in self.faces[face].e:
                if ((vertex_1,vertex_0) not in self.edges):
                    self.edges[(vertex_0,vertex_1)] = (Edge(np.array(locs_xyz[vertex_0]),
                                                            np.array(locs_xyz[vertex_1])))

    def normalize(self,bp):
        min_length = min([self.edges[i].rawlen for i in self.edges])
        for i in self.edges:
            self.edges[i].normalize(min_length,bp)

    def apply_connections(self):
        #add connections
        connection_index = 0
        for face in self.faces:
            for A_index, B_index in self.faces[face].connections:
                
                A_dir,B_dir = [index in self.edges for index in [A_index,B_index]] 

                self.connections[connection_index] = (
                        Connection( self.edges[A_index] if A_dir else self.edges[A_index[::-1]],
                                    self.edges[B_index] if B_dir else self.edges[B_index[::-1]],
                                    A_dir,
                                    B_dir))
                                    
                connection_index += 1
    
    def make_mrdna_helices(self):
        for edge in self.edges:
            self.edges[edge].create_helix()
    def generate_mrdna_connections(self,ss_length):
        for connection in self.connections:
            self.connections[connection].apply_connection(ss_length)
            
    def get_mrdna_helices(self):
        helices = []
        for e in self.edges:
            helices.append(self.edges[e].helix)
        return helices
    def get_ss_dna_helices(self):
        helices = []
        for e in self.connections:
            helices.append(self.connections[e].ss)
        return helices

class Connection():
    def __init__(self,A,B,forward_A,forward_B):
        self.A = A # we connect the 3' side of A
        self.B = B # to the 5' side of B
        self.forward_A = forward_A # is this at the front location of A?
        self.forward_B = forward_B # is this at the back location of B?
    def apply_connection(self,ss_length):
        if ss_length == 0:
            self.A.helix.connect_end3(self.B.helix.start5 if self.forward_B else self.B.helix.end5, type_ = "terminal_crossover") if self.forward_A else self.A.helix.connect_start3(self.B.helix.start5 if self.forward_B else self.B.helix.end5, type_ ="terminal_crossover")
        else:
            self.ss = mrdna.SingleStrandedSegment(str(np.random.rand()),
                        start_position = self.A.nStop_c if self.forward_A else self.A.nStart_c,
                        stop_position =  self.B.nStart_c if self.forward_B else self.B.nStop_c ,
                        num_nt=ss_length)

            self.A.helix.connect_end3(self.ss) if self.forward_A else self.A.helix.connect_start3(self.ss)
            self.B.helix.connect_start5(self.ss) if self.forward_B else self.B.helix.connect_end5(self.ss)

class Face():
    def __init__(self,vertices,index):
        self.v = vertices
        self.e = list(zip(self.v,np.roll(self.v,-1)))

        self.index = index
        self.connections = list(zip(self.e,
                                    map(tuple,np.roll(self.e,-1,axis = 0))))
    
class Edge(): 

    def __init__(self,start,stop):

        self.start = start
        self.stop = stop
        self.rawlen = np.sqrt(np.sum((np.array(stop)-np.array(start))**2))  

    def normalize(self,min_length,bp):

        compress_factor = 0.05
        self.nStart = self.start / min_length * bp * 3.4
        self.nStop = self.stop / min_length * bp * 3.4

        vector = self.nStop-self.nStart
        self.nStart_c = self.nStart + vector * compress_factor 
        self.nStop_c  = self.nStop - vector * compress_factor
        self.nLen = np.sqrt(np.sum((np.array(self.nStop)-np.array(self.nStart))**2))
        self.nBp = int(bp*self.rawlen / min_length)

    def create_helix(self):
        self.helix = mrdna.DoubleStrandedSegment(
                                str(np.random.rand()),
                                self.nBp, 
                                start_position= self.nStart_c,
                                end_position =  self.nStop_c)
    
def get_segments(FNAME, LENGTH_OF_SMALLEST, SPACERS,nicks=1,overhangs = None,turberfieldnicks = 0,ldmoverhangs = None):
    locs_xyz,faces_index = read_ply(FNAME)

    poly = Polyhedron(locs_xyz,faces_index)
    poly.normalize(LENGTH_OF_SMALLEST)
    #here we add the overhangs...
    poly.apply_connections()
    
    #here we do the addition

    poly.make_mrdna_helices()
    poly.generate_mrdna_connections(SPACERS)

    helices = poly.get_mrdna_helices()
    ss = poly.get_ss_dna_helices()

    return helices + ss

        