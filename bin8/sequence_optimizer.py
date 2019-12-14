import numpy as np


class domain:
    def __init__(self,seq,index,num_nt):
        self.seq = seq
        self.index = index
        self.num_nt = num_nt

def write_NUPACK_files(domains):
    with open('NP_design.np','w+') as f:
        f.write('material = dna\ntemperature = 25.0\nseed = 93\n')
        #define domains 
        
        for d in domains:
            seq = d.seq if type(d.seq) != type(None) else 'N%s'%(d.num_nt,)
            f.write("domain d_%s = %s\n"%(d.index,seq))

        long_domains = list(filter(lambda x: x.num_nt > 8,domains))

        #define strands from domains
        for d in long_domains:
            f.write("strand s_%s = d_%s\n"%(d.index,d.index,))
            f.write("strand s_%s_comp = d_%s*\n"%(d.index,d.index,))

        #define complexes from strands
        for d in long_domains:
            f.write("complex c_%s = s_%s s_%s_comp\n"%(d.index,d.index,d.index,))
            f.write("complex c_lonely_%s = s_%s\n"%(d.index,d.index,))

        #define target structures for each complex
        for d in long_domains:
            f.write("c_%s.structure = D%s+\n"%(d.index,d.num_nt))
            f.write("c_lonely_%s.structure = U%s\n"%(d.index,d.num_nt,))
        # define list of all dsDNA
        all_complexes = ''
        for d in long_domains:
            all_complexes += "c_%s "%(d.index,)

        #time to define the tubes 

        #CROSSTALK TUBE
        f.write('tube Crosstalk = %s\n'%(all_complexes,))
        for d in long_domains:
            f.write("Crosstalk.c_%s.conc[M] = 1e-5\n"%(d.index,))
        f.write("Crosstalk.offtargets = {maxsize = 2}\n")


        '''
        #EVERY INDIVIDUAL TUBE... this is our previous way where we considered mostly orthogonal domains.
        for d in long_domains:
            f.write('tube Lonely_%s = %s\n'%(d.index,"c_lonely_%s"%(d.index)))
            f.write('Lonely_%s.offtargets = {maxsize = 2}\n'%(d.index,))
            #giving concentrations here leads to a segfault... I don't know why...
        '''
        #PREVENT
        #DOMAINS!

        f.write("prevent = AAAA, CCC, GGG, UUUU, KKKKKK, MMMMMM, RRRRRR, SSSSSS, WWWWWW, YYYYYY\n")

        #objective
        f.write("stop[%] = 5\n")

def read_nupack():

    with open('NP_design_0.npo','r') as f:
        data = f.readlines()

    start = np.arange(len(data))[np.array(data) == '    domains:\n'][0]
    stop = np.arange(len(data))[np.array(data) == '    strands:\n'][0]

    domains = data[start+1:stop]
    stripped_domains = [i.replace(' ','').replace('\n','').split(':') for i in domains]
    domains_dict = {}
    for i in stripped_domains:
        domains_dict[i[0]] = i[1]

    return domains_dict


def optimize_sequence(segs_list):
    import os

    #TODO: ignore really short sequences...

    domains = []

    for index,seg in enumerate(segs_list):
        domains.append(
            domain(seg.sequence,index,seg.num_nt))

    write_NUPACK_files(domains)
    os.system('multitubedesign NP_design.np')
    
    nupack_domain_dict = read_nupack()

    for domain_name in nupack_domain_dict:
        if '*' not in domain_name:
            index = int(domain_name.split('_')[1])
            segs_list[index].sequence = nupack_domain_dict[domain_name]

    '''
    for i in segs_list: #there's one domain for every segment
        i.sequence = 'T'*i.num_nt
    '''
    
    breakpoint()
    


















