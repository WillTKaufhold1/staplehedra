import numpy as np
with open('NP_design_0.npo','r') as f:

    data = f.readlines()

start = np.arange(len(data))[np.array(data) == '    domains:\n'][0]
stop = np.arange(len(data))[np.array(data) == '    strands:\n'][0]

domains = data[start+1:stop]
stripped_domains = [i.replace(' ','').replace('\n','').split(':') for i in domains]
domains_dict = {}
for i in stripped_domains:
    domains_dict[i[0]] = i[1]
#'sequences:\n'... '    domains:\n'...'    strands:\n'
