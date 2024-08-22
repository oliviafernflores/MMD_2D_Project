'''
Mu = mutation rate * # syn sites / total sites
L = total # sites syn and nsyn * (# nsyn sites / # syn sites) 

mutation rate: 5.4e-9 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4509997/
'''

import pickle

dd_syn = pickle.load(open('/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/data_dictionaries/mus_all_pops.dd_syn_with_ancestral.bpkl', 'rb'))
dd_nsyn = pickle.load(open('/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/data_dictionaries/mus_all_pops.dd_nsyn_with_ancestral.bpkl', 'rb'))

syn_sites = len(dd_syn)
nsyn_sites = len(dd_nsyn)

total_sites = syn_sites + nsyn_sites

syn_divided_by_total = syn_sites / total_sites

nsyn_divided_by_syn = nsyn_sites / syn_sites

L = total_sites * nsyn_divided_by_syn

mutation_rate = 5.4e-9

mu = mutation_rate * syn_divided_by_total

f = open('calculation_mu_L_critical_values.txt', 'w')

f.write('Number of syn sites: ' + str(syn_sites) + '\n')
f.write('Number of nsyn sites: ' + str(nsyn_sites) + '\n')
f.write('Number of total sites: ' + str(total_sites) + '\n')
f.write('Number of syn sites divided by number of total sites: ' + str(syn_divided_by_total) + '\n')
f.write('Mutation rate: ' + str(mutation_rate) + '\n')
f.write('L: number of total sites * nsyn/syn ratio: ' + str(L) + '\n')
f.write('Mu: mutation rate times fraction of syn sites (number of syn sites divided by number of total sites): ' + str(mu) + '\n')


out = open('mu_L.txt', 'w')
out.write('L' + '\t' + 'Mu' + '\n')
out.write(str(L) + '\t' + str(mu))