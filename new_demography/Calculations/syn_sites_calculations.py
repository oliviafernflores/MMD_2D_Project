file = open('/Users/oliviafernflores/Documents/annovar/new_gene_annotation.exonic_variant_function', 'r')
syn = 0
nsyn = 0
count = 0
for line in file.readlines():
    count += 1
    if line.split('\t')[1] == 'synonymous SNV':
        syn += 1
    if line.split('\t')[1] == 'nonsynonymous SNV':
        nsyn += 1
print('total exonic sites: ' + str(count))
print('syn sites: ' + str(syn))
print('nsyn sites: ' + str(nsyn))