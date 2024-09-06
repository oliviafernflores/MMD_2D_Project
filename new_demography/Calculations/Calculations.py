file = open('/Users/oliviafernflores/Documents/annovar/new_gene_annotation.variant_function', 'r')
types = []
count = 0
exons = 0
for line in file.readlines():
    count += 1
    if line.split('\t')[0] not in types:
        types.append(line.split('\t')[1])
    if line.split('\t')[0] == 'exonic':
        exons += 1
print('in Annovar file: ' + str(types))
print('total lines: ' + str(count))
print('exonic sites: ' + str(exons))
