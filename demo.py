# -*- coding: utf-8 -*-
import os


def list_dir(rootDir):
    file_list = []

    for lists in os.listdir(rootDir):
        path = os.path.join(rootDir, lists)
        if os.path.isdir(path):
            file_list += list_dir(path)
        else:
            file_list.append(path)
    return file_list


f = open('merge.test.vcf')
nf = open('small.merge.test.vcf','w')
count = 0
for line in f:
    if line[0] == '#':
        nf.write(line)
    else:
        count+=1
        nf.write(line)
        if 'chr2\t92307833' in line:
            break
print line
print count
