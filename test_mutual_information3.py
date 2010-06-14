#!/usr/bin/env python
import sys, os, string, re, random, math
from data_mining_class import DataProperties

try:
    infilename  = sys.argv[1]
    #outfilename = sys.argv[2]
except:
    print "Usage:",sys.argv[0], "in_all_data"
    print "Example:",sys.argv[0], "mutual_info_testset.tab"
    sys.exit(1)
# python test_data_mining.py cytokines.tab > output.txt

#print 'Output created by python ', sys.argv[0], ' ', infilename
full_data = DataProperties(infilename)
#print 'instances = ',        full_data.num_instances
#print 'attributes = ',       full_data.num_attributes
#print 'status key = ',       full_data.status_key

(mutual_info_dict, MI_sorted) = full_data.mutual_information()
#print '\t mutual info \t chi-square \n'
#for i in range(full_data.num_attributes):
#    print i,'\t',MI_sorted[i],'\t',chi2_sorted[i] 

#print mutual_info_dict

# print matrix pair-wise class interaction informations
var_name_list = []
for key in mutual_info_dict:
    var_name_list.append(key)

# define a interaction information cutoff, below which is considered non-interacting
cutoff=0.0
# upper triangle and diagonal only
num_vars=len(var_name_list)
name1_interactions=[]
for i in range(num_vars):
    name1=var_name_list[i]
    for j in range(num_vars):
        name2=var_name_list[j]
        if j<i:
            inter_info_score = 0.0
        elif j==i:
            inter_info_score=-1*(full_data.interaction_information(name1,name2))
        else:
            inter_info_score=full_data.interaction_information(name1,name2)
        if abs(inter_info_score) < cutoff:
            inter_info_score=0.0
        name1_interactions.append(inter_info_score)

# print matrix to console
full_data.printMatrix(var_name_list, name1_interactions)

# directed
#for name1 in var_name_list:
#    name1_interactions=[]
#    for name2 in var_name_list:
#        inter_info_score=full_data.interaction_information(name1,name2)
#        if abs(inter_info_score) < cutoff:
#            inter_info_score = 0.0
#        if name1!=name2:
#            name1_interactions.append(str(inter_info_score))
#        else:
#            # define class interaction information between variable and itself to be 0
#            name1_interactions.append('0.0')
#    print '\t'.join(name1_interactions)





