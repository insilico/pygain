#!/usr/bin/env python
import sys, os, string, math
from data_mining_class import DataProperties

try:
	infilename  = sys.argv[1]
except:
	print "Usage: ",sys.argv[0], " in_all_data"
	print "Example: ",sys.argv[0], " mutual_info_testset.tab"
	sys.exit(1)

# Create DataProperties object from data file
full_data = DataProperties(infilename)

# compute GAIN
mat = full_data.calculate_gain()

# Sort attributes by main effect
MI_sorted = full_data.mutual_information()

# print matrix to console
full_data.print_matrix(MI_sorted, mat)
