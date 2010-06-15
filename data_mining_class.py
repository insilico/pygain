#!/usr/bin/env python
import sys, os, string, re, random, math

class DataProperties(object):
	"""
	Data Mining Class
	Brett McKinney March 06, 2006
	usage:
			from data_mining_class import DataProperties
			data = DataProperties('tab_dlm_file.tab')
			num_attributes = data.num_attributes
			print num_attributes
	(ave_score, sorted_attributes)=data.relieff_analysis(nearest_nbs = 10)
			
	"""
	def __init__(self, infilename):
		# define infilename as an attribute for run_relief
		self.infilename = infilename
		# open file for reading		
		tab_infile = open( infilename,'r')
		# reads all lines at once
		self.data = []
		self.data.extend([line.strip() for line in tab_infile])
		tab_infile.close()
		# creates list of lists of self.data
		self.row_lists = [line.split('\t') for line in self.data]
		# pair up parallel items in lists
		# * has a transposing effect
		# list of attributes (i.e., first row of data except for last element)
		self.data_transpose = zip(*self.row_lists)
		self.attribute_list = self.data[0].split('\t')[0:-1]
		# status symbol
		self.status_key = string.strip(self.data[0].split('\t')[-1])
		# number of attributes and instances
		self.num_instances  = len(self.data)-1
		self.num_attributes = len(self.attribute_list)
		# Create attribute_name -> data dictionary
		self.attribute_dictionary = {}
		for row in self.data_transpose:
			attribute_key = string.strip(row[0])
			data_str = row[1:]
			#print attribute_key
			#print data_str
			self.attribute_dictionary[attribute_key] = data_str
			# remember when constructing data set with sampled attributes
			# to include status_key at the end

	def print_matrix(self, header, mat, colspace = 2):
		"""Prints the GAIN matrix in a nicely formatted fashion."""
		for name in header:
			# formatted numbers take up 7 characters
			# determine spacing based on header column lengths
			# if negative, spaces will be 0
			n_spaces = 7 - len(name)
			sys.stdout.write(name + " " * n_spaces)
			# always add column spaces in between
			sys.stdout.write(" " * colspace)
		print
		for row in range(len(header)):
			for col in range(len(header)):
				n_spaces = len(header[col]) - 7
				# Format matrix values using float formatting rules
				fltfmt = str("%.5f"%mat[row * len(header) + col])
				sys.stdout.write(fltfmt + " " * n_spaces)
				sys.stdout.write(" " * colspace)
			print
   
	def sort_value(self, d, reverse=False):
		""" Returns the keys of dictionary d sorted by their values, default is low to high,
	if parameter reverse = True, function sorts high to low """
		items=d.items()
		backitems=[ [v[1],v[0]] for v in items]
		backitems.sort()
		if reverse:
		   backitems.reverse() # sort from high to low
		return [ backitems[i][1] for i in range(0,len(backitems))]

	def entropy(self,attribute_key):
		# calculate entropy of attribute given its key name
		data_str = self.attribute_dictionary[attribute_key]
		# frequency table
		freq_dict = {}
		for val in data_str:
			if freq_dict.has_key(val):
				freq_dict[val] = freq_dict[val] + 1
			else:
				freq_dict[val] = 1
		# probability table
		prob_dict = {}
		for key in freq_dict:
			prob_dict[key] = float(freq_dict[key])/self.num_instances
		entropy = 0
		for key in prob_dict:
			p = prob_dict[key]
			entropy = entropy - p*math.log(p,2)
		return entropy

	def joint_entropy(self,attribute_key1,attribute_key2):
		# calculate joint entropy given attr key names
		data_str1 = self.attribute_dictionary[attribute_key1]
		data_str2 = self.attribute_dictionary[attribute_key2]
  
		# frequency table
		combo_freq_dict = {}
		for x in zip(data_str1,data_str2):
			combination = ",".join(x)
			if combo_freq_dict.has_key(combination):
				combo_freq_dict[combination] += 1
			else:
				combo_freq_dict[combination] = 1
		# probability table
		combo_prob_dict = {}
		for key in combo_freq_dict:
			combo_prob_dict[key] = float(combo_freq_dict[key])/self.num_instances
		joint_entropy = 0.0
		for key in combo_prob_dict:
			p = combo_prob_dict[key]
			joint_entropy -= p*math.log(p,2)
		return joint_entropy
		
	def interaction_information(self,attribute_key1,attribute_key2):
		# I(A;B;C)=I(A;B|C)-I(A;B)
		# I(A;B;C)=H(AB)+H(BC)+H(AC)-H(A)-H(B)-H(C)-H(ABC)
		# where C is the class
		data_str1  = self.attribute_dictionary[attribute_key1]
		data_str2  = self.attribute_dictionary[attribute_key2]
		data_class = self.attribute_dictionary[self.status_key]

		## Compute H(ABC)
		# frequency table
		combo_freq_dict = {}
		for x in zip(data_str1,data_str2,data_class):
			combination = ",".join(x)
			if combo_freq_dict.has_key(combination):
				combo_freq_dict[combination] = combo_freq_dict[combination] + 1
			else:
				combo_freq_dict[combination] = 1
		# probability table
		combo_prob_dict = {}
		for key in combo_freq_dict:
			combo_prob_dict[key] = float(combo_freq_dict[key])/self.num_instances
		H_ABC = 0.0
		for key in combo_prob_dict:
			p = combo_prob_dict[key]
			H_ABC = H_ABC - p*math.log(p,2)
		H_AB=self.joint_entropy(attribute_key1,attribute_key2)
		H_AC=self.joint_entropy(attribute_key1,self.status_key)
		H_BC=self.joint_entropy(attribute_key2,self.status_key)
		H_A =self.entropy(attribute_key1)
		H_B =self.entropy(attribute_key2)
		H_C =self.entropy(self.status_key)
		return H_AB+H_BC+H_AC-H_A-H_B-H_C-H_ABC
		
	def mutual_information(self):
		#print 'class entropy: ', self.entropy(self.status_key)
		attrs_minus_class_keys = []
		for key in self.attribute_dictionary:
			if key != self.status_key:
				attrs_minus_class_keys.append(key)	
		self.entropy_dict = {}
		self.mutual_info_dict = {}
		#norm = 0  # if you want to normalize I's
		for key in attrs_minus_class_keys:
			self.entropy_dict[key]	 = self.entropy(key)
			self.mutual_info_dict[key] = self.entropy_dict[key] + self.entropy(self.status_key) - self.joint_entropy(key,self.status_key)
		#	norm = norm + self.mutual_info_dict[key]
		#for key in attrs_minus_class_keys:
		#	self.mutual_info_dict[key] = self.mutual_info_dict[key]/float(norm)
		#print 'entropies: ', self.entropy_dict
		MI_sorted_attrs = self.sort_value(self.mutual_info_dict,True)
		return (self.mutual_info_dict, MI_sorted_attrs)
		#print 'sorted mutual informations: ', sorted_MI
