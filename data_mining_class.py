#!/usr/bin/env python
from __future__ import division
import sys, os, string, re, random, math
import csv

def transpose(a):
	"""Transpose a list of lists"""
	# This is mid-level magic. A function can take an arbitrary list of
	# arguments like so:
	#	def function(*args):
	#		...
	# On the calling end, we can prepend an asterisk to an iterable object
	# to indicate we want its contents, rather than itself, to be the arguments
	# to a function.
	# 
	# zip() returns a list where the nth element is a tuple of the nth
	# elements of each argument.
	#
	# Thus, the first "row" of output from zip(*a) is the first element of
	# each list in a.
	return zip(*a)

def count(xs):
	"""Return a dictionary with a count for each object in xs"""
	ret = {}
	for x in xs:
		if x in ret:
			ret[x] += 1
		else:
			ret[x] = 1
	return ret

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
		tsv = open( infilename,'r')
		reader = csv.reader(tsv, delimiter="\t")

		# CSV readers are iterable, so just pull all of our data in at once.
		data = [row for row in reader]

		tsv.close()

		data_transpose = transpose(data)

		# status symbol
		self.status_key = string.strip(data[0][-1])

		self.num_instances  = len(data) - 1

		# Create attribute_name -> data dictionary
		self.attribute_dictionary = {}
		for row in data_transpose:
			attribute_key = row[0]
			data_str = row[1:]
			#print attribute_key
			#print data_str
			self.attribute_dictionary[attribute_key] = data_str
			# remember when constructing data set with sampled attributes
			# to include status_key at the end

	def print_matrix(self, header, mat, colspace = 2, ndigits = 5):
		"""Prints the GAIN matrix in a nicely formatted fashion."""
		for name in header:
			# formatted numbers take up ndigits + 2 characters (for 0. in the number )
			# determine spacing based on header column lengths
			# if negative, spaces will be 0
			n_spaces = (ndigits + 2) - len(name)
			sys.stdout.write(name + " " * n_spaces)
			# always add column spaces in between
			sys.stdout.write(" " * colspace)
		print
		for row in range(len(header)):
			for col in range(len(header)):
				n_spaces = len(header[col]) - (ndigits + 2) 
				# Format matrix values using float formatting rules
				fltstr = "%." + str(ndigits) + "f"
				fltfmt = str(fltstr%mat[row * len(header) + col])
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
		"""Calculate entropy of attribute, given its key name"""

		# frequency table
		freq_dict = count(self.attribute_dictionary[attribute_key])

		# probabilities
		probs = [freq / self.num_instances for freq in freq_dict.itervalues()]

		entropy = sum([-p * math.log(p,2) for p in probs])

		return entropy

	def joint_entropy(self,attribute_key1,attribute_key2):
		# calculate joint entropy given attr key names
		data_str1 = self.attribute_dictionary[attribute_key1]
		data_str2 = self.attribute_dictionary[attribute_key2]
  
		# frequency table
		combo_freq_dict = count(zip(data_str1,data_str2))

		probs = [freq / self.num_instances for freq in combo_freq_dict.itervalues()]

		entropy = sum([-p * math.log(p,2) for p in probs])

		return entropy
		
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
			combo_prob_dict[key] = combo_freq_dict[key] / self.num_instances
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
		#	self.mutual_info_dict[key] = self.mutual_info_dict[key] / norm
		#print 'entropies: ', self.entropy_dict
		MI_sorted_attrs = self.sort_value(self.mutual_info_dict,True)
		return (self.mutual_info_dict, MI_sorted_attrs)
		#print 'sorted mutual informations: ', sorted_MI
