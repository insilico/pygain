#!/usr/bin/env python
from __future__ import division
import sys, os, string, math, csv

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

def cartesian(L,*lists):
	if not lists:
		for x in L:
			yield (x,)
	else:
		for x in L:
			for y in cartesian(lists[0],*lists[1:]):
				yield (x,) + y

def autoproduct(l,n):
	return cartesian(*[l for _ in range(n)])

class Entropy:
	"""An object simulating a memoizing entropy function"""
	def __init__(self, attributes, num_instances, precompute=0):
		"""Initialize memoized entropy function.

		Can specify to precompute to a certain depth. For example,
		precompute = 2 will cache the entropies for all keys, as well
		as all ordered pairs of keys.
		"""
		self.attributes = attributes
		self.num_instances = num_instances
		self.cache = {}

		keys = attributes.keys()
		while precompute > 0:
			for key in autoproduct(keys,precompute):
				self.cache[key] = self.entropy(key)
			precompute -= 1

	def __call__(self, *keys):
		"""Behaves identically to entropy, but returns cached value, if available."""
		try:
			return self.cache[keys]
		except KeyError:
			self.cache[keys] = self.entropy(keys)
			return self.cache[keys]
	
	def entropy(self,keys):
		"""Calculate entropy of attribute(s), given key name(s)"""
		attrs = [self.attributes[key] for key in keys]

		freqs = count(zip(*attrs))

		probs = [freq / self.num_instances for freq in freqs.itervalues()]

		return sum([-p * math.log(p,2) for p in probs])

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
		self.attributes = {}
		for row in data_transpose:
			self.attributes[row[0]] = row[1:]
			# remember when constructing data set with sampled attributes
			# to include status_key at the end

		self.entropy = Entropy(self.attributes,self.num_instances,precompute=2)

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

	def calculate_gain(self, header, cutoff = 0.0):
		"""Computes and returns matrix of GAIN scores.  Scores below the cutoff are set to 0"""
		# only computes upper triangle and diagonal
		num_cols = len(header)
		gain_mat = [] 
		for i in range(num_cols):
			name1 = header[i]
			for j in range(num_cols):
				name2 = header[j]
				if j < i:
					inter_info_score = 0.0
				elif j == i:
					inter_info_score = -1 * self.interaction_information(name1,name2)
				else:
					inter_info_score = self.interaction_information(name1,name2)
				if abs(inter_info_score) < cutoff:
					inter_info_score = 0.0
				gain_mat.append(inter_info_score)
		return gain_mat
   
	def sort_value(self, d, reverse=False):
		""" Returns the keys of dictionary d sorted by their values, default is low to high,
	if parameter reverse = True, function sorts high to low """
		items=d.items()
		backitems=[ [v[1],v[0]] for v in items]
		backitems.sort()
		if reverse:
		   backitems.reverse() # sort from high to low
		return [ backitems[i][1] for i in range(0,len(backitems))]

	def interaction_information(self, attrA, attrB):
		# I(A;B;C)=I(A;B|C)-I(A;B)
		# I(A;B;C)=H(AB)+H(BC)+H(AC)-H(A)-H(B)-H(C)-H(ABC)
		# where C is the class

		H_ABC	= self.entropy(attrA, attrB, self.status_key)
		H_AB	= self.entropy(attrA, attrB)
		H_AC	= self.entropy(attrA, self.status_key)
		H_BC	= self.entropy(attrB, self.status_key)
		H_A	= self.entropy(attrA)
		H_B	= self.entropy(attrB)
		H_C	= self.entropy(self.status_key)
		return H_AB+H_BC+H_AC-H_A-H_B-H_C-H_ABC
		
	def mutual_information(self):
		#print 'class entropy: ', self.entropy(self.status_key)
		attrs_minus_class_keys = [key for key in self.attributes if key != self.status_key]
		self.entropy_dict = {}
		self.mutual_info_dict = {}
		#norm = 0  # if you want to normalize I's
		for key in attrs_minus_class_keys:
			self.entropy_dict[key]	 = self.entropy(key)
			self.mutual_info_dict[key] = self.entropy_dict[key] + self.entropy(self.status_key) - self.entropy(key,self.status_key)
		#	norm = norm + self.mutual_info_dict[key]
		#for key in attrs_minus_class_keys:
		#	self.mutual_info_dict[key] = self.mutual_info_dict[key] / norm
		#print 'entropies: ', self.entropy_dict
		MI_sorted_attrs = self.sort_value(self.mutual_info_dict,reverse=True)
		return (self.mutual_info_dict, MI_sorted_attrs)
		#print 'sorted mutual informations: ', sorted_MI
