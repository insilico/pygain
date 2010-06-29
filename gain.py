#!/usr/bin/env python
from __future__ import division
import sys, math, csv
import getopt
from itertools import izip

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

def cmp2((a,x),(b,y)):
	return cmp(x,y)

class Entropy:
	"""An object simulating a memoizing entropy function"""
	def __init__(self, data, class_idx, precompute=0):
		"""Initialize memoized entropy function.

		Can specify to precompute to a certain depth. For example,
		precompute = 2 will cache the entropies for all keys, as well
		as all ordered pairs of keys.
		"""
		self.data = data

		self.num_instances = len(data[0])
		
		# Cache plogp values
		self.plogp = [0] + [-p * math.log(p,2)
				for p in (x / self.num_instances
					 for x in xrange(1, self.num_instances + 1)
				)]

		self.cache = {}

		keys = (lambda x: x[:class_idx] + x[class_idx + 1:])(range(len(data)))
		self.cache[(class_idx,)] = self.entropy((class_idx,))
		while precompute > 0:
			for key in self.uppertrin(keys,precompute):
				self.cache[key] = self.entropy(key)
				self.cache[key + (class_idx,)] = self.entropy(key + (class_idx,))
			precompute -= 1

	def __call__(self, *keys):
		"""Behaves identically to entropy, but returns cached value, if available."""

		if keys in self.cache:
			return self.cache[keys]
		else:
			self.cache[keys] = self.entropy(keys)
			return self.cache[keys]
	
	def entropy(self,keys):
		"""Calculate entropy of attribute(s), given key name(s)"""
		attrs = (self.data[key] for key in keys)

		freqs = self.count(izip(*attrs))

		return sum(self.plogp[freq] for freq in freqs)

	def count(self,xs):
		"""Return a list of counts of objects in a list"""
		ret = {}
		for x in xs:
			if x in ret:
				ret[x] += 1
			else:
				ret[x] = 1

		return ret.itervalues()
	
	def uppertrin(self,L,n):
		"""Return the upper-triangle of an n-dimensional matrix L^n
	
		Returns all tuples of n elements chosen out of L, with
		replacement, but without permutations. In 2-dimensions,
		this is the upper-triangular matrix, with the diagonal.
		"""
		def tails(tail):
		        while tail:
		                head, tail = tail[0], tail[1:]
		                yield (head, tail)
	        if n == 1:
	                for l in L:
	                        yield (l,)
	        else:
	                for a,t in tails(L):
	                        for b in self.uppertrin(t, n - 1):
	                                yield (a,) + b

class GAIN:
	def __init__(self, infile):
		reader = csv.reader(infile, delimiter=" ")

		# CSV readers are iterable, so just pull all of our data in at once.
		self.attributes = reader.next()[5:]
		rows = [map(self.translate,row[5:]) for row in reader]

		self.data = transpose(rows)

		self.class_idx = 0

		self.entropy = Entropy(self.data,self.class_idx,precompute=2)

	def translate(self,x):
		if x == "NA":
			return None
		return int(x)

	def print_tsv(self, outfile, idcs, mat, ndigits=5):
		"""Print the GAIN matrix in tab-separated value format"""
		format = lambda f: "%%.%df" % ndigits % f

		reordered = [[format(mat[i][j]) for i in idcs] for j in idcs]
		reordered.insert(0,[self.attributes[i] for i in idcs])
		reordered.append([])

		outfile.write('\n'.join('\t'.join(row) for row in reordered))
		return

	def print_matrix(self, outfile, idcs, mat, colspace = 2, ndigits = 5):
		"""Prints the GAIN matrix in a nicely formatted fashion."""
		for i in idcs:
			name = self.attributes[i]
			# formatted numbers take up ndigits + 2 characters (for 0. in the number )
			# determine spacing based on header column lengths
			# if negative, spaces will be 0
			n_spaces = (ndigits + 2) - len(name)
			outfile.write(name + " " * n_spaces)
			# always add column spaces in between
			outfile.write(" " * colspace)
		outfile.write('\n')
		for i in idcs:
			for j in idcs:
				name = self.attributes[j]
				n_spaces = len(name) - (ndigits + 2) 
				# Format matrix values using float formatting rules
				fltstr = "%." + str(ndigits) + "f"
				fltfmt = fltstr % mat[i][j]
				outfile.write(fltfmt + " " * n_spaces)
				outfile.write(" " * colspace)
			outfile.write('\n')

	def calculate_gain(self, cutoff = 0.0):
		"""Computes and returns matrix of GAIN scores.  Scores below the cutoff are set to 0"""
		attrs = len(self.attributes)

		gain_mat = [[0.0 for j in range(attrs)] for i in range(attrs)]
		for i in range(attrs):
			if i == self.class_idx:
				continue

			gain_mat[i][i] = -self.autointeraction(i)

			for j in range(i+1,attrs):
				if j == self.class_idx:
					continue
				gain_mat[i][j] = self.interaction_information(i,j)
				gain_mat[j][i] = gain_mat[i][j]

		return gain_mat

	def sort_value(self, d, reverse=False):
		"""Returns the keys of dictionary d sorted by their values, default is low to high,
	if parameter reverse = True, function sorts high to low """
		items = sorted(d,cmp2)
		if reverse:
		   items.reverse() # sort from high to low
		return [item[0] for item in items]

	def interaction_information(self, attrA, attrB):
		# I(A;B;C)=I(A;B|C)-I(A;B)
		# I(A;B;C)=H(AB)+H(BC)+H(AC)-H(A)-H(B)-H(C)-H(ABC)
		# where C is the class

		H_ABC	= self.entropy(attrA, attrB, self.class_idx)
		H_AB	= self.entropy(attrA, attrB)
		H_AC	= self.entropy(attrA, self.class_idx)
		H_BC	= self.entropy(attrB, self.class_idx)
		H_A	= self.entropy(attrA)
		H_B	= self.entropy(attrB)
		H_C	= self.entropy(self.class_idx)
		return H_AB+H_BC+H_AC-H_A-H_B-H_C-H_ABC

	def autointeraction(self, attr):
		return self.entropy(attr, self.class_idx) - self.entropy(attr) - self.entropy(self.class_idx)

	def mutual_information(self):

		idcs = xrange(len(self.attributes))

		return self.sort_value(
			((i,-self.autointeraction(i))
				for i in idcs
					if i != self.class_idx),
			reverse=True)

def main(argv):
	help = """Usage: %s [OPTIONS]

Construct GAIN matrix from PLINK RAW file

Options:
    --input	-i	Input file (default: stdin)
    --output	-o	Output file (default: stdout)
    --help		display this help and exit
	""" % argv.pop(0).split('/')[-1]

	try:
		opts, args = getopt.getopt(argv, "i:o:h",
			["input=","output=","help"])
	except getopt.error, msg:
		print msg
		return 0

	infile = sys.stdin
	outfile = sys.stdout

	# Parse arguments
	for opt, arg in opts:
		if opt in ('-i', '--input'):
			infile = open(arg)
		if opt in ('-o', '--output'):
			output = open(arg, 'w')
		if opt in ('-h','--help'):
			print help
			return 0

	gain = GAIN(infile)
	gmatrix = gain.calculate_gain()
	ranked_attrs = gain.mutual_information()

	gain.print_tsv(outfile,ranked_attrs,gmatrix)

if __name__ == '__main__':
	sys.exit(main(sys.argv))
