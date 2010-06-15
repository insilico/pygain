import sys, os, string, re, random, math
from data_mining_class import DataProperties

	def relieff_analysis(data_prop, **kwds):
		"""
		kwd defaults: nearest_nbs = 10 (./relieff nearest neighbors)
		usage: (score, sorted_attributes)=data.relieff_analysis(nearest_nbs= 5)
		
		./relieff creates output file called infilename.relieff, which
		we open after running in order to read in the relief scores
		and zip with attributes and sort. Returns total relief score
		and relieff-sorted attributes from data_prop.data.
		"""	  
		if kwds.has_key('nearest_nbs') == False:
			kwds['nearest_nbs'] = 10
		relieff_cmd = './relieff ' + data_prop.infilename + ' ' + \
			  str(data_prop.num_instances) + ' ' + str(kwds['nearest_nbs']) +  \
			  ' 0 > /dev/null'
		#print relieff_cmd
		failure = os.system(relieff_cmd)
		if failure:
			print '%s: running %s failed' % (sys.argv[0], relieff_cmd)
			sys.exit(1)
		# open relief file and get attribute scores
		relieff_resultsfile = open(data_prop.infilename + '.relieff','r') 
		relieff_results	= relieff_resultsfile.readlines()
		relieff_resultsfile.close()
		result_list = []
		# collect scores and convert from string to float
		for nums in relieff_results:
			result_list.append(float(nums))
		max_score = max(result_list)
		# convert score to something like an energy, where the worst (lowest)
		# relief score has the largest energy: subtract the largest
		# score and reflect about 0
		energies = [-(i-max_score) for i in result_list]
		total_energy = sum(energies)
		ave_energy = total_energy/data_prop.num_attributes
		# normalize to make a probability field
		#energies = map(lambda x: (float(x)/total_energy), energies)
		# create an attribute -> relieff score dictionary,
		# so we can return sorted attributes
		energy_dictionary = {}
		for (k,v) in zip(data_prop.attribute_list,energies):
			energy_dictionary[k] = v
		# sort by score
		sorted_energies = data_prop.sort_value(energy_dictionary) 
		return (ave_energy, max(energies), energy_dictionary, sorted_energies)

	def irelieff_sort(data_prop, **kwds):
		irelieff_cmd = './irelieff ' + data_prop.infilename + ' ' + \
			  str(data_prop.num_instances) + ' ' + str(kwds['nearest_nbs']) +  \
			  ' 0 100 > /dev/null'
		failure = os.system(irelieff_cmd)
		if failure:
			print '%s: running %s failed' % (sys.argv[0], irelieff_cmd)
			sys.exit(1)
		# open irelief file and get attribute scores
		irelieff_resultsfile = open(data_prop.infilename + '.irelieff','r') 
		irelieff_results	= irelieff_resultsfile.readlines()
		irelieff_resultsfile.close()
		iresult_list = []
		# collect scores and convert from string to float
		for nums in irelieff_results:
			iresult_list.append(float(nums))
		iresults_dictionary = {}
		for (k,v) in zip(data_prop.attribute_list,iresult_list):
			iresults_dictionary[k] = v
		# sort by score
		sorted_attributes = data_prop.sort_value(iresults_dictionary,True) 
		return (sorted_attributes)

if __name__ == "__main__":
	try:
		infilename  = sys.argv[1]
		n_nbs = sys.argv[2]
	except:
		print "Usage:",sys.argv[0], "in_all_data n_neighbors"
		print "Example:",sys.argv[0], "mutual_info_testset.tab n_neighbors"
		sys.exit(1)
	dp = DataProperties(infilename)
	relieff_analysis(dp, n_nbs)
