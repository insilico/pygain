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
        for line in tab_infile: 
            self.data.append(line.replace('\n', ''))
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

    def sort_value(self, d, n):
        """ Returns the keys of dictionary d sorted by their values, default is low to high,
	if parameter n = 1, function sorts high to low """
        items=d.items()
        backitems=[ [v[1],v[0]] for v in items]
        backitems.sort()
        if n == 1:
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
            self.entropy_dict[key]     = self.entropy(key)
            self.mutual_info_dict[key] = self.entropy_dict[key] + self.entropy(self.status_key) - self.joint_entropy(key,self.status_key)
        #    norm = norm + self.mutual_info_dict[key]
        #for key in attrs_minus_class_keys:
        #    self.mutual_info_dict[key] = self.mutual_info_dict[key]/float(norm)
        #print 'entropies: ', self.entropy_dict
        MI_sorted_attrs = self.sort_value(self.mutual_info_dict,1)
        return (self.mutual_info_dict, MI_sorted_attrs)
        #print 'sorted mutual informations: ', sorted_MI

    def individual_chi_square(self,attribute_key1,attribute_key2):
        # calculate joint entropy given attr key names
        data_str1 = self.attribute_dictionary[attribute_key1]
        data_str2 = self.attribute_dictionary[attribute_key2]
        data_str1 = string.replace(str(data_str1),'(','')
        data_str1 = string.replace(str(data_str1),')','')
        data_str1 = string.replace(str(data_str1),' ','')        
        data_str2 = string.replace(str(data_str2),'(','')
        data_str2 = string.replace(str(data_str2),')','')
        data_str2 = string.replace(str(data_str2),' ','')
        data_str1 = data_str1.split(',')
        data_str2 = data_str2.split(',')
        # count the number of instances in each class state
        # it would speed up things to move this to a separate method
        # because you don't have to calculate everytime
        class_totals_dict = {}
        for x in data_str2: # attribute_key2 is usually the class variable
            if class_totals_dict.has_key(x):
                class_totals_dict[x] = class_totals_dict[x] + 1
            else:
                class_totals_dict[x] = 1
        class_fracs_dict = {}
        for class_state in class_totals_dict:
            class_fracs_dict[class_state] = float(class_totals_dict[class_state])/self.num_instances
        attr1_totals_dict = {}
        for x in data_str1: # attribute_key2 is usually the class variable
            if attr1_totals_dict.has_key(x):
                attr1_totals_dict[x] = attr1_totals_dict[x] + 1
            else:
                attr1_totals_dict[x] = 1
        # create expected frequencies
        expected_freq_dict = {}
        for attr in attr1_totals_dict:
            for class_state in class_fracs_dict:
                key = ",".join([attr,class_state])
                expected_freq_dict[key] =  class_fracs_dict[class_state]*attr1_totals_dict[attr]
        #print class_fracs_dict
        #print attr1_totals_dict
        #print expected_freq_dict
        # observed frequency dictionary
        observed_freq_dict = {}
        for x in zip(data_str1,data_str2):
            #print x[0], x[1]
            combination = ",".join(x)
            if observed_freq_dict.has_key(combination):
                observed_freq_dict[combination] = observed_freq_dict[combination] + 1
            else:
                observed_freq_dict[combination] = 1
        #print observed_freq_dict
        chi2 = 0
        for key in expected_freq_dict:
            key_split = key.split(',')
            if observed_freq_dict.has_key(key):
                chi2 = chi2 + math.pow((observed_freq_dict[key] - expected_freq_dict[key]),2)/expected_freq_dict[key]
            else:
                chi2 = chi2 + math.pow((expected_freq_dict[key]),2)/expected_freq_dict[key]
        return chi2 

    def total_chi_squares(self):
        # total chi-square for each attribute
        attrs_minus_class_keys = []
        for key in self.attribute_dictionary:
            if key != self.status_key:
                attrs_minus_class_keys.append(key)    
        self.chi2_dict = {}
        for key in attrs_minus_class_keys:
            self.chi2_dict[key] = self.individual_chi_square(key,self.status_key)
        chi2_sorted_attrs = self.sort_value(self.chi2_dict, 1)
        return (self.chi2_dict, chi2_sorted_attrs)
        #print 'sorted mutual informations: ', sorted_MI
            
    def relieff_analysis(self, **kwds):
       """
       kwd defaults: nearest_nbs = 10 (./relieff nearest neighbors)
       usage: (score, sorted_attributes)=data.relieff_analysis(nearest_nbs= 5)
       
        ./relieff creates output file called infilename.relieff, which
        we open after running in order to read in the relief scores
        and zip with attributes and sort. Returns total relief score
        and relieff-sorted attributes from self.data.
        """      
       if kwds.has_key('nearest_nbs') == False:
           kwds['nearest_nbs'] = 10
       relieff_cmd = './relieff ' + self.infilename + ' ' + \
              str(self.num_instances) + ' ' + str(kwds['nearest_nbs']) +  \
              ' 0 > /dev/null'
       #print relieff_cmd
       failure = os.system(relieff_cmd)
       if failure:
           print '%s: running %s failed' % (sys.argv[0], relieff_cmd)
           sys.exit(1)
       # open relief file and get attribute scores
       relieff_resultsfile = open(self.infilename + '.relieff','r') 
       relieff_results    = relieff_resultsfile.readlines()
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
       ave_energy = total_energy/self.num_attributes
       # normalize to make a probability field
       #energies = map(lambda x: (float(x)/total_energy), energies)
       # create an attribute -> relieff score dictionary,
       # so we can return sorted attributes
       energy_dictionary = {}
       for (k,v) in zip(self.attribute_list,energies):
           energy_dictionary[k] = v
       # sort by score
       sorted_energies = sort_value(energy_dictionary) 
       return (ave_energy, max(energies), energy_dictionary, sorted_energies)

    def irelieff_sort(self, **kwds):
       irelieff_cmd = './irelieff ' + self.infilename + ' ' + \
              str(self.num_instances) + ' ' + str(kwds['nearest_nbs']) +  \
              ' 0 100 > /dev/null'
       failure = os.system(irelieff_cmd)
       if failure:
           print '%s: running %s failed' % (sys.argv[0], irelieff_cmd)
           sys.exit(1)
       # open irelief file and get attribute scores
       irelieff_resultsfile = open(self.infilename + '.irelieff','r') 
       irelieff_results    = irelieff_resultsfile.readlines()
       irelieff_resultsfile.close()
       iresult_list = []
       # collect scores and convert from string to float
       for nums in irelieff_results:
           iresult_list.append(float(nums))
       iresults_dictionary = {}
       for (k,v) in zip(self.attribute_list,iresult_list):
           iresults_dictionary[k] = v
       # sort by score
       sorted_attributes = sort_value(iresults_dictionary,1) 
       return (sorted_attributes)
