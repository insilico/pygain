pygain
========

#### Genetic Association Interaction Network (GAIN) tool ####

### Description ###
GAIN[1] is based on interaction information between three attributes; in this case,
between two single nucleotide polymporphisms (SNPs) and a class or phenotype 
attribute. Interaction information is the gain in phenotype information obtained
by considering SNP A and SNP B jointly beyond the phenotype information that 
would be gained by considering SNPs A and B independently.

This software was created as a bioinformatics tool for usage by our research 
group, [In Silico](http://insilico.utulsa.edu), as well as other researchers 
and interested parties.  

### Dependencies and Usage ###
pygain is developed and tested on 64-bit Linux (Ubuntu), but should work on any 
platform supported by Python (Python version 2.6.5 tested).

To run pysnprank from command-line:

    ./gain.py -i plink-data.raw -o gain-matrix.txt

Additional parameters:
	Usage: gain.py [OPTIONS]

	Construct GAIN matrix from PLINK .raw or tab-delimited file

	Options:
		--version       display program version
		--help      -h	display this help and exit
		--input	-i		Input file (default: stdin)
		--output	-o	Output file (default: stdout)
		--export-sif -e	Export Cytoscape .sif file

### Contributors ###
See AUTHORS file.

### References ###
[1]B.A. McKinney, J.Guo, J.E. Crowe, Jr., and D. Tian. Capturing the spectrum of 
interaction effects in genetic association studies by simulated evaporative 
cooling network analysis. PLoS Genetics 2009, 5(3): e1000432. 
doi:10.1371/journal.pgen.1000432. [open access](http://www.plosgenetics.org/article/info:doi/10.1371/journal.pgen.1000432)
