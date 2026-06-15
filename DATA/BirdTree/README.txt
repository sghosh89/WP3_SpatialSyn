Phylogeny related data 

"BLIOCPhyloMasterTax.csv" is the file containing all bird species by far, I downloaded this taxonomy file on 18th Apr, 2023 from 
https://birdtree.org/subsets/
Species name might differ in their BirdTree database with my given species, I resolved those differing names in 
"species_0_250km_nbin_4_filledin_min32yr.csv file" (see last three columns). Related code is given here: "R/get_birdspecies_phylotree_with_Nyrs_threshold.R".

The text files are named as follows:
"unique_speciesnameBirdTree_0_250km_nbin4_tailsig95.txt" is the text file containing the species that showed significant tail-dependent synchrony in abundance 
(based on 95% CI) within 0-250Km between-sites distance category and at least 40 years sampled as used in the main analysis.
I used these species names to download 1000 trees from https://birdtree.org/subsets/, these trees are given in the corresponding folder as output.nex file.

We also provide other files (text files and nexus files) set of species for:
	1) "unique_speciesnameBirdTree_0_250km_nbin4_tailsig95.txt" is the text file containing the species that showed significant tail-dependent synchrony in abundance 
	   (based on 95% CI) within 0-250Km between-sites distance category and at least 40 years sampled as used in the main analysis.
	2) For 0-250km distance category we also provide those files based on 95% CI significant tail-dependent synchrony considering routes that are sampled for a     
           minimum of 32 and 36 years, respectively.
	3) For default (a minimum of 40 years) sampling routes (as used in the main analysis) we also provide files for other two distance categories: 0-100km, and 100-250km.
	4) filenames ended with "...sitepairthr_45.txt" mean species were further filtered based on if they were sampled at least at 10 sites (i.e., 45 unique site pairs).  