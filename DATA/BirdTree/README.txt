Phylogeny related data 

"BLIOCPhyloMasterTax.csv" is the file containing all bird species by far, I downloaded this taxonomy file on 18th Apr, 2023 from 
https://birdtree.org/subsets/
Species name might differ in their BirdTree database with my given species, I resolved those differing names in 
"species_0_250km_nbin_4_filledin.csv file" (see last three columns). Related code is given here: "R/get_birdspecies_phylotree.R".

"unique_speciesnameBirdTree_0_250km_nbin4_tailsig75.txt" is the text file containing the species name compatible with BirdTree database, 
I used these species names to download 1000 trees from https://birdtree.org/subsets/, these trees are given in the folder "sig75_0_250km_..." as output.nex file.
We also provide other files (text files and nexus files) set of species for 0-100 km, and 100-250 km distance categories.