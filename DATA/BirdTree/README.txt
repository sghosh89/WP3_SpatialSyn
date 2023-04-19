Phylogeny related data 

"BLIOCPhyloMasterTax.csv" is the file containing all bird species by far, I downloaded this taxonomy file on 18th Apr, 2023 from 
https://birdtree.org/subsets/
Species name might differ in their BirdTree database with my given species, I resolved those differing names in 
"species_0_250km_filledin.csv file" (see last three columns). Related code is given here: "R/get_birdspecies_phylotree.R".
"unique_speciesnameBirdTree_0_250km.txt" is the text file containing the species name compatible with BirdTree database, I used these species names 
to download 1000 trees from https://birdtree.org/subsets/, these trees are given in "tree-pruner-f1d9b817-3739-4e7d-bbaf-1227c85c4a2c" folder.
