Results folder

We have provided 4 csv files here.

1) df_abund_climate_spatsyn_0_250km.csv = Output from R/summarize_res.R, contains info for tail-dependent spatial synchrony 
					in abundance, and in climates within 0-250 Km between-sites distance, 
					and for 262 species (1 species AOU=7470 was omitted as no tail-dep for precipitation)

2) df_abund_climate_spatsyn_0_250km_with_optimal_biovar.csv = Output from R/get_bio_opt_for_species.R, contains same info
								as of the previous file with added info about the bio1, 
								bio12 variables (optimal T, P niche), but for 253 species
								that were considered for phylogenetic relationship 
								(after matching names from BirdTree).

3) species263_127LT_136UT.csv = AOU code and the scientific names for all 263 bird species considered for tail-dependent spatial synchrony
				within 0-250 Km between-sites distance. First 127 species showed lower tail-dependence (synchronously 
				rare across sites), and the rest of 136 species showed upper tail-dependence (synchronously 
				common across sites).

4) species_dietcat_edited.csv = Diet category (5 types) for all 373 species that were sampled at least at two sites for a minimum of 20 years 
				from BBS dataset (Pardieck et al. 2020). The output csv file from R/diet_cat.R was filled in manually as
				diet info for some species were not available in EltonTraits databse (Wilman et al. 2014). And we filled 
				in those missing info from Birds of the World (Billerman et al. 2022). 