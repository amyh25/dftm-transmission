
GENERAL INFORMATION
This README.txt file was updated on March 14, 2025

A. Paper associated with this archive
Citation: XXXXXX
Brief abstract:
The theory of host-pathogen interactions has successfully shown that persistent pathogen virulence may be explained through tradeoffs between different pathogen fitness components, but classical theory cannot explain pathogen coexistence.  More recent theory invokes both tradeoffs and environmental heterogeneity, but resembles classical theory in focusing on a limited range of possible tradeoffs, and therefore has seen few applications.  To better understand the usefulness of tradeoff theory for explaining pathogen coexistence in nature, we measured components of pathogen fitness in two distantly related morphotypes of a baculovirus that infects larvae of the Douglas-fir tussock moth (Orgyia pseudotsugata). We show that the two morphotypes vary in multiple components of fitness, including the probability of infection given exposure to the pathogen, the incubation time of the pathogen, variability in the incubation time of the pathogen, and the detectability of the pathogen.  Moreover, because the baculovirus is transmitted when host larvae accidentally consume infectious virus particles while feeding on foliage of the insect's host trees, the strength and direction of the differences in fitness components of the two morphotypes depends on the host-tree species on which host larvae consume the virus.  Through simulations of a model parameterized using our experimental data, we demonstrate how several varying fitness components can work in concert to promote strain coexistence, particularly highlighting the role of variability in incubation time.  Our results suggest that the two morphotypes may coexist because of variation in forest tree-species composition, providing important empirical evidence that tradeoffs and environmental heterogeneity can together modulate pathogen competition. 

B. Originators
- Ari S. Freedman, University of Chicago Department of Ecology and Evolution (now at Princeton University Department of Ecology and Evolutionary Biology)
- Amy Huang, University of Chicago Department of Ecology and Evolution (now at Massachusetts Institute of Technology Computational Systems Biology)
- Katherine Dixon, University of Chicago Department of Ecology and Evolution
- Carlos Polivka, Wenatchee Forestry Sciences Laboratory
- Greg Dwyer, University of Chicago Department of Ecology and Evolution

C. Corresponding author contact information
Greg Dwyer, gdwyer@uchicago.edu

D. Dates of data collection
Summer 2019

E. Geographic Location(s) of data collection
Wenatchee, Washington and surrounding area

F. Funding Sources
- NSF EEID grant DEB-2109774
- NSF OPUS grant 285 2043796
- University of Chicago Ecology and Evolution Fellowship
- University of Chicago Biological Sciences Collegiate Division Fellowship


ACCESS INFORMATION
1. Licenses/restrictions placed on the data or code
NA
2. Data derived from other sources
NA
3. Recommended citation for this data/code archive
XXXXXX

DATA & CODE FILE OVERVIEW
This data repository consist of 2 data files, 4 R code scripts, 11 Stan model files, and this README document, with the following data and code filenames and variables

Data files and variables
    1. data/SOK_data.csv, data from infection/speed-of-kill experiments
	capsid - virus morphotype used for infecting larvae (or none for control plates)
	strain - isolate of virus used for infecting larvae (or CTRL for control plates)
	density - low, medium, or high, corresponding to 1, 2, or 3 drops of ~1800 occlusion bodies/microliter virus solution
	tree_sp - tree species that larvae were fed before infection, either grand fir or Douglas-fir
	numeric_dose - 1, 2, or 3 corresponding to the number of drops or ~1800 occlusion bodies/microliter virus solution
	total_n - number of larvae infected in this treatment
	total_virus - number of larvae in this treatment killed by the virus
	dose_response - proportion of larvae in this treatment killed by the virus, equal to total_virus / total_n
	day_6, day_7, ..., day_28 - number of larvae in this treatment killed by the virus 6, 7, ..., 28 days after infection
	ob_per_drop - exact dose in occlusion bodies/microliter of the virus solution used for this treatment
	ob_count - number of occlusion bodies used for infection in this treatment, equal to ob_per_drop * numeric_dose
    2. data/manual_avoid_effect.csv, data from avoidance experiments
	isolate - isolate of virus used for infecting larvae (or CTRL if control plate)
	morphotype - morphotype of virus used for infecting larvae (SNPV for control plates)
	tree_sp - tree species on which choice test was run, either grand fir or Douglas-fir
	rep - plate number within that treatment
	control - percent of leaves eaten from the control side of the plate
	treatment - percent of leaves eaten from the treatment side of the plate

Code scripts and workflow
    1. mortality/mortality_models.R, runs the Stan models mortality/morphotype_and_tree.stan, mortality/morphotype_only.stan, mortality/tree_only.stan, mortality/neither_morphotype_nor_tree.stan, and mortality/no_hierarchy.stan, analyzes their goodness of fits and produces related figures
    2. SOK/SOK_models.R, runs the Stan models SOK/morphotype_and_tree.stan, SOK/morphotype_only.stan, SOK/tree_only.stan, and SOK/neither_morphotype_nor_tree.stan, analyzes their goodness of fits and produces related figures
    3. avoidance/avoidance_models.R, runs the Stan models avoidance/model.stan and avoidance/model_intercept.stan, analyzes their goodness of tis and produces related figures
    4. SOK/SOK_variance.R (to be run after SOK/SOK_models.R to use its Stan model fits), runs an SEIR model parameterized by the speed-of-kill distribution models' fits
    5. results_summary.R (to be run after mortality/mortality_models.R, SOK/SOK_models.R, and avoidance/avoidance_models.R), plots a results summary figure

SOFTWARE VERSIONS
R version 4.0.2, Stan version 2.32.2
R packages: rstan, tidyverse, loo, ggpubr, deSolve, pracma, binom, boot, bayesplot
