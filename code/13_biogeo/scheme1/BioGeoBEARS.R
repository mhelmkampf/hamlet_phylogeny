#---- Ancestral range estimations (BioGeoBEARS) ------#
#--------------- Hamlets (Hypoplectrus)----------------#

#----- Loading packages -----#
#install.packages("rexpokit")
#install.packages("cladoRcpp")
#library(devtools)
#devtools::install_github(repo="nmatzke/BioGeoBEARS", dependencies=FALSE)

library(optimx)
library(ade4)
library(FD)     
library(snow)     
library(parallel)
library(BioGeoBEARS)
library(cladoRcpp)
library(rexpokit)
library(roxygen2)
library(phytools)

#----- Setting WD -----#
#wd = np("") # Set your working directory
#setwd(wd)
setwd("~/scheme1")

#----- Setting up extension data directory -----#
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))

#----- Loading tree -----#
tree.ham.data = "scheme1.tre"
phylo.ham = read.tree(tree.ham.data)
phylo.ham = ladderize(phylo.ham)


#----- Loading geographical data -----#
geog.ham.data = "Geo_data_scheme1.txt"
moref(geog.ham.data)

tipranges.ham = getranges_from_LagrangePHYLIP(lgdata_fn=geog.ham.data)

####################################
#----- Biogeographical models -----#
####################################
#--------- DEC AND DEC+J ----------#
####################################

#--------- DEC ----------#
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$force_sparse=FALSE   
BioGeoBEARS_run_object$speedup=TRUE          
BioGeoBEARS_run_object$calc_ancprobs=TRUE    
BioGeoBEARS_run_object$use_optimx = TRUE

# Input the maximum range size
max_range_size = 3
BioGeoBEARS_run_object$max_range_size = max_range_size

# Multicore processing (only if needed)
BioGeoBEARS_run_object$num_cores_to_use=20

# Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+) - Not the present case
BioGeoBEARS_run_object$force_sparse=FALSE

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geog.ham.data

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = tree.ham.data

# Loading the dispersal multiplier matrix and running some checks
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE

# Set up DEC model
# (nothing to do; defaults)

# Final check
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

# RUN
runslow = TRUE
resfn = "Hamlets_DEC_scheme1.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDEC = res
} else {
  # Loads to "res"
  load(resfn)
  resDEC = res
}

#--------- DEC+J ----------#

BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$force_sparse=FALSE   
BioGeoBEARS_run_object$speedup=TRUE          
BioGeoBEARS_run_object$calc_ancprobs=TRUE    
BioGeoBEARS_run_object$use_optimx = TRUE


# Input the maximum range size
max_range_size = 3
BioGeoBEARS_run_object$max_range_size = max_range_size

# Multicore processing (only if needed)
BioGeoBEARS_run_object$num_cores_to_use=20

# Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+) - Not the present case
BioGeoBEARS_run_object$force_sparse=FALSE

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geog.ham.data

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = tree.ham.data

# Loading the dispersal multiplier matrix and running some checks
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE

# Set up DEC+J model
dstart = resDEC$outputs@params_table["d","est"]
estart = resDEC$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Final check
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "Hamlets_DEC+J_scheme1.Rdata"
runslow = TRUE
if (runslow)
{
  
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resDECj = res
} else {
  # Loads to "res"
  load(resfn)
  resDECj = res
}


####################################
#----- DIVALIKE AND DIVALIKE+J ----#
####################################

#--------- DIVALIKE ----------#
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$force_sparse=FALSE   
BioGeoBEARS_run_object$speedup=TRUE          
BioGeoBEARS_run_object$calc_ancprobs=TRUE    
BioGeoBEARS_run_object$use_optimx = TRUE

# Input the maximum range size
max_range_size = 3
BioGeoBEARS_run_object$max_range_size = max_range_size

# Multicore processing (only if needed)
BioGeoBEARS_run_object$num_cores_to_use=20

# Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+) - Not the present case
BioGeoBEARS_run_object$force_sparse=FALSE

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geog.ham.data

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = tree.ham.data

# Loading the dispersal multiplier matrix and running some checks
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE

# Set up DIVALIKE model
# Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# Final check
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

# RUN
runslow = TRUE
resfn = "Hamlets_DIVALIKE_scheme1.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDIVALIKE = res
} else {
  # Loads to "res"
  load(resfn)
  resDIVALIKE = res
}


#--------- DIVALIKE+J ----------#

BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$force_sparse=FALSE   
BioGeoBEARS_run_object$speedup=TRUE          
BioGeoBEARS_run_object$calc_ancprobs=TRUE    
BioGeoBEARS_run_object$use_optimx = TRUE

# Input the maximum range size
max_range_size = 3
BioGeoBEARS_run_object$max_range_size = max_range_size

# Multicore processing (only if needed)
BioGeoBEARS_run_object$num_cores_to_use=20

# Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+) - Not the present case
BioGeoBEARS_run_object$force_sparse=FALSE

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geog.ham.data

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = tree.ham.data

# Loading the dispersal multiplier matrix and running some checks
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE

# Set up DIVALIKE+J model
dstart = resDIVALIKE$outputs@params_table["d","est"]
estart = resDIVALIKE$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# Add jump dispersal/founder-event speciation
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in DEC+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999

# Final check
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

# RUN
resfn = "Hamlets_DIVALIKE+J_scheme1.Rdata"
runslow = TRUE
if (runslow)
{
  
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resDIVALIKEj = res
} else {
  # Loads to "res"
  load(resfn)
  resDIVALIKEj = res
}


####################################
#- BAYAREALIKE AND BAYAREALIKE+J--#
####################################

#--------- BAYAREALIKE ----------#
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$force_sparse=FALSE   
BioGeoBEARS_run_object$speedup=TRUE          
BioGeoBEARS_run_object$calc_ancprobs=TRUE    
BioGeoBEARS_run_object$use_optimx = TRUE

# Input the maximum range size
max_range_size = 3
BioGeoBEARS_run_object$max_range_size = max_range_size

# Multicore processing (only if needed)
BioGeoBEARS_run_object$num_cores_to_use=20

# Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+) - Not the present case
BioGeoBEARS_run_object$force_sparse=FALSE

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geog.ham.data

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = tree.ham.data

# Loading the dispersal multiplier matrix and running some checks
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE

# Set up BAYAREALIKE model
# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# Final check
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

# RUN
runslow = TRUE
resfn = "Hamlets_BAYAREALIKE_scheme1.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resBAYAREALIKE = res
} else {
  # Loads to "res"
  load(resfn)
  resBAYAREALIKE = res
}

#--------- BAYAREALIKE+J----------#

BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$force_sparse=FALSE   
BioGeoBEARS_run_object$speedup=TRUE          
BioGeoBEARS_run_object$calc_ancprobs=TRUE    
BioGeoBEARS_run_object$use_optimx = TRUE

# Input the maximum range size
max_range_size = 3
BioGeoBEARS_run_object$max_range_size = max_range_size

# Multicore processing (only if needed)
BioGeoBEARS_run_object$num_cores_to_use=20

# Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+) - Not the present case
BioGeoBEARS_run_object$force_sparse=FALSE

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geog.ham.data

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = tree.ham.data

# Loading the dispersal multiplier matrix and running some checks
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE

# Set up DIVALIKE+W model
dstart = resBAYAREALIKE$outputs@params_table["d","est"]
estart = resBAYAREALIKE$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# *DO* allow jump dispersal/founder-event speciation (set the starting value close to 0)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under DIVALIKE+W, the max of "j" should be 1, not 3 (as is default in DEC+J) or 2 (as in DIVALIKE+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# setting min and max for d and e to avoid crashing:
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

# Final check
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

# RUN
resfn = "Hamlets_BAYAREALIKE+J_scheme1.Rdata"
runslow = TRUE
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resBAYAREALIKEj = res
} else {
  # Loads to "res"
  load(resfn)
  resBAYAREALIKEj = res
}


####################################
#------- SUMMARY STATISTICS -------#
####################################

# Set up empty table to hold the statistical results

restable = NULL

####################################
#----- EXTRACTING PARAMETERS ------#
####################################

#----- DEC ------#

# DEC
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDEC, returnwhat="table", addl_params=c("j","w"), paramsstr_digits=4)
# +J
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDECj, returnwhat="table", addl_params=c("j","w"), paramsstr_digits=4)

restable = rbind(restable, res1, res2)

#----- DIVALIKE ------#

# DIVALIKE
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKE, returnwhat="table", addl_params=c("j","w"), paramsstr_digits=4)
# +J
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKEj, returnwhat="table", addl_params=c("j","w"), paramsstr_digits=4)

restable = rbind(restable, res1, res2)

#----- BAYAREALIKE ------#

# BAYAREALIKE
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKE, returnwhat="table", addl_params=c("j","w"), paramsstr_digits=4)
# +J
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKEj, returnwhat="table", addl_params=c("j","w"), paramsstr_digits=4)

restable = rbind(restable, res1, res2)


####################################
#-------- Assembling table --------#
####################################
row.names(restable) = c("DEC", "DEC+J","DIVALIKE", "DIVALIKE+J",
                        "BAYAREALIKE", "BAYAREALIKE+J")
restable = put_jcol_after_ecol(restable)

restable2 = restable

# With AICcs -- factors in sample size
samplesize = length(phylo.ham$tip.label)
AICtable = calc_AICc_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams, samplesize=samplesize)
restable2 = cbind(restable2, AICtable)
restable_AICc_rellike = AkaikeWeights_on_summary_table(restable=restable2, colname_to_use="AICc")
restable_AICc_rellike = put_jcol_after_ecol(restable_AICc_rellike)
restable_AICc_rellike


# Save with nice conditional formatting
write.table(conditional_format_table(restable_AICc_rellike), file="Hamlets_restable_AICc_formatted_scheme1.txt", quote=FALSE, sep="\t")

####################################
#---- Plot results best model -----#
####################################

# Loading packages
library(ape)
library(ggtree)
library(ggplot2)
library(phytools)
library(BioGeoBEARS)
library(cladoRcpp)
library(treeio)
library(ggimage)


#----- Setting WD -----#
# Set your working directory
setwd("~/scheme1")

#----- Setting up extension data directory -----#
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))

#----- Loading tree -----#
tree.ham.data = "scheme1.tre"
phylo.ham = read.tree(tree.ham.data)
phylo.ham = ladderize(phylo.ham)

#----- Loading geographical data -----#
geog.ham.data = "Geo_data_scheme1.txt"
moref(geog.ham.data)

tipranges.ham = getranges_from_LagrangePHYLIP(lgdata_fn=geog.ham.data)

areas = getareas_from_tipranges_object(tipranges.ham)

numareas = length(areas)

tips = 1:length(phylo.ham$tip.label)
nodes = (length(phylo.ham$tip.label)+1):(length(phylo.ham$tip.label)+phylo.ham$Nnode)

load("Hamlets_BAYAREALIKE+J_scheme1.Rdata") #the best supported model

if (!is.na(res$inputs$max_range_size))
{
  max_range_size = res$inputs$max_range_size
} else {
  max_range_size = length(areas)
}


if (is.null(res$inputs$states_list))
{
  numstates = numstates_from_numareas(numareas=length(areas), maxareas=max_range_size, include_null_range=res$inputs$include_null_range)
  numstates
  states_list_areaLetters = areas_list_to_states_list_new(areas, maxareas=max_range_size, include_null_range=res$inputs$include_null_range)
  states_list_0based_index = rcpp_areas_list_to_states_list(areas, maxareas=max_range_size, include_null_range=res$inputs$include_null_range)
} else {
  states_list_0based_index = res$inputs$states_list
}


relprobs_matrix = res$ML_marginal_prob_each_state_at_branch_top_AT_node

if (length(nodes) > 1)
{
  relprobs_matrix_for_internal_states = relprobs_matrix[nodes,]	# subset to just internal nodes
} else {
  relprobs_matrix_for_internal_states = relprobs_matrix[nodes,]	# subset to just internal nodes
  # Convert back to a matrix
  relprobs_matrix_for_internal_states = matrix(data=relprobs_matrix_for_internal_states, nrow=1, ncol=ncol(relprobs_matrix))
}


if (is.null(states_list_0based_index))
{
  statenames = areas_list_to_states_list_new(areas, maxareas=max_range_size, include_null_range=res$inputs$include_null_range, split_ABC=FALSE)
  ranges_list = as.list(statenames)
  statenames
} else {
  ranges_list = states_list_0based_to_ranges_txt_list(state_indices_0based=states_list_0based_index, areanames=areas)
  ranges_list
  statenames = unlist(ranges_list)
  statenames
}


MLprobs = get_ML_probs(relprobs_matrix)

MLstates = get_ML_states_from_relprobs(relprobs_matrix, statenames, returnwhat="states", if_ties="takefirst")


#############
# CHRONOGRAM Plotting
#############

chr = ggtree(phylo.ham,size = 0.4, lineend = "round", ladderize = FALSE)
chr = revts(chr)

biogeo = read.csv("biogeonames.txt", sep = ",", header = T, stringsAsFactors = F)
for (i in 1:ncol(biogeo)){
  biogeo[which(biogeo[,i] == ""),i] <- NA
}
rownames(biogeo) <- biogeo$X
biogeo <- biogeo[,-1]

relprobs_matrix1<-as.data.frame(relprobs_matrix)
relprobs_matrix1[,9]<-1:337
colnames(relprobs_matrix1)<-c("non","C","A","G","X","X1","X2","X3","node")
#write.csv(relprobs_matrix1,file="relprobs_matrix1.csv")
only_internal_nodes<-relprobs_matrix1[-c(1:169),]
#write.csv(only_internal_nodes,file="only_internal_nodes.csv")

pies<-nodepie(only_internal_nodes, cols=1:8)
colors_list_for_states = c("A"="#698B22","C"="#4169E1","G"="#EE6A50","X"="#008080","X1"="#800080","X2"="#808000","X3"="#FFFFFF")
pies <- lapply(pies, function(g) g+scale_fill_manual(values = colors_list_for_states))

plt_pies<- chr + theme_tree2() + scale_x_continuous(breaks=seq(-3,0,0.5), labels=seq(3,0,-0.5),limits=c(-3,2)) +
  scale_y_continuous(limits=c(0,180)) + geom_inset(pies, width = .1, height = .1) +
   geom_tiplab(size=1.5, offset=0.5)

fig_pies <- gheatmap(plt_pies, biogeo, offset=0.1, width=0.1, font.size=1.5, 
                hjust=0,colnames_position = "top") + scale_fill_manual(values = colors_list_for_states)


pdf(file="Hamlets_BioGeo_BAYAREALIKE+J_scheme1_pies.pdf", height = 9, width = 8, useDingbats = FALSE)
fig_pies
dev.off()

####### color by species at the tip ####### 
#add tips and species
only_tips<-relprobs_matrix1[c(1:169),]
only_tips[,10]<-phylo.ham$tip.label
#write.csv(only_tips,file="only_tips.csv")
relprobs_matrix1_table<-read.csv(file="relprobs_matrix1.csv")
color_list_for_species = c("aberrans"="#DFDF8D","affinis"="#FFB1D8","atlahua"="#5C7A1E","castroaguirrei"="#E3280D","chlorurus"="#8B4513","ecosur"="#FFBF70","espinozaperezei"="#E93E95","floridae"="#A155E7","gemma"="#1E90FF","gummigutta"="#F99729","guttavarius"="#FFEA00","indigo"="#22198E","liberte"="#C3A8DC","maya"="#7196AE","nigricans"="#333333","providencianus"="#8BCF4C","puella"="#E48175","randallorum"="#8AC3BA","tan"="#D2B48C","unicolor"="#B3B3B3")

fig2<- chr %<+% relprobs_matrix1_table + geom_tippoint(aes(color=Species),size=1)+
  scale_color_manual(values=color_list_for_species)+ theme_tree2() + scale_x_continuous(breaks=seq(-3,0,0.5), labels=seq(3,0,-0.5),limits=c(-3,2)) +
  scale_y_continuous(limits=c(0,180))
pdf(file="Hamlets_scheme1_species.pdf", height = 9, width = 8, useDingbats = FALSE)
fig2
dev.off()

#plt <- chr + theme_tree2() + scale_x_continuous(breaks=seq(-3,0,0.5), labels=seq(3,0,-0.5),limits=c(-3,2)) +
#  scale_y_continuous(limits=c(0,180)) +
#  geom_label2(aes(subset=!isTip, label=MLstates, fill = MLstates), size=1.1, colour = "Black", alpha = 0.85,
#              fontface = "bold", label.padding = unit(0.13, "lines"), label.size = 0) + geom_tiplab(size=1.5, offset=0.5) +scale_fill_manual(values = colors_list_for_states)
#
#fig <- gheatmap(plt, biogeo, offset=0.1, width=0.1, font.size=1.5,
#                hjust=0,colnames_position = "top") + scale_fill_manual(values = colors_list_for_states)
#
#
#pdf(file="Hamlets_BioGeo_BAYAREALIKE+J_scheme1_mostprobablestate.pdf", height = 9, width = 8, useDingbats = FALSE)
#fig
#dev.off()