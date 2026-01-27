############################# SIMULATION ################################

# Libraries
library(dplyr) # To format tables
library(ggplot2)
library(lme4)
library(pbapply)
library(parallel)
library(tidyverse)
library(nanoparquet)

# Functions
source("ssd_null.R")
#source("ssd_inf.R")
source("helpers_simulation.R")

# common factors --------------
var_u0 <- c(0, 0.25, 0.5, 1)
eff_size_beta <- c(0.5, 1, 1.5)
eff_size_prob <- 1 / (1 + exp(-eff_size_beta))
BF_thres <- c(1, 3, 5)
b_fract <- 3
eta <- 0.8
ndatasets <- 5000
Max <- 1000
batch_size <- 500

########## Hypothesis Set 1: Equality vs. informative ##########################


## Find the number of clusters ====
path <- "~/"
results_folder <- "results_set1"
dir.create(results_folder)

###====== Design matrix ======
n1 <- c(5, 10, 40)
design_matrixN2 <- expand.grid(var_u0, eff_size_prob, BF_thres, eta, n1, fixed <- "n1")
nrow_designN2 <- nrow(design_matrixN2)
colnames(design_matrixN2) <- c("var_u0", "p_int", "BF_threshold", "eta", "n1", "fixed")
design_matrixN2 <- mutate(design_matrixN2, seed = as.integer(sample(2^32 /
                                                                        2, n())))
write_parquet(design_matrixN2, "design_matrix_findN2_set1")

run_null_wrapper <- function(Row) {
    run_null(
        Row = Row,
        design_matrix = design_matrixN2,
        ndatasets = ndatasets,
        Max = Max,
        batch_size = batch_size,
        results_folder = results_folder,
        b = b_fract
    )
}

###======== Running the simulation =======
clusters <- makeForkCluster(detectCores() * 0.75)

#clusters <- makeCluster(detectCores() * 0.75) #For Windows
clusterEvalQ(clusters, {
    library(tidyverse)
    library(lme4)
    library(dplyr)
    library(bain)
})
# clusterExport(
#   clusters,
#   c(
#     "design_matrixN2",
#     "ndatasets",
#     "b_fract",
#     "Max",
#     "batch_size",
#     "results_folder",
#     "binary_search_eq",
#     "collect_results",
#     "collect_times",
#     "extract_res_bain",
#     "filter_underpowered",
#     "binary_search_eq",
#     "fit_glmer",
#     "get_variance",
#     "marker_func",
#     "print_results",
#     "reached_condition",
#     "round2",
#     "run_null",
#     "SSD_crt_null_binary",
#     "varcov",
#     "run_null_wrapper"
#   ), envir = environment()
# )

output <- parallel::parLapply(
    cl = clusters,
    X = 4:nrow_designN2,
    fun = run_null_wrapper
)

stopCluster(clusters)

###======= Collect results ========
#pair = 1 refering to the pair of hypotheses H0 vs H1
res_findN2_set1 <- collect_results(
    design_matrix = design_matrixN2,
    results_folder = results_folder,
    finding = "N2",
    pair = 1
)
res_time_findN2_set1 <- collect_times(
    design_matrix = design_matrixN2,
    pair = 1,
    finding = "N2",
    results_folder = results_folder
)

##========= Find cluster size ==========

results_folder

###======= Design matrix ============
n2 <- c(30, 60, 90)
design_matrixN1 <- expand.grid(var_u0, eff_size_prob, BF_thres, eta, n2, fixed <- "n2")
colnames(design_matrixN1) <- c("var_u0", "eff_size", "BF_threshold", "eta", "n2", "fixed")
nrow_designN1 <- nrow(design_matrixN1)
design_matrixN1 <- mutate(design_matrixN1, seed = as.integer(sample(2^32 /
                                                                        2, n())))
write_parquet(design_matrixN1, "design_matrix_findN1_set1")

###========= Running the simulation =========
run_null_wrapper <- function(Row) {
    run_null(
        Row = Row,
        design_matrix = design_matrixN1,
        ndatasets = ndatasets,
        Max = Max,
        batch_size = batch_size,
        results_folder = results_folder,
        b = b_fract
    )
}
clusters <- makeForkCluster(detectCores() * 0.75) #For Linux

output <- parallel::parLapply(
    cl = clusters,
    X = 1:nrow_designN1,
    fun = run_null_wrapper
)
stopCluster(clusters)

### Collect results -----
#pair = 1 refering to the pair of hypotheses H0 vs H1
res_findN1_set1 <- collect_results(
    design_matrix = design_matrixN1,
    results_folder = results_folder,
    finding = "N1",
    pair = 1
)
res_time_findN1_set1 <- collect_times(
    design_matrix = design_matrixN1,
    pair = 1,
    finding = "N1",
    results_folder = results_folder
)

##================ Plots =================




############## Hypothesis Set 2: Informative vs. informative ###################
var_u0 <- c(0, 0.25, 0.5, 1)
eff_size_beta <- c(0.5, 1, 1.5)
eff_size_prob <- 1 / (1 + exp(-eff_size_beta))
BF_thres <- c(1, 3, 5)
eta <- 0.8
ndatasets <- 5000
Max <- 1000
batch_size <- 500


##========== Find the number of clusters =========
path <- "~/"
results_folder <- "results_set2"
dir.create(results_folder)


###========= Design matrix ===========
n1 <- c(5, 10, 40)
design_matrixN2 <- expand.grid(var_u0, eff_size_prob, BF_thres, eta, n1, fixed <- "n1")
nrow_designN2 <- nrow(design_matrixN2)
colnames(design_matrixN2) <- c("var_u0", "p_int", "BF_threshold", "eta", "n1", "fixed")
design_matrixN2 <- mutate(design_matrixN2, seed = as.integer(sample(2^32 /
                                                                        2, n())))
write_parquet(design_matrixN2, "design_matrix_findN2_set2")

###======= Running the simulation =========
run_inf_wrapper <- function(Row) {
    run_inf(
        Row = Row,
        design_matrix = design_matrixN2,
        ndatasets = ndatasets,
        Max = Max,
        batch_size = batch_size,
        results_folder = results_folder
    )
}

clusters <- makeForkCluster(detectCores() * 0.75)

#clusters <- makeCluster(detectCores() * 0.75) #For Windows
clusterEvalQ(clusters, {
    library(tidyverse)
    library(lme4)
    library(dplyr)
    library(bain)
})

output <- parallel::parLapply(
    cl = clusters,
    X = 1:nrow_designN2,
    fun = run_inf_wrapper
)

stopCluster(clusters)

###========== Collect results =============
res_findN2_set2 <- collect_results(
    design_matrix = design_matrixN2,
    results_folder = results_folder,
    finding = "N2",
    pair = 2
)
res_time_findN2_set2 <- collect_times(
    design_matrix = design_matrixN2,
    pair = 2,
    finding = "N2",
    results_folder = results_folder
)


##============ Find cluster size =============

###=========== Design matrix ==============
n2 <- c(30, 60, 90)
design_matrixN1 <- expand.grid(var_u0, eff_size_prob, BF_thres, eta, n2, fixed <- "n2")
colnames(design_matrixN1) <- c("var_u0", "eff_size", "BF_threshold", "eta", "n2", "fixed")
nrow_designN1 <- nrow(design_matrixN1)
design_matrixN1 <- mutate(design_matrixN1, seed = as.integer(sample(2^32 /
                                                                        2, n())))
write_parquet(design_matrixN1, "design_matrix_findN1_set2")

###========== Running the simulation ============
run_inf_wrapper <- function(Row) {
    run_inf(
        Row = Row,
        design_matrix = design_matrixN1,
        ndatasets = ndatasets,
        Max = Max,
        batch_size = batch_size,
        results_folder = results_folder
    )
}
clusters <- makeForkCluster(detectCores() * 0.75) #For Linux

output <- parallel::parLapply(
    cl = clusters,
    X = 1:nrow_designN1,
    fun = run_null_wrapper
)
stopCluster(clusters)

###========== Collect results ============
res_findN1_set2 <- collect_results(
    design_matrix = design_matrixN1,
    results_folder = results_folder,
    finding = "N1",
    pair = 2
)
res_time_findN1_set2 <- collect_times(
    design_matrix = design_matrixN1,
    pair = 2,
    finding = "N1",
    results_folder = results_folder
)

##================= Plots ====================