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
source("ssd_inf.R")
source("helpers_simulation.R")

# common factors --------------
calc_var_u0 <- function(x){
    num <- (pi^2) / 3 * x
    v_u0 <- num / (1 - x)
    return(v_u0)
}
icc <- c(0.01, 0.05, 0.075, 0.1)
var_u0 <- calc_var_u0(icc)
eff_size_OR <- c(1.5, 2, 3)
logit_beta <- log(eff_size_OR)
p_ctrl <- c(0.05, 0.1, 0.2, 0.4)
logit_ctrl <- vector(mode = "numeric")
for (i in 1:length(p_ctrl)) {
    logit_ <- log(p_ctrl[i]/(1 - p_ctrl[i]))
    logit_ctrl <- c(logit_ctrl, logit_)
}
p_int <- vector(mode = "numeric")
eff_sizes <- expand.grid(logit_beta = logit_beta, logit_ctrl = logit_ctrl)
p_int <- 1/(1 + exp(-(eff_sizes[, 1] + eff_sizes[, 2])))
p_ctr <-  exp(eff_sizes[, "logit_ctrl"])/(1 + exp(eff_sizes[, "logit_ctrl"]))
eff_sizes <- cbind(eff_sizes, p_int, p_ctr)
BF_thres <- c(1, 3, 5)
b_fract <- 3
eta <- 0.8
ndatasets <- 5000
Max <- 500
batch_size <- 500

########## Hypothesis Set 1: Equality vs. informative ##########################
## Find the number of clusters ====
path <- "~/"
results_folder <- "results_set1v2"
dir.create(results_folder)

###====== Design matrix ======
n1 <- c(5, 10, 40)
design_matrixN2 <- expand.grid(var_u0, BF_thres, eta, n1, fixed <- "n1")
colnames(design_matrixN2) <- c("var_u0", "BF_threshold", "eta", "n1", "fixed")
design_matrixN2 <- cross_join(eff_sizes, design_matrixN2)
design_matrixN2 <- mutate(design_matrixN2, seed = as.integer(sample(2^32 /
                                                                        2, n())))
nrow_designN2 <- nrow(design_matrixN2)
write_parquet(design_matrixN2, "design_matrix_findN2_set1_v2")

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
clusters <- makeForkCluster(detectCores() * 0.8)

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
    fun = run_null_wrapper
)

stopCluster(clusters)

###======= Collect results ========
results_folder <- "results_set1.2/results_set1.3"
#pair = 1 refering to the pair of hypotheses H0 vs H1
res_findN2_set1 <- collect_results(
    design_matrix = design_matrixN2,
    results_folder = results_folder,
    finding = "N2",
    pair = 1,
    name_results = "FindN2Row", b = 3
)
res_time_findN2_set1 <- collect_times(
    design_matrix = design_matrixN2,
    pair = 1,
    finding = "N2",
    results_folder = results_folder,
    times_name = "timeN2Row"
)

res_findN2_set1 <- rbind(final_results_findN2_set1, final_results_findN2_set1_extra) 
saveRDS(res_findN2_set1, file = "res_findN2_set1.RDS")
### ===========Plots ====================
# Bayes Factor H0 vs H1, n1, var_u0, p_int
ggplot(res_findN2_set1[res_findN2_set1$b == 1, ], aes(y = log(median.BF01), x = n2.final, color = as.factor(n1.final), shape = as.factor(n1.final))) +
    geom_point() + geom_line() +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    facet_grid(rows = vars(var_u0), cols = vars(p_int)) +
    labs(title = "Bayes Factor H0 vs H1", color = "Cluster \nsizes", shape = "Cluster \nsizes") +
    xlab("Number of clusters") + ylab("Bayes Factor") + theme(legend.position = "bottom")

# n2, n1, BF_thresh, p_int, var_u0
ggplot(res_findN2_set1[res_findN2_set1$b == 1, ], 
       aes(x = n1.final, y = n2.final, color = as.factor(BF_threshold), 
           shape = as.factor(BF_threshold))) + geom_point() + geom_line() +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    facet_grid(rows = vars(var_u0), cols = vars(p_int)) +
    labs(title = "Determining Number of Clusters in Function of Cluster Sizes, Bayes Factor Thresholds, \nresponse rate,and between-cluster variance",  
         color = "Bayes Factor \nThresholds", shape = "Bayes Factor \nThresholds") +
    xlab("Cluster sizes") + ylab("Number of clusters") + theme(legend.position = "bottom") +
    ylim(c(0, 100))

# n2, n1, BF_thresh, p_int, var_u0 including b
ggplot(res_findN2_set1, 
       aes(x = n1.final, y = n2.final, color = as.factor(BF_threshold), 
           shape = as.factor(b))) + geom_point() + geom_line() +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    facet_grid(rows = vars(var_u0), cols = vars(p_int)) +
    labs(title = "Determining Number of Clusters in Function of Cluster Sizes, Bayes Factor Thresholds, \nresponse rate,and between-cluster variance",  
         color = "Bayes Factor \nThresholds", shape = "Bayes Factor \nThresholds") +
    xlab("Cluster sizes") + ylab("Number of clusters") + theme(legend.position = "bottom")

# Table
library(xtable)
table_results <- c("var_u0", "p_int", "BF_threshold", "n1.final", "n2.final", "eta.BF10", "eta.BF01")
round_col <- c("BF_threshold", "n1", "n2.final")
res_findN2_set1[which(res_findN2_set1$b == 2 ),table_results]



##========= Find cluster size ==========

results_folder <- "results_set1.2"

###======= Design matrix ============
n2 <- c(30, 60, 90)
design_matrixN1 <- expand.grid(var_u0, BF_thres, eta, n2, fixed <- "n2")
colnames(design_matrixN1) <- c("var_u0", "BF_threshold", "eta", "n2", "fixed")
design_matrixN1 <- cross_join(eff_sizes, design_matrixN1)
design_matrixN1 <- mutate(design_matrixN1, seed = as.integer(sample(2^32 /
                                                                        2, n())))
nrow_designN1 <- nrow(design_matrixN1)
write_parquet(design_matrixN1, "design_matrix_findN1_set1_v2")

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
clusters <- makeForkCluster(detectCores() * 0.8) #For Linux

output <- parallel::parLapply(
    cl = clusters,
    X =  1:nrow_designN1,
    fun = run_null_wrapper
)
stopCluster(clusters)

### Collect results -----
#pair = 1 refering to the pair of hypotheses H0 vs H1
res_findN1_set1 <- collect_results(
    design_matrix = design_matrixN1,
    results_folder = results_folder,
    finding = "N1",
    pair = 1,
    name_results = "FindN1Row",
    b = 3
    
)
res_time_findN1_set1 <- collect_times(
    design_matrix = design_matrixN1,
    pair = 1,
    finding = "N1",
    results_folder = results_folder,
    times_name = "timeN1Row"
)

##================ Plots =================

# Bayes Factor H0 vs H1, n1, var_u0, p_int
ggplot(res_findN1_set1[res_findN1_set1$b == 1, ], aes(y = log(median.BF01), x = n1.final, color = as.factor(n2), shape = as.factor(n2))) +
    geom_point() + geom_line() +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    facet_grid(rows = vars(var_u0), cols = vars(p_int)) +
    labs(title = "Bayes Factor H0 vs H1", color = "Cluster \nsizes", shape = "Cluster \nsizes") +
    xlab("Number of clusters") + ylab("Bayes Factor") + theme(legend.position = "bottom")

# n2, n1, BF_thresh, p_int, var_u0
ggplot(res_findN1_set1[res_findN1_set1$b == 1, ], 
       aes(x = n2.final, y = n1.final, color = as.factor(BF_threshold), 
           shape = as.factor(BF_threshold))) + geom_point() + geom_line() +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    facet_grid(rows = vars(var_u0), cols = vars(p_int)) +
    labs(title = "Determining Number of Clusters in Function of Cluster Sizes, Bayes Factor Thresholds, \nresponse rate,and between-cluster variance",  
         color = "Bayes Factor \nThresholds", shape = "Bayes Factor \nThresholds") +
    xlab("Number of clusters") + ylab("Cluster sizes") + theme(legend.position = "bottom")

# n2, n1, BF_thresh, p_int, var_u0 including b
ggplot(res_findN1_set1, 
       aes(x = n2.final, y = n1.final, color = as.factor(BF_threshold), 
           shape = as.factor(b))) + geom_point() + geom_line() +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    facet_grid(rows = vars(var_u0), cols = vars(p_int)) +
    labs(title = "Determining Number of Clusters in Function of Cluster Sizes, Bayes Factor Thresholds, \nresponse rate,and between-cluster variance",  
         color = "Bayes Factor \nThresholds", shape = "Bayes Factor \nThresholds") +
    xlab("Number of clusters") + ylab("Cluster sizes") + theme(legend.position = "bottom")

#Table
table_results <- c("var_u0", "p_int", "BF_threshold", "n1.final", "n2.final", "eta.BF10", "eta.BF01")
results_colnames <- c("var_u0", "$BF_{thresh}$", "N1", "N2", "$P(BF_{10}>BF_{thresh})$", "$P(BF_{01}>BF_{thresh})$")
res_findN1_set1[which(res_findN1_set1$b == 2 & res_findN1_set1$var_u0 != 0.25 & res_findN1_set1$p_int > 0.7), table_results]

############## Hypothesis Set 2: Informative vs. informative ###################
icc <- c(0.01, 0.05, 0.075, 0.1)
var_u0 <- calc_var_u0(icc)
eff_size_OR <- c(1.5, 2, 3)
logit_beta <- log(eff_size_OR)
p_ctrl <- c(0.05, 0.1, 0.2, 0.4)
logit_ctrl <- vector(mode = "numeric")
for (i in 1:length(p_ctrl)) {
    logit_ <- log(p_ctrl[i]/(1 - p_ctrl[i]))
    logit_ctrl <- c(logit_ctrl, logit_)
}
p_int <- vector(mode = "numeric")
eff_sizes <- expand.grid(logit_beta = logit_beta, logit_ctrl = logit_ctrl)
p_int <- 1/(1 + exp(-(eff_sizes[, 1] + eff_sizes[, 2])))
p_ctr <-  exp(eff_sizes[, "logit_ctrl"])/(1 + exp(eff_sizes[, "logit_ctrl"]))
eff_sizes <- cbind(eff_sizes, p_int, p_ctr)
BF_thres <- c(1, 3, 5)
eta <- 0.8
ndatasets <- 5000
Max <- 500
batch_size <- 500


##========== Find the number of clusters =========
path <- "~/"
results_folder <- "results_set2v2"
dir.create(results_folder)

###========= Design matrix ===========
n1 <- c(5, 10, 40)
design_matrixN2 <- expand.grid(var_u0, BF_thres, eta, n1, fixed <- "n1")
colnames(design_matrixN2) <- c("var_u0", "BF_threshold", "eta", "n1", "fixed")
design_matrixN2 <- cross_join(eff_sizes, design_matrixN2)
design_matrixN2 <- mutate(design_matrixN2, seed = as.integer(sample(2^32 /
                                                                        2, n())))
nrow_designN2 <- nrow(design_matrixN2)
write_parquet(design_matrixN2, "design_matrix_findN2_set2_v2")
#design_matrixN2 <- read_parquet("design_matrix_findN2_set2")

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

clusters <- makeForkCluster(detectCores() * 0.85)


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
    pair = 2,
    name_results = "ResultsN2Row"
)
res_time_findN2_set2 <- collect_times(
    design_matrix = design_matrixN2,
    pair = 2,
    finding = "N2",
    results_folder = results_folder,
    times_name = "timeN2Row"
)

res_findN2_set2 <- rbind(final_results_findN2_set2 , final_results_findN2_set2_extra)
save(res_findN2_set2, file = "res_findN2_set2.RDS")

### ========== Plots ==================
# Bayes Factor H0 vs H1, n1, var_u0, p_int
ggplot(res_findN2_set2, aes(y = log(median.BF12), x = n2.final, color = as.factor(n1.final), shape = as.factor(n1.final))) +
    geom_point() + geom_line() +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    facet_grid(rows = vars(var_u0), cols = vars(p_int)) +
    labs(title = "Bayes Factor H0 vs H1", color = "Cluster \nsizes", shape = "Cluster \nsizes") +
    xlab("Number of clusters") + ylab("Bayes Factor") + theme(legend.position = "bottom")

# n2, n1, BF_thresh, p_int, var_u0
ggplot(res_findN2_set2, 
       aes(x = n1.final, y = n2.final, color = as.factor(BF_threshold), 
           shape = as.factor(BF_threshold))) + geom_point() + geom_line() +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    facet_grid(rows = vars(var_u0), cols = vars(p_int)) +
    labs(title = "Determining Number of Clusters in Function of Cluster Sizes, Bayes Factor Thresholds, \nresponse rate,and between-cluster variance",  
         color = "Bayes Factor \nThresholds", shape = "Bayes Factor \nThresholds") +
    xlab("Cluster sizes") + ylab("Number of clusters") + theme(legend.position = "bottom") + 
    ylim(c(0,110))

#Table
table_results <- c("var_u0", "p_int", "BF_threshold", "n1.final", "n2.final", "eta.BF12")
res_findN2_set2[which(res_findN2_set2$p_int < 0.8 & res_findN2_set2$var_u0 != 0.25), table_results]

##============ Find cluster size =============

###=========== Design matrix ==============
n2 <- c(30, 60, 90)
design_matrixN1 <- expand.grid(var_u0, BF_thres, eta, n2, fixed <- "n2")
colnames(design_matrixN1) <- c("var_u0", "BF_threshold", "eta", "n2", "fixed")
design_matrixN1 <- cross_join(eff_sizes, design_matrixN1)
nrow_designN1 <- nrow(design_matrixN1)
design_matrixN1 <- mutate(design_matrixN1, seed = as.integer(sample(2^32 /
                                                                        2, n())))
write_parquet(design_matrixN1, "design_matrix_findN1_set2_v2")

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
clusters <- makeForkCluster(detectCores() * 0.8) #For Linux

output <- parallel::parLapply(
    cl = clusters,
    X = c(16, 28),
    fun = run_inf_wrapper
)
stopCluster(clusters)

run_inf_wrapper(1)

###========== Collect results ============
res_findN1_set2 <- collect_results(
    design_matrix = design_matrixN1,
    results_folder = results_folder,
    finding = "N1",
    pair = 2,
    name_results = "ResultsN1Row"
)
res_time_findN1_set2 <- collect_times(
    design_matrix = design_matrixN1,
    pair = 2,
    finding = "N1",
    results_folder = results_folder,
    times_name = "timeN1Row"
)

##================= Plots ====================
# Bayes Factor H0 vs H1, n1, var_u0, p_int
ggplot(res_findN1_set2, aes(y = log(median.BF12), x = n1.final, color = as.factor(n2.final), shape = as.factor(n2.final))) +
    geom_point() + geom_line() +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    facet_grid(rows = vars(var_u0), cols = vars(p_int)) +
    labs(title = "Bayes Factor H0 vs H1", color = "Number of \nclusters", shape = "Number of \nclusters") +
    xlab("Cluster sizes") + ylab("Bayes Factor") + theme(legend.position = "bottom")

# n2, n1, BF_thresh, p_int, var_u0
ggplot(res_findN1_set2, 
       aes(x = n2.final, y = n1.final, color = as.factor(BF_threshold), 
           shape = as.factor(BF_threshold))) + geom_point() + geom_line() +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    facet_grid(rows = vars(var_u0), cols = vars(p_int)) +
    labs(title = "Determining Number of Clusters in Function of Cluster Sizes, Bayes Factor Thresholds, \nresponse rate,and between-cluster variance",  
         color = "Bayes Factor \nThresholds", shape = "Bayes Factor \nThresholds") +
    xlab("Number of clusters") + ylab("Cluster sizes") + theme(legend.position = "bottom") + 
    ylim(c(0,80))

#Table
table_results <- c("var_u0", "p_int", "BF_threshold", "n2.final", "n1.final", "eta.BF12")
res_findN1_set2[which(res_findN1_set2$p_int < 0.8 & res_findN1_set2$var_u0 != 0.25), table_results]
