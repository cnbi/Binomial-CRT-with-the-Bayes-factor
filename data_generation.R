################## BINARY DATA GENERATION ############################

#"@title Generation of data sets for two treatment-condition cluster randomized trial
#"@description
#"@arguments
## ndatasets: Numeric. Number of data sets that the user wants to generate to determine the sample size
## n1: Numeric. Cluster size
## n2: Numeric. Total number of clusters.
## var_u0: Numeric. Between-cluster variance. Variance at the cluster level.
## mean_interv: Numeric. Equivalent to the effect size.
## batch_size: This parameter determines the size of batches used during the fitting of the multilevel model.
## p_intv: Probability in the intervention condition of obtaining the desired outcome
## p_ctrl: Probability in the control condition of the ocurrence of the desired outcome.



gen_CRT_binarydata <- function(ndatasets = ndatasets,
                               n1 = n1,
                               n2 = n2,
                               var_u0 = var_u0,
                               p_intv,
                               p_ctrl,
                               batch_size,
                               seed) {
    
    interv_logit <- log(p_intv / (1 - p_intv))
    control_logit <- log(p_ctrl / (1 - p_ctrl))
    # Create variables id  of the cluster and condition
    id <- rep(1:n2, each = n1)
    if (n2 %% 2 == 0) {
        condition <- rep(c(0, 1), each = n1 * n2 / 2)
    } else {
        # Odd number of clusters
        # the extra cluster goes to control condition
        half <- floor(n2 / 2)
        condition <- c(rep(0, n1 * half), rep(1, n1 * half), rep(0, n1))
    }
    # Dummy variables for no intercept model
    intervention <- condition
    control <- 1 - intervention
    marker <- 0
    
    # Tables for results
    output_glmer <- vector(mode = "list", length = ndatasets)
    data_list <- vector(mode = "list", length = ndatasets)
    if (missing(seed)) {
        seeds <- sample(2^32 / 2, ndatasets)
    } else {
        set.seed(seed)
        seeds <- sample(2^32 / 2, ndatasets)
    }

    # Data generation ----------------------------------------------------------
    data_list <- lapply(seeds, function(s) {
        set.seed(s)
        u0 <- rnorm(n2, 0, sqrt(var_u0))
        u0 <- rep(u0, each = n1)
        resp_logit <- control_logit * control + interv_logit * intervention + u0
        prob <- exp(resp_logit) / (1 + exp(resp_logit))
        resp <- rbinom(n1 * n2, size = 1, prob = prob)
        data_matrix <- cbind(id,
                             u0,
                             control_logit,
                             control,
                             interv_logit,
                             intervention,
                             resp_logit,
                             resp)
        data_matrix
    })
    
    data_list <- lapply(data_list, as.data.frame) # Necessary to use glmer
    
    # Multilevel analysis --------------------------------------------------------
    # Batches
    batch_size <- batch_size
    ifelse((ndatasets / batch_size) %% 1 == 0,
           num_batches <- ndatasets / batch_size,
           num_batches <- (ndatasets / batch_size) + 1
    )
    for (batch in seq(num_batches)) {
        #Indexes
        start_index <- (batch_size * (batch - 1)) + 1
        end_index <- min(batch * batch_size, ndatasets)
        #Multilevel fitting
        output_glmer[start_index:end_index] <- lapply(data_list[start_index:end_index], fit_glmer)
    }
    marker <- lapply(output_glmer, marker_func) # Mark singularity
    singular_datasets <- Reduce("+", marker) # How many are singular?
    
    estimates <- lapply(output_glmer, fixef)                 # Means
    cov_intervention <- lapply(output_glmer, varcov, "intervention")      # Covariance
    cov_control <- lapply(output_glmer, varcov, "control")
    cov_list <- Map(list, "intervention" = cov_intervention, "control" = cov_control)
    var_u0_data <- unlist(lapply(output_glmer, get_variance))
    total_var_data <- (pi^2) / 3 + var_u0_data
    rho_data <- var_u0_data / total_var_data
    
    rm(
        id,
        condition,
        intervention,
        control,
        control_logit,
        interv_logit,
        data_list,
        output_glmer,
        batch_size,
        cov_intervention,
        cov_control,
        var_u0_data,
        total_var_data
    )
    
    return(
        output <- list(
            "rho_data" = rho_data,
            "estimates" = estimates, # These are in logit scale!!
            "cov_list" = cov_list,
            "singularity" = singular_datasets
        )
    )
}
