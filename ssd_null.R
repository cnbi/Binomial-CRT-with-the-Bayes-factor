########## FIND SAMPLE SIZE FOR CRT WITH BINARY OUTCOME ##############



SSD_crt_null_binary <- function(p_intv, p_ctrl = 0, n1 = 15, n2 = 30, ndatasets = 1000, 
                                var_u0, BF_thresh1, BF_thresh0, eta1 = 0.8, eta0 = 0.8,
                                fixed = "n2", b_fract = 3, max = 1000, batch_size = 100) {
    # Libraries ----
    library(lme4)
    library(dplyr)
    
    # Warnings
    if (is.numeric(c(n1, n2, ndatasets, BF_thresh1, BF_thresh0, var_u0,
                     eta1, eta0, b_fract, max, batch_size)) == FALSE) 
        stop("All arguments, except 'fixed', must be numeric")
    if (eta1 > 1 | eta0 > 1) stop("The probability of exceeding Bayes Factor threshold cannot be larger than 1")
    if (eta1 < 0 | eta0 < 0) stop("The probability of exceeding Bayes Factor threshold must be a positive value")
    if (is.character(fixed) == FALSE) stop("Fixed can only be a character indicating n1 or n2.")
    if (fixed %in% c("n1", "n2") == FALSE) stop("Fixed can only be character indicating n1 or n2.")
    if ((b_fract == round(b_fract)) == FALSE) stop("The fraction of information (b) must be an integer")
    
    #Functions ----------------
    source("data_generation.R")
    source("helper_functions.R")
    source("print_results.R")
    source("bayes_factor.R")
    
    # Starting values ----------------------------------------------------------
    conditions_met <- FALSE          #Indication we met the power criteria.
    ultimate_sample_sizes <- FALSE  #Indication that we found the required sample size.
    results_H0 <- matrix(NA, nrow = ndatasets, ncol = 4)
    results_H1 <- matrix(NA, nrow = ndatasets, ncol = 4)
    
    # Binary search start ------------------------------------------------------
    if (fixed == "n1") {
        min_sample <- 6                     # Minimum number of clusters
        low <- min_sample                   #lower bound
    } else if (fixed == "n2") {
        min_sample <- 5                     # Minimum cluster size
        low <- min_sample                   #lower bound
    }
    high <- max                    #higher bound
    
    # Hypotheses ----------------------------------------------------------------
    hypothesis1 <- "Intervention>Control"
    null <- "Intervention=Control"
    final_SSD <- vector(mode = "list", length = b_fract)
    type <- "Equality"
    b <- 1
    previous_high <- 0
    previous_eta <- 0
    current_eta <- 0
    singular_warn <- 0
    
    # Simulation of data and evaluation of condition  ----------------------------------
    while (ultimate_sample_sizes == FALSE) {
        # If H1 is true
        data_H1 <- do.call(gen_CRT_binarydata, list(ndatasets, n1, n2, var_u0,
                                                    p_intv, p_ctrl = 0,
                                                    batch_size = batch_size))
        
        # If H0 is true
        data_H0 <- do.call(gen_CRT_binarydata, list(ndatasets, n1, n2, var_u0,
                                                    p_intv = 0, p_ctrl = 0,
                                                    batch_size = batch_size))
        
        while (b < (b_fract + 1)) {
            #Approximated adjusted fractional Bayes factors------------------------------
            n_eff_H1 <- ((n1 * n2) / (1 + (n1 - 1) * data_H1$rho_data)) / 2
            output_AAFBF_H1 <- Map(calc_aafbf, type, data_H1$estimates, data_H1$cov_list, list(b), n_eff_H1)
            
            n_eff_H0 <- ((n1 * n2) / (1 + (n1 - 1) * data_H0$rho_data)) / 2
            output_AAFBF_H0 <- Map(calc_aafbf, type, data_H0$estimates, data_H0$cov_list, list(b), n_eff_H0)
            
            # Results ---------------------------------------------------------------------
            results_H1[, 1] <- unlist(lapply(output_AAFBF_H1, extract_res, 1)) # Bayes factor H1vsH0
            results_H1[, 2] <- unlist(lapply(output_AAFBF_H1, extract_res, 4)) #posterior model probabilities of H1
            results_H1[, 3] <- unlist(lapply(output_AAFBF_H1, extract_res, 2)) # Bayes factor H0vsH1
            results_H1[, 4] <- unlist(lapply(output_AAFBF_H1, extract_res, 3)) #posterior model probabilities of H0
            
            colnames(results_H1) <- c("BF.10", "PMP.1", "BF.01", "PMP.0")
            
            results_H0[, 1] <- unlist(lapply(output_AAFBF_H0, extract_res, 1)) # Bayes factor H1vsH0
            results_H0[, 2] <- unlist(lapply(output_AAFBF_H0, extract_res, 4)) #posterior model probabilities of H1
            results_H0[, 3] <- unlist(lapply(output_AAFBF_H0, extract_res, 2)) # Bayes factor H0vsH1
            results_H0[, 4] <- unlist(lapply(output_AAFBF_H0, extract_res, 3)) #posterior model probabilities of H0
            
            colnames(results_H0) <- c("BF.10", "PMP.1", "BF.01", "PMP.0")
            
            #Evaluation of condition -------------------------------------------
            # Proportion
            prop_BF10 <- length(which(results_H1[, "BF.10"] > BF_thresh1)) / ndatasets
            prop_BF01 <- length(which(results_H0[, "BF.01"] > BF_thresh0)) / ndatasets
            
            # Evaluation
            ifelse(prop_BF01 > eta0 & prop_BF10 > eta1, conditions_met <- TRUE, conditions_met <- FALSE)
            previous_eta <- current_eta
            ## Finding the lowest proportion
            if (prop_BF01 < eta0 & prop_BF10 < eta1) {
                current_eta <- min(prop_BF10, prop_BF01)
                if (abs(prop_BF01 - eta0) < abs(prop_BF10 - eta1)) {
                    eta <- eta0
                } else {
                    eta <- eta1
                }
            } else if (prop_BF01 < eta0 | prop_BF10 < eta1 ) {
                if (prop_BF01 < eta0) {
                    current_eta <- prop_BF01
                    eta <- eta0
                } else if (prop_BF10 < eta1) {
                    current_eta <- prop_BF10
                    eta <- eta1
                }
            } else if (conditions_met) {
                diff1 <- eta1 - prop_BF10
                diff0 <- eta0 - prop_BF01
                current_eta <- ifelse(diff1 < diff0, prop_BF10, prop_BF01)
                eta <- ifelse(current_eta == prop_BF10, eta1, eta0)
            }
            # print("Bayes factor check!")
            
            # Binary search algorithm ------------------------------------------
            updated_sample <- binary_search_eq(condition_met = conditions_met,
                                               fixed = fixed, n1 = n1, n2 = n2,
                                               low = low, high = high, max = max,
                                               eta = eta,
                                               current_eta = current_eta,
                                               previous_eta = previous_eta,
                                               previous_high = previous_high,
                                               min_sample = min_sample, b = b,
                                               prop_BF01 = prop_BF01,
                                               prop_BF10 = prop_BF10,
                                               results_H0 = results_H0,
                                               results_H1 = results_H1,
                                               data_H0 = data_H0,
                                               data_H1 = data_H1,
                                               final_SSD = final_SSD)
            list2env(updated_sample, environment())
        } # Finishes b loop
        # Break loop
        if (b == b_fract + 1) {
            ultimate_sample_sizes <- TRUE
        }
        rm(data_H0, data_H1)
    } # Finishes final sample size loop
    # Final output
    final_SSD[[b + 1]] <- c(BF_thresh0, BF_thresh1)
    final_SSD[[b + 2]] <- c(eta0, eta1)
    
    print_results(list_results = final_SSD, type_hyp = "eq", b = b)
    if (any(singular_warn > 0)) warning("At least one of the fitted models is singular. For more information about singularity see help('isSingular').
                               The number of models that are singular can be found in the output object.")
    invisible(final_SSD)
}