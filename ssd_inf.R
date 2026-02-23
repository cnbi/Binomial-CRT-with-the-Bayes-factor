################ FIND SAMPLE SIZE FOR INFORMATIVE HYPOTHESES ####################

SSD_crt_inf_binary <- function(p_intv,
                               p_ctrl,
                               n1 = 15,
                               n2 = 30,
                               ndatasets = 1000,
                               var_u0,
                               BF_thresh1,
                               eta1 = 0.8,
                               fixed = "n2",
                               max = 1000,
                               batch_size = 100,
                               seed) {
    # Libraries ----
    library(lme4)
    library(dplyr)
    library(stargazer)
    library(bain)
    
    #Functions ----------------
    source("data_generation.R")
    source("helper_functions.R")
    source("print_results.R")
    
    # Starting values ----------------------------------------------------------
    conditions_met <- FALSE          #Indication we met the power criteria.
    ultimate_sample_sizes <- FALSE  #Indication that we found the required sample size.
    results_H1 <- matrix(NA, nrow = ndatasets, ncol = 4)

    # Warnings
    if (is.numeric(
        c(
            n1,
            n2,
            ndatasets,
            BF_thresh1,
            var_u0,
            eta1,
            max,
            batch_size
        )
    ) == FALSE)
    stop("All arguments, except 'fixed', must be numeric")
    if (eta1 > 1)
        stop("The probability of exceeding Bayes Factor threshold cannot be larger than 1")
    if (is.character(fixed) == FALSE)
        stop("Fixed can only be a character indicating n1 or n2.")
    if (fixed %in% c("n1", "n2") == FALSE)
        stop("Fixed can only be character indicating n1 or n2.")
    if (p_ctrl > 1 | p_intv > 1)
        stop("Proportion with the outcome of interest cannot be larger than 1")
    if (p_ctrl == 0 )
        stop("Proportion with the outcome of interest in control condition cannot be equal to 0")
    
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
    hypothesis2 <- "Intervention<Control"
    final_SSD <- vector(mode = "list", length = b_fract)
    previous_high <- 0
    previous_eta <- 0
    current_eta <- 0
    singular_warn <- 0
    
    # Simulation of data and evaluation of condition  ----------------------------------
    while (ultimate_sample_sizes == FALSE) {
        # If H1 is true
        data_H1 <- do.call(
            gen_CRT_binarydata,
            list(ndatasets, n1, n2, var_u0, p_intv, p_ctrl, batch_size = batch_size, seed)
        )
        
        #Approximated adjusted fractional Bayes factors------------------------------
        n_eff_H1 <- ((n1 * n2) / (1 + (n1 - 1) * data_H1$rho_data)) / 2
        n_eff_H1_list <- lapply(n_eff_H1, function(x) c(x, x))

        output_bain_H1 <- Map(
            bain,
            x = data_H1$estimates,
            hypothesis = "intervention>control",
            n = n_eff_H1_list,
            Sigma = data_H1$cov_list,
            group_parameters = 1
        )
        
        # Results ---------------------------------------------------------------------
        #browser()
        results_H1[, 1] <- sapply(output_bain_H1, extract_res_bain, "BF.u", 1) #H1vsHu
        results_H1[, 2] <- sapply(output_bain_H1, extract_res_bain, "PMPc", 1) #PMP1
        results_H1[, 3] <- sapply(output_bain_H1, extract_res_bain, "BF.c", 1) #H1vs2
        results_H1[, 4] <- sapply(output_bain_H1, extract_res_bain, "PMPc", 3) #PMP2
        colnames(results_H1) <- c("BF.1u", "PMP.1", "BF.12", "PMP.2")
        
        #Evaluation of condition -------------------------------------------
        # Proportion
        prop_BF12 <- length(which(results_H1[, "BF.12"] > BF_thresh1)) / ndatasets

        # Evaluation
        ifelse(prop_BF12 > eta1 |
                        prop_BF12 == eta1,
               conditions_met <- TRUE,
               conditions_met <- FALSE
        )
        current_eta <- prop_BF12
        # print("Bayes factor check!")
        
        # Binary search algorithm ------------------------------------------
        updated_sample <- binary_search_ineq(
            condition_met = conditions_met,
            fixed = fixed,
            n1 = n1,
            n2 = n2,
            low = low,
            high = high,
            max = max,
            eta = eta1,
            current_eta = current_eta,
            previous_eta = previous_eta,
            previous_high = previous_high,
            min_sample = min_sample,
            prop_BF12 = prop_BF12,
            results_H1 = results_H1,
            data_H1 = data_H1
        )
        list2env(updated_sample, environment())
        rm(data_H1)
    } # Finishes final sample size loop
    # Final output
    final_SSD[[b_fract + 1]] <- BF_thresh1
    final_SSD[[b_fract + 2]] <- eta1
    
    print_results(list_results = final_SSD,
                  type_hyp = "ineq")
    if (any(singular_warn > 0))
        warning(
            "At least one of the fitted models is singular. For more information about singularity see help('isSingular').
                               The number of models that are singular can be found in the output object."
        )
    invisible(final_SSD)
}