########################## SMALL FUNCTIONS ###########################

# Model fitting
fit_glmer <- function(x) {
    fitted_model <- tryCatch({
        suppressMessages(suppressWarnings({
            glmer(resp ~ intervention + control - 1 + (1 | id), 
                  data = x, family = binomial)
        }))
    }, error = function(e) return(NULL))
    
    # Check for Convergence Issues
    # If the model didn't converge (max|grad| > 0.002), try another optimiser
    if (!is.null(fitted_model) && !is.null(fitted_model@optinfo$conv$lme4$messages)) {
        
        # Check if any message contains "failed to converge"
        if (any(grepl("failed to converge", fitted_model@optinfo$conv$lme4$messages))) {
            
            message("Initial model failed to converge. Re-running with nloptwrap optimiser")
            
            fitted_model <- suppressMessages(suppressWarnings({
                glmer(resp ~ intervention + control - 1 + (1 | id), 
                      data = x, 
                      family = binomial,
                      control = glmerControl(optimizer = "nloptwrap", 
                                             optCtrl = list(maxfun = 2e5)))
            }))
        }
    }
     # Try another optimiser
    if (is.null(fitted_model) || !is.null(fitted_model@optinfo$conv$lme4$messages)) {
        fitted_model <- tryCatch({suppressMessages(suppressWarnings({
            glmer(resp ~ intervention + control - 1 + (1 | id), 
                  data = x, 
                  family = binomial,
                  control = glmerControl(optimizer = "bobyqa", 
                                         optCtrl = list(maxfun = 2e5)))
        }))
        }, error = function(e) return(NULL))
    }
    #browser()
    # Check for non-positive definite matrix
    if (!is.null(fitted_model)) {
        hessian_ <- fitted_model@optinfo$derivs$Hessian
        if (!is.null(hessian_)) {
            eigen_val <- eigen(hessian_)$values
            if (any(eigen_val <= 0)) {
                attr(fitted_model, "NonPositDef_warning") <- TRUE
            } else {
                attr(fitted_model, "NonPositDef_warning") <- FALSE
            }
        } else {
            attr(fitted_model, "NonPositDef_warning") <- isSingular(fitted_model)
        }
        
    }
    
    return(fitted_model)
}

fit_glmer2 <- function(x, use_glm) {
    fitted_model <-  suppressMessages(suppressWarnings({
        glmer(resp ~ intervention + control - 1 + (1 | id), 
              data = x, family = binomial)
    }))
    
    # Retry with nlotpwrap optimiser if singular
    if (!is.null(fitted_model) && isSingular(fitted_model)) {
        fitted_model <- tryCatch(suppressMessages(suppressWarnings(
            glmer(resp ~ intervention + control - 1 + (1 | id), 
                  data = x, 
                  family = binomial,
                  control = glmerControl(optimizer = "nloptwrap", 
                                         optCtrl = list(maxfun = 2e5)))
        )),
        error = function(e) NULL)
    }
    
    # Use GLM instead when everything fails
    # if (is.null(fitted_model) || isSingular(fitted_model)) {
    #     if (use_glm) {
    #         fitted_model <- tryCatch({
    #             glm(resp ~ intervention + control - 1, data = x, family = binomial) },
    #             error = function(e) {message("GLM also failed: ", e$message)
    #                 return(NULL)
    #             })
    #     } else if (is.null(fitted_model)) {
    #         message("The glmer optimiser failed.")
    #     } else {
    #         message("Singular fit is detected (var_u0 ≈ 0). Consider use_glm = TRUE.")
    #     }
    # }
    
    # Check for non-positive definite matrix
    if (!is.null(fitted_model)) {
        attr(fitted_model, "use_glm") <- use_glm
        # Check Hessian
        if (!use_glm) {
            nonpos <- tryCatch({
                hessian_ <- fitted_model@optinfo$derivs$Hessian
                if (!is.null(hessian_)) {
                    eigen_val <- eigen(hessian_)$values
                    if (any(eigen_val <= 0)) {
                        attr(fitted_model, "NonPositDef_warning") <- TRUE
                    }
                } else {
                    isSingular(fitted_model)
                }
            }, error = function(e) TRUE)
            
            attr(fitted_model, "NonPositDef_warning") <- nonpos
        }
    }
    return(fitted_model)
}


# glm was used instead
glm_used_func <- function(output.model) {
    if (is.null(output.model)) return(NA)
    isTRUE(attr(output.model, "glm_used"))
}

# Obtain the variance covariance matrix for fixed effects
varcov <- function(output.glmer, name) {
    if (is.null(output.glmer)) return(matrix(NA_real_, 1, 1))
    results <- tryCatch(
        suppressWarnings(vcov(output.glmer)[name, name]),
        error = function(e) {
            message("vcov() failed (non-positive definite matrix): returning NA")
            NA_real_
        }
    )
    return(matrix(results, nrow = 1, ncol = 1))
}

# Obtain variance for random intercept
get_variance <- function(output.glmer) {
    if (is.null(output.glmer) || glm_used_func(output.glmer)) return(NA)
    variance <- as.data.frame(VarCorr(output.glmer))
    value <- variance[, 4]
    return(value)
}

# Marking singular matrices
marker_func <- function(output.glmer) {
    if (is.null(output.glmer)) return(NA)
    if (glm_used_func(output.glmer)) return(NA)
    ifelse(isSingular(output.glmer), marker <- 1, marker <- 0)
}

#Marking non-positive definite
npd_func <- function(output.glmer) {
    if (is.null(output.glmer)) return(NA)
    attr(output.glmer, "NonPositDef_warning")
}

# Check data before bain
check_data <- function(index, data_) {
    estimates <- data_$estimates[[index]]
    covar_var <- data_$cov_list[[index]]
    rho <- data_$rho_data
    
    !is.null(estimates) && !any(is.na(estimates)) && !any(is.infinite(estimates)) &&
        !is.null(covar_var) && !any(is.na(unlist(covar_var))) && !any(is.infinite(unlist(covar_var))) &&
        !is.null(rho) && !any(is.na(rho)) && !any(is.infinite(rho))
}

# Extract results
extract_res_bain <- function(bain_object, element, number) {
    result <- bain_object$fit[[element]][number]
    if (is.null(result)) {
        result <- NaN
    }
    return(result)
}

# Rounding half away from zero
round2 <- function(number, decimals = 0) {
    sign_number <- sign(number)
    number <- abs(number) * 10^decimals
    number <- number + 0.5 + sqrt(.Machine$double.eps)
    number <- trunc(number)
    number <- number / 10 ^ decimals
    number * sign_number
}
# Source: https://stackoverflow.com/questions/66600344/commercial-rounding-in-r-i-e-always-round-up-from-5/66600470#66600470


# Binary Search Algorithm with equality constraint-----------------------------
binary_search_eq <- function(condition_met, fixed, n1, n2, low, high, max, eta,
                             current_eta, previous_eta, previous_high, min_sample,
                             b, prop_BF01, prop_BF10, results_H0, results_H1,
                             data_H0, data_H1, final_SSD){
    if (!condition_met) {
        print(c("Using cluster size:", n1,
                "and number of clusters:", n2,
                "prop_BF01: ", prop_BF01, "prop_BF10: ", prop_BF10,
                "low:", low, "high:", high, "b:", b,  "current_eta:", current_eta,
                "previous_eta:", previous_eta))
        message("Increasing sample size")
        if (fixed == "n1") {
            if (n2 == max)    { # If the sample size reaches the maximum
                final_SSD[[b]] <- list("n1" = n1,
                                       "n2" = n2,
                                       "Proportion.BF01" = prop_BF01,
                                       "Proportion.BF10" = prop_BF10,
                                       "b.frac" = b,
                                       "data_H0" = results_H0,
                                       "data_H1" = results_H1,
                                       "singularity" = cbind(H0 = data_H0$singularity,
                                                             H1 = data_H1$singularity))
                b <- b + 1
                return(list(low = min_sample,
                            high = max,
                            n1 = n1,
                            n2 = n2,
                            previous_eta = 0,
                            previous_high = 0,
                            b = b,
                            final_SSD = final_SSD
                ))
            } else {
                # Increase the number of clusters since eta is too small
                low <- n2                         #lower bound
                high <- high                      #higher bound
                n2 <- round2((low + high) / 2)     #point in the middle
                if (!n2 %% 2 == 0) n2 <- n2 + 1 # To ensure number of clusters is even
                
                # Adjust higher bound when there is a ceiling effect
                if (low + n2 == high * 2) {
                    low <- n2                         #lower bound
                    if (previous_high > 0) {
                        high <- previous_high
                    } else {
                        high <- max                       #higher bound
                    }
                    n2 <- round2((low + high) / 2)     #point in the middle
                }
                return(list(low = low,
                            high = high,
                            n1 = n1,
                            n2 = n2,
                            b = b,
                            previous_eta = current_eta))
            }
        } else if (fixed == "n2") {
            if (n1 == max)    {# If the sample size reaches the maximum
                final_SSD[[b]] <- list("n1" = n1,
                                       "n2" = n2,
                                       "Proportion.BF01" = prop_BF01,
                                       "Proportion.BF10" = prop_BF10,
                                       "b.frac" = b,
                                       "data_H0" = results_H0,
                                       "data_H1" = results_H1,
                                       "singularity" = cbind(H0 = data_H0$singularity,
                                                             H1 = data_H1$singularity))
                b <- b + 1
                return(list(low = min_sample,
                            high = max,
                            n1 = n1,
                            n2 = n2,
                            previous_eta = 0,
                            previous_high = 0,
                            b = b,
                            final_SSD = final_SSD))
            } else {
                # Increase the cluster sizes since eta is too small
                low <- n1                        #lower bound
                high <- high                     #higher bound
                n1 <- round2((low + high) / 2)    #point in the middle
                
                # Adjust higher bound when there is a ceiling effect or no increase of power
                if ((low + n1 == high * 2) | (current_eta == previous_eta)) {
                    low <- n1                        #lower bound
                    #Set the higher bound based on the previous high or the maximum
                    if (previous_high > 0 ) {
                        if (high == previous_high) {
                            high <- max
                        } else {
                            high <- previous_high
                        }
                    } else {
                        high <- max
                    }
                    n1 <- round2((low + high) / 2)    # point in the middle
                }
                return(list(low = low,
                            high = high,
                            n1 = n1,
                            n2 = n2,
                            b = b,
                            previous_eta = current_eta))
            }
        }
    } else if (condition_met) {
        print(c("previous:", previous_eta))
        previous_high <- high
        SSD_object <- list("n1" = n1,
                           "n2" = n2,
                           "Proportion.BF01" = prop_BF01,
                           "Proportion.BF10" = prop_BF10,
                           "b.frac" = b,
                           "data_H0" = results_H0,
                           "data_H1" = results_H1,
                           "singularity" = cbind(H0 = data_H0$singularity,
                                                 H1 = data_H1$singularity))
        print(c("Using cluster size:", n1,
                "and number of clusters:", n2,
                "prop_BF01: ", prop_BF01, "prop_BF10: ", prop_BF10,
                "low:", low, "high:", high, "b:", b))
        if (fixed == "n1") {
            # Eta is close enough to the desired eta or
            # there is no change in eta and the lower bound is close to the middle point
            if ((current_eta - eta < 0.1 & n2 - low == 2) | 
                (previous_eta == current_eta & n2 - low == 2)) {
                final_SSD[[b]] <- SSD_object
                b <- b + 1
                low <- min_sample
                high <- max
                
                return(list(
                    low = low,
                    high = max,
                    n1 = n1,
                    n2 = n2,
                    b = b,
                    final_SSD = final_SSD,
                    previous_eta = 0,
                    previous_high = 0))
                
            } else if (low + n2 == high * 2) {
                # Adjust higher bound when there is a ceiling effect
                low <- n2                         #lower bound
                high <- n2 * 2                       #higher bound
                n2 <- round2((low + high) / 2)     #point in the middle
                return(list(low = low,
                            high = high,
                            n1 = n1,
                            n2 = n2,
                            b = b,
                            previous_eta = current_eta))  
            } else {
                # Decreasing to find the ultimate number of clusters
                message("Lowering number of clusters")
                low <- low                         #lower bound
                high <- n2                         #higher bound
                n2 <- round2((low + high) / 2)      #point in the middle
                if (!n2 %% 2 == 0) n2 <- n2 + 1
                
                return(list(
                    low = low,
                    high = high,
                    n1 = n1,
                    n2 = n2,
                    previous_high = high, #The higher bound when the power criterion was met
                    previous_eta = current_eta
                ))
            }
        } else if (fixed == "n2") {
            # Eta is close enough to the desired eta or
            # there is no change in eta and the lower bound is close to the middle point or
            # reached the minimum number that meets the Bayesian power condition
            if ((current_eta - eta < 0.1 && current_eta - eta > 0 && n1 - low == 1) ||
                (current_eta == previous_eta && n1 - low == 1) ||
                (current_eta == previous_eta && low + n1 == high * 2)) {
                final_SSD[[b]] <- SSD_object
                b <- b + 1
                low <- min_sample
                high <- max
                
                return(list(low = min_sample,
                            high = max,
                            n1 = n1,
                            n2 = n2,
                            b = b,
                            final_SSD = final_SSD,
                            previous_eta = 0,
                            previous_high = 0))
                
            } else if (low + n1 == high * 2) {
                # Adjust higher bound when there is a ceiling effect
                low <- n1                         #lower bound
                high <- n1*2                       #higher bound
                n1 <- round2((low + high) / 2)     #point in the middle
                return(list(low = low,
                            high = high,
                            n1 = n1,
                            n2 = n2,
                            b = b,
                            previous_eta = current_eta))
            } else {
                message("Lowering cluster size")
                # Decreasing the cluster size to find the ultimate sample size
                low <- low                         #lower bound
                high <- n1                         #higher bound
                n1 <- round2((low + high) / 2)      #point in the middle
                return(list(low = low,
                            high = high,
                            n1 = n1,
                            n2 = n2,
                            previous_high = high,
                            previous_eta = current_eta))
            }
        } # Finish if n2
    } # Finish condition met
}


# Binary Search Algorithm with Inequality Constraints------------------------
binary_search_ineq <- function(condition_met, fixed, n1, n2, low, high, max, eta,
                             current_eta, previous_eta, previous_high, min_sample,
                             prop_BF12, results_H1, data_H1){
    if (!condition_met) {
        print(c("Using cluster size:", n1,
                "and number of clusters:", n2,
                "prop_BF12: ", prop_BF12,
                "low:", low, "high:", high, "current_eta:", current_eta,
                "previous_eta:", previous_eta))
        message("Increasing sample size")
        if (fixed == "n1") {
            if (n2 == max)    { # If the sample size reaches the maximum
                final_SSD <- list("n1" = n1,
                                       "n2" = n2,
                                       "Proportion.BF12" = prop_BF12,
                                       "data_H1" = results_H1,
                                       "singularity" = data_H1$singularity)
                low <- min_sample
                previous_eta <- 0
                previous_high <- 0
                high <- max
                return(list(low = low,
                            high = high,
                            n2 = n2,
                            previous_eta = current_eta,
                            previous_high = previous_high,
                            final_SSD = final_SSD,
                            ultimate_sample_sizes = TRUE
                ))
            } else {
                # Increase the number of clusters since eta is too small
                low <- n2                         #lower bound
                high <- high                      #higher bound
                n2 <- round2((low + high) / 2)     #point in the middle
                if (!n2 %% 2 == 0) n2 <- n2 + 1 # To ensure number of clusters is even
                
                # Adjust higher bound when there is a ceiling effect
                if (low + n2 == high * 2) {
                    low <- n2                         #lower bound
                    if (previous_high > 0) {
                        high <- n2 * 2                       #higher bound
                    } else {
                        high <- max                       #higher bound
                    }
                    n2 <- round2((low + high) / 2)     #point in the middle
                }
                return(list(low = low,
                            high = high,
                            n1 = n1,
                            n2 = n2,
                            previous_eta = current_eta))
            }
        } else if (fixed == "n2") {
            if (n1 == max)    {# If the sample size reaches the maximum
                final_SSD <- list("n1" = n1,
                                       "n2" = n2,
                                       "Proportion.BF12" = prop_BF12,
                                       "data_H1" = results_H1,
                                       "singularity" = data_H1$singularity)
                low <- min_sample
                previous_eta <- 0
                previous_high <- 0
                high <- max
                return(list(low = low,
                            high = high,
                            n1 = n1,
                            previous_eta = current_eta,
                            previous_high = previous_high,
                            final_SSD = final_SSD,
                            ultimate_sample_sizes = TRUE))
            } else {
                # Increase the cluster sizes since eta is too small
                low <- n1                        #lower bound
                high <- high                     #higher bound
                n1 <- round2((low + high) / 2)    #point in the middle
                
                # Adjust higher bound when there is a ceiling effect or no increase of power
                if (low + n1 == high * 2) {
                    low <- n1                        #lower bound
                    #Set the higher bound based on the previous high or the maximum
                    if (previous_high > 0 ) {
                        high <- n1 * 2                       #higher bound
                    } else {
                        high <- max
                    }
                    n1 <- round2((low + high) / 2)    # point in the middle
                }
                return(list(low = low,
                            high = high,
                            n1 = n1,
                            previous_eta = current_eta))
            }
        }
    } else if (condition_met) {
        print(c("previous:", previous_eta))
        previous_high <- high
        final_SSD <- list("n1" = n1,
                           "n2" = n2,
                           "Proportion.BF12" = prop_BF12,
                           "data_H1" = results_H1,
                           "singularity" = data_H1$singularity)
        print(c("Using cluster size:", n1,
                "and number of clusters:", n2, "prop_BF12: ", prop_BF12,
                "low:", low, "high:", high))
        if (fixed == "n1") {
            # Finish search because eta is close enough to the desired eta or
            # there is no change in eta and the lower bound is close to the middle point
            if ((current_eta - eta < 0.1 && n2 - low == 2) ||
                (previous_eta == current_eta && n2 - low == 2)) {
                low <- min_sample
                high <- max
                return(list(
                    low = low,
                    high = max,
                    n1 = n1,
                    n2 = n2,
                    final_SSD = final_SSD,
                    ultimate_sample_sizes = TRUE))
                
            } else if (low + n2 == high * 2) {
                # Adjust higher bound when there is a ceiling effect
                low <- n2                         #lower bound
                high <- n2 * 2                       #higher bound
                n2 <- round2((low + high) / 2)     #point in the middle
                return(list(low = low,
                            high = high,
                            n1 = n1,
                            n2 = n2,
                            previous_eta = current_eta))  
            } else {
                # Decreasing to find the ultimate number of clusters
                message("Lowering number of clusters")
                low <- low                         #lower bound
                high <- n2                         #higher bound
                n2 <- round2((low + high) / 2)      #point in the middle
                if (!n2 %% 2 == 0) n2 <- n2 + 1
                
                
                return(list(
                    low = low,
                    high = high,
                    n1 = n1,
                    n2 = n2,
                    previous_high = high, #The higher bound when the power criterion was met
                    previous_eta = current_eta
                ))
            }
        } else if (fixed == "n2") {
            # Finish search because Eta is close enough to the desired eta or
            # there is no change in eta and the lower bound is close to the middle point or
            # reached the minimum number that meets the Bayesian power condition
            if ((current_eta - eta < 0.1 && n1 - low == 1) ||
                (current_eta == previous_eta && n1 - low == 1) ||
                (current_eta == previous_eta && low + n1 == high * 2)) {
                low <- min_sample
                high <- max
                
                return(list(low = min_sample,
                            high = max,
                            n1 = n1,
                            n2 = n2,
                            final_SSD = final_SSD,
                            ultimate_sample_sizes = TRUE))
                
            } else if (low + n1 == high * 2) {
                # Adjust higher bound when there is a ceiling effect
                low <- n1                         #lower bound
                high <- n1 * 2                       #higher bound
                n1 <- round2((low + high) / 2)     #point in the middle
                return(list(low = low,
                            high = high,
                            n1 = n1,
                            n2 = n2,
                            previous_eta = current_eta))
            } else {
                message("Lowering cluster size")
                # Decreasing the cluster size to find the ultimate sample size
                low <- low                         #lower bound
                high <- n1                         #higher bound
                n1 <- round2((low + high) / 2)      #point in the middle
                return(list(low = low,
                            high = high,
                            n1 = n1,
                            n2 = n2,
                            previous_high = high,
                            previous_eta = current_eta))
            }
        } # Finish if n2
    } # Finish condition met
}
