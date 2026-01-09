########################## SMALL FUNCTIONS ###########################

# Model fitting
fit_glmer <- function(x) {
    suppressMessages({
        fitted_model <- glmer(resp ~ intervention + control - 1 + (1 | id), data = x, 
                              family = binomial)})
    return(fitted_model)
}

# Obtain the variance covariance matrix for fixed effects
varcov <- function(output.glmer, name) {
    varcov <- matrix(vcov(output.glmer)[name, name], nrow = 1, ncol = 1)
    return(varcov)
}

# Obtain variance for random intercept
get_variance <- function(output.glmer) {
    variance <- as.data.frame(VarCorr(output.glmer))
    value <- variance[, 4]
    return(value)
}

# Marking singular matrices
marker_func <- function(output.glmer) {
    ifelse(isSingular(output.glmer), marker <- 1, marker <- 0)
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
                low <- min_sample
                previous_eta <- 0
                previous_high <- 0
                high <- max
                return(list(low = low,
                            high = high,
                            n1 = n1,
                            n2 = n2,
                            previous_eta = current_eta,
                            previous_high = previous_high,
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
                low <- min_sample
                previous_eta <- 0
                previous_high <- 0
                high <- max
                return(list(low = low,
                            high = high,
                            n1 = n1,
                            n2 = n2,
                            previous_eta = current_eta,
                            previous_high = previous_high,
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
                        high <- previous_high
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
            if ((current_eta - eta < 0.1 && n2 - low == 2) ||
                (previous_eta == current_eta && n2 - low == 2)) {
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
                    final_SSD = final_SSD))
                
            } else {
                # Decreasing to find the ultimate number of clusters
                message("Lowering number of clusters")
                low <- low                         #lower bound
                high <- n2                         #higher bound
                n2 <- round2((low + high) / 2)      #point in the middle
                if (!n2 %% 2 == 0) n2 <- n2 + 1
                if (n2 < 30) warning("The number of groups is less than 30.
                                             This may cause problems in convergence and singularity.")
                
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
            if ((current_eta - eta < 0.1 && n1 - low == 1) ||
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
                            final_SSD = final_SSD))
                
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