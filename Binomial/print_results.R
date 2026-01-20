########################### PRINT RESULT ###################################

print_results <- function(list_results, type_hyp, b) {
    # Title
    title <- "Final sample size"
    cat(paste("\n", title, "\n", sep = ""))
    row <- paste(rep("=", nchar(title)), collapse = "")
    cat(row, "\n")
    
    if (type_hyp == "eq") {
        n_object <- length(list_results)
        results_matrix <- matrix(NA, nrow = b, ncol = 5)
        results_matrix[, 1] <- seq(b)
        
        object_result_b <- list_results[1:b]
        results_matrix[, 2] <- sapply(object_result_b, `[[`, "n1") #n1
        results_matrix[, 3] <- sapply(object_result_b, `[[`, "n2") #n2
        results_matrix[, 4] <- sapply(object_result_b, `[[`, "Proportion.BF01") #BF_01
        results_matrix[, 5] <- sapply(object_result_b, `[[`, "Proportion.BF10") #BF_10
        colnames(results_matrix) <- c("b", "n1", "n2", paste0("P(BF.01 >", list_results[[(n_object - 1)]][1], "| H0) > ", 
                                                             list_results[[n_object]][1]), 
                                      paste0("P(BF.10 >", list_results[[(n_object - 1)]][2], "| H1) > ", 
                                            list_results[[n_object]][2]))
        
        cat("Hypotheses:", "\n")
        cat("    H1: Intervention > Control", "\n")
        cat("    H2: Intervention = Control", "\n")
        
        stargazer::stargazer(results_matrix, type = "text", summary = FALSE)
        cat("n1: Cluster sizes", "\n")
        cat("n2: Number of clusters", "\n")
        }
}