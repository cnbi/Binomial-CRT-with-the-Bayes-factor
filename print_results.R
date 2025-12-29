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
        results_matrix[, 2] <- round(as.numeric(head(unlist(lapply(object_result_b, `[[`, 1)), b_number))) #n1
        results_matrix[, 3] <- round(as.numeric(head(unlist(lapply(object_result_b, `[[`, 2)), b_number))) #n2
        results_matrix[, 4] <- round(as.numeric(head(unlist(lapply(object_result_b, `[[`, 3)), b_number)), 3) #BF_01
        results_matrix[, 5] <- round(as.numeric(head(unlist(lapply(object_result_b, `[[`, 4)), b_number)), 3) #BF_10
        colnames(results_matrix) <- c("b", "n1", "n2", paste("P(BF.01 >", list_results[[(n_object - 1)]][1], "| H0) > ", 
                                                             list_results[[n_object - 1]][2], sep = " "), 
                                      paste("P(BF.10 >", list_results[[(n_object)]][1], "| H1) > ", 
                                            list_results[[n_object]][2], sep = " "))
        
        cat("Hypotheses:", "\n")
        cat("    H1: Intervention > Control", "\n")
        cat("    H2: Intervention = Control", "\n")
        
        cat("***********************************************************************", "\n")
        print(format(results_matrix, justify = "centre"))
        cat("***********************************************************************", "\n")
        cat("n1: Cluster sizes", "\n")
        cat("n2: Number of clusters", "\n")
        }
}