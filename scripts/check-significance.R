# Compute means of replicates.
summarise <- function(data, method) {
  
  splitted <- lapply(split(data[, names(data) %in% c("recall", "precision", "adjustedrand")], data$t), colMeans)
  
  t <- unique(data$t)
  method <- rep(method, length(t))
  
  res <- data.frame(do.call("rbind", splitted))
  res <- cbind(method, t, res)
  
}


# Perform a paired t-test.
compare_diff <- function(data, method1, method2, thresholds, metric, alt) {

  data_method1 <- data[(data$method == method1) & (data$t %in% thresholds),]
  data_method2 <- data[(data$method == method2) & (data$t %in% thresholds),]
  
  summary_method1 <- summarise(data_method1, method1)
  summary_method2 <- summarise(data_method2, method2)
  
  tryCatch({

    t.test(summary_method1[, metric], summary_method2[, metric], var.equal = F, paired = T, alternative = alt)
    
  }, error = function(e) {

    n <- c("statistic", "parameter", "p.value", "conf.int", "estimate", "null.value", "alternative", "method", "data.name")
    vals <- rep("n/a", length(n))
    names(vals) <- n
    as.list(vals)
    
  })
  
}


# Compare metric values of two methods over the specified thresholds.
fill_row <- function(data, method1, method2, thresholds, metrics = c("precision", "recall", "adjustedrand"), alts = c("greater", "less", "less")) {
  
  mat <- matrix(nrow = 1, ncol = length(metrics))
  colnames(mat) <- metrics
  rownames(mat) <- "p"
  
  for (j in 1:length(metrics)) {
    mat[1,j] <- compare_diff(data, method1, method2, thresholds, metrics[j], alts[j])$p.value
  }
  
  mat
  
}


# Run a set of paired t-tests using the data from metrics_file. 
# The i-th t-test compares the metric values between the i-th method 
# in method1 and the i-th in method2 over the specified thresholds.
# The kind of t-test is specified in alternatives.
#
# Returns a matrix of p-values with one row per t-test and one column 
# per metric.
run_test <- function(metrics_file, thresholds, method1, method2, alternatives) {

  data <- read.csv(metrics_file, header = F, sep = ",")
  names(data) <- c("method", "shuffling", "t", "recall", "precision", "nmi", "rand", "adjustedrand")
  
  mat <- matrix(nrow = length(method1), ncol = 3)
  colnames(mat) <- c("precision", "recall", "adjustedrand")
  rownames(mat) <- paste(method1, method2, sep = " vs. ")
  
  for (i in 1:length(method1)) {
    mat[i,] <- fill_row(data, method1[i], method2[i], thresholds, alts = alternatives)
  }
  
  cat(">>> p-values of t-tests:\n", 
      "[alternatives (per column): ", paste(alternatives, collapse = ", "), "]\n",
      "[thresholds considered: ", paste(thresholds, collapse = ", "), "]\n", 
      sep = ""
  )
  
  mat
  
}

# === Tests ===

# non-fastidious vs. fastidious and fastidious (t + 1) vs. fastidious (2 * t)
run_test("uneven.metrics", 1:6, 
         c("gefast-e", "gefast-e", "gefast-ef1", "gefast-s", "gefast-s", "gefast-sf1"), 
         c("gefast-ef1", "gefast-e2f", "gefast-e2f", "gefast-sf1", "gefast-s2f", "gefast-s2f"), 
         c("greater", "less", "less")
)

# edit-distance mode vs. scoring-function mode
run_test("uneven.metrics", 1:10, 
         c("gefast-e", "gefast-ef1", "gefast-e2f"), 
         c("gefast-s", "gefast-sf1", "gefast-s2f"), 
         rep("two.sided", 3)
)