# Based on:
#
# Mahé F, Rognes T, Quince C, de Vargas C, Dunthorn M. (2014) 
# Swarm: robust and fast clustering method for amplicon-based studies. 
# PeerJ 2:e593 https://doi.org/10.7717/peerj.593
#
# Supplement 8: https://doi.org/10.7717/peerj.593/supp-8
#
# and
#
# vsearch-eval
# https://github.com/torognes/vsearch-eval
# 
# File: https://github.com/torognes/vsearch-eval/blob/master/cluster/scripts/plot.R

library(ggplot2)
library(reshape2)



# Facet grid of boxplots of all three metrics for all methods (facet per method-metric pair)
plot_metrics <- function(metrics_file, plot_file, thresholds, col_names = c("method", "shuffling", "threshold", "metric", "value")) {
  
  d <- read.table(metrics_file, dec = ".", sep =",")
  colnames(d) <- col_names
  
  # Threshold values treated as class. Force numerical sorting
  d$threshold <- as.character(d$threshold)
  d$threshold <- factor(d$threshold,
                        levels=as.character(thresholds),
                        ordered=TRUE)
  
  # Subset and rename
  d2 <- subset.data.frame(d, d$metric != "NMI" & d$metric != "rand")
  
  method_levels = c("Swarm (v1.2.3)", "Swarm (v2)", "GeFaST (edit distance, w/o f.)", "GeFaST (edit distance, t + 1)", "GeFaST (edit distance, 2 * t)", "GeFaST (scoring function, w/o f.)", "GeFaST (scoring function, t + 1)", "GeFaST (scoring function, 2 * t)")
  d2$method <- gsub("swarm-old", method_levels[1], d2$method)
  d2$method <- gsub("swarm", method_levels[2], d2$method)
  d2$method <- gsub("gefast-ef1", method_levels[4], d2$method)
  d2$method <- gsub("gefast-e2f", method_levels[5], d2$method)
  d2$method <- gsub("gefast-e", method_levels[3], d2$method)
  d2$method <- gsub("gefast-sf1", method_levels[7], d2$method)
  d2$method <- gsub("gefast-s2f", method_levels[8], d2$method)
  d2$method <- gsub("gefast-s", method_levels[6], d2$method)
  d2$method <- factor(d2$method, levels = method_levels)
  
  levels(d2$metric)[1] <- "Adjusted Rand index"
  levels(d2$metric)[3] <- "Precision"
  levels(d2$metric)[4] <- "Rand index"
  levels(d2$metric)[5] <- "Recall"
  
  # Faceted plot
  plot_labels=as.character(thresholds)
  ggplot(d2, aes(x = threshold, y = value)) +
    geom_boxplot(outlier.size = 1.5) +
    scale_x_discrete(breaks = as.character(thresholds), labels=plot_labels) +
    xlab("Threshold t") +
    ylab("Metric values") +
    facet_grid(metric ~ method) + 
    labs(colour="Method", shape = "Fastidious")
  
  ggsave(file = plot_file, width=16, height=5, device = "pdf")
  
}



# Facet grid of median values of all three metrics for all methods (facet per method-metric pair)
plot_medians <- function(metrics_file, plot_file, thresholds, medians_file = NA, col_names = c("method", "shuffling", "threshold", "metric", "value")) {
  
  d <- read.table(metrics_file, dec = ".", sep =",")
  colnames(d) <- col_names
  
  # Threshold values treated as class. Force numerical sorting
  d$threshold <- as.character(d$threshold)
  d$threshold <- factor(d$threshold,
                        levels=as.character(thresholds),
                        ordered=TRUE)
  
  # Compute median values
  d4 <- melt(acast(d, metric~threshold~method, median))
  colnames(d4) <- c("metric", "threshold", "method", "value")
  d4$threshold <- as.character(d4$threshold)
  d4$threshold <- factor(d4$threshold,
                         levels = as.character(thresholds),
                         ordered=TRUE)
  d5 <- subset.data.frame(d4, d4$metric != "NMI" & d4$metric != "rand")
  
  levels(d5$metric)[1] <- "Adjusted Rand index"
  levels(d5$metric)[3] <- "Precision"
  levels(d5$metric)[4] <- "Rand index"
  levels(d5$metric)[5] <- "Recall"

  method_levels = c("Swarm (v1.2.3)", "Swarm (v2)", "GeFaST (edit distance, w/o f.)", "GeFaST (edit distance, t + 1)", "GeFaST (edit distance, 2 * t)", "GeFaST (scoring function, w/o f.)", "GeFaST (scoring function, t + 1)", "GeFaST (scoring function, 2 * t)")
  d5$method <- gsub("swarm-old", method_levels[1], d5$method)
  d5$method <- gsub("swarm", method_levels[2], d5$method)
  d5$method <- gsub("gefast-ef1", method_levels[4], d5$method)
  d5$method <- gsub("gefast-e2f", method_levels[5], d5$method)
  d5$method <- gsub("gefast-e", method_levels[3], d5$method)
  d5$method <- gsub("gefast-sf1", method_levels[7], d5$method)
  d5$method <- gsub("gefast-s2f", method_levels[8], d5$method)
  d5$method <- gsub("gefast-s", method_levels[6], d5$method)
  d5$method <- factor(d5$method, levels = method_levels)

  if (!is.na(medians_file)) {
    
    # Output median values to a file
    d7 <- acast(d5, metric ~ method ~ threshold)
    write.csv(d7, file = medians_file, quote = TRUE, eol = "\n", na = "NA")
    
  }
  
  # Faceted plot using median values
  ggplot(d5, aes(x = threshold, y = value)) +
    geom_boxplot() +
    xlab("Threshold t") +
    ylab("Metric values") +
    facet_grid(metric ~ method) + 
    labs(colour="Method", shape = "Fastidious")
  
  ggsave(file = plot_file, width=16, height=5, device = "pdf")
  
}



# Facet grid of median values of all three metrics for all methods (facet per metric, jittering)
plot_medians_jitter <- function(metrics_file, plot_file, thresholds, medians_file = NA, col_names = c("method", "shuffling", "threshold", "metric", "value")) {
  
  d <- read.table(metrics_file, dec = ".", sep =",")
  colnames(d) <- col_names
  
  # Threshold values treated as class. Force numerical sorting
  d$threshold <- as.character(d$threshold)
  d$threshold <- factor(d$threshold,
                        levels=as.character(thresholds),
                        ordered=TRUE)
  

  # Compute median values
  d4 <- melt(acast(d, metric~threshold~method, median))
  colnames(d4) <- c("metric", "threshold", "method", "value")
  d4$threshold <- as.character(d4$threshold)
  d4$threshold <- factor(d4$threshold,
                         levels = as.character(thresholds),
                         ordered=TRUE)
  d5 <- subset.data.frame(d4, d4$metric != "NMI" & d4$metric != "rand")
  
  levels(d5$metric)[1] <- "Adjusted Rand index"
  levels(d5$metric)[3] <- "Precision"
  levels(d5$metric)[4] <- "Rand index"
  levels(d5$metric)[5] <- "Recall"

  if (!is.na(medians_file)) {
    
    # Output median values to a file
    d7 <- acast(d5, metric ~ method ~ threshold)
    write.csv(d7, file = medians_file, quote = TRUE, eol = "\n", na = "NA")
    
  }
  
  # Faceted plot using median values
  fastidious <- rep("deactivated", length(d5$method))
  fastidious[grepl("f1", d5$method)] <- "t + 1" 
  fastidious[grepl("2f", d5$method)] <- "2 * t"
  fastidious <- factor(fastidious, levels = c("deactivated", "t + 1", "2 * t"))
  
  method_levels = c("Swarm (v1.2.3)", "Swarm (v2)", "GeFaST (edit distance)", "GeFaST (scoring function)")
  d5$method <- gsub("2f", "", gsub("f1", "", d5$method))
  d5$method <- gsub("swarm-old", method_levels[1], d5$method)
  d5$method <- gsub("swarm", method_levels[2], d5$method)
  d5$method <- gsub("gefast-e", method_levels[3], d5$method)
  d5$method <- gsub("gefast-s", method_levels[4], d5$method)
  d5$method <- factor(d5$method, levels = method_levels)
  
  d5 <- cbind(d5, fastidious)
  
  ggplot(d5, aes(x = threshold, y = value)) + scale_shape_manual(values=c(1, 0, 5)) +
    geom_jitter(aes(colour = method, shape = fastidious), width = 0.25) +
    xlab("Threshold t") +
    ylab("Metric values") + 
    facet_grid(. ~ metric) + 
    labs(colour="Method", shape = "Fastidious") + 
    theme(legend.position = c(.5, .25), legend.direction = "vertical", legend.box = "horizontal")
  
  ggsave(file = plot_file, width=16, height=5, device = "pdf")
  
}



# Facet grid of median values of all three metrics for all methods (facet per metric, line plots)
plot_medians_line <- function(metrics_file, plot_file, thresholds, medians_file = NA, col_names = c("method", "shuffling", "threshold", "metric", "value")) {
  
  d <- read.table(metrics_file, dec = ".", sep =",")
  colnames(d) <- col_names
  
  # Threshold values treated as class. Force numerical sorting
  d$threshold <- as.character(d$threshold)
  d$threshold <- factor(d$threshold,
                        levels=as.character(thresholds),
                        ordered=TRUE)
  
  
  # Compute median values
  d4 <- melt(acast(d, metric~threshold~method, median))
  colnames(d4) <- c("metric", "threshold", "method", "value")
  d4$threshold <- as.character(d4$threshold)
  d4$threshold <- factor(d4$threshold,
                         levels = as.character(thresholds),
                         ordered=TRUE)
  d5 <- subset.data.frame(d4, d4$metric != "NMI" & d4$metric != "rand")
  
  levels(d5$metric)[1] <- "Adjusted Rand index"
  levels(d5$metric)[3] <- "Precision"
  levels(d5$metric)[4] <- "Rand index"
  levels(d5$metric)[5] <- "Recall"
  
  if (!is.na(medians_file)) {
    
    # Output median values to a file
    d7 <- acast(d5, metric ~ method ~ threshold)
    write.csv(d7, file = medians_file, quote = TRUE, eol = "\n", na = "NA")
    
  }
  
  # Faceted plot using median values
  fastidious <- rep("deactivated", length(d5$method))
  fastidious[grepl("f1", d5$method)] <- "t + 1" 
  fastidious[grepl("2f", d5$method)] <- "2 * t"
  fastidious <- factor(fastidious, levels = c("deactivated", "t + 1", "2 * t"))
  
  method_levels = c("Swarm (v1.2.3)", "Swarm (v2)", "GeFaST (edit distance)", "GeFaST (scoring function)")
  d5$method <- gsub("2f", "", gsub("f1", "", d5$method))
  d5$method <- gsub("swarm-old", method_levels[1], d5$method)
  d5$method <- gsub("swarm", method_levels[2], d5$method)
  d5$method <- gsub("gefast-e", method_levels[3], d5$method)
  d5$method <- gsub("gefast-s", method_levels[4], d5$method)
  d5$method <- factor(d5$method, levels = method_levels)
  
  d5 <- cbind(d5, fastidious)
  
  ggplot(d5) +  
    aes(x = threshold, y = value, colour = method, shape = fastidious, group=interaction(method, fastidious)) + scale_shape_manual(values=c(1, 0, 5)) + 
    geom_point() + geom_line() +
    xlab("Threshold t") +
    ylab("Metric values") + 
    facet_grid(. ~ metric) + 
    labs(colour="Method", shape = "Fastidious") + 
    theme(legend.position = c(.5, .25), legend.direction = "vertical", legend.box = "horizontal")
  
  ggsave(file = plot_file, width=16, height=5, device = "pdf")
  
}



# Facet grid of metric value of all three metrics for all methods (facet per metric, jittering)
plot_single_metric_jitter <- function(metrics_file, plot_file, thresholds, metrics_out_file = NA, col_names = c("method", "shuffling", "threshold", "metric", "value")) {
  
  d <- read.table(metrics_file, dec = ".", sep =",")
  colnames(d) <- col_names
  
  # Threshold values treated as class. Force numerical sorting
  d$threshold <- as.character(d$threshold)
  d$threshold <- factor(d$threshold,
                        levels=as.character(thresholds),
                        ordered=TRUE)
  
  d5 <- subset.data.frame(d[,c("metric", "threshold", "method", "value")], d$metric != "NMI" & d$metric != "rand")
  
  levels(d5$metric)[1] <- "Adjusted Rand index"
  levels(d5$metric)[3] <- "Precision"
  levels(d5$metric)[4] <- "Rand index"
  levels(d5$metric)[5] <- "Recall"
  
  if (!is.na(metrics_out_file)) {
    
    # Output median values to a file
    d7 <- acast(d5, metric ~ method ~ threshold)
    write.csv(d7, file = metrics_out_file, quote = TRUE, eol = "\n", na = "NA")
    
  } 
  
  # Faceted plot using median values
  fastidious <- rep("deactivated", length(d5$method))
  fastidious[grepl("f1", d5$method)] <- "t + 1" 
  fastidious[grepl("2f", d5$method)] <- "2 * t"
  fastidious <- factor(fastidious, levels = c("deactivated", "t + 1", "2 * t"))
  
  method_levels = c("Swarm (v1.2.3)", "Swarm (v2)", "GeFaST (edit distance)", "GeFaST (scoring function)")
  d5$method <- gsub("2f", "", gsub("f1", "", d5$method))
  d5$method <- gsub("swarm-old", method_levels[1], d5$method)
  d5$method <- gsub("swarm", method_levels[2], d5$method)
  d5$method <- gsub("gefast-e", method_levels[3], d5$method)
  d5$method <- gsub("gefast-s", method_levels[4], d5$method)
  d5$method <- factor(d5$method, levels = method_levels)
  
  d5 <- cbind(d5, fastidious)
  
  ggplot(d5, aes(x = threshold, y = value)) +
    geom_jitter(aes(colour = method, shape = fastidious), width = 0.25) + scale_shape_manual(values=c(1, 0, 5)) +
    xlab("Threshold t") +
    ylab("Metric values") + 
    facet_grid(. ~ metric) + 
    labs(colour="Method", shape = "Fastidious") + 
    theme(legend.position = c(.5, .25), legend.direction = "vertical", legend.box = "horizontal")
  
  ggsave(file = plot_file, width=16, height=5, device = "pdf")
  
}



# Facet grid of metric value of all three metrics for all methods (facet per metric, line plot)
plot_single_metric_line <- function(metrics_file, plot_file, thresholds, metrics_out_file = NA, col_names = c("method", "shuffling", "threshold", "metric", "value")) {
  
  d <- read.table(metrics_file, dec = ".", sep =",")
  colnames(d) <- col_names
  
  # Threshold values treated as class. Force numerical sorting
  d$threshold <- as.character(d$threshold)
  d$threshold <- factor(d$threshold,
                        levels=as.character(thresholds),
                        ordered=TRUE)
  
  d5 <- subset.data.frame(d[,c("metric", "threshold", "method", "value")], d$metric != "NMI" & d$metric != "rand")
  
  levels(d5$metric)[1] <- "Adjusted Rand index"
  levels(d5$metric)[3] <- "Precision"
  levels(d5$metric)[4] <- "Rand index"
  levels(d5$metric)[5] <- "Recall"
  
  if (!is.na(metrics_out_file)) {
    
    # Output median values to a file
    d7 <- acast(d5, metric ~ method ~ threshold)
    write.csv(d7, file = metrics_out_file, quote = TRUE, eol = "\n", na = "NA")
    
  }
  
  # Faceted plot using median values
  fastidious <- rep("deactivated", length(d5$method))
  fastidious[grepl("f1", d5$method)] <- "t + 1" 
  fastidious[grepl("2f", d5$method)] <- "2 * t"
  fastidious <- factor(fastidious, levels = c("deactivated", "t + 1", "2 * t"))
  
  method_levels = c("Swarm (v1.2.3)", "Swarm (v2)", "GeFaST (edit distance)", "GeFaST (scoring function)")
  d5$method <- gsub("2f", "", gsub("f1", "", d5$method))
  d5$method <- gsub("swarm-old", method_levels[1], d5$method)
  d5$method <- gsub("swarm", method_levels[2], d5$method)
  d5$method <- gsub("gefast-e", method_levels[3], d5$method)
  d5$method <- gsub("gefast-s", method_levels[4], d5$method)
  d5$method <- factor(d5$method, levels = method_levels) 
  
  d5 <- cbind(d5, fastidious)
   
  ggplot(d5) +  
    aes(x = threshold, y = value, colour = method, shape = fastidious, group=interaction(method, fastidious)) + scale_shape_manual(values=c(1, 0, 5)) +
    geom_point() + geom_line() +
    xlab("Threshold t") +
    ylab("Metric values") + 
    facet_grid(. ~ metric) + 
    labs(colour="Method", shape = "Fastidious") + 
    theme(legend.position = c(.5, .25), legend.direction = "vertical", legend.box = "horizontal")
  
  ggsave(file = plot_file, width=16, height=5, device = "pdf")
  
}
