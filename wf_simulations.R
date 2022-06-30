library(vegan)
library(cowplot)

setwd("~/Documents/EBI/Papers/transmission_blp/")

Barcode_Frequencies <- read_excel("Barcode_Frequencies.xlsx", 
                                  col_names = FALSE, skip = 2)
inoculum <- Barcode_Frequencies$...2
n_mutants <- length(inoculum)
barcodes <- seq_len(length(inoculum))

run_sim <- function(Ne, n_generations, n_repeats) {

  population <- array(data = 0, dim = c(n_mutants, n_generations, n_repeats))
  
  hills_n1 <- array(dim = c(n_generations, n_repeats))
  max_abun <- array(dim = c(n_generations, n_repeats))
  
  # Set the initial condition at the first step
  for (repeat_idx in 1:n_repeats) {
    initial_sample <- table(sample(barcodes, Ne, replace = TRUE, prob = inoculum))
    population[as.numeric(names(initial_sample)), 1, repeat_idx] <- initial_sample
    
    hills_n1[1, repeat_idx] <- renyi(population[, 1, repeat_idx], scales = 1, hill = TRUE)
    max_abun[1, repeat_idx] <- max(population[, 1, repeat_idx]) / Ne
  }
  
  # Run forward
  for (repeat_idx in 1:n_repeats) {
    for (generation in 2:n_generations) {
      population[, generation, repeat_idx] <-
        rmultinom(1, Ne, prob = population[, generation - 1, repeat_idx] / Ne)
      
      # Could do this at the end with apply etc, but looping anyway
      hills_n1[generation, repeat_idx] <- renyi(population[, generation, repeat_idx], scales = 1, hill = TRUE)
      max_abun[generation, repeat_idx] <- max(population[, generation, repeat_idx]) / Ne
    }
  }
  
  df <- data.frame(
    Ne = as.character(Ne),
    mean_h1 = mean(hills_n1[,n_generations]),
    lower_h1 = quantile(hills_n1[,n_generations], probs = 0.025),
    upper_h1 = quantile(hills_n1[,n_generations], probs = 0.975),
    mean_maxa = mean(max_abun[,n_generations]),
    lower_maxa = quantile(max_abun[,n_generations], probs = 0.025),
    upper_maxa = quantile(max_abun[,n_generations], probs = 0.975)
  )
  df
}

time <- 24 #hours
generation_time <- 1.1 # hours https://pubmed.ncbi.nlm.nih.gov/3699893/
# NB https://journals.asm.org/doi/full/10.1128/IAI.00527-13
# has generation time 161/60 = 2.68
n_repeats <- 100
n_generations <- round(time / generation_time)

wt_h1 <- apply(Barcode_Frequencies[,c(3, 4, 5)], 2, renyi, scales = 1, hill = TRUE)
wt_max <- apply(Barcode_Frequencies[,c(3, 4, 5)], 2, max) / 100

df <- data.frame(
  Ne = "WT",
  mean_h1 = wt_h1,
  lower_h1 = NA,
  upper_h1 = NA,
  mean_maxa = wt_max,
  lower_maxa = NA,
  upper_maxa = NA
)

for (Ne in c(67, 133, 266, 532, 1064)) {
  df <- rbind(df, run_sim(Ne, n_generations, n_repeats))
}

p1 <- ggplot(df) +
  geom_point(aes(x = Ne, y = mean_h1), stroke=1, shape=1) +
  geom_errorbar(aes(x = Ne, ymin = lower_h1, ymax = upper_h1)) +
  theme_bw() +
  scale_x_discrete(limits = c("WT", "67", "133", "266", "532", "1064")) +
  xlab("Effective population size (Ne)") +
  ylab("Hills N1")

p2 <- ggplot(df) +
  geom_point(aes(x = Ne, y = mean_maxa), stroke=1, shape=1) +
  geom_errorbar(aes(x = Ne, ymin = lower_maxa, ymax = upper_maxa)) +
  theme_bw() +
  scale_x_discrete(limits = c("WT", "67", "133", "266", "532", "1064")) +
  xlab("Effective population size (Ne)") +
  ylab("Abundance of most dominant clone")

pdf("fig2b.pdf", width = 12, height = 6)
plot_grid(p1, p2, labels = c('A', 'B'), ncol = 2, label_size = 12)
dev.off()

