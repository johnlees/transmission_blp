library(vegan)
library(cowplot)
library(EasyABC)
library(readxl)
library(ggplot2)

setwd("~/Documents/EBI/Papers/transmission_blp/")

run_sim <- function(inoculum, Ne, n_generations, n_repeats) {
  n_mutants <- length(inoculum)
  barcodes <- seq_len(length(inoculum))

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
    mean_h1 = mean(hills_n1[n_generations,]),
    lower_h1 = quantile(hills_n1[n_generations,], probs = 0.025),
    upper_h1 = quantile(hills_n1[n_generations,], probs = 0.975),
    mean_maxa = mean(max_abun[n_generations,]),
    lower_maxa = quantile(max_abun[n_generations,], probs = 0.025),
    upper_maxa = quantile(max_abun[n_generations,], probs = 0.975)
  )
  df
}

time <- 24 #hours
generation_time <- 3 # hours
n_repeats <- 300
n_generations <- round(time / generation_time)

Barcode_Frequencies <- read_excel("Barcode_Frequencies_07062022.xlsx", 
                                  col_names = FALSE, skip = 2)
inoculum <- Barcode_Frequencies$...2

wt_h1 <- apply(Barcode_Frequencies[,c(3, 4, 5, 6, 7)], 2, renyi, scales = 1, hill = TRUE)
wt_max <- apply(Barcode_Frequencies[,c(3, 4, 5, 6, 7)], 2, max) / 100

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
  df <- rbind(df, run_sim(inoculum, Ne, n_generations, n_repeats))
}

p1 <- ggplot(df) +
  geom_point(aes(x = Ne, y = mean_h1), stroke=1, shape=1) +
  geom_errorbar(aes(x = Ne, ymin = lower_h1, ymax = upper_h1)) +
  theme_bw() +
  theme(axis.title = element_text(family = "Arial", size = 12),
        axis.text = element_text(family = "Arial", size = 11)) +
  scale_x_discrete(limits = c("WT", "67", "133", "266", "532", "1064")) +
  xlab("Effective population size (Ne)") +
  ylab("Hills N1")

ggsave(p1, filename = "fig2b1.pdf", device = cairo_pdf, 
       width = 6, height = 6, units = "in")

p2 <- ggplot(df) +
  geom_point(aes(x = Ne, y = mean_maxa), stroke=1, shape=1) +
  geom_errorbar(aes(x = Ne, ymin = lower_maxa, ymax = upper_maxa)) +
  theme_bw() +
  theme(axis.title = element_text(family = "Arial", size = 12),
        axis.text = element_text(family = "Arial", size = 11)) +
  scale_x_discrete(limits = c("WT", "67", "133", "266", "532", "1064")) +
  xlab("Effective population size (Ne)") +
  ylab("Abundance of most dominant clone")

ggsave(p2, filename = "fig2b2.pdf", device = cairo_pdf, 
       width = 6, height = 6, units = "in")

# ABC
n <- 20
tolerance <- c(1.25,0.75)
model <- function(Ne) {
  Ne <- round(Ne)
  if(Ne < 1) {
    sum <- c(NA, NA)
  } else {
    df_ret <- run_sim(inoculum, Ne, n_generations, 1)
    sum <- c(df_ret$mean_h1, df_ret$mean_maxa)
  }
  sum
}
prior <- list(c("unif", 1, 2000))
sum_stat_obs <- c(mean(wt_h1), mean(wt_max))

ABC_wt <- ABC_sequential(method="Beaumont",
                         model=model,
                         prior=prior,
                         prior_test="X1 > 1",
                         progress_bar = TRUE,
                         nb_simul = n,
                         summary_stat_target=sum_stat_obs,
                         tolerance_tab=tolerance)
posterior <- sample(ABC_wt$param, size=100, replace=TRUE, prob = ABC_wt$weights)
mean(posterior)
quantile(posterior, c(0.025, 0.975))
# Ne: 33 (6-76 95% CrI)

BlpC_BarcodeFrequencies <- read_excel("BlpC_BarcodeFrequencies.xlsx",
                                      col_names = FALSE, skip = 2)
blpc_h1 <- apply(BlpC_BarcodeFrequencies[,c(3, 4, 5, 6, 7, 8)], 2, renyi, scales = 1, hill = TRUE)
blpc_max <- apply(BlpC_BarcodeFrequencies[,c(3, 4, 5, 6, 7, 8)], 2, max) / 100

df <- data.frame(
              Ne = "blpC",
              mean_h1 = blpc_h1,
              lower_h1 = NA,
              upper_h1 = NA,
              mean_maxa = blpc_max,
              lower_maxa = NA,
              upper_maxa = NA)


for (Ne in c(67, 133, 266, 532, 1064, 10000)) {
  df <- rbind(df, run_sim(inoculum, Ne, n_generations, n_repeats))
}

p3 <- ggplot(df) +
  geom_point(aes(x = Ne, y = mean_h1), stroke=1, shape=1) +
  geom_errorbar(aes(x = Ne, ymin = lower_h1, ymax = upper_h1)) +
  theme_bw() +
  theme(axis.title = element_text(family = "Arial", size = 12),
        axis.text = element_text(family = "Arial", size = 11)) +
  scale_x_discrete(limits = c("blpC", "67", "133", "266", "532", "1064", "10000"),
                   labels = c("blpC" = "\u0394blpC")) +
  xlab("Effective population size (Ne)") +
  ylab("Hills N1")

ggsave(p3, filename = "fig4h1.pdf", device = cairo_pdf, 
       width = 6, height = 6, units = "in")

p4 <- ggplot(df) +
  geom_point(aes(x = Ne, y = mean_maxa), stroke=1, shape=1) +
  geom_errorbar(aes(x = Ne, ymin = lower_maxa, ymax = upper_maxa)) +
  theme_bw() +
  theme(axis.title = element_text(family = "Arial", size = 12),
        axis.text = element_text(family = "Arial", size = 11)) +
  scale_x_discrete(limits = c("blpC", "67", "133", "266", "532", "1064", "10000"),
                   labels = c("blpC" = "\u0394blpC")) +
  xlab("Effective population size (Ne)") +
  ylab("Abundance of most dominant clone")

ggsave(p4, filename = "fig4h2.pdf", device = cairo_pdf, 
       width = 6, height = 6, units = "in")

# Combined plot
df <- rbind(df, 
            data.frame(
              Ne = "WT",
              mean_h1 = wt_h1,
              lower_h1 = NA,
              upper_h1 = NA,
              mean_maxa = wt_max,
              lower_maxa = NA,
              upper_maxa = NA
            ))

p5 <- ggplot(df) +
  geom_point(aes(x = Ne, y = mean_h1), stroke=1, shape=1) +
  geom_errorbar(aes(x = Ne, ymin = lower_h1, ymax = upper_h1)) +
  theme_bw() +
  theme(axis.title = element_text(family = "Arial", size = 12),
        axis.text = element_text(family = "Arial", size = 11)) +
  scale_x_discrete(limits = c("WT", "blpC", "67", "133", "266", "532", "1064", "10000"),
                   labels = c("blpC" = "\u0394blpC")) +
  xlab("Effective population size (Ne)") +
  ylab("Hills N1")

ggsave(p5, filename = "combined2_4_H1.pdf", device = cairo_pdf, 
       width = 6, height = 6, units = "in")

p6 <- ggplot(df) +
  geom_point(aes(x = Ne, y = mean_maxa), stroke=1, shape=1) +
  geom_errorbar(aes(x = Ne, ymin = lower_maxa, ymax = upper_maxa)) +
  theme_bw() +
  theme(axis.title = element_text(family = "Arial", size = 12),
        axis.text = element_text(family = "Arial", size = 11)) +
  scale_x_discrete(limits = c("WT", "blpC", "67", "133", "266", "532", "1064", "10000"),
                   labels = c("blpC" = "\u0394blpC")) +
  xlab("Effective population size (Ne)") +
  ylab("Abundance of most dominant clone")

ggsave(p6, filename = "combined2_4_abun.pdf", device = cairo_pdf, 
       width = 6, height = 6, units = "in")

n <- 20
tolerance <- c(1.25,0.75)
model <- function(Ne) {
  Ne <- round(Ne)
  if(Ne < 1) {
    sum <- c(NA, NA)
  } else {
    df_ret <- run_sim(inoculum, Ne, n_generations, 1)
    sum <- c(df_ret$mean_h1, df_ret$mean_maxa)
  }
  sum
}
prior <- list(c("unif", 1, 2000))
sum_stat_obs <- c(mean(blpc_h1), mean(blpc_max))

ABC_blpc <- ABC_sequential(method="Beaumont",
                         model=model,
                         prior=prior,
                         prior_test="X1 > 1",
                         nb_simul = n,
                         summary_stat_target=sum_stat_obs,
                         tolerance_tab=tolerance)
posterior <- sample(ABC_blpc$param, size=100, replace=TRUE, prob = ABC_blpc$weights)
mean(posterior)
quantile(posterior, c(0.025, 0.975))
# Ne = 792 (620-904 95% CrI)

# As above, but with a bottleneck
run_sim_trans <- function(inoculum, Ne, bottleneck_size, bottleneck_gen,
                          n_generations, n_repeats) {
  n_mutants <- length(inoculum)
  barcodes <- seq_len(length(inoculum))
  
  population <- array(data = 0, dim = c(n_mutants, n_generations, n_repeats))
  
  hills_n1 <- array(dim = c(n_generations, n_repeats))
  n_clones <- array(dim = c(n_generations, n_repeats))
  
  # Set the initial condition at the first step
  for (repeat_idx in 1:n_repeats) {
    initial_sample <- table(sample(barcodes, Ne, replace = TRUE, prob = inoculum))
    population[as.numeric(names(initial_sample)), 1, repeat_idx] <- initial_sample
    
    hills_n1[1, repeat_idx] <- renyi(population[, 1, repeat_idx], scales = 1, hill = TRUE)
    n_clones[1, repeat_idx] <- sum(population[, 1, repeat_idx] > 0)
  }
  
  # Run forward
  for (repeat_idx in 1:n_repeats) {
    for (generation in 2:n_generations) {
      if (generation == bottleneck_gen) {
        bottleneck_pop <- rmultinom(1, bottleneck_size, prob = population[, generation - 1, repeat_idx] / Ne)
        population[, generation, repeat_idx] <-
          rmultinom(1, Ne, prob = bottleneck_pop)
      } else {
        population[, generation, repeat_idx] <-
          rmultinom(1, Ne, prob = population[, generation - 1, repeat_idx] / Ne)    
      }

      
      # Could do this at the end with apply etc, but looping anyway
      hills_n1[generation, repeat_idx] <- renyi(population[, generation, repeat_idx], scales = 1, hill = TRUE)
      n_clones[generation, repeat_idx] <- sum(population[, generation, repeat_idx] > 0)
    }
  }
  
  df <- data.frame(
    Ne = as.character(Ne),
    mean_h1 = mean(hills_n1[n_generations,]),
    lower_h1 = quantile(hills_n1[n_generations,], probs = 0.025),
    upper_h1 = quantile(hills_n1[n_generations,], probs = 0.975),
    mean_clones = mean(n_clones[n_generations,]),
    lower_clones = quantile(n_clones[n_generations,], probs = 0.025),
    upper_clones = quantile(n_clones[n_generations,], probs = 0.975)
  )
  df
}

time <- 240 #hours
generation_time <- 3 # hours
n_generations <- round(time / generation_time)
Ne <- 33
Barcode_Frequencies <- read_excel("Barcode_Frequencies_07062022.xlsx", 
                                  col_names = FALSE, skip = 2)
inoculum <- Barcode_Frequencies$...2

n <- 20
tolerance <- c(1.25,0.75)
model <- function(bottleneck) {
  bottleneck_size <- round(bottleneck)
  bottleneck_gen <- round(n_generations / 2)
  
  df_ret <- run_sim_trans(inoculum, Ne, bottleneck_size, bottleneck_gen,
                                      n_generations, 1)
  sum <- c(df_ret$mean_h1, df_ret$mean_clones)
  sum
}
prior <- list(c("unif", 1, Ne))
sum_stat_obs <- c(mean(c(1.118,
                         2.43,
                         3.413,
                         1.136)),
                  mean(c(13.774,
                         30.621,
                         31.543,
                         31.826)))

ABC_trans <- ABC_sequential(method="Beaumont",
                           model=model,
                           prior=prior,
                           nb_simul = n,
                           summary_stat_target=sum_stat_obs,
                           tolerance_tab=tolerance)
posterior <- sample(ABC_trans$param, size=100, replace=TRUE, prob = ABC_trans$weights)
mean(posterior)
quantile(posterior, c(0.025, 0.975))
