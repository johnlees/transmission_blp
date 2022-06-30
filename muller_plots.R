library(ggplot2)
library(readxl)
library(reshape2)
library(cowplot)
library(ggmuller)
library(EvoFreq)

setwd("~/Documents/EBI/Papers/transmission_blp/")

Barcode_Frequencies <- read_excel("Barcode_Frequencies.xlsx", 
                                  col_names = FALSE, skip = 2)

# Two mixed
inoculum <- sample(Barcode_Frequencies$...2, prob = Barcode_Frequencies$...2, size = 20)
inoculum <- 100 * inoculum/sum(inoculum)
ex1 <- cbind(inoculum,
             tail(sort(Barcode_Frequencies$...3), n = 20),
             tail(sort(Barcode_Frequencies$...6), n = 20),
             tail(sort(Barcode_Frequencies$...13), n = 20),
             deparse.level = 0)

# Single dominant
ex2 <- cbind(inoculum[1:10],
             tail(sort(Barcode_Frequencies$...5), n = 10),
             tail(sort(Barcode_Frequencies$...6), n = 10),
             tail(sort(Barcode_Frequencies$...12), n = 10),
             deparse.level = 0)

# Apply average colonisation at each time point
Colonization_Levels <- read_excel("Colonization Levels.xlsx")

colonisation <- c(mean(Colonization_Levels$`2 hours`, na.rm = TRUE),
                  mean(Colonization_Levels$`1 day`, na.rm = TRUE),
                  mean(Colonization_Levels$`2 weeks`, na.rm = TRUE),
                  mean(Colonization_Levels$`4 weeks`, na.rm = TRUE))

abundance_ex1 <- sweep(ex1, 2, colonisation, `*`) / 100
abundance_ex2 <- sweep(ex2, 2, colonisation, `*`) / 100

m_plot <- function(abundance, method = "ggmuller") {
  adjacency_df <- data.frame(Parent=1, Identity=seq(2, 21))
  populations_df <- reshape2::melt(abundance)
  populations_df$Var2[populations_df$Var2 == 1] <- 0
  populations_df$Var2[populations_df$Var2 == 2] <- 1
  populations_df$Var2[populations_df$Var2 == 3] <- 14
  populations_df$Var2[populations_df$Var2 == 4] <- 28
  populations_df$Var1 <- populations_df$Var1 + 1
  colnames(populations_df) <- c("Identity", "Time", "Population")
  
  dummy_df <- data.frame(Identity=1, Time=c(0, 1, 14, 28), Population=0)
  populations_df <- rbind(populations_df, dummy_df)
  
  if (method == "ggmuller") {
    Muller_df <- get_Muller_df(adjacency_df, populations_df)
    Muller_df <- Muller_df[Muller_df$Identity != 1,]
    
    Muller_df_pop <- add_empty_pop(Muller_df)
    id_list <- sort(unique(Muller_df_pop$Identity)) # list of legend entries, omitting NA
    p <- ggplot(Muller_df_pop, aes_string(x = "Time", y = "Population", group = "Group_id", fill = "Identity", colour = "Identity")) + 
      geom_area() +
      guides(linetype = FALSE, color = FALSE) + 
      scale_fill_manual(name = "Identity", values = my_palette, breaks = id_list, na.value="transparent") +
      scale_color_manual(values = my_palette, na.value="transparent") +
      theme_classic() +
      theme(legend.position = "none") +
      xlab("Time (days)") +
      ylab("Population (CFUs)")
  } else if (method == "evofreq") {
    wide_df <- long_to_wide_freq_ready(
      populations_df,
      time_col_name = "Time",
      clone_col_name = "Identity",
      parent_col_name = "Parent",
      size_col_name = "Population",
      edges_df = adjacency_df)
    clones <- wide_df$clones
    parents <- wide_df$parents
    size_df <- wide_df$wide_size_df
    #rgb_clone_colors <- sapply(seq(1, length(clones)), function(x){paste(sample(0:255,size=3,replace=TRUE),collapse=",")})
    freq_frame <- get_evofreq(size_df, clones, parents, clone_cmap = "inferno")
    p <- plot_evofreq(freq_frame)
    p <- p + xlab("Time (days)") + ylab("Population (CFUs)")
  }
  p
}

p1 <- m_plot(abundance_ex1)
p2 <- m_plot(abundance_ex2)

## with EvoFreq (smoothed)
p1 <- m_plot(abundance_ex1, method = "evofreq")
p2 <- m_plot(abundance_ex2, method = "evofreq")

pdf("fig2.pdf", width = 12, height = 6)
plot_grid(p1, p2, labels = c('A', 'B'), ncol = 2, align = "v", label_size = 12)
dev.off()

