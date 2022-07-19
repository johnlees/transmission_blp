library(ggplot2)
library(readxl)
library(reshape2)
library(cowplot)

setwd("~/Documents/EBI/Papers/transmission_blp/")

HillsN1 <- read_excel("HillsN1.xlsx")
HillsN1 <- reshape2::melt(HillsN1)
colnames(HillsN1) <- c("Time", "HillsN1")
HillsN1$Time <- as.character(HillsN1$Time)

HillsN1 <- HillsN1[!is.na(HillsN1$HillsN1),]
HillsN1$Time[HillsN1$Time == "Inoculum"] <- 0
HillsN1$Time[HillsN1$Time == "2 hours"] <- 2
HillsN1$Time[HillsN1$Time == "6 hours"] <- 6
HillsN1$Time[HillsN1$Time == "1 day"] <- 24
HillsN1$Time[HillsN1$Time == "2 weeks"] <- 24*2*7
HillsN1$Time[HillsN1$Time == "4 weeks"] <- 24*4*7
HillsN1$Time[HillsN1$Time == "6 weeks"] <- 24*6*7
HillsN1$Time <- as.numeric(HillsN1$Time)

p1 <- ggplot(HillsN1, aes(x=Time, y=HillsN1)) +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = FALSE) +
  geom_jitter(width = 0.2, stroke=1, shape=16) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90),
        axis.title = element_text(family = "Arial", size = 12),
        axis.text = element_text(family = "Arial", size = 11)) +
  scale_x_continuous(name="Time", 
                     breaks=c(0, 2, 6, 24, 24*2*7, 24*4*7, 24*6*7),
                     minor_breaks = NULL,
                     labels=c("Inoculum", "2 hours", "6 hours", "1 day", "2 weeks", "4 weeks", "6 weeks"), 
                     limits=c(-20, 24*6*7*1.05),
                     expand=expansion(),
                     trans="sqrt") +
  scale_y_continuous(name="Hills N1", 
                     limits=c(0.1, 1200),
                     trans="sqrt")

ggsave(p1, filename = "fig1F.pdf", device = cairo_pdf, 
       width = 6, height = 6, units = "in")

summary(glm(HillsN1 ~ Time, HillsN1, family = "poisson"))


UniqueClonesNumber <- read_excel("UniqueClonesNumber.xlsx")
UniqueClonesNumber <- reshape2::melt(UniqueClonesNumber)
colnames(UniqueClonesNumber) <- c("Time", "Unique")
UniqueClonesNumber$Time <- as.character(UniqueClonesNumber$Time)

UniqueClonesNumber <- UniqueClonesNumber[!is.na(UniqueClonesNumber$Unique),]
UniqueClonesNumber$Time[UniqueClonesNumber$Unique == "Inoculum"] <- 0
UniqueClonesNumber$Time[UniqueClonesNumber$Unique == "2 hours"] <- 2
UniqueClonesNumber$Time[UniqueClonesNumber$Unique == "6 hours"] <- 6
UniqueClonesNumber$Time[UniqueClonesNumber$Unique == "1 day"] <- 24
UniqueClonesNumber$Time[UniqueClonesNumber$Unique == "2 weeks"] <- 24*2*7
UniqueClonesNumber$Time[UniqueClonesNumber$Unique == "4 weeks"] <- 24*4*7
UniqueClonesNumber$Time[UniqueClonesNumber$Unique == "6 weeks"] <- 24*6*7
UniqueClonesNumber$Time <- as.numeric(HillsN1$Time)

p2 <- ggplot(UniqueClonesNumber, aes(x=Time, y=Unique)) +
  geom_smooth(method = "glm", 
              method.args = list(family = "quasipoisson"), 
              se = FALSE) +
  geom_jitter(width = 0.4, stroke=1, shape=16) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90),
        axis.title = element_text(family = "Arial", size = 12),
        axis.text = element_text(family = "Arial", size = 11)) +
  scale_x_continuous(name="Time", 
                     breaks=c(0, 2, 6, 24, 24*2*7, 24*4*7, 24*6*7),
                     minor_breaks = NULL,
                     labels=c("Inoculum", "2 hours", "6 hours", "1 day", "2 weeks", "4 weeks", "6 weeks"), 
                     limits=c(-1, 24*6*7*1.05),
                     expand=expansion(),
                     trans="sqrt") +
  scale_y_continuous(name="Number of Unique Clones", 
                     limits=c(0.1, 1750),
                     trans="sqrt")

ggsave(p2, filename = "fig1D.pdf", device = cairo_pdf, 
       width = 6, height = 6, units = "in")

summary(glm(Unique ~ Time, UniqueClonesNumber, family = "quasipoisson"))


PercAbun <- read_excel("PercentAbundance.xlsx")
PercAbun <- reshape2::melt(PercAbun)
colnames(PercAbun) <- c("Time", "Abundance")
PercAbun$Time <- as.character(PercAbun$Time)
PercAbun$Abundance <- PercAbun$Abundance / 100

PercAbun <- PercAbun[!is.na(PercAbun$Abundance),]
PercAbun$Time[PercAbun$Time == "Inoculum"] <- 0
PercAbun$Time[PercAbun$Time == "2 hours"] <- 2
PercAbun$Time[PercAbun$Time == "6 hours"] <- 6
PercAbun$Time[PercAbun$Time == "1 day"] <- 24
PercAbun$Time[PercAbun$Time == "2 weeks"] <- 24*2*7
PercAbun$Time[PercAbun$Time == "4 weeks"] <- 24*4*7
PercAbun$Time[PercAbun$Time == "6 weeks"] <- 24*6*7
PercAbun$Time <- as.numeric(PercAbun$Time)

p3 <- ggplot(PercAbun, aes(x=Time, y=Abundance)) +
  geom_smooth(method = "glm", 
              method.args = list(family = "quasibinomial"), 
              se = FALSE) +
  geom_jitter(width = 0.2, stroke=1, shape=16) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90),
        axis.title = element_text(family = "Arial", size = 12),
        axis.text = element_text(family = "Arial", size = 11)) +
  scale_x_continuous(name="Time", 
                     breaks=c(0, 2, 6, 24, 24*2*7, 24*4*7, 24*6*7),
                     minor_breaks = NULL,
                     labels=c("Inoculum", "2 hours", "6 hours", "1 day", "2 weeks", "4 weeks", "6 weeks"), 
                     limits=c(-1, 24*6*7*1.05),
                     expand=expansion(),
                     trans="sqrt") +
  scale_y_continuous(name="Proportion of the Most Abundant Clone", 
                     limits=c(0, 1),
                     labels=scales::percent,
                     trans="sqrt")

ggsave(p3, filename = "fig1E.pdf", device = cairo_pdf, 
       width = 6, height = 6, units = "in")


summary(glm(Abundance ~ Time, PercAbun, family = "quasibinomial"))

pdf("fig1.pdf", width = 16, height = 6)
plot_grid(p2, p3, p1, labels = c('D', 'E', 'F'), ncol = 3, align = "v", label_size = 12)
dev.off()


