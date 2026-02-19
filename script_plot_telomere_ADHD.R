library(ggplot2)
library(dplyr)
library(forcats)
library(haven)
library(readxl)
library(viridis)  # for color scale


scatter_reg_lines <- read_excel("scatter_reg_lines_telomere.xlsx")
scatter_reg_lines$intercept <- as.numeric(scatter_reg_lines$intercept)


ggplot(dat, aes(x = beta.exposure, y = beta.outcome)) +
  geom_point() + 
  geom_errorbarh(aes(xmin = beta.exposure - se.exposure, xmax = beta.exposure + se.exposure), color = "black", width = 0) + 
  geom_errorbar(aes(ymin = beta.outcome - se.outcome, ymax = beta.outcome + se.outcome), color = "black", height = 0) + 
  labs(x = "SNP effect on ADHD", y = "SNP effect on TL") + 
  theme_minimal() + 
  geom_abline(data = scatter_reg_lines, aes(intercept = intercept, slope = beta, color = Method), linetype = "solid") +
  xlim(-0.25, 0.254) +   # Define your desired x-axis range
  ylim(-0.15, 0.15) +  # Define your desired y-axis range
  theme(axis.title = element_text(size = 10)) + # Adjust text size for axis labels
  theme(axis.text = element_text(size = 8))  # Adjust text size for both x and y axes
