rm(list=ls())

library(tidyverse)
library(ggplot2)
library(ggthemes)

filtdf <- read_csv("Compiled-data-for-R.csv") %>%
  filter(Sample == 'BG-PEG-NH2')
myplot <- ggplot(data = filtdf, mapping = aes(x = `Concentration (uM)`, y = `BG-PEG-NH2 (9.6 - 10.0 min)`)) +
  geom_point()
myplot +
  geom_smooth(method = lm, se = FALSE)
lm1 <- lm(data = filtdf, `BG-PEG-NH2 (9.6 - 10.0 min)`~`Concentration (uM)`); lm1
anova(lm1)
summary(lm1)

italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2
format(coef(lm1)[1], digits = 2)

lm_eqn <- function(df){
  lmx <- lm(`BG-PEG-NH2 (9.6 - 10.0 min)`~`Concentration (uM)`, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(lmx)[1], digits = 4), 
                        b = format(coef(lmx)[2], digits = 4), 
                        r2 = format(summary(lmx)$r.squared, digits = 4)))
  as.character(as.expression(eq));                 
}
myplot + 
  geom_smooth(method = lm, se = FALSE) +
  labs(title = "BG-PEG-NH2 Standard Curve", x = "Concentration (uM)", y = "Peak area (mAU*s)") +
  theme_classic() +
  theme(plot.title = element_text(size=20), axis.title = element_text(size=15))

# But this won't give anything sensible at the lower concentrations
# (which we're interested in), and we believe this is largely due to the 5 uM
# measurement being spuriously high (it 'bleeds into' adjacent peaks). So
# exclude this and analyse again.

filtdf2 <- read_csv("Compiled-data-for-R.csv") %>%
  filter(Sample == 'BG-PEG-NH2', `Concentration (uM)` != 5)
myplot <- ggplot(data = filtdf2, mapping = aes(x = `Concentration (uM)`, y = `BG-PEG-NH2 (9.6 - 10.0 min)`)) +
  geom_point()
myplot +
  geom_smooth(method = lm, se = FALSE)
lm2 <- lm(data = filtdf2, `BG-PEG-NH2 (9.6 - 10.0 min)`~`Concentration (uM)`); lm2
anova(lm2)


lm_eqn <- function(df){
  lmx <- lm(`BG-PEG-NH2 (9.6 - 10.0 min)`~`Concentration (uM)`, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(lmx)[1], digits = 4), 
                        b = format(coef(lmx)[2], digits = 4), 
                        r2 = format(summary(lmx)$r.squared, digits = 6)))
  as.character(as.expression(eq));                 
}
dev.new()
myplot + 
  geom_smooth(method = lm, se = FALSE) +
  labs(title = "BG-PEG-NH2 Standard Curve", x = "Concentration (uM)", y = "Peak area (mAU*s)") +
  theme_classic() +
  theme(plot.title = element_text(size=20), axis.title = element_text(size=15))
ggsave("Plot.png")
dev.off()


# The 'coefficients' in summary() of the linear model for this standard 'curve'
summary(lm2)$coefficients
# give the intercept (c):
intcpt <- summary(lm2)$coefficients[1]
# and the gradient (m):
grad <- summary(lm2)$coefficients[2]

# so y = mx + c, with y as peak area and x as concentration of BG-PEG-NH2
# Use this linear model to estimate the remaining amount of BG-PEG-NH2 
# (substrate) remaining after reaction to form BG-PEG-MA
# (use x = (y - c)/m):

ma_bg_area <- read_csv("Compiled-data-for-R.csv") %>%
  filter(Sample == 'MA-BG') %>%
  select(`BG-PEG-NH2 (9.6 - 10.0 min)`) %>%
  as.numeric()

ma_bg_substrate_remaining <- (ma_bg_area - intcpt)/grad

# Input amount of substrate:
ma_bg_input_substrate <- read_csv("Compiled-data-for-R.csv") %>%
  filter(Sample == 'MA-BG') %>%
  select(`Concentration (uM)`) %>%
  as.numeric()

formation_efficiency_estimate <- 100 - (ma_bg_substrate_remaining/
                                    ma_bg_input_substrate)*100

formation_efficiency_estimate
# Rounding --> 99.9 % efficient. 

