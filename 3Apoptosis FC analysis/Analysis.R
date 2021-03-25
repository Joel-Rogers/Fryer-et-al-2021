rm(list=ls())

library(tidyverse)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(minpack.lm)
library(outliers)
library(scales)

# Load and somewhat tidy data

dataparent <- getwd()

maincsvs <- list.files(dataparent, full.names = TRUE, pattern = "^statistics_plate") %>%
  lapply(function(x) {
    read.csv(x, encoding="UTF-8", na.strings = "NA") %>%
      as_tibble()
})

maindf <- NULL
for(i in 1:length(maincsvs)){
  maindf <- bind_rows(maindf, maincsvs[[i]])
}


maindf2 <- maindf %>%
  select(-Experiment, -Group, -Filename, -Workspace, -Compensation.Source,
         -Plot.Title, -Volume, -Concentration, -Y.Mean, -Y.Median, -Y.SD,
         -Y.Peak, -Y..CV, -Y.rSD, -Y..rCV, -X..CV, -X.rSD, -X..rCV) %>%
  na_if("N/A")


View(maindf2)

# Condense the data down into a slightly more sensible representation

maindf3 <- maindf2 %>%
  mutate(`Total events` = ifelse(Gate == "All Events", Count, NA),
    `Events in cell gate` = as.double(ifelse(Gate == "R1" & X.Parameter == "FSC-A",
                                   Count, NA)),
    `Cell Mean BluFL1` = as.double(ifelse(Gate == "R1" & X.Parameter == "BL1-A",
                                       X.Mean, NA)),
    `Cell Median BluFL1` = as.double(ifelse(Gate == "R1" & X.Parameter == "BL1-A",
                                  X.Median, NA)),
    `Cell BluFL1 SD` = as.double(ifelse(Gate == "R1" & X.Parameter == "BL1-A",
                              X.SD, NA)),
    `Events in positive gate` = as.double(ifelse(Gate == "R2", Count, NA)),
    `Percentage apoptotic` = as.double(ifelse(Gate == "R2", X.Gated, NA)),
    Plate = ifelse(Plate == "2020-12-16 - Plate scFv2", 2, 
                        ifelse(Plate == "2020-12-17 - plate scFv3", 3, 1))) %>%
  select(-X.Mean, -X.Median, -X.SD, -X.Peak, -X.Total, -X.Gated,
         -X.Parameter, -Y.Parameter, -Gate, -Count)

maindf4 <- maindf3 %>%
  replace(is.na(.),0) %>%
  group_by(Plate, Sample) %>%
  summarise(across(`Total events`:`Percentage apoptotic`, sum))

View(maindf4)


# Convert all the samples into their construct valency and 
# concentrations using the originally generated (randomised) plate layouts.
# Samples in column 8 are staining controls - their layout was not randomised.

lookupcsvs <- list.files(dataparent, full.names = TRUE,
                         pattern = "^shuffle_output") %>%
  lapply(function(x) {
    read.csv(x, encoding="UTF-8", na.strings = "NA") %>%
      as_tibble()
})

# Define a function to add row numbers to a dataframe
# From: https://rdrr.io/github/skgrange/threadr/src/R/add_row_numbers.R
add_row_numbers <- function(df, name = "row_number", zero_based = FALSE) {
  # Drop variable if exists
  if (name %in% names(df)) df[, name] <- NULL
  # Create sequence of integers
  sequence <- seq.int(1, to = nrow(df))
  if (zero_based) sequence <- sequence - 1L
  # Add sequence to data frame
  df[, name] <- sequence
  # Move variable to the first column position
  df <- select(df, !!name, everything())
  return(df)
}


lookupdf <- NULL
for(i in 1:length(lookupcsvs)){
  for(c in 1:length(lookupcsvs[[i]])){
    columnname = colnames(lookupcsvs[[i]][c])
    lookupdf_transient <- lookupcsvs[[i]][c] %>%
      add_row_numbers(name = "Plate row") %>%
      mutate(Plate = i,
             Valency = substr(columnname,2,nchar(columnname)-1),
             `Plate column` = c, `Plate row` = chartr("12345678", "ABCDEFGH",
                                                      `Plate row`),
             Sample = as.character(paste(`Plate row`,`Plate column`, sep='')),
             `Concentration (nM)` =
               as.double(substr(
                 as.character(.[[2]]),1,nchar(as.character(.[[2]]))-3))) %>%
      select(Plate, Sample, Valency, `Concentration (nM)`)
    lookupdf <- bind_rows(lookupdf, lookupdf_transient)
  }
}

View(lookupdf)


# Now correct the concentrations to reflect the change in the highest dosage 
# and the serial dilution that I made:
# (a trial experiment revealed that the original dosages were too high and 
# covered too narrow a dynamic range, so 
# we halved the highest dosage and switched from a 2-fold serial dilution to 
# four-fold - this was done after the randomised layouts were already 
# generated).
original_dosages = c(230,115,57.5,28.75,14.375,7.1875,3.59375,1.796875)
new_dosages = c(115,28.75,7.1875,1.796875,0.44921875,0.112304688,0.028076172,
                0.007019043)

lookupdf2 <- lookupdf %>%
  mutate(Corrected.concentration.nM =
           new_dosages[match(`Concentration (nM)`, original_dosages)]) %>%
  select(-`Concentration (nM)`)

View(lookupdf2)


# And finally, use these lookup tables to un-randomise the data in the maindfs

maindf5 <- maindf4 %>%
  mutate(Sample = as.character(Sample)) %>%
  full_join(lookupdf2)

View(maindf5)


# ~~ Some sanity-checking code/ analysis is omitted here ~~ #


# Plot the data, and fit a set of Hill eqns to it 
# to determine the EC50s and Hill slopes of each.

maindf7 <- maindf5 %>%
  filter(!is.na(Valency)) %>%
  group_by(Valency, Corrected.concentration.nM) %>%
  summarise(Mean.percentage.apoptotic = mean(`Percentage apoptotic`),
            Standard.deviation.apoptotic = sd(`Percentage apoptotic`)) %>%
  mutate(Valency = as.factor(Valency)) %>%
  group_by(Valency)

View(maindf7)


# Initialise some colourblind-friendly palettes:

# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# Gonna drop this stupid pale yellow if possible, and the
# light blue as well for this one specifically
cbPalette2 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# As a connected scatter plot, without jitter
dev.new()
ggplot(data = maindf7,
       mapping = aes(x = Corrected.concentration.nM,
                     y = Mean.percentage.apoptotic,
                     color = Valency)) +
  geom_point(size = 3) +
  geom_errorbar(aes(
    ymin=ifelse(Mean.percentage.apoptotic-Standard.deviation.apoptotic < 0, 0,
                Mean.percentage.apoptotic-Standard.deviation.apoptotic),
    ymax=Mean.percentage.apoptotic+Standard.deviation.apoptotic),
    width = 0.1, size=0.7) +
  geom_line(size=1) +
  geom_hline(yintercept=50, linetype="dashed", color = "black") + 
  scale_x_log10(labels = c(0.01, 0.1, 1, 10, 100), 
                breaks = c(0.01, 0.1, 1, 10, 100)) +
  #scale_x_log10(labels = number_format(accuracy = 0.01)) +
  scale_colour_manual(values=cbPalette2) +
  # Can consider dropping the line fitting
  #geom_smooth(se = FALSE) + 
  #annotate(geom = "text", -Inf, Inf, hjust = -0.2, vjust = 5, label = lm_eqnhillphosph(fitpredictphosphstd), parse = TRUE, size = 5, colour = "blue")+
  #geom_smooth(method = "lm", se = FALSE) +
  labs(title = paste("Induction of apoptosis by anti-TRAIL-R1 scFv displayed\n",
       " at various valencies"), color = "scFv valency") +
  xlab("Effective scFv concentration (nM)") +
  ylab("Proportion of apoptotic cells (%)") +
  theme_classic() +
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.title = element_text(size=15), legend.text = element_text(size=15),
        axis.title.x = element_text(vjust=-1),
        legend.title = element_text(size=15),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        aspect.ratio = 1/1.5)
ggsave("Raw-data-plot-all.png")

# Looks good. Time to fit some useful curves to each.


# Hill eqn that I'll use: 
# y = bottom + (top-bottom)/(1+10^((log(EC50)-x)*Hillslope))
# Or, expressed more simply and for this particular situation:
# y = 100/(1+10^((log(E)-x)*H))


# nsl doesn't handle space well, so mutate the column name
maindf8 <- maindf5 %>%
  filter(!is.na(Valency)) %>%
  mutate(Percentage.apoptotic = `Percentage apoptotic`) %>%
  select(-`Percentage apoptotic`)

View(maindf8)


# Allowing fits to bottom part of the sigmoid as well...
# (Because the monovalent construct doesn't get near 100% apoptosis, allowing 
# the upper part of the sigmoid to be fit doesn't give a sensible result)


# Fitting the replicates altogether, rather than just the means.
# (Yes, I could have done the below as a function, but I didn't want to mess 
# around with converting between strings and variable names)

E0 = 10; H0 = 1; B0 = 0; Top = 100; fiteqnall1tb <- nlsLM(
  Percentage.apoptotic ~ B +
    (Top-B)/(1+(E/Corrected.concentration.nM)^H),
  start = list(E = E0, H = H0, B = B0),
  lower = c(0,0.1,0), upper = c(Inf, 10,50),
  data = filter(maindf8, Valency == 1),
  control = nls.control(maxiter = 1000, minFactor = 1/(4096**9),
                        warnOnly = TRUE, tol = 1e-08)); fiteqnall1tb

E0 = 10; H0 = 1; B0 = 0; Top = 100; fiteqnall2tb <- nlsLM(
  Percentage.apoptotic ~ B +
    (Top-B)/(1+(E/Corrected.concentration.nM)^H),
  start = list(E = E0, H = H0, B = B0),
  lower = c(0,0.1,0), upper = c(Inf, 10,50),
  data = filter(maindf8, Valency == 2),
  control = nls.control(maxiter = 1000, minFactor = 1/(4096**9),
                        warnOnly = TRUE, tol = 1e-08)); fiteqnall2tb

E0 = 10; H0 = 1; B0 = 0; Top = 100; fiteqnall3tb <- nlsLM(
  Percentage.apoptotic ~ B +
    (Top-B)/(1+(E/Corrected.concentration.nM)^H),
  start = list(E = E0, H = H0, B = B0),
  lower = c(0,0.1,0), upper = c(Inf, 10,50),
  data = filter(maindf8, Valency == 3),
  control = nls.control(maxiter = 1000, minFactor = 1/(4096**9),
                        warnOnly = TRUE, tol = 1e-08)); fiteqnall3tb

E0 = 10; H0 = 1; B0 = 0; Top = 100; fiteqnall4tb <- nlsLM(
  Percentage.apoptotic ~ B +
    (Top-B)/(1+(E/Corrected.concentration.nM)^H),
  start = list(E = E0, H = H0, B = B0),
  lower = c(0,0.1,0), upper = c(Inf, 10,50),
  data = filter(maindf8, Valency == 4),
  control = nls.control(maxiter = 1000, minFactor = 1/(4096**9),
                        warnOnly = TRUE, tol = 1e-08)); fiteqnall4tb

E0 = 10; H0 = 1; B0 = 0; Top = 100; fiteqnall5tb <- nlsLM(
  Percentage.apoptotic ~ B +
    (Top-B)/(1+(E/Corrected.concentration.nM)^H),
  start = list(E = E0, H = H0, B = B0),
  lower = c(0,0.1,0), upper = c(Inf, 10,50),
  data = filter(maindf8, Valency == 5),
  control = nls.control(maxiter = 1000, minFactor = 1/(4096**9),
                        warnOnly = TRUE, tol = 1e-08)); fiteqnall5tb

E0 = 10; H0 = 1; B0 = 0; Top = 100; fiteqnall6tb <- nlsLM(
  Percentage.apoptotic ~ B +
    (Top-B)/(1+(E/Corrected.concentration.nM)^H),
  start = list(E = E0, H = H0, B = B0),
  lower = c(0,0.1,0), upper = c(Inf, 10,50),
  data = filter(maindf8, Valency == 6),
  control = nls.control(maxiter = 1000, minFactor = 1/(4096**9),
                        warnOnly = TRUE, tol = 1e-08)); fiteqnall6tb

coefficients12 = tibble(Valency = NA, EC50 = NA, EC50.sd = NA,
                        Hill.coefficient = NA,
                        Hill.coefficient.sd = NA, Bottom = NA,
                        Bottom.sd = NA) %>%
  add_row(Valency = "1", EC50 =summary(fiteqnall1tb)$coefficients[1],
          EC50.sd = as.double(
            as_tibble(summary(fiteqnall1tb)$coefficients)[2][[1]][1]),
          Hill.coefficient = summary(fiteqnall1tb)$coefficients[2],
          Hill.coefficient.sd = as.double(
            as_tibble(summary(fiteqnall1tb)$coefficients)[2][[1]][2]),
          Bottom = summary(fiteqnall1tb)$coefficients[3],
          Bottom.sd = as.double(
            as_tibble(summary(fiteqnall1tb)$coefficients)[2][[1]][3])) %>%
  add_row(Valency = "2", EC50 =summary(fiteqnall2tb)$coefficients[1],
          EC50.sd = as.double(
            as_tibble(summary(fiteqnall2tb)$coefficients)[2][[1]][1]),
          Hill.coefficient = summary(fiteqnall2tb)$coefficients[2],
          Hill.coefficient.sd = as.double(
            as_tibble(summary(fiteqnall2tb)$coefficients)[2][[1]][2]),
          Bottom = summary(fiteqnall2tb)$coefficients[3],
          Bottom.sd = as.double(
            as_tibble(summary(fiteqnall2tb)$coefficients)[2][[1]][3])) %>%
  add_row(Valency = "3", EC50 =summary(fiteqnall3tb)$coefficients[1],
          EC50.sd = as.double(
            as_tibble(summary(fiteqnall3tb)$coefficients)[2][[1]][1]),
          Hill.coefficient = summary(fiteqnall3tb)$coefficients[2],
          Hill.coefficient.sd = as.double(
            as_tibble(summary(fiteqnall3tb)$coefficients)[2][[1]][2]),
          Bottom = summary(fiteqnall3tb)$coefficients[3],
          Bottom.sd = as.double(
            as_tibble(summary(fiteqnall3tb)$coefficients)[2][[1]][3])) %>%
  add_row(Valency = "4", EC50 =summary(fiteqnall4tb)$coefficients[1],
          EC50.sd = as.double(
            as_tibble(summary(fiteqnall4tb)$coefficients)[2][[1]][1]),
          Hill.coefficient = summary(fiteqnall4tb)$coefficients[2],
          Hill.coefficient.sd = as.double(
            as_tibble(summary(fiteqnall4tb)$coefficients)[2][[1]][2]),
          Bottom = summary(fiteqnall4tb)$coefficients[3],
          Bottom.sd = as.double(
            as_tibble(summary(fiteqnall4tb)$coefficients)[2][[1]][3])) %>%
  add_row(Valency = "5", EC50 =summary(fiteqnall5tb)$coefficients[1],
          EC50.sd = as.double(
            as_tibble(summary(fiteqnall5tb)$coefficients)[2][[1]][1]),
          Hill.coefficient = summary(fiteqnall5tb)$coefficients[2],
          Hill.coefficient.sd = as.double(
            as_tibble(summary(fiteqnall5tb)$coefficients)[2][[1]][2]),
          Bottom = summary(fiteqnall5tb)$coefficients[3],
          Bottom.sd = as.double(
            as_tibble(summary(fiteqnall5tb)$coefficients)[2][[1]][3])) %>%
  add_row(Valency = "6", EC50 =summary(fiteqnall6tb)$coefficients[1],
          EC50.sd = as.double(
            as_tibble(summary(fiteqnall6tb)$coefficients)[2][[1]][1]),
          Hill.coefficient = summary(fiteqnall6tb)$coefficients[2],
          Hill.coefficient.sd = as.double(
            as_tibble(summary(fiteqnall6tb)$coefficients)[2][[1]][2]),
          Bottom = summary(fiteqnall6tb)$coefficients[3],
          Bottom.sd = as.double(
            as_tibble(summary(fiteqnall6tb)$coefficients)[2][[1]][3])) %>%
  mutate(Set = "coefficients12") %>%
  filter(!is.na(Valency)); coefficients12


outdf <- tibble(coefficients12, date = Sys.Date()); outdf
write_csv(outdf, "Analysis.csv", na = "NA", append = TRUE, quote_escape = "double")

hilleqnall1tb = function(x){
  as.double(coefficients12$Bottom[coefficients12$Valency == 1]) +
    (100 - as.double(coefficients12$Bottom[coefficients12$Valency == 1]))/
    (1+(as.double(coefficients12$EC50[coefficients12$Valency == 1])/
          x)^as.double(coefficients12$Hill.coefficient[coefficients12$Valency == 1]))
}
hilleqnall2tb = function(x){
  as.double(coefficients12$Bottom[coefficients12$Valency == 2]) +
    (100 - as.double(coefficients12$Bottom[coefficients12$Valency == 2]))/
    (1+(as.double(coefficients12$EC50[coefficients12$Valency == 2])/
          x)^as.double(coefficients12$Hill.coefficient[coefficients12$Valency == 2]))
}
hilleqnall3tb = function(x){
  as.double(coefficients12$Bottom[coefficients12$Valency == 3]) +
    (100 - as.double(coefficients12$Bottom[coefficients12$Valency == 3]))/
    (1+(as.double(coefficients12$EC50[coefficients12$Valency == 3])/
          x)^as.double(coefficients12$Hill.coefficient[coefficients12$Valency == 3]))
}
hilleqnall4tb = function(x){
  as.double(coefficients12$Bottom[coefficients12$Valency == 4]) +
    (100 - as.double(coefficients12$Bottom[coefficients12$Valency == 4]))/
    (1+(as.double(coefficients12$EC50[coefficients12$Valency == 4])/
          x)^as.double(coefficients12$Hill.coefficient[coefficients12$Valency == 4]))
}
hilleqnall5tb = function(x){
  as.double(coefficients12$Bottom[coefficients12$Valency == 5]) +
    (100 - as.double(coefficients12$Bottom[coefficients12$Valency == 5]))/
    (1+(as.double(coefficients12$EC50[coefficients12$Valency == 5])/
          x)^as.double(coefficients12$Hill.coefficient[coefficients12$Valency == 5]))
}
hilleqnall6tb = function(x){
  as.double(coefficients12$Bottom[coefficients12$Valency == 6]) +
    (100 - as.double(coefficients12$Bottom[coefficients12$Valency == 6]))/
    (1+(as.double(coefficients12$EC50[coefficients12$Valency == 6])/
          x)^as.double(coefficients12$Hill.coefficient[coefficients12$Valency == 6]))
}

ggplot(data = maindf7,
       mapping = aes(x = Corrected.concentration.nM,
                     y = Mean.percentage.apoptotic,
                     color = Valency)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(
    ymin=ifelse(Mean.percentage.apoptotic-Standard.deviation.apoptotic < 0, 0,
                Mean.percentage.apoptotic-Standard.deviation.apoptotic),
    ymax=Mean.percentage.apoptotic+Standard.deviation.apoptotic),
    width = 0.1, size=0.7) +
  geom_hline(yintercept=50, linetype="dashed", color = "black") + 
  scale_x_log10(labels = c(0.01, 0.1, 1, 10, 100), 
                breaks = c(0.01, 0.1, 1, 10, 100)) +
  scale_colour_manual(values=cbPalette2) +
  stat_function(fun = hilleqnall1tb, colour = "#999999", size = 1) +
  stat_function(fun = hilleqnall2tb, colour = "#E69F00", size = 1) +
  stat_function(fun = hilleqnall3tb, colour = "#56B4E9", size = 1) +
  stat_function(fun = hilleqnall4tb, colour = "#009E73", size = 1) +
  stat_function(fun = hilleqnall5tb, colour = "#D55E00", size = 1) +
  stat_function(fun = hilleqnall6tb, colour = "#CC79A7", size = 1) +
  labs(title = paste("Induction of apoptosis by anti-TRAIL-R1 scFv\ndisplayed",
                     "at various valencies"), color = "scFv valency") +
  xlab("Effective scFv concentration (nM)") +
  ylab("Proportion of apoptotic cells (%)") +
  theme_classic() +
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.title = element_text(size=15), legend.text = element_text(size=15),
        axis.title.x = element_text(vjust=-1),
        legend.title = element_text(size=15),
        axis.text.x = element_text(size = 12, colour = "black",
                                   margin = 
                                     unit(c(t = 2.5, r = 0, b = 0, l = 0), 
                                          "mm")),
        axis.text.y = element_text(size = 12, colour = "black",
                                   margin = 
                                     unit(c(t = 0, r = 2.5, b = 0, l = 0), 
                                          "mm")),
        axis.ticks.length = unit(-1.4, "mm"),
        aspect.ratio = 1/1.5)
ggsave("scatter-Hilleqn-tb-all-reps.png")

ggplot(data = coefficients12, mapping = aes(
  x = as.factor(Valency),
  y = EC50,
  fill = as.factor(Valency))) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(
    ymin=ifelse(EC50 - EC50.sd <0, 0, EC50 - EC50.sd),
    ymax=EC50 + EC50.sd),
    width=.2, size = 1) +
  scale_y_log10() +
  scale_fill_manual(values=cbPalette2) +
  labs(title = paste("Induction of apoptosis by anti-TRAIL-R1 scFv\ndisplayed",
                     "at various valencies"), fill = "scFv valency") +
  xlab("scFv presentation valency") +
  ylab("EC50 (nM scFv)") +
  theme_classic() +
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.title = element_text(size=15), legend.text = element_text(size=15),
        axis.title.x = element_text(vjust=-1),
        legend.title = element_text(size=15),
        axis.text.x = element_text(size = 12, colour = "black",
                                   margin = 
                                     unit(c(t = 2.5, r = 0, b = 0, l = 0), 
                                          "mm")),
        axis.text.y = element_text(size = 12, colour = "black",
                                   margin = 
                                     unit(c(t = 0, r = 2.5, b = 0, l = 0), 
                                          "mm")),
        axis.ticks.length = unit(-1.4, "mm"),
        aspect.ratio = 1/1.5)
ggsave("Bar-Hilleqn-tb-all-reps.png")


# Let's try a statistical test.
# It looks like a Welch's two-tailed t-test would work well for this. There are 
# 21 "degrees of freedom" (8*3 datapoints, minus 3 parameters allowed to vary)
coefficients12


t_tested <- coefficients12 %>%
  mutate(t.test.1x = sqrt(((coefficients12$EC50[Valency == 1] - EC50)/
           (sqrt((coefficients12$EC50.sd[Valency == 1]^2)/3 + 
                   (EC50.sd^2)/3)))^2),
         t.test.5x = sqrt(((coefficients12$EC50[Valency == 5] - EC50)/
           (sqrt((coefficients12$EC50.sd[Valency == 5]^2)/3 + 
                   (EC50.sd^2)/3)))^2))

t_tested
# So using the lookup table for 21 df, two-tailed, the t-value is greater 
# than the critical value in the table for p = 0.001. So p < 0.001 for all
# valencies (except 1x of course, which is undefined). Same for the 5x
# comparison, except that it's only p < 0.005 for 4x.
