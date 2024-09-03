library(tidyverse)
library(haven)

D<-"shortdata/HN2016.dta" |>
  read_dta() 

###############################################################################
# Verify some things
###############################################################################

# 1. prob1 is not used ever

sum(D$prob1L>0)
sum(D$prob1R>0)

# 2. Prizes are all the same:

mean(D$prize1L==D$prize1R)
mean(D$prize2L==D$prize2R)
mean(D$prize3L==D$prize3R)
mean(D$prize4L==D$prize4R)

# 3. Probabilities add up to one:

min(D$prob2L+D$prob3L+D$prob4L)
min(D$prob2R+D$prob3R+D$prob4R)

###############################################################################
# Extract the just parts of the data I need
###############################################################################


D <- D |>
  select(id,choice,prob2L,prob3L,prob4L,prob2R,prob3R,prob4R) |>
  rename(
    probL1 = prob2L,
    probL2 = prob3L,
    probL3 = prob4L,
    probR1 = prob2R,
    probR2 = prob3R,
    probR3 = prob4R
  ) |>
  mutate(choiceLeft = 1-choice) |>
  select(-choice)
