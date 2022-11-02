# Packages
library(tidyverse)

# Read in results
res1 <- read_csv("./3_Aim2/2_pipeline/2_aim2Power_1.csv")

res1 %>%
  group_by(l, mort2yr_dom, mort2yr_notDom, recur2yr_dom, recur2yr_notDom, amp2yr_dom, amp2yr_notDom) %>%
  summarise(value = mean(value, na.rm = TRUE)) %>%
  filter(!is.na(value)) %>%
  mutate(pctDiffOracle = abs(value - 549.8051)/549.8051) %>% 
  ungroup() %>%
  add_count() %>% 
  ungroup() %>%
  summarise(power = sum(pctDiffOracle <= 0.1)/n)
