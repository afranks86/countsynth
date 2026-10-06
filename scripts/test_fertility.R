library(countsynth)
library(tidyverse)

cfg <- read_countsynth_config("configs/fertility_variational.yml")
res <- countsynth_run(cfg)
