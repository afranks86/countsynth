library(bpnmf)
library(tidyverse)

cfg <- read_bpnmf_config("configs/fertility_variational.yml")
res <- bpnmf_run(cfg)
