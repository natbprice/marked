library(PurchaseSims)
library(tidyverse)
# install.packages("marked")
# library(marked)

# Calculate super population size 
superPopSize <- function(N, nocc, phi) {
  nocc * (1- phi) * N + N
}

# Calculate entrance probability
entranceProb <- function(N, Ns, nocc) {
  (Ns - N) / (nocc * Ns)
}

N <- 10e3
phi <- 0.5
p <- 0.5
nocc <- 20
Ns <- superPopSize(N, nocc, phi) 
pent <- entranceProb(N, Ns, nocc)
pBirth <- 1

output <- simul_popan(
  phi = rep(phi, nocc),
  p = rep(p, nocc + 1),
  pent = rep(pent, nocc),
  Ns = Ns,
  pBirth = rep(pBirth, nocc + 1)
)

CH <- unite(as_tibble(output$CH[, -1]), "ch", sep = "") %>%
  mutate(ID = row_number()) %>%
  filter(ch != paste(rep(0, nocc), collapse = ""))

captureData <- CH %>% 
  group_by(ch) %>% 
  summarize(freq = n()) 

p.ch <- process.ch(captureData %>% pull(ch), captureData %>% pull(freq))

proc = process.data(captureData %>% as.data.frame(), model = "js", begin.time = 1)
ddl = make.design.data(proc)
mod = crm(proc, ddl, hessian = TRUE)

# Calculate model predictions of retention and purchase probabilities
modPred <-
  predict(mod,
          ddl = ddl,
          se = TRUE)

Ns.pred <- sum(captureData$freq) + modPred$N[1]
pent.pred <- modPred$pent$estimate[1]
B <- Ns * pent
B.pred <- Ns.pred * pent.pred
(B.pred - B) / B

Ns.pred/Ns
pent.pred/pent
