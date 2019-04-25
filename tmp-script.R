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

N <- 1e3
phi <- 0.7
p <- 0.2
nocc <-10
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

p.ch <- process.ch(captureData %>% pull(ch), captureData %>% pull(freq), all = T)

proc = process.data(captureData %>% as.data.frame(), model = "js", begin.time = 1)
ddl = make.design.data(proc)
# initial <- list(Phi = c(`(Intercept)` = 2.17221527553774), p = c(`(Intercept)` = -2.15371289211996),
#      pent = c(`(Intercept)` = -2.31309999137221), N = c(`(Intercept)` = 8.94595618338772))
initial <- list(Phi = c(`(Intercept)` = log(phi)), p = c(`(Intercept)` = log(p)),
                pent = c(`(Intercept)` = log(pent)), N = c(`(Intercept)` = log(Ns - sum(captureData$freq))))
mod = crm(proc, ddl, hessian = TRUE, initial = initial)

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

(Ns.pred - Ns) / Ns
(pent.pred - pent) / pent

e <- c((B.pred - B) / B,
       (Ns.pred - Ns) / Ns,
       (pent.pred - pent) / pent)
