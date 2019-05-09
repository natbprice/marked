library(PurchaseSims)
library(tidyverse)
# install.packages("marked")
# library(marked)


# Helper functions --------------------------------------------------------

# Calculate super population size 
superPopSize <- function(N, nocc, phi) {
  # nocc * (1 - phi) * N + N
  (nocc - 1) * (1 - phi) * N + N
}

# Calculate entrance probability
entranceProb <- function(N, Ns, nocc) {
  (Ns - N) / ((nocc - 1) * Ns)
}


# Simulate data -----------------------------------------------------------

N.true <- 1e4
phi.true <- 0.85
p.true <- 0.70
nocc <- 8
Ns.true <- superPopSize(N.true, nocc, phi.true) 
pent.true <- entranceProb(N.true, Ns.true, nocc)
pBirth <- 1

output <- simul_popan(
  phi = rep(phi.true, nocc),
  p = rep(p.true, nocc),
  pent = rep(pent.true, nocc - 1),
  Ns = Ns.true,
  pBirth = rep(pBirth, nocc)
)

alive <- colSums(output$state)
plot(alive)


CH <- unite(as_tibble(output$CH), "ch", sep = ",") %>%
  mutate(ID = row_number()) %>%
  filter(ch != paste(rep(0, nocc), collapse = ","))

captureData <- CH %>% 
  group_by(ch) %>% 
  summarize(freq = n())

p.ch <- process.ch(captureData %>% pull(ch), captureData %>% pull(freq), all = T)

# Calulate number of unmarked
if (is.null(p.ch$group)) {
  u <- tapply(p.ch$freq, list(p.ch$first), sum)
} else {
  u <- tapply(p.ch$freq,
              list(p.ch$first, model_data$group),
              sum)
}
udot <- sum(u)

# Run CJS model
proc.cjs = process.data(captureData %>% as.data.frame(), model = "cjs", begin.time = 1)
ddl.cjs = make.design.data(proc.cjs)
mod.cjs = crm(proc.cjs, ddl.cjs, hessian = TRUE)

proc = process.data(captureData %>% as.data.frame(), model = "js", begin.time = 1)
ddl = make.design.data(proc)


initial.cjs <- list(Phi = mod.cjs$results$beta$Phi, 
                     p  = mod.cjs$results$beta$p,
                     pent  = log(mean(u/udot)),
                     N  = log(udot))
mod = crm(proc, ddl, hessian = TRUE, 
          # initial = initial.cjs,
          method = "CMAES",
          itnmax = 10e3)



# Calculate model predictions of retention and purchase probabilities
modPred <-
  predict(mod,
          ddl = ddl,
          se = TRUE)

Ns.pred <- sum(captureData$freq) + modPred$N[1]
pent.pred <- modPred$pent$estimate[1]
B <- Ns.true * pent.true
B.pred <- Ns.pred * pent.pred
(B.pred - B) / B

(Ns.pred - Ns.true) / Ns.true
(pent.pred - pent.true) / pent.true


# Groups ------------------------------------------------------------------

N.true2 <- 1e4
phi.true2 <- phi.true/2
p.true2 <- p.true/2
Ns.true2 <- superPopSize(N.true2, nocc, phi.true2) 
pent.true2 <- entranceProb(N.true2, Ns.true2, nocc)

output2 <- simul_popan(
  phi = rep(phi.true2, nocc),
  p = rep(p.true2, nocc),
  pent = rep(pent.true2, nocc - 1),
  Ns = Ns.true2,
  pBirth = rep(pBirth, nocc)
)

CH2 <- unite(as_tibble(output2$CH), "ch", sep = ",") %>%
  mutate(ID = row_number()) %>%
  filter(ch != paste(rep(0, nocc), collapse = ","))

captureData <- full_join(CH %>% mutate(groupID = "A"),
                         CH2 %>% mutate(groupID = "B")) %>% 
  mutate(groupID = as.factor(groupID)) %>% 
  group_by(groupID, ch) %>% 
  summarize(freq = n())

p.ch <- process.ch(captureData %>% pull(ch), captureData %>% pull(freq), all = T)

# Calulate number of unmarked
if (is.null(captureData$groupID)) {
  u <- tapply(p.ch$freq, list(p.ch$first), sum)
} else {
  u <- tapply(p.ch$freq,
              list(p.ch$first, pull(captureData, groupID)),
              sum)
}
udot <- colSums(u)

# Run CJS model
proc.cjs = process.data(captureData %>% as.data.frame(),
                        model = "cjs", 
                        begin.time = 1,
                        groups = "groupID")
ddl.cjs = make.design.data(proc.cjs)
mod.cjs = crm(
  proc.cjs,
  ddl.cjs,
  model.parameters = list(Phi = list(formula = ~ groupID),
                          p = list(formula = ~ groupID)),
  hessian = TRUE
)

modPred.cjs <-
  predict(mod.cjs,
          ddl = ddl.cjs,
          se = TRUE)

proc = process.data(captureData %>% as.data.frame(), 
                    model = "js", 
                    begin.time = 1,
                    groups = "groupID")
ddl = make.design.data(proc)


pent.ini1 <- u[nocc,1]/udot[[1]]
pent.ini2 <- u[nocc,2]/udot[[2]]
initial.cjs <- list(Phi = mod.cjs$results$beta$Phi,
                    p = mod.cjs$results$beta$p,
                    pent  = c(`(Intercept)` = log(pent.ini1),
                              groupIDB = log(pent.ini2) - log(pent.ini1)),
                    N  = c(`(Intercept)` = log(udot[[1]]), 
                           groupIDB = log(udot[[2]]) - log(udot[[1]])))


mod = crm(proc, 
          ddl, 
          model.parameters = list(Phi = list(formula = ~ groupID),
                                  p = list(formula = ~ groupID),
                                  pent = list(formula = ~ groupID),
                                  N = list(formula = ~ groupID)),
          initial = initial.cjs,
          method = "CMAES",
          detectAllBirths = T,
          hessian = TRUE
          # control = list(maxit = 1e3)
          )



# Calculate model predictions of retention and purchase probabilities
modPred <-
  predict(mod,
          ddl = ddl,
          se = TRUE)

Ns.pred1 <- udot[1] + modPred$N[1,2]
Ns.pred2 <- udot[2] + modPred$N[2,2]

(Ns.pred1 - Ns.true) / Ns.true
(Ns.pred2 - Ns.true2) / Ns.true2

pent.pred1 <- modPred$pent$estimate[1]
pent.pred2 <- modPred$pent$estimate[9]

(pent.pred1 - pent.true) / pent.true
(pent.pred2 - pent.true2) / pent.true2

B <- Ns.true * pent.true
B2 <- Ns.true2 * pent.true2
B.pred1 <- Ns.pred1 * pent.pred1
B.pred2 <- Ns.pred2 * pent.pred2

(B.pred1 - B) / B
(B.pred2 - B2) / B2



