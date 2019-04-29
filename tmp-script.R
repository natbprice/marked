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

N.true <- 400e3
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


CH <- unite(as_tibble(output$CH), "ch", sep = "") %>%
  mutate(ID = row_number()) %>%
  filter(ch != paste(rep(0, nocc), collapse = ""))

captureData <- CH %>% 
  group_by(ch) %>% 
  summarize(freq = n()) 

ch <- t(captureData %>%
          pull(ch) %>%
          sapply(function(s)
            substring(s, 1:nchar(s), 1:nchar(s))) %>% 
          apply(2, as.numeric))

freq <- captureData %>% 
  pull(freq)

first <- apply(ch, 1, function(x) {first(which(x == 1))})
firstMat <- t(sapply(first, function(x, nocc) {
  return(c(rep(0, x[1] - 1), 1, rep(0, nocc - x[1])))
}, nocc = nocc))
u <- colSums(freq * firstMat)
udot = sum(u)

# Coefficients
beta.N <- log(Ns.true - udot)
beta.phi <- log(phi.true)
beta.p <- log(p.true)
beta.pent <- log(pent.true)

# Parameters
Ns <- exp(beta.N) + udot
phi <- rep(exp(beta.phi), nocc)
p <- rep(exp(beta.p), nocc)
beta <- rep(exp(beta.pent), nocc - 1)
beta <- c(1 - sum(beta), beta)

# psi <- numeric(nocc)
# psi[1] <- beta[1]
# for(i in 1:(nocc - 1)) {
#   psi[i + 1] <- beta[i + 1] + beta[1] * prod((1 - p[1:i]) * phi[1:i])
#   # psi[i + 1] <- psi[i] * (1 - p[i]) * phi[i] + beta[i + 1]
# }
p.entry <- numeric(nocc)
p.entry[1] <- beta[1] * p[1]
for(i in 1:(nocc - 1)) {
  p.entry[i + 1] <- beta[i + 1] + p[i+1] * beta[1] * prod((1 - p[1:i]) * phi[1:i])
  # psi[i + 1] <- psi[i] * (1 - p[i]) * phi[i] + beta[i + 1]
}

# Multinomial coefficient
lmultinomial <- function (x) {
  lfactorial(sum(x)) - sum(lfactorial(x))
}

Ns.vec <- seq(udot, Ns.true * 1.2, length.out = 20)
lnl1a <- numeric(length(Ns.vec))
l1a <- numeric(length(Ns.vec))
for(i in 1:length(Ns.vec)) {
  
  # psi is multiplied by p, but all births are detected...
  lnl1a[i] <- lchoose(Ns.vec[i], udot) +
    udot * log(sum(p.entry)) +
    (Ns.vec[i] - udot) * log(1 - sum(p.entry))
  
}

plot(Ns.vec, lnl1a)
lines(c(Ns,Ns), c(min(lnl1a), max(lnl1a)))

lnl1b <- lmultinomial(u) + 
  sum(u * log((psi * p)/sum(psi * p)))

if(f_eval %% 100 == 0 | f_eval == 1){browser()}

lnl <- lnl - (lnl1a + lnl1b)


p.ch <- process.ch(captureData %>% pull(ch), captureData %>% pull(freq), all = T)

# Run CJS model
proc.cjs = process.data(captureData %>% as.data.frame(), model = "cjs", begin.time = 1)
ddl.cjs = make.design.data(proc.cjs)
mod.cjs = crm(proc.cjs, ddl.cjs, hessian = TRUE)

proc = process.data(captureData %>% as.data.frame(), model = "js", begin.time = 1)
ddl = make.design.data(proc)

initial.true <- list(Phi = c(`(Intercept)` = log(phi.true)), 
                     p = c(`(Intercept)` = log(p.true)),
                     pent = c(`(Intercept)` = log(pent.true)),
                     N = c(`(Intercept)` = log(Ns.true - sum(captureData$freq))))
initial.cjs <- list(Phi = c(`(Intercept)` = mod.cjs$results$beta$Phi), 
                     p = c(`(Intercept)` = mod.cjs$results$beta$p),
                     pent = c(`(Intercept)` = log(mean(u/udot))),
                     N = c(`(Intercept)` = log(udot)))
initial.rand <- lapply(initial.true, function(x) x*runif(1, 0.8, 1.2))
mod = crm(proc, ddl, hessian = TRUE, 
          initial = initial.cjs)



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

