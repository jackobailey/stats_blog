set.seed(43)
n.sims <- 100
library(mirt)
library(mitools)
library(tidyverse)
alpha <- 1
beta_1 <- 1
beta_2 <- 1
n <- 2000
n.cats <- 2
beta_cat <- 0.5
x1.x2.beta <- 0.5

include.pv <- TRUE
sim.out <- as.list(rep(NA, n.sims)) 

n.years <- 5
x1.items.per.year <- 4
total.items <- 10
n.per.survey <- 2000
max.x <- 1.5
max.x1 <- max.x
max.x2 <- max.x


x1.means <- runif(n= (n.years+1), min = -max.x1, max = max.x1)
x2.means <- runif(n= (n.years+1), min = -max.x2, max = max.x2)
item.subsets <- t(replicate(n.years, (1:total.items) %in% 
                              sample(1:total.items, size = x1.items.per.year, replace = F)))
item.subsets <- rbind(item.subsets, TRUE)

beta1.by.year <- runif(n= (n.years+1), min = -1, max = 1)
beta2.by.year <- runif(n= (n.years+1), min = -1, max = 1)

surveys <- c(LETTERS[1:n.years], "rosetta")

names(beta1.by.year) <- names(beta2.by.year) <- rownames(item.subsets) <- 
  names(x1.means) <- names(x2.means) <- 
  surveys
# generate x1 by year


survey.ids <- inverse.rle(list(values = surveys, lengths = rep(n.per.survey, length(x1.means))))
 
x1 <- rnorm(mean = x1.means[survey.ids], sd = 1, n = length(survey.ids))
x2 <- rnorm(mean = x2.means[survey.ids], sd = 1, n = length(survey.ids)) + x1.x2.beta * x1
x1 <- scale(x1)[, ]
x2 <- scale(x2)[, ]

y <- (beta1.by.year[survey.ids] * x1)  + (beta2.by.year[survey.ids] * x2) + rnorm(length(x1))

item.parameters <- cbind(a = runif(min = 0.5, max = 2, total.items), b = runif(min = -2, max = 2, total.items))

generateItemResponses <- function(item) {
  rbinom(prob = plogis(item.parameters[item, "a"] * (x1 - item.parameters[item, "b"])), size = 1, n = length(x1))
}
all.responses <- sapply(1:total.items, generateItemResponses)
all.responses[!item.subsets[survey.ids, ]] <- NA
all.responses <- data.frame(all.responses)

library(mirt)
x1.mod <- mirt(data = data.frame(all.responses), model = 1)
x1.hat <- fscores(x1.mod, full.scores = TRUE, full.scores.SE = TRUE)

marginal_rxx2 <- function (mod, density = dnorm, var_theta = 1, which.items = 1:extract.mirt(mod, "nitems"), ...) 
{
  stopifnot(extract.mirt(mod, "nfact") == 1L)
  stopifnot(is(mod, "SingleGroupClass"))
  fn <- function(theta, mod, den, which.items, ...) {
    TI <- testinfo(mod, matrix(theta), which.items = which.items)
    TI/(TI + var_theta) * den(theta, ...)
  }
  integrate(fn, lower = -Inf, upper = Inf, mod = mod, den = density, which.items,
            ...)$value
}

# sd(x1[surveys=="rosetta"])
sd(x1)

marginal.reliabilities <- sapply(surveys, function(x) marginal_rxx2(mod = x1.mod, which.items = which(item.subsets[surveys==x, ])) )
empirical.reliabilities <- sapply(surveys, function(x) empirical_rxx(x1.hat[survey.ids==x, ]))
true.reliabilities <- sapply(surveys, function(x) cor(cbind(x1[survey.ids==x], x1.hat[survey.ids==x, ]))[1,2])^2
names(marginal.reliabilities) <- names(true.reliabilities) <- names(empirical.reliabilities) <- surveys
library(eivtools)

data <-  data.frame(survey = survey.ids,
                    x1 = x1.hat[, 1], x2 = x2, y = y)

rhoToS <- function(rho) {
  s <- sqrt((1 / rho)^2 - 1)  
  return(s)
}

eivIRTForSurvey <- function(current.survey, reliabilities) {
  reliability <- reliabilities[current.survey]
  data.temp <- data[data$survey==current.survey, ]  
  data.temp$x1 <- data.temp$x1 / sd(data.temp$x1) * sqrt((1+rhoToS(sqrt(reliability))^2))
  eirt.mod.temp <- eivreg(data = data.temp, 
                               formula = y~x1+x2, 
                               reliability = c(x1 = as.vector(reliability)))  
  eirt.mod.sum <- summary(eirt.mod.temp)
  coefs.to.extract <- c("x1", "x2")
  out <- data.frame(estimator = "EIV",
                    measurement = "IRT",
                    var = coefs.to.extract, 
             coef.est = eirt.mod.sum$coefficients[coefs.to.extract, c("Estimate")], 
             se = eirt.mod.sum$coefficients[coefs.to.extract, c("Std. Error")], 
             coef.true = c(beta1.by.year[current.survey],  beta2.by.year[current.survey]),
             reliability.est = c(reliability, 1), 
             reliablity.true = c((cor(data.temp$x1, x1[data$survey==current.survey]))^2, 1), 
             survey = current.survey,
             x1x2cor = cor(x2[data$survey==current.survey], x1[data$survey==current.survey]))
  return(out)
}

all.eirt.for.rxx <- function(reliabilities) {
  all.eirt.results <- lapply(surveys, eivIRTForSurvey, reliabilities = reliabilities)
  eirt.comb <- do.call(rbind, all.eirt.results)
  eirt.comb$uci <- eirt.comb$coef.est + 1.96 * eirt.comb$se
  eirt.comb$lci <- eirt.comb$coef.est - 1.96 * eirt.comb$se
  
  eirt.comb$coef.true.mag <- (eirt.comb$coef.true * sign(eirt.comb$coef.true))
  eirt.comb$coef.est.mag <- (eirt.comb$coef.est * sign(eirt.comb$coef.true))
  
  eirt.comb$resid.coef.mag <- eirt.comb$coef.est.mag - eirt.comb$coef.true.mag
  eirt.comb$resid.coef <- eirt.comb$coef.est - eirt.comb$coef.true
  
  return(eirt.comb)
}

eirt.comb.true.rxx <- all.eirt.for.rxx(reliabilities = true.reliabilities)
eirt.comb.emp.rxx <-  all.eirt.for.rxx(reliabilities = empirical.reliabilities)
eirt.comb.marg.rxx <-  all.eirt.for.rxx(reliabilities = marginal.reliabilities)
eirt.comb.true.rxx$reliability.method <- 'true'
eirt.comb.emp.rxx$reliability.method <- 'empirical'
eirt.comb.marg.rxx$reliability.method <- 'marginal'
library(ggplot2)

eirt.comb<- rbind(eirt.comb.true.rxx, eirt.comb.emp.rxx, eirt.comb.marg.rxx)

eirt.comb$in95ci <- eirt.comb$coef.est >eirt.comb$lci & eirt.comb$coef.est <eirt.comb$uci
eirt.comb$num.indicators <- factor(ifelse(eirt.comb$survey=="rosetta", 10, 4))
# eirt.comb.marg.rxx

library(gridExtra)
library(mellonMisc)

eirt.validation.plot <- ggplot(eirt.comb, aes(x = coef.true, y = coef.est, colour = var,
                                                            shape = num.indicators)) +
  geom_point() + xlab("True coefficient") +
  facet_grid(cols = vars(reliability.method)) + ylab("Estimated coefficient") + 
  ylim(c(-1,1)) + xlim(c(-1,1)) + 
  geom_errorbar(aes(ymin = coef.true,ymax = coef.est), size=0.5) + 
  geom_abline(size = 0.25, linetype = 2) + 
  theme_bw()

eirt.validation.plot

within.sd <- mean(tapply(x1, survey.ids, sd)) 
between.sd <- sd(tapply(x1, survey.ids, mean))
# within.sd.prop = 
(within.sd / (within.sd + between.sd))