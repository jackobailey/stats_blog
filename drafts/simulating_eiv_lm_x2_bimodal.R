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


for(sim in which(is.na(sim.out))) {
  print(sim)
  cat <- as.numeric(sample(1:n.cats, n, replace = T)==1)
  x1 <- rnorm(n) 
  x1 <- scale(x1)[, ]
  x2 <- rnorm(n)  + x1.x2.beta*x1 + beta_cat * cat
  x2 <- scale(x2)[, ]
  
  
  y <- alpha + beta_1 * x1 + beta_2 * x2  + rnorm(n, sd = 1)
  
  x1a <- rbinom(prob = plogis(x1 + rnorm(n, sd=1)+1), size = 1, n=n)
  x1b <- rbinom(prob = plogis(x1 *0.5+ rnorm(n, sd=2)-1), size = 1, n=n)
  x1c <- rbinom(prob = plogis(x1 *2), size = 1, n=n)
  x1d <- rbinom(prob = plogis(x1 * 1.5 + rnorm(n, sd=2)+0.25), size = 1, n=n)
  x1e <- rbinom(prob = plogis(x1 *0.1 + rnorm(n, sd=1)), size = 1, n=n)
  x1f <- rbinom(prob = plogis(x1 *1 + rnorm(n, sd=0.5)), size = 1, n=n)
  x1g <- rbinom(prob = plogis(x1 *1 + rnorm(n, sd=0.5)), size = 1, n=n)
  x1h <- rbinom(prob = plogis(x1 *2 + rnorm(n, sd=0.5)), size = 1, n=n)
  
  x2a <- rbinom(prob = plogis(x2 + rnorm(n, sd=1)+1), size = 1, n=n)
  x2b <- rbinom(prob = plogis(x2 *0.5+ rnorm(n, sd=2)-1), size = 1, n=n)
  x2c <- rbinom(prob = plogis(x2 *2), size = 1, n=n)
  x1dat <- data.frame(x1a,x1b, x1c, x1d, x1e, x1f, x1g, x1h)
  x2dat <- data.frame(x2a,x2b, x2c)
  x1.mod <- mirt(data = x1dat,
                 model = 1, verbose = FALSE)
  x2.mod <- mirt(data = x2dat, 
                 model = 1, verbose = FALSE)
  x1.hat <- fscores(x1.mod, full.scores= TRUE, full.scores.SE = TRUE, method = "EAP")
  x2.hat <- fscores(x2.mod, full.scores= TRUE, full.scores.SE = TRUE, method = "EAP")
  reliabilities <- c(x1 = as.vector(empirical_rxx(x1.hat)),
                     x2 = as.vector(empirical_rxx(x2.hat)))
  
  reliabilities2 <- c(x1 = as.vector(marginal_rxx(x1.mod)),
                      x2 = as.vector(marginal_rxx(x2.mod)))
 
  
  rhoToS <- function(rho) {
    s <- sqrt((1 / rho)^2 - 1)  
    return(s)
  }
  
  
  rho.x1.est <- sqrt(empirical_rxx(x1.hat))
  rho.x2.est <- sqrt(empirical_rxx(x2.hat))
  
  rho.x1.est.m <- sqrt(marginal_rxx(x1.mod))
  rho.x2.est.m <- sqrt(marginal_rxx(x2.mod))
  
  rho.x1.true <- cor(cbind(x1.hat[, 1], x1))[1,2]
  rho.x2.true <- cor(cbind(x2.hat[, 1], x2))[1,2]
  
  data.hat.empirical <- data.frame(y, x1 = (x1.hat[, 1]/ sd(x1.hat[, 1])) * sqrt((1+rhoToS(rho.x1.est)^2)),
                                   x2 = (x2.hat[, 1]/ sd(x2.hat[, 1])) *  sqrt (1+rhoToS(rho.x2.est)^2), 
                                   rxx1 = rho.x1.est^2,
                                   rxx2 = rho.x2.est^2,
                                   r.type = "empirical")
  
  data.hat.marginal <- data.frame(y, x1 = (x1.hat[, 1]/ sd(x1.hat[, 1])) * sqrt((1+rhoToS(rho.x1.est.m)^2)),
                                  x2 = (x2.hat[, 1]/ sd(x2.hat[, 1])) *  sqrt (1+rhoToS(rho.x2.est.m)^2), 
                                  rxx1 = rho.x1.est.m^2,
                                  rxx2 = rho.x2.est.m^2,
                                  r.type = "marginal")
  
  data.hat.real <- data.frame(y, x1 = (x1.hat[, 1]/ sd(x1.hat[, 1])) * sqrt((1+rhoToS(rho.x1.true)^2)),
                              x2 = (x2.hat[, 1]/ sd(x2.hat[, 1])) *  sqrt (1+rhoToS(rho.x2.true)^2), 
                              rxx1 = rho.x1.true^2,
                              rxx2 = rho.x2.true^2,
                              r.type = "real")
  data.hat.irt <- data.frame(y, x1 = x1.hat[, 1],
                             x2 = x2.hat[, 1], 
                             rxx1 = rho.x1.est.m^2,
                             rxx2 = rho.x2.est.m^2,
                             r.type = "marginal")
  data.classical <- data.frame(y, 
                               x1 = x1 + rnorm(n, sd = rhoToS(rho.x1.est)) ,
                               x2 = x2 + rnorm(n, sd = rhoToS(rho.x2.est)), 
                               rxx1 = rho.x1.est^2,
                               rxx2 = rho.x2.est^2,
                               r.type = "simulated")
  data.classical2 <- data.classical
  
  data.classical2$rxx1 <- cor(cbind(data.classical$x1, x1))[1,2]
  data.classical2$rxx2 <- cor(cbind(data.classical$x2, x2))[1,2]
  data.classical2$r.type <- "real"
  
  library(eivtools)
  
  datasets <- list(classical = data.classical,
                   classical.real = data.classical2,
                   irt = data.hat.irt, 
                   irt.cool.empirical = data.hat.empirical,
                   irt.cool.marginal = data.hat.marginal,
                   irt.cool.real = data.hat.real)
  
  temp.sim.res <- list()
  for(dat.name in names(datasets)) {
    lm.mod.temp <- lm(data = datasets[[dat.name]], 
                      formula = y~x1+x2)
    # using empirical rxx 
    eirt.mod.temp <- eivreg(data = datasets[[dat.name]], 
                            formula = y~x1+x2, 
                            reliability = c(x1 = datasets[[dat.name]]$rxx1[1], 
                                            x2 = datasets[[dat.name]]$rxx2[2]))
    
    ols.temp.out <- data.frame(var = names(lm.mod.temp$coefficients),
                               coefficient = lm.mod.temp$coefficients, 
                               se = summary(lm.mod.temp)$coefficients[, "Std. Error"] , 
                               error = dat.name, 
                               estimator = "OLS", 
                               r.method = "none",
                               reliability = NA)
    
    eiv.temp.out <- data.frame(var = names(lm.mod.temp$coefficients),
                               coefficient = eirt.mod.temp$coefficients, 
                               se = summary(eirt.mod.temp)$coefficients[, "Std. Error"], 
                               error = dat.name, 
                               estimator = "EIV", 
                               r.method = datasets[[dat.name]]$r.type[1],
                               reliability = eirt.mod.temp$reliability[names(lm.mod.temp$coefficients)])
    
    dat.temp.out <- rbind(ols.temp.out, eiv.temp.out)
    true.cors <- c(x1 = cor(cbind(x1, datasets[[dat.name]]$x1))[1,2], 
                   x2 = cor(cbind(x2, datasets[[dat.name]]$x2))[1,2])
    dat.temp.out$true.cor.latent <- true.cors[dat.temp.out$var]
    temp.sim.res[[dat.name]] <-   dat.temp.out
  }
  
  
  
  if(include.pv) {
    bothdat <- cbind(x1dat, x2dat)
    sim_data <- cbind(y, bothdat)
    mod.text <- paste0("F1 = 1-", ncol(x1dat), "\n", 
                       "F2 = ", (ncol(x1dat)+1), "-",
                       ncol(bothdat), "\n",
                       "COV = F1*F2")
    model <- mirt.model(mod.text)
    comb.mod <- mirt(data = bothdat, 
                     model = model, 
                     covdata = data.frame(y), 
                     formula = ~.)
    
    pv_ycov <- fscores(comb.mod,  plausible.draws = 50)
    
    
    draw_merger_ycov <- function(draw_num = 1){
      sim_data_draw <- sim_data %>% 
        mutate(x1 = pv_ycov[[draw_num]][,1] / sd(pv_ycov[[draw_num]][,1]), 
               x2 = pv_ycov[[draw_num]][,2] / sd(pv_ycov[[draw_num]][,2]))
      return(sim_data_draw)
    }
    pv_imputes_ycov <- map(1:length(pv_ycov), draw_merger_ycov) %>%
      imputationList()
    mi_mods_ycov <- with(pv_imputes_ycov,  lm(y ~ x1 + x2))
    pv.mod <- MIcombine(mi_mods_ycov)
    
    pv.temp.out <- data.frame(var = names(lm.mod.temp$coefficients),
                               coefficient = pv.mod$coefficients, 
                               se = diag(sqrt(pv.mod$variance)), 
                               error = "PV: covariates + covariance", 
                               estimator = "MI", 
                               r.method = "None",
                               reliability = NA)
    pv.temp.out$true.cor.latent <- true.cors[pv.temp.out$var]
    temp.sim.res[["pvfull"]] <- pv.temp.out
  }
  
  
  temp.sim.res <- do.call(rbind, temp.sim.res)
  temp.sim.res <- temp.sim.res[temp.sim.res$var %in% c("x1", "x2"), ]
  temp.sim.res$uci <- temp.sim.res$coefficient + 1.96 * temp.sim.res$se 
  temp.sim.res$lci <- temp.sim.res$coefficient - 1.96 * temp.sim.res$se 
  
  temp.sim.res$true <- NA
  temp.sim.res$true[temp.sim.res$var=="x1"] <- beta_1
  temp.sim.res$true[temp.sim.res$var=="x2"] <- beta_2
  
  temp.sim.res$in95ci <- temp.sim.res$uci >temp.sim.res$true & temp.sim.res$lci <temp.sim.res$true
  temp.sim.res$resid  <- temp.sim.res$coefficient - temp.sim.res$true
  temp.sim.res$rx1x2 <- cor(cbind(x1,x2))[1,2]
  temp.sim.res$sim <- sim
  
  temp.sim.res <- cbind(beta_1 = beta_1, beta_2 = beta_2, 
        x1.x2.beta = x1.x2.beta, 
        beta_cat = beta_cat,
        n.cats = n.cats, 
        n = n, 
        alpha = alpha, temp.sim.res)
 
  sim.out[[sim]] <- temp.sim.res
}


sim.all <- do.call(rbind, sim.out[!is.na(sim.out)])

sim.all$abs.resid <- abs(sim.all$resid)
sim.all$resid2 <- sim.all$resid^2
sim.all$true.reliability <- sim.all$true.cor.latent^2

sim.all$resid.reliability <- sim.all$reliability - sim.all$true.reliability
sim.all$abs.resid.reliability <- abs(sim.all$reliability - sim.all$true.reliability)
sim.all$sig <- sim.all$lci>0
sim.all$ci.width <- sim.all$uci - sim.all$lci

sim.summary <- aggregate(sim.all[, c("in95ci", "resid", "abs.resid", "coefficient", 
                                     "reliability", "resid.reliability", 
                                     "abs.resid.reliability",
                                     "resid2", "ci.width", "sig")], sim.all[, c("var", "error", "estimator", "r.method")], 
                         FUN = mean)
sim.summary$rmse <- sqrt(sim.summary$resid2)
sim.summary$resid2 <- NULL

sim.summary$ratio <- NA
sim.summary$ratio[sim.summary$var=="x1"] <- sim.summary$coefficient[sim.summary$var=="x1"] / sim.summary$coefficient[sim.summary$var=="x2"]
sim.summary <- sim.summary[!(sim.summary$error=="classical" & sim.summary$r.method %in% c("marginal", "empirical")), ]
sim.summary$n.sims <- sum(!is.na(sim.out))
library(dplyr)
sim.summary <- sim.summary %>% arrange(estimator, desc(error), r.method, var)
sim.summary 
save(sim.summary, file = "resources/irt_eiv_simsx2bimodal100.rda")
