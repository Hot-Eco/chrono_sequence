dat.predict <- cbind(
  DAT.predict,
  flag.site = 0,
  flag.ms = 1,
  fSiteId = levels(DAT$fSiteId)[1]
)

Nc <- ncol(dat.predict)
Nsim <- 10000
Nplot <- 25

sim <- postSim(Nsim, mtotal.site, dat.predict, join = FALSE)
sim <- exp(sim)

# standard interval
fit <- do_predict(mtotal.site, dat.predict, "percentC", exp)


dat.simlines <- cbind(dat.predict, sim) %>% 
  select(c(1:Nc, sample((Nc+1):(Nsim+Nc), Nplot))) %>%
  gather(rep, percentC, -(1:ncol(dat.predict)))



# estimated coverage of standard point-wise interval
allIn <- function(lwr, upr, isim) all(isim > lwr & isim < upr)
coverage <- apply(sim, 2, function(isim) allIn(fit$lwr, fit$upr, isim))
cat("Coverage of standard interval =", sum(coverage) / ncol(sim))


# simulataneous interval
Cg <- predict(mtotal.site, dat.predict, type = "lpmatrix")
Vb <- vcov(mtotal.site)
BUdiff <- rmvn(10000, mu = rep(0, nrow(Vb)), vc = Vb)
simDev <- Cg %*% t(BUdiff)
se.fit <- predict(mtotal.site, dat.predict, se.fit = TRUE)$se.fit
absDev <- abs(sweep(simDev, 1, se.fit, "/"))
masd <- apply(absDev, 2, max)
crit <- quantile(masd, prob = 0.95, type = 8)

fit.sim <- cbind(fit, se.fit) %>%
  mutate(lwr = exp(log(percentC) - crit * se.fit),
         upr = exp(log(percentC) + crit * se.fit))


# estimated coverage of simultaneous interval
coverage <- apply(sim, 2, function(isim) allIn(fit.sim$lwr, fit.sim$upr, isim))
cat("Coverage of simultaneous interval =", sum(coverage) / ncol(sim))

ggplot(data = DAT, aes(x = tsf.years)) +
  geom_point(aes(y = percentC, colour = fireType)) + 
  FireTypeColour + 
  facet_grid(depth ~ microsite) + 
  theme_bw() + theme(legend.position = "bottom") +

  # standard interval
  geom_ribbon(data = fit,
              aes(ymin = lwr, ymax = upr, fill = fireType), alpha = 0.5) +
  
  # simultaneous interval
  geom_ribbon(data = fit.sim, 
              aes(ymin = lwr, ymax = upr, fill = fireType), alpha = 0.3) +

  # posterior simulations
  geom_line(data = dat.simlines,
            aes(y = percentC, group = interaction(rep, depth, microsite, fireType)), 
            alpha = 0.2)
