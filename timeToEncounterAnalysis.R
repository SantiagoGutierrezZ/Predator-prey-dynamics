library(glmmTMB)
library(DHARMa)
library(effects)
library(purrr)

#########
load("timeToEconuterAnalysis.RData")

##################
#
# comparation rabbit vs predators
# 
###################################

rabbit_badger$time_diff <- as.numeric(rabbit_badger$time_diff)
mod1.rabbit_badger <- glmmTMB(time_diff ~ interval, data = rabbit_badger, family = lognormal)
simulateResiduals(mod1.rabbit_badger, plot = TRUE)
# plot(allEffects(mod1.rabbit_badger))
summary(mod1.rabbit_badger)

rabbit_fox$time_diff <- as.numeric(rabbit_fox$time_diff)
mod1.rabbit_fox <- glmmTMB(time_diff ~ interval, data = rabbit_fox, family = lognormal)
simulateResiduals(mod1.rabbit_fox, plot = TRUE)
# plot(allEffects(mod1.rabbit_fox))
summary(mod1.rabbit_fox)

rabbit_genet$time_diff <- as.numeric(rabbit_genet$time_diff)
mod1.rabbit_genet <- glmmTMB(time_diff ~ interval, data = rabbit_genet, family = lognormal)
simulateResiduals(mod1.rabbit_genet, plot = TRUE)
# plot(allEffects(mod1.rabbit_genet))
summary(mod1.rabbit_genet)
mod2.rabbit_genet <- glmmTMB(time_diff ~ interval, data = rabbit_genet, family = Gamma(link = "log"))
simulateResiduals(mod2.rabbit_genet)
summary(mod2.rabbit_genet)

rabbit_mongoose$time_diff <- as.numeric(rabbit_mongoose$time_diff)
mod1.rabbit_mongoose <- glmmTMB(time_diff ~ interval, data = rabbit_mongoose, family = lognormal)
simulateResiduals(mod1.rabbit_mongoose, plot = TRUE)
# plot(allEffects(mod1.rabbit_mongoose))
summary(mod1.rabbit_mongoose)

rabbit_lynx$time_diff <- as.numeric(rabbit_lynx$time_diff)
mod1.rabbit_lynx <- glmmTMB(time_diff ~ interval, data = rabbit_lynx, family = lognormal)
simulateResiduals(mod1.rabbit_lynx, plot = TRUE)
# plot(allEffects(mod1.rabbit_lynx))
summary(mod1.rabbit_lynx)

##################
#
# comparation hare vs predators
# 
###################################

hare_badger$time_diff <- as.numeric(hare_badger$time_diff)
mod1.hare_badger <- glmmTMB(time_diff ~ interval, data = hare_badger, family = lognormal)
simulateResiduals(mod1.hare_badger, plot = TRUE)
# plot(allEffects(mod1.hare_badger))
summary(mod1.hare_badger)

hare_fox$time_diff <- as.numeric(hare_fox$time_diff)
hare_fox$time_diff <- ifelse(hare_fox$time_diff == 0, 1, hare_fox$time_diff)
mod1.hare_fox <- glmmTMB(time_diff ~ interval, data = hare_fox, family = lognormal)
simulateResiduals(mod1.hare_fox, plot = TRUE)
# plot(allEffects(mod1.hare_fox))
summary(mod1.hare_fox)
mod2.hare_fox <- glmmTMB(time_diff ~ interval, data = hare_fox, family = Gamma(link = "log"))
simulateResiduals(mod2.hare_fox, plot = TRUE)
summary(hare_fox$time_diff)


hare_genet$time_diff <- as.numeric(hare_genet$time_diff)
mod1.hare_genet <- glmmTMB(time_diff ~ interval, data = hare_genet, family = lognormal)
simulateResiduals(mod1.hare_genet, plot = TRUE)
# plot(allEffects(mod1.hare_genet))
summary(mod1.hare_genet)

hare_mongoose$time_diff <- as.numeric(hare_mongoose$time_diff)
mod1.hare_mongoose <- glmmTMB(time_diff ~ interval, data = hare_mongoose, family = lognormal)
simulateResiduals(mod1.hare_mongoose, plot = TRUE)
# plot(allEffects(mod1.hare_mongoose))
summary(mod1.hare_mongoose)

hare_lynx$time_diff <- as.numeric(hare_lynx$time_diff)
mod1.hare_lynx <- glmmTMB(time_diff ~ interval, data = hare_lynx, family = lognormal)
simulateResiduals(mod1.hare_lynx, plot = TRUE)
# plot(allEffects(mod1.hare_lynx))
summary(mod1.hare_lynx)

