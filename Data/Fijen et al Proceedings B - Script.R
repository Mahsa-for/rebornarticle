
#Libraries#######################################
library(stringr)
library(multcomp)
library(lme4)
library(lubridate)
library(ggplot2)
library(visreg)
library(car)
library(gridExtra)
library(lemon) #for extracting the legend in ggplot
library(effects)
library(ggsignif)

se <- function(x) sd(x, na.rm=TRUE) /  sqrt(length(x[!is.na(x)])) 
Sys.setlocale("LC_TIME", "English_United States")

#Load Files#################################
flow.cov <- read.csv('flower.cover.period.csv', header = T)
flow.cov.av <- read.csv('flower.cover.average.csv', header = T)
sp.boot <- read.csv("bootstrapped.species.richness.csv", header = T)
tot.dens <- read.csv("total.densities.csv", header = T)
mean.dens <- read.csv("mean.densities.csv", header = T)

### Landscape species vs leek species #####
land.species <- aggregate(species ~ field, sp.boot[sp.boot$leek.land == 'land' &
                                                     sp.boot$period == 'both',], sum)
leek.species <- sp.boot[sp.boot$leek.land == 'leek',]
leek.species$land.species <- land.species$species[match(leek.species$field, land.species$field)]

### First the general relationship
general.leek.species <- aggregate(species ~ field, leek.species, sum) # for the general relationship
general.leek.species$land.species <- leek.species$land.species[match(general.leek.species$field, leek.species$field)]
m1 <- lm(species ~ land.species, general.leek.species)
summary(m1)

### Then for each functional group
m4 <- lmer(species ~land.species*dominant + (1|field), leek.species) 
summary(m4)
m4b <- lmer(species ~land.species + dominant + (1|field), leek.species)
summary(m4b)
anova(m4,m4b)
m4c <- lmer(species ~ dominant + (1|field), leek.species)
anova(m4c,m4b)
m4d <- lmer(species ~ land.species + (1|field), leek.species)
anova(m4d,m4b)

# Partial residuals
visreg.m4b <- visreg(m4, xvar = 'land.species', by = 'dominant', band = F, overlay = T, overlay = T, 
                     plot = F, rug = F)
# Figure 2a
p1 <- ggplot(visreg.m4b$fit, aes(land.species, visregFit, colour = factor(dominant), fill=factor(dominant))) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.2, 
              colour="grey50", linetype=1, size=0.2) +
  geom_line(size = 1.5, aes(linetype = dominant)) + 
  geom_point(data=visreg.m4b$res, aes(land.species, visregRes, colour = dominant, shape = dominant), size=3) +
  scale_fill_grey(start=0.5, end=0.5) +
  labs(linetype="Weight", fill="Weight", tag = 'A') +
  theme_bw() + labs(x = "Local pollinator species richness", y = "Crop pollinator species richness") +
  theme(axis.title = element_text(size=18),
        title = element_text(size = 18),
        axis.text = element_text(size=16, colour="#000000"),
        legend.text=element_text(size=14), legend.title=element_text(size=14),
        strip.background =element_rect(fill="white"), 
        strip.text = element_text(size=14),
        legend.position = "none") +
  scale_color_manual(name = 'Functional group', values=c('#354ea7', '#CF0000')) +
  scale_shape_manual(name = 'Functional group',
                     values = c(16,17,15,8,7)) +
  scale_linetype_manual(name = 'Functional group',
                        values=c("twodash", "dashed", "solid"))
p1


### leek density ~  landscape species
leek.species$dens.leek <- tot.dens$dens.leek[match(paste(leek.species$field, leek.species$dominant),
                                                   paste(tot.dens$field, tot.dens$dominant))]
leek.species$flow.cov <- flow.cov.av$flow.cov[match(leek.species$field, flow.cov.av$field)]

# transform data
leek.species$dens.leek.log <- log10(leek.species$dens.leek)

### First the general relationship
general.leek.dens <- aggregate(dens.leek ~ field, leek.species, sum)
general.leek.dens$dens.leek.log <- log10(general.leek.dens$dens.leek)
general.leek.dens$land.species <- leek.species$land.species[match(general.leek.dens$field, leek.species$field)]
m1 <- lm(dens.leek.log ~ land.species, general.leek.dens)
summary(m1)

### Then per functional group
m5 <- lmer(dens.leek.log ~land.species*dominant + (1|field), leek.species)
summary(m5)
m5b <- lmer(dens.leek.log ~land.species+dominant+ (1|field), leek.species)
summary(m5b)
anova(m5b,m5)

m5c <- lmer(dens.leek.log ~ land.species + (1|field), leek.species)
summary(m5c)
m5d <- lmer(dens.leek.log ~dominant+ (1|field), leek.species)
summary(m5d)
anova(m5b,m5d )
anova(m5b,m5c)

effects.m5b <- as.data.frame(Effect(c("land.species", "dominant"), m5b))

# partial residuals
visreg.m5b <- visreg(m5b, xvar = 'land.species', by = "dominant", band = F, overlay = T, overlay = T, 
                     plot = F, rug = F, trans = function(x) 10^x)

# Figure 2b
p2 <- ggplot(visreg.m5b$fit, aes(land.species, visregFit, colour = factor(dominant), fill=factor(dominant))) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.2, 
              colour="grey50", linetype=1, size=0.2) +
  geom_line(size = 1.5, aes(linetype = dominant)) + 
  geom_point(data=visreg.m5b$res, aes(land.species, visregRes, colour = dominant, shape = dominant), size=3) +
  scale_fill_grey(start=0.5, end=0.5) +
  labs(linetype="Weight", fill="Weight", tag = 'B') +
  theme_bw() + labs(x = "Local pollinator species richness", y = "Crop pollinator abundance") +
  theme(axis.title = element_text(size=18),
        title = element_text(size = 18),
        axis.text = element_text(size=16, colour="#000000"),
        legend.text=element_text(size=14), legend.title=element_text(size=14),
        strip.background =element_rect(fill="white"), 
        strip.text = element_text(size=14),
        legend.position = "none") +
  scale_color_manual(name = 'Functional group', values=c('#354ea7', '#CF0000'))  +
  scale_shape_manual(name = 'Functional group',
                     values = c(16,17,15,8,7)) +
  scale_linetype_manual(name = 'Functional group',
                        values=c("twodash", "dashed", "solid"))
p2

# Legend for figure 2
p2b <- ggplot(visreg.m5b$fit, aes(land.species, visregFit, colour = factor(dominant), fill=factor(dominant))) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.2, 
              colour="grey50", linetype=1, size=0.2) +
  geom_line(size = 1.5, aes(linetype = dominant)) + 
  geom_point(data=visreg.m5b$res, aes(land.species, visregRes, colour = dominant, shape = dominant), size=4.5) +
  scale_fill_grey(start=0.5, end=0.5) +
  labs(linetype="Weight", fill="Weight") +
  theme_bw() + labs(x = "Local pollinator species richness", y = "Crop pollinator density") +
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=16, colour="#000000"),
        legend.text=element_text(size=14), legend.title=element_text(size=14),
        legend.position= 'bottom',
        strip.background =element_rect(fill="white"), 
        strip.text = element_text(size=14)) +
  scale_color_manual(name = 'Functional group', values=c('#354ea7', '#CF0000'),
                     labels = c("Dominant crop pollinators", "Opportunistic crop pollinators"),
                     guide = guide_legend(keywidth = unit(2, 'cm')))  +
  scale_fill_manual(name = "Functional group", values = c('grey50','grey50'),
                    labels = c("Dominant crop pollinators", "Opportunistic crop pollinators")) +
  scale_shape_manual(name = 'Functional group', values = c(16,17,15,8,7),
                     labels = c("Dominant crop pollinators", "Opportunistic crop pollinators")) +
  scale_linetype_manual(name = 'Functional group',
                        values=c("twodash", "dashed", "solid"),
                        labels = c("Dominant crop pollinators", "Opportunistic crop pollinators"))

legend <- g_legend(p2b)

# Figure 2
grid.arrange(p1,p2, legend, ncol = 2,
             layout_matrix = rbind(c(1,2), c(3)),
             heights = c(1,0.15))


### Densities before and during in landscape, and in the crop ####
### statistics on densities landscape 
m1 <- lmer(log10(dens.land) ~ period*dominant + (1|field), mean.dens)
summary(m1)

mean.dens$interaction <- interaction(mean.dens$period, mean.dens$dominant)
m1b <- lmer(log10(dens.land) ~ interaction + (1|field), mean.dens)
summary(m1b)
summary(glht(m1b, linfct = mcp(interaction = "Tukey")))

# backtransform and tidying
effects.m1b <- as.data.frame(Effect("interaction", m1b))
effects.m1b <- effects.m1b[c(1,2,5,6,3,4),]
mean.land.se <- data.frame(dominant = c('dominant', 'dominant', 'opportunistic', 'opportunistic',
                                        'non.crop', 'non.crop'),
                           period = c('before', 'during'))
mean.land.se$fit <- c(10^(effects.m1b$fit))
mean.land.se$lower <- c(10^(effects.m1b$lower))
mean.land.se$upper <- c(10^(effects.m1b$upper))

mean.land.se$dominant <- factor(mean.land.se$dominant, levels = c('dominant', 'opportunistic', "non.crop"),
                                labels = c('Dominant crop pollinators', 'Opportunistic crop pollinators', 'Non-crop pollinators'))

mean.land.se$period <- factor(mean.land.se$period, levels = c('before', 'during'),
                              labels = c('Before', 'During'))

# Figure 5a
annotation_df <- data.frame(dominant=c('Dominant crop pollinators', 'Opportunistic crop pollinators', 'Non-crop pollinators'), 
                            start=c("Before", "Before", "Before"), 
                            end=c("During","During","During"),
                            y=c(45),
                            label= c("n.s.", "p = 0.001", "n.s."))

limits <- aes(ymax = mean.land.se$upper,
              ymin = mean.land.se$lower)
d1 <- ggplot(mean.land.se, aes(x = period, y = fit)) +
  geom_bar(stat = "identity" ) + geom_errorbar(limits, size=0.5, width=.25) +
  facet_wrap(~dominant) + theme_bw() + geom_hline(yintercept = 0) +
  labs(x = "Period", y = "Abundance per \n semi-natural habitat transect", tag = "A") +
  scale_y_continuous(limits = c(0,48))+  
  geom_signif(data=annotation_df, aes(xmin=start, xmax=end, annotations=label, y_position=y),
              tip_length = 0.00, vjust = -0.2, textsize = 4,
              manual=TRUE) +
  theme(axis.title.x = element_text(size=14),
        title = element_text(size = 18),
        axis.title.y = element_text(size=12),
        axis.text = element_text(size=12, colour="#000000"),
        legend.text=element_text(size=14),
        strip.background =element_rect(fill="white"),
        strip.text = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 
d1

###  average density all transects
#statistics leek land densities
t <- aggregate(dens.land ~ field + dominant, mean.dens, mean)
t$leek.land <- 'land'
names(t) <- c('field', 'dominant', 'dens', 'leek.land')
t2 <- aggregate(dens.leek ~ field + dominant, mean.dens, mean)
t2$leek.land <- 'leek'
names(t2) <- c('field', 'dominant', 'dens', 'leek.land')
temporary <- rbind(t, t2)

temporary$interaction <- interaction(temporary$dominant, temporary$leek.land)
temporary$flow.cov <- 0
temporary$flow.cov[temporary$leek.land == 'land'] <- flow.cov.av$flow.cov[match(temporary$field[temporary$leek.land == 'land'], flow.cov.av$field)]
temporary$dens <-log10(temporary$dens+1)

m1 <- lmer(dens ~ leek.land*dominant + (1|field), temporary)
summary(m1)

m1b <- lmer(dens ~ interaction + (1|field), data = temporary, REML = FALSE)
summary(m1b)
summary(glht(m1b, linfct = mcp(interaction = "Tukey")))

#Backtransform and tidying
effects.m1b <- as.data.frame(Effect("interaction", m1b))
effects.m1b <- effects.m1b[c(1,4,3,5,2),]

dens.leek.land.all <- data.frame(dominant = c('dominant', 'dominant', 'opportunistic', 'opportunistic',
                                              'non.crop', 'non.crop'),
                                 leek.land = c('land', 'leek'))
dens.leek.land.all$fit <- c(10^(effects.m1b$fit)-1,NA)
dens.leek.land.all$lower <- c(10^(effects.m1b$lower)-1,NA)
dens.leek.land.all$upper <- c(10^(effects.m1b$upper)-1,NA)

dens.leek.land.all$dominant <- factor(dens.leek.land.all$dominant, levels = c('dominant', 'opportunistic', "non.crop"),
                                      labels = c('Dominant crop pollinators', 'Opportunistic crop pollinators', 'Non-crop pollinators'))

dens.leek.land.all$leek.land <- factor(dens.leek.land.all$leek.land, levels = c('land', 'leek'),
                                       labels = c('Semi-natural \n habitat', 'Crop'))

# Figure 5b
annotation_df <- data.frame(dominant=c('Dominant crop pollinators', 'Opportunistic crop pollinators'), 
                            start=c("Semi-natural \n habitat", "Semi-natural \n habitat"), 
                            end=c("Crop","Crop"),
                            y=c(80),
                            label= c("p < 0.001", "n.s."))
limits <- aes(ymax = dens.leek.land.all$upper,
              ymin = dens.leek.land.all$lower)
d2 <- ggplot(dens.leek.land.all, aes(x = leek.land, y = fit)) +
  geom_bar(stat = "identity" ) + geom_errorbar(limits, size=0.5, width=.25) +
  facet_wrap(~dominant) + theme_bw() + geom_hline(yintercept = 0) +
  labs(x = "Transect location", y = " \n Abundance per transect", tag = "B") +
  scale_y_continuous(limits = c(0,85))+  
  geom_signif(data=annotation_df, aes(xmin=start, xmax=end, annotations=label, y_position=y),
              tip_length = 0.00, vjust = -0.2, textsize = 4,
              manual=TRUE) +
  theme(axis.title = element_text(size=14),
        title = element_text(size = 18),
        axis.title.y = element_text(size=12),
        axis.text = element_text(size=12, colour="#000000"),
        legend.text=element_text(size=14),
        strip.background =element_rect(fill="white"),
        strip.text = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 
d2
# Figure 5
grid.arrange(d1,d2)


### Pollinator densities and snh cover ####
dens.landscape <- mean.dens

# First general relation of density with snh
general.dens <- aggregate(dens.land ~ field, dens.landscape, sum)
t <- aggregate(dens.leek ~ field, dens.landscape, sum)
general.dens$dens.leek <- t$dens.leek
general.dens$snh <- dens.landscape$snh[match(general.dens$field, dens.landscape$field)]
general.dens$flow.cov <- dens.landscape$flow.cov[match(general.dens$field, dens.landscape$field)]

m1 <- lm(log10(dens.land) ~ snh + flow.cov, general.dens)
summary(m1)

#transform densities
dens.landscape$dens.land <- log10(dens.landscape$dens.land)
dens.landscape$dens.leek <- log10(dens.landscape$dens.leek)

m1 <- lm(dens.land ~ flow.cov + snh, 
         data = dens.landscape[dens.landscape$dominant == 'dominant' & 
                                 dens.landscape$period == 'before',])
m2 <- lm(dens.land ~ flow.cov + snh, 
         data = dens.landscape[dens.landscape$dominant == 'opportunistic' & 
                                 dens.landscape$period == 'before',])
m3 <- lm(dens.land ~ flow.cov + snh, 
         data = dens.landscape[dens.landscape$dominant == 'non.crop' & 
                                 dens.landscape$period == 'before',])
m4 <- lm(dens.land ~ flow.cov + snh, 
         data = dens.landscape[dens.landscape$dominant == 'dominant' & 
                                 dens.landscape$period == 'during',])
m5 <- lm(dens.land ~ flow.cov + snh, 
         data = dens.landscape[dens.landscape$dominant == 'opportunistic' & 
                                 dens.landscape$period == 'during',])
m6 <- lm(dens.land ~ flow.cov + snh, 
         data = dens.landscape[dens.landscape$dominant == 'non.crop' & 
                                 dens.landscape$period == 'during',])
m7 <- lm(dens.leek ~  snh , 
         data = dens.landscape[dens.landscape$dominant == 'dominant' & 
                                 dens.landscape$period == 'during',])
m8 <- lm(dens.leek ~  snh , 
         data = dens.landscape[dens.landscape$dominant == 'opportunistic' & 
                                 dens.landscape$period == 'during',])

summary(m1)
summary(m2)
summary(m3)
summary(m4)
summary(m5)
summary(m6)
summary(m7)
summary(m8)

### backtransform everything and plot partial residuals. 10^x for log10 transformation of Y
### for snh-cover x^2 as was sqrt transformed. 
par(mfrow = c(3,3), mar = c(4,6,2,2), oma = c(1,1,1,1), xpd = NA)
v1 <- visreg(m1, xvar = 'snh', 
             trans = function(x) 10^(x),
             partial = T, line.par = list(lwd =2, col = 'black'),
             points.par = list(cex = 1, col = 'black'), overlay = F,
             strip.names = c('Before crop flowering', 'During crop flowering'),
             xlab = 'Semi-natural habitat cover (%)',
             ylab ='Abundance landscape \n dominant before',
             main = 'Dominant crop pollinators',
             cex.lab = 1.2, cex.axis = 1.2,
             xtrans = function(x) x^2)
text(x = -10, y = 45, labels =  c('A'), cex = 2)

v2 <- visreg(m2, xvar = 'snh', 
             trans = function(x) 10^(x),
             partial = T, line.par = list(lwd =2, col = 'black'),
             points.par = list(cex = 1, col = 'black'), overlay = F,
             strip.names = c('Before crop flowering', 'During crop flowering'),
             xlab = 'Semi-natural habitat cover (%)',
             ylab ='Abundance landscape \n opportunistic before',
             main = 'Opportunistic crop pollinators',
             cex.lab = 1.2, cex.axis = 1.2,
             xtrans = function(x) x^2)
text(x = -10, y = 147, labels =  c('B'), cex = 2)

v3 <- visreg(m3, xvar = 'snh', 
             trans = function(x) 10^(x),
             partial = T, line.par = list(lwd =2, col = 'black'),
             points.par = list(cex = 1, col = 'black'), overlay = F,
             strip.names = c('Before crop flowering', 'During crop flowering'),
             xlab = 'Semi-natural habitat cover (%)',
             ylab ='Abundance landscape \n non-crop before',
             main = 'Non-crop pollinators',
             cex.lab = 1.2, cex.axis = 1.2,
             xtrans = function(x) x^2)
text(x = -10, y = 50, labels =  c('C'), cex = 2)

v4 <- visreg(m4, xvar = 'snh', 
             trans = function(x) 10^(x),
             partial = T, line.par = list(lwd =2, col = 'black'),
             points.par = list(cex = 1, col = 'black'), overlay = F,
             strip.names = c('Before crop flowering', 'During crop flowering'),
             xlab = 'Semi-natural habitat cover (%)',
             ylab ='Abundance landscape \n dominant during',
             cex.lab = 1.2, cex.axis = 1.2,
             xtrans = function(x) x^2)
text(x = -10, y = 18, labels =  c('D'), cex = 2)

v5 <- visreg(m5, xvar = 'snh', 
             trans = function(x) 10^(x),
             partial = T, line.par = list(lwd =2, col = 'black'),
             points.par = list(cex = 1, col = 'black'), overlay = F,
             strip.names = c('Before crop flowering', 'During crop flowering'),
             xlab = 'Semi-natural habitat cover (%)',
             ylab ='Abundance landscape \n opportunistic during',
             cex.lab = 1.2, cex.axis = 1.2,
             xtrans = function(x) x^2)
text(x = -10, y = 26, labels =  c('E'), cex = 2)

v6 <- visreg(m6, xvar = 'snh', 
             trans = function(x) 10^(x),
             partial = T, line.par = list(lwd =2, col = 'black'),
             points.par = list(cex = 1, col = 'black'), overlay = F,
             strip.names = c('Before crop flowering', 'During crop flowering'),
             xlab = 'Semi-natural habitat cover (%)',
             ylab ='Abundance landscape \n non-crop during',
             cex.lab = 1.2, cex.axis = 1.2,
             xtrans = function(x) x^2)
text(x = -10, y = 32.5, labels =  c('F'), cex = 2)

v7 <- visreg(m7, xvar = 'snh', 
             trans = function(x) 10^(x),
             partial = T, line.par = list(lwd =2, col = 'black'),
             points.par = list(cex = 1, col = 'black'), overlay = F,
             strip.names = c('Before crop flowering', 'During crop flowering'),
             xlab = 'Semi-natural habitat cover (%)',
             ylab ='Abundance crop dominant',
             cex.lab = 1.2, cex.axis = 1.2,
             xtrans = function(x) x^2)
text(x = -10, y = 320, labels =  c('G'), cex = 2)
text(x = 50, y = 270, labels = c('**'), cex = 4)

v8 <- visreg(m8, xvar = 'snh', 
             trans = function(x) 10^(x),
             partial = T, line.par = list(lwd =2, col = 'black'),
             points.par = list(cex = 1, col = 'black'), overlay = F,
             strip.names = c('Before crop flowering', 'During crop flowering'),
             xlab = 'Semi-natural habitat cover (%)',
             ylab ='Abundance crop opportunistic',
             cex.lab = 1.2, cex.axis = 1.2,
             xtrans = function(x) x^2)
text(x = -10, y = 125, labels =  c('H'), cex = 2)

par(mfrow = c(1,1), mar = c(4,4,2,2))

### Species vs snh per period (separate for wider landscape and crop) #####
sp.boot.land <- sp.boot[sp.boot$leek.land == 'land',]
sp.boot.leek <- sp.boot[sp.boot$leek.land == 'leek',]

### SNH relationship on species richness in wider landscape
sp.boot.land$snh <- dens.landscape$snh[match(sp.boot.land$field, dens.landscape$field)]
sp.boot.land$flow.cov <- flow.cov$flow.cov[match(paste(sp.boot.land$field, sp.boot.land$period), 
                                                 paste(flow.cov$field, flow.cov$period))]
data <- sp.boot.land[sp.boot.land$period != 'both',]

m1 <- lmer(species ~ snh*dominant*period + flow.cov + (1|field), data)
summary(m1)

m1b <- lmer(species ~ snh*dominant + snh*period + dominant*period + flow.cov + (1|field), data)
summary(m1b)
anova(m1,m1b) # m1b is marginally better and simpler.

m1bb <- lmer(species ~ snh*dominant + dominant*period + flow.cov + (1|field), data)
summary(m1bb)
anova(m1bb,m1b) #m1bb is simpler and not worse

m1bbb <- lmer(species ~ snh + dominant*period + flow.cov + (1|field), data)
summary(m1bbb)
anova(m1bb,m1bbb) #m1bb is significantly better

m1bc <- lmer(species ~ snh*dominant + period + flow.cov + (1|field), data)
summary(m1bc)
anova(m1bb,m1bc) #m1bb is signficantly better

data$interaction <- interaction(data$dominant, data$period)

m1bd <- lmer(species ~ snh +  flow.cov + interaction + (1|field), data)
summary(m1bd)
summary(glht(m1bd, linfct = mcp(interaction = "Tukey"))) # to test significant interactions

# Figure 3a
data$predict_species <- predict(m1bb) 

fit.m1 <- as.data.frame(Effect(c('dominant', "snh", "period"), m1bb, xlevels = 50))
labels <- c(before = "Before crop flowering", during = "During crop flowering")

t1 <- ggplot(fit.m1, aes(snh^2, fit, colour = dominant, linetype = dominant)) + 
  theme_bw() + 
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70", alpha = 0.5, linetype = 1)+ 
  geom_line(size = 2) +
  geom_point(data = data, aes(y= predict_species, colour = dominant, 
                              shape = dominant)) + 
  labs(x = "Semi-natural habitat cover (%)", y = "Number of pollinator species", tag = "A") +
  #  scale_y_continuous(limits = c(-2,20)) + 
  facet_grid(~period, labeller=labeller(period = labels)) +
  theme(axis.title.x = element_text(size=14),
        title = element_text(size = 18),
        axis.title.y = element_text(size=14),
        axis.text = element_text(size=14, colour="#000000"),
        legend.text=element_text(size=12),
        legend.title = element_text(size = 14),
        legend.position = 'none',
        strip.background =element_rect(fill="white"),
        strip.text = element_text(size=12))+
  scale_color_manual(name = 'Functional group', values=c('#354ea7', '#000000','#CF0000'),
                     labels = c("Dominant crop pollinators", "Non-crop pollinators", "Opportunistic crop pollinators"))  +
  scale_fill_manual(name = "Functional group", values = c('grey50','grey50'),
                    labels = c("Dominant crop pollinators", "Non-crop pollinators", "Opportunistic crop pollinators")) +
  scale_shape_manual(name = 'Functional group', values = c(16,15,17,8,7),
                     labels = c("Dominant crop pollinators", "Non-crop pollinators", "Opportunistic crop pollinators")) +
  scale_size_manual(values=c(10, 1.5)) +
  scale_linetype_manual(name = 'Functional group',
                        values=c("twodash", "solid", "dashed"),
                        labels = c("Dominant crop pollinators", "Non-crop pollinators", "Opportunistic crop pollinators"))

t1

#legend
t1b <- ggplot(fit.m1, aes(snh^2, fit, colour = dominant, linetype = dominant)) + 
  theme_bw() + 
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70", alpha = 0.5, linetype = 1)+ 
  geom_line(size = 2) +
  geom_point(data = data, aes(y= predict_species, colour = dominant, 
                              shape = dominant), size = 4) + 
  labs(x = "Semi-natural habitat cover (%)", y = "Number of pollinator species") +
  #  scale_y_continuous(limits = c(-2,20)) + 
  facet_grid(~period, labeller=labeller(period = labels)) +
  theme(axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text = element_text(size=14, colour="#000000"),
        legend.text=element_text(size=12),
        legend.title = element_text(size = 14),
        strip.background =element_rect(fill="white"),
        strip.text = element_text(size=12))+
  scale_color_manual(name = 'Functional group', values=c('#354ea7', '#000000','#CF0000'),
                     labels = c("Dominant crop pollinators", "Non-crop pollinators", "Opportunistic crop pollinators"),
                     guide = guide_legend(keywidth = unit(3, 'cm'), size = 4))  +
  scale_fill_manual(name = "Functional group", values = c('grey50','grey50'),
                    labels = c("Dominant crop pollinators", "Non-crop pollinators", "Opportunistic crop pollinators")) +
  scale_shape_manual(name = 'Functional group', values = c(16,15,17,8,7),
                     labels = c("Dominant crop pollinators", "Non-crop pollinators", "Opportunistic crop pollinators")) +
  scale_linetype_manual(name = 'Functional group',
                        values=c("twodash", "solid", "dashed"),
                        guide = guide_legend(keywidth = unit(3, 'cm')),
                        labels = c("Dominant crop pollinators", "Non-crop pollinators", "Opportunistic crop pollinators"))
t1b
legend <- g_legend(t1b)

# Post-hoc tests on single groups per period
m1 <- lm(species ~ snh + flow.cov, data[data$period == 'before' &
                                          data$dominant == 'dominant',])
summary(m1)
m2 <- lm(species ~ snh + flow.cov, data[data$period == 'before' &
                                          data$dominant == 'opportunistic',])
summary(m2)
m3 <- lm(species ~ snh + flow.cov, data[data$period == 'before' &
                                          data$dominant == 'non.crop',])
summary(m3)
m4 <- lm(species ~ snh + flow.cov, data[data$period == 'during' &
                                          data$dominant == 'dominant',])
summary(m4)
m5 <- lm(species ~ snh + flow.cov, data[data$period == 'during' &
                                          data$dominant == 'opportunistic',])
summary(m5)
m6 <- lm(species ~ snh + flow.cov, data[data$period == 'during' &
                                          data$dominant == 'non.crop',])
summary(m6)

# General relationship of snh on species richness in wider landscape
data2 <- aggregate(species ~ field, sp.boot.land[sp.boot.land$period == 'both',], sum)
data2$snh <- data$snh[match(data2$field, data$field)]
data2$flow.cov <- flow.cov.av$flow.cov[match(data2$field, flow.cov.av$field)]
plot(species ~ snh, data2)
model1 <- lm(species ~ snh + flow.cov, data2)
summary(model1)

### SNH relationship on species richness in crop fields
sp.boot.leek$snh <- dens.landscape$snh[match(sp.boot.leek$field, dens.landscape$field)]
data <- sp.boot.leek

m1 <- lmer(species ~ snh*dominant + (1|field), data)
summary(m1)
m1b <- lmer(species ~ dominant + snh + (1|field), data)
anova(m1,m1b) # model with interaction (m1) significantly better
m1c <- lmer(species ~ dominant  + (1|field), data)
anova(m1,m1c)

# Figure 3b
data$predict_species <- predict(m1) 
t <- as.data.frame(Effect(c('dominant', "snh"), m1, partial.residuals = T))

fit.m1 <- as.data.frame(Effect(c('dominant', "snh"), m1, xlevels = 50))

labels <- c(before = "Before crop flowering", during = "During crop flowering")

t2 <- ggplot(fit.m1, aes(snh^2, fit, colour = dominant, linetype = dominant)) + 
  theme_bw() + 
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70", alpha = 0.5, linetype = 1)+ 
  geom_line(size = 2) +
  geom_point(data = data, aes(y= predict_species, colour = dominant, 
                              shape = dominant)) + 
  labs(x = "Semi-natural habitat cover (%)", y = "Number of pollinator species", tag = "B") +
  theme(axis.title.x = element_text(size=14),
        title = element_text(size = 18),
        axis.title.y = element_text(size=14),
        axis.text = element_text(size=14, colour="#000000"),
        legend.text=element_text(size=12),
        legend.title = element_text(size = 14),
        legend.position = 'none',
        strip.background =element_rect(fill="white"),
        strip.text = element_text(size=12))+
  scale_color_manual(name = 'Functional group', 
                     values=c('#354ea7', '#CF0000')) +
  guides(shape=FALSE) +
  scale_linetype_manual(name = 'Functional group',
                        values=c("twodash","dashed" )) +
  scale_shape_manual(name = 'Functional group',
                     values = c(16,17,15,8,7))
t2

# Post-hoc tests on single groups
m1 <- lm(species ~ snh, data[data$dominant== 'dominant',])
summary(m1)
m2 <- lm(species ~ snh, data[data$dominant== 'opportunistic',])
summary(m2)

# General relationship of snh on species richness in crop field
data <- aggregate(species ~ field, sp.boot.leek, sum)
data$snh <- sp.boot.leek$snh[match(data$field, sp.boot.leek$field)]
m1 <- lm(species ~ snh, data)
summary(m1)

# Figure 3
grid.arrange(t1,t2,legend, nrow = 2, widths = c(0.66,0.33),
             heights = c(0.8,0.2))
