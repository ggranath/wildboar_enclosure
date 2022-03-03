################################################
# Rooting by wild boar (Sus scrofa) 
# has minor effects on soil nutrient and 
# carbon dynamics, but reduces tree growth
#
# Lundgren et al.
#
# Under review
#
################################################

# Load required packages
library(lme4)
library(lmerTest)
library(tidyr)
library(dplyr)
library(car)
library(vegan)
library(readxl)
library(cowplot)
library(egg)
library(ggfortify)


# Load soil data ####
df <- read_excel("Lundgren_etal_soil_data.xlsx")
df$Treatment <- relevel(factor(df$Treatment), ref="Reference")
df$Environment <- factor(df$Environment)

#function for standard error
se = function (x) sd(x)/sqrt(length(x))

# Soil Density ####
  # mean density across soil layers
df$Density <- (df$Top_Density_g_cm3 + df$Lower_Density_g_cm3)/2
Dens.test <- lm(Density ~ Treatment * Environment, df)
summary(Dens.test)
plot(Dens.test)
Anova(Dens.test, type = 2)

# Complete model
df.dens.long <- gather(df, key = "Layer", value = "Density",
             Top_Density_g_cm3, Lower_Density_g_cm3)
df.dens.long$Layer <- factor(df.dens.long$Layer)
df.dens.long$Environment <- factor(df.dens.long$Environment)
Dens.Layer <- lmer(Density ~ Treatment * Environment*Layer + (1|Plot_number), df.dens.long)
summary(Dens.Layer)
anova(Dens.Layer, type=2)
plot(Dens.Layer)
ls_means(Dens.Layer)

mean.dens <- aggregate(Density ~ Treatment + Environment + Layer, df.dens.long, mean)
mean.dens$se <- aggregate(Density ~ Treatment + Environment + Layer, df.dens.long, se)$Dens
mean.dens$Layer <- rep(c("Lower", "Top"), each=4)

den <- ggplot(mean.dens, aes(x = Layer, y = Density, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge", color="black") +
  geom_errorbar(aes(ymax = Density + se, ymin = Density - se),
                position = position_dodge(0.9), width = 0.25, color = "Gray25") +
  ylab(expression(Soil~density~'('~g~cm^-3~')')) +
  xlab(c("Soil layer")) +
  theme_cowplot() +
  panel_border() +
  theme(legend.position = c(0.05, 0.85),
        legend.title = element_blank()) +
  scale_fill_manual(values = c("white", "gray")) +
  scale_y_continuous(limits = c(0,1.5), expand = c(0, 0)) +
  facet_grid(.~Environment) 
den <- tag_facet(den, open="", close="", tag_pool = c("(a) PINE", "(b) SPRUCE"))
den

# Organic matter content ####

df$OMC <- (df$Top_Soil_organic_matter + df$Lower_Soil_organic_matter)/2
OMC.test <- lm(OMC ~ Treatment, df)
Anova(OMC.test, type = 2)

OMC.Env <- lm(OMC ~ Treatment * Environment, df)
Anova(OMC.Env, type = 2)

# mean OMC across soil layers
df$OMC <- (df$Top_Soil_organic_matter + df$Lower_Soil_organic_matter)/2
omc.test <- lm(OMC ~ Treatment * Environment, df)
summary(omc.test)
plot(omc.test)
Anova(omc.test, type = 2)

# Complete model
df.omc.long <- gather(df, key = "Layer", value = "OMC",
                      Top_Soil_organic_matter, Lower_Soil_organic_matter)
df.omc.long$Layer <- factor(df.omc.long$Layer)
df.omc.long$Environment <- factor(df.omc.long$Environment)
omc.Layer <- lmer(log(OMC) ~ Treatment * Environment*Layer + (1|Plot_number), df.omc.long)
summary(omc.Layer)
anova(omc.Layer, type=2)
plot(omc.Layer)
ls_means(omc.Layer)

mean.omc <- aggregate(OMC*1000 ~ Treatment + Environment + Layer, df.omc.long, mean)
mean.omc$se <- aggregate(OMC*1000 ~ Treatment + Environment + Layer, df.omc.long, se)$OMC
mean.omc$Layer <- rep(c("Lower", "Top"), each=4)
colnames(mean.omc)[4] <- "OMC"

om <- ggplot(mean.omc, aes(x = Layer, y = OMC, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge", color="black") +
  geom_errorbar(aes(ymax = OMC + se, ymin = OMC - se),
                position = position_dodge(0.9), width = 0.25, color = "Gray25") +
  ylab(expression(Soil~organic~matter~'('~mg~g^-1~')')) +
  xlab(c("Soil layer")) +
  theme_cowplot() +
  panel_border() +
  theme(legend.position = c(0.05, 0.85),
        legend.title = element_blank()) +
  scale_fill_manual(values = c("white", "gray")) +
  scale_y_continuous(limits = c(0,225), expand = c(0, 0)) +
  facet_grid(.~Environment) 
om <- tag_facet(om, open="", close="", tag_pool = c("(c) PINE", "(d) SPRUCE"))
om

fig2 <- plot_grid(den, om, ncol = 1, align = "v")
save_plot("fig2.png", fig2, ncol=2, base_height = 6.71, base_width = 4)

# Soil temperature and moisture ####
df.resp <- read_excel("Lundgren_etal_soil_resp_data.xlsx")
df.temp <- gather(df.resp, key = "Month", value = "Temperature",
             Soil_temperature_June, Soil_temperature_July, Soil_temperature_August,
             Soil_temperature_September, Soil_temperature_November) %>%
  group_by(Month, Treatment, Environment, Plot_ID) %>%
  summarise(Temperature = mean(Temperature))
Temp.test <- lmer(Temperature ~ Treatment*Environment*Month + (1|Plot_ID), df.temp)
summary(Temp.test)
anova(Temp.test, type = 2)

df.mois <- gather(df.resp, key = "Month", value = "Moisture",
             Soil_moisture_June, Soil_moisture_July, Soil_moisture_August,
             Soil_moisture_September)%>%
  group_by(Month, Treatment, Environment, Plot_ID) %>%
  summarise(Moisture = mean(Moisture))
Moist.test <- lmerTest::lmer(Moisture ~ Treatment*Environment*Month + (1|Plot_ID), df.mois)
summary(Moist.test)
anova(Moist.test, type = 2)

# PLOT temp + moisture
# temp
mean.temp <- aggregate(Temperature~Treatment+Environment+Month, df.temp, mean)
mean.temp$se <- aggregate(Temperature~Treatment+Environment+Month, df.temp, se)$Temperature
mean.temp$month <- as.Date(rep(c("2020-08-15","2020-07-15", "2020-06-15", "2020-11-15", "2020-09-15"), each=4)) 
ggplot(mean.temp, aes(x=month, y=Temperature, color=Treatment, shape=Environment))+
  geom_line() +
  geom_point()
# similar patterns across environments so lets show month * treatment interaction
mean.temp <- aggregate(Temperature~Treatment+Month, df.temp, mean)
mean.temp$se <- aggregate(Temperature~Treatment+Month, df.temp, se)$Temperature
mean.temp$month <- as.Date(rep(c("2020-08-15","2020-07-15", "2020-06-15", "2020-11-15", "2020-09-15"), each=2)) 
mean.temp$Treatment <- factor(mean.temp$Treatment, levels=c("Reference", "Enclosure"))
tem <- ggplot(mean.temp, aes(x=month, y=Temperature, shape=Treatment))+
  geom_line(aes(linetype = Treatment)) +
  geom_point()+
  geom_errorbar(aes(ymax = Temperature + se, ymin = Temperature - se), width = 0.1, color = "Gray25") +
  scale_shape_manual(values=c(21,19)) +
  scale_y_continuous(limits = c(0,20),  breaks=seq(0,20,2), expand = c(0, 0)) +
  xlab("Month") +
  ylab("Temperature (°C)") +
  theme_cowplot() +
  theme(legend.position = c(0.05, 0.2),
        legend.title = element_blank()) + draw_plot_label("(a)", x= as.Date("2020-06-05"), y=20 )
tem

# moist
mean.mois <- aggregate(Moisture~Treatment+Environment+Month, df.mois, mean)
mean.mois$se <- aggregate(Moisture~Treatment+Environment+Month, df.mois, se)$Moisture
mean.mois$month <- as.Date(rep(c("2020-08-15","2020-07-15", "2020-06-15", "2020-09-15"), each=4)) 
mean.mois$Treatment <- factor(mean.mois$Treatment, levels=c("Reference", "Enclosure"))
moi <- ggplot(mean.mois, aes(x=month, y=Moisture, fill=Treatment, shape=Environment))+
  geom_line(aes(linetype = Treatment)) +
  geom_point(aes(shape= Environment), position = position_dodge(2)) + 
  geom_errorbar(aes(ymax = Moisture + se, ymin = Moisture - se), position = position_dodge(2),
                width = 0.1, color = "Gray25") +
  scale_shape_manual(values=c(22,24)) +
  scale_fill_manual(values=c("white", "black"), guide=guide_legend(override.aes=list(shape=22) )) +
  scale_y_continuous(limits = c(0,0.7),  breaks=seq(0,0.7,0.1), expand = c(0, 0)) +
  xlab("Month") +
  ylab(expression(Moisture~"("~dm^3~dm^-3~")")) +
  theme_cowplot() +
  theme(legend.position = c(0.2, 0.15),
        legend.direction = "horizontal",
        legend.spacing.x = unit(x = 0.1, 'cm'),
        legend.spacing.y = unit(x = 0, 'cm'),
        legend.title = element_blank()) + draw_plot_label("(b)", x= as.Date("2020-06-05"), y=0.7)
moi

# C, N & P ####
df$Carbon <- (df$Top_Carbon_mg_g + df$Lower_Carbon_mg_g)/2
Carbon.test <- lm(Carbon ~ Treatment, df)
Anova(Carbon.test, type = 2)
Carbon.Env <- lm(Carbon ~ Treatment * Environment, df)
Anova(Carbon.Env, type = 2)

df$Nitrogen <- (df$Top_Nitrogen_mg_g + df$Lower_Nitrogen_mg_g)/2
Nitrogen.test <- lm(Nitrogen ~ Treatment, df)
Anova(Nitrogen.test, type = 2)
Nitrogen.Env <- lm(Nitrogen ~ Treatment * Environment, df)
Anova(Nitrogen.Env, type = 2)


# Complete model C, N, P
     # C conc
df.c.long <- gather(df, key = "Layer", value = "C",
                      Top_Carbon_mg_g, Lower_Carbon_mg_g)
df.c.long$Layer <- factor(df.c.long$Layer)
df.c.long$Environment <- factor(df.c.long$Environment)
c.Layer <- lmer(C ~ Treatment * Environment*Layer + (1|Plot_number), df.c.long)
summary(c.Layer)
anova(c.Layer, type=2)
plot(c.Layer)
ls_means(c.Layer)

mean.cconc <- aggregate(C ~ Treatment + Environment + Layer, df.c.long, mean)
mean.cconc$se <- aggregate(C ~ Treatment + Environment + Layer, df.c.long, se)$C
mean.cconc$Layer <- rep(c("Lower", "Top"), each=4)

cconc <- ggplot(mean.cconc, aes(x = Layer, y = C, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge", color="black") +
  geom_errorbar(aes(ymax = C + se, ymin = C - se),
                position = position_dodge(0.9), width = 0.25, color = "Gray25") +
  ylab(expression(Carbon~concentration~'('~mg~g^-1~')')) +
  xlab(c("Soil layer")) +
  theme_cowplot() +
  panel_border() +
  theme(legend.position = c(0.05, 0.85),
        legend.title = element_blank()) +
  scale_fill_manual(values = c("white", "gray")) +
  scale_y_continuous(limits = c(0,90), expand = c(0, 0)) +
  facet_grid(.~Environment) 
cconc <- tag_facet(cconc, open="", close="", tag_pool = c("(a) PINE", "(b) SPRUCE"))
cconc

    # C stock
df$C_m2 <- (df$Top_Carbon_mg_g*df$Top_Density_g_cm3*2 + 
              df$Lower_Carbon_mg_g*df$Lower_Density_g_cm3*3)*10000/1000
c.tot <- lm(C_m2 ~ Treatment * Environment, df)
summary(c.tot)
Anova(c.tot, type=2)
plot(c.tot)

mean.cm2 <- aggregate(C_m2 ~ Treatment + Environment, df, mean)
mean.cm2$se <- aggregate(C_m2 ~ Treatment + Environment, df, se)$C_m2

cstock <- ggplot(mean.cm2, aes(x = Environment, y = C_m2, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge", color="black") +
  geom_errorbar(aes(ymax = C_m2 + se, ymin = C_m2 - se),
                position = position_dodge(0.9), width = 0.25, color = "Gray25") +
  ylab(expression(Carbon~'('~g~m^-2~')')) +
  xlab(c("Environment")) +
  theme_cowplot() +
  panel_border() +
  theme(legend.position = c(0.05, 0.85),
        legend.title = element_blank()) +
  scale_fill_manual(values = c("white", "gray")) +
  scale_y_continuous(limits = c(0,4000), expand = c(0, 0))
cstock <- tag_facet(cstock, open="", close="", tag_pool = c("(a)"))
cstock

    # N conc
df.n.long <- gather(df, key = "Layer", value = "N",
                    Top_Nitrogen_mg_g, Lower_Nitrogen_mg_g)
df.n.long$Layer <- factor(df.n.long$Layer)
df.n.long$Environment <- factor(df.n.long$Environment)
n.Layer <- lmer(N ~ Treatment * Environment*Layer + (1|Plot_number), df.n.long)
summary(n.Layer)
anova(n.Layer, type=2)
plot(n.Layer)
ls_means(n.Layer)

mean.nconc <- aggregate(N ~ Treatment + Environment + Layer, df.n.long, mean)
mean.nconc$se <- aggregate(N ~ Treatment + Environment + Layer, df.n.long, se)$N
mean.nconc$Layer <- rep(c("Lower", "Top"), each=4)

nconc <- ggplot(mean.nconc, aes(x = Layer, y = N, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge", color="black") +
  geom_errorbar(aes(ymax = N + se, ymin = N - se),
                position = position_dodge(0.9), width = 0.25, color = "Gray25") +
  ylab(expression(Nitrogen~concentration~'('~mg~g^-1~')')) +
  xlab(c("Soil layer")) +
  theme_cowplot() +
  panel_border() +
  theme(legend.position = c(0.05, 0.85),
        legend.title = element_blank()) +
  scale_fill_manual(values = c("white", "gray")) +
  scale_y_continuous(limits = c(0,6), expand = c(0, 0)) +
  facet_grid(.~Environment) 
nconc <- tag_facet(nconc, open="", close="", tag_pool = c("(c) PINE", "(d) SPRUCE"))
nconc

    # N stock
df$N_m2 <- (df$Top_Nitrogen_mg_g*df$Top_Density_g_cm3*2 + 
              df$Lower_Nitrogen_mg_g*df$Lower_Density_g_cm3*3)*10000/1000
n.tot <- lm(N_m2 ~ Treatment * Environment, df)
summary(n.tot)
Anova(n.tot, type=2)
plot(n.tot)

mean.nm2 <- aggregate(N_m2 ~ Treatment + Environment, df, mean)
mean.nm2$se <- aggregate(N_m2 ~ Treatment + Environment, df, se)$N_m2

nstock <- ggplot(mean.nm2, aes(x = Environment, y = N_m2, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge", color="black") +
  geom_errorbar(aes(ymax = N_m2 + se, ymin = N_m2 - se),
                position = position_dodge(0.9), width = 0.25, color = "Gray25") +
  ylab(expression(Nitrogen~'('~g~m^-2~')')) +
  xlab(c("Environment")) +
  theme_cowplot() +
  panel_border() +
  theme(legend.position = c(0.05, 0.85),
        legend.title = element_blank()) +
  scale_fill_manual(values = c("white", "gray")) +
  scale_y_continuous(limits = c(0,300), expand = c(0, 0))
nstock <- tag_facet(nstock, open="", close="", tag_pool = c("(b)"))
nstock

    # P conc
df.p.long <- gather(df, key = "Layer", value = "P",
                    Top_Phosphorous_mg_g, LowerPhosphorous_mg_g)
df.p.long$Layer <- factor(df.p.long$Layer)
df.p.long$Environment <- factor(df.p.long$Environment)
p.Layer <- lmer(P ~ Treatment * Environment*Layer + (1|Plot_number), df.p.long)
summary(p.Layer)
anova(p.Layer, type=2)
plot(p.Layer)
ls_means(p.Layer)

mean.pconc <- aggregate(P ~ Treatment + Environment + Layer, df.p.long, mean)
mean.pconc$se <- aggregate(P ~ Treatment + Environment + Layer, df.p.long, se)$P
mean.pconc$Layer <- rep(c("Lower", "Top"), each=4)

pconc <- ggplot(mean.pconc, aes(x = Layer, y = P, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge", color="black") +
  geom_errorbar(aes(ymax = P + se, ymin = P - se),
                position = position_dodge(0.9), width = 0.25, color = "Gray25") +
  ylab(expression(Phosphorous~concentration~'('~mg~g^-1~')')) +
  xlab(c("Soil layer")) +
  theme_cowplot() +
  panel_border() +
  theme(legend.position = c(0.05, 0.85),
        legend.title = element_blank()) +
  scale_fill_manual(values = c("white", "gray")) +
  scale_y_continuous(limits = c(0,0.7), expand = c(0, 0)) +
  facet_grid(.~Environment) 
pconc <- tag_facet(pconc, open="", close="", tag_pool = c("(e) PINE", "(f) SPRUCE"))
pconc

    # P stock
df$P_m2 <- (df$Top_Phosphorous_mg_g*df$Top_Density_g_cm3*2 + 
              df$LowerPhosphorous_mg_g*df$Lower_Density_g_cm3*3)*10000/1000
p.tot <- lm(P_m2 ~ Treatment * Environment, df)
summary(p.tot)
Anova(p.tot, type=2)
plot(p.tot)

p.mean <- aggregate(P_m2 ~ Treatment, df, mean)
p.mean$se <- aggregate(P_m2 ~ Treatment, df, FUN=se)$P_m2
t.test(P_m2 ~ Treatment, df, var.equal = FALSE)

mean.pm2 <- aggregate(P_m2 ~ Treatment + Environment, df, mean)
mean.pm2$se <- aggregate(P_m2 ~ Treatment + Environment, df, se)$P_m2

pstock <- ggplot(mean.pm2, aes(x = Environment, y = P_m2, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge", color="black") +
  geom_errorbar(aes(ymax = P_m2 + se, ymin = P_m2 - se),
                position = position_dodge(0.9), width = 0.25, color = "Gray25") +
  ylab(expression(Phosphorous~'('~g~m^-2~')')) +
  xlab(c("Environment")) +
  theme_cowplot() +
  panel_border() +
  theme(legend.position = c(0.05, 0.85),
        legend.title = element_blank()) +
  scale_fill_manual(values = c("white", "gray")) +
  scale_y_continuous(limits = c(0,35), expand = c(0, 0))
pstock <- tag_facet(pstock, open="", close="", tag_pool = c("(c)"))
pstock

fig4 <- plot_grid(cconc, nconc, pconc, ncol = 1, align = "v")
save_plot("fig4.png", fig4, ncol=2, base_height = 10.71, base_width = 4)

fig5 <- plot_grid(cstock, nstock, pstock, ncol = 1, align = "v")
save_plot("fig5.png", fig5, ncol=2, base_height = 9.71, base_width = 3)

        # NO3- and NH4+ accum
NO3.test <- lm(log(Nitrate_ug_cm2_day) ~ Treatment * Environment, df)
summary(NO3.test)
Anova(NO3.test)
plot(NO3.test)
    # test for enclosure effect in spruce forest
NO3.test <- lm(log(Nitrate_ug_cm2_day) ~ Treatment * relevel(Environment, ref="Spruce"), df)
summary(NO3.test)

NH4.test <- lm(log(Ammonium_ug_cm2_day) ~ Treatment * Environment, df)
Anova(NH4.test)
summary(NH4.test)
plot(NH4.test)

mean.NO3 <- aggregate(Nitrate_ug_cm2_day ~ Treatment + Environment, df, mean)
mean.NO3$se <- aggregate(Nitrate_ug_cm2_day ~ Treatment + Environment, df, se)$Nitrate_ug_cm2_day

NO3 <- ggplot(mean.NO3, aes(x = Environment, y = Nitrate_ug_cm2_day, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge", color="black") +
  geom_errorbar(aes(ymax = Nitrate_ug_cm2_day + se, ymin = Nitrate_ug_cm2_day - se),
                position = position_dodge(0.9), width = 0.25, color = "Gray25") +
  ylab(expression(NO[3]~'('~µg~cm^-2~day^-1~')')) +
  xlab(c("Environment")) +
  theme_cowplot() +
  panel_border() +
  theme(legend.position = c(0.05, 0.85),
        legend.title = element_blank()) +
  scale_fill_manual(values = c("white", "gray")) +
  scale_y_continuous(limits = c(0,1.2), expand = c(0, 0))
NO3 <- tag_facet(NO3, open="", close="", tag_pool = c("(a)"))
NO3

mean.NH4 <- aggregate(Ammonium_ug_cm2_day ~ Treatment + Environment, df, mean)
mean.NH4$se <- aggregate(Ammonium_ug_cm2_day ~ Treatment + Environment, df, se)$Ammonium_ug_cm2_day

NH4 <- ggplot(mean.NH4, aes(x = Environment, y = Ammonium_ug_cm2_day, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge", color="black") +
  geom_errorbar(aes(ymax = Ammonium_ug_cm2_day + se, ymin = Ammonium_ug_cm2_day - se),
                position = position_dodge(0.9), width = 0.25, color = "Gray25") +
  ylab(expression(NH[4]~'('~µg~cm^-2~day^-1~')')) +
  xlab(c("Environment")) +
  theme_cowplot() +
  panel_border() +
  theme(legend.position = c(0.05, 0.85),
        legend.title = element_blank()) +
  scale_fill_manual(values = c("white", "gray")) +
  scale_y_continuous(limits = c(0,2.7), expand = c(0, 0))
NH4 <- tag_facet(NH4, open="", close="", tag_pool = c("(b)"))
NH4

fig6 <- plot_grid(NO3, NH4, ncol = 1, align = "v")
save_plot("fig6.png", fig6, ncol=2, base_height = 6.71, base_width = 2)

# PLFA ####
    # PERMANOVA # Lower_cy17:0 not exists for top
cols <- unlist(lapply(strsplit(colnames(df[,17:30]), "_"),'[', 2))
df.plfa.long1 <- df[,17:30]
df.plfa.long2 <- df[,34:47]
colnames(df.plfa.long1) <- cols ; colnames(df.plfa.long2) <- cols
df.plfa.long <- rbind(df.plfa.long1, df.plfa.long2)
df.plfa.long <- cbind(df[,1:4], df.plfa.long)

#df.plfa.long <- gather(df, key = "Layer", value = "Density",
#                       'Top_14:0': 'Top_19:0')

df.plfa.long$Layer <- factor(rep(c("Top", "Lower"), each=20))
df.plfa.long$Environment <- factor(df.plfa.long$Environment)

# PCA plot

df.PCA <- df.plfa.long %>% group_by(Plot_number, Treatment, Environment) %>%
  summarise(across(.cols = 2:15, list(mean)))

For.PCA <- prcomp(df.PCA[,c(4:17)], center=TRUE, scale=TRUE)
summary(For.PCA)

PCA <- autoplot(For.PCA, data = df.PCA, fill = "Treatment", shape = "Environment",
                  colour = "Treatment", frame = T,
                loadings = F, loadings.label = F, size = 3, scale = 1) +
  scale_shape_manual(values=c(22,24)) +
  scale_color_manual(values=c("light grey", "dark grey")) +
  scale_fill_manual(values=c("light grey", "black"), guide=guide_legend(override.aes=list(shape=22) )) +
  theme_cowplot() +
  theme(legend.position = c(0.02, 0.1),
        legend.direction = "horizontal",
        legend.spacing.x = unit(x = 0.1, 'cm'),
        legend.spacing.y = unit(x = 0, 'cm'),
        legend.title = element_blank()) + draw_plot_label("", x=-0.3, y=-0.7) +
  ylim(-0.7, 0.5) +
  xlim(-0.5, 0.6)
PCA <- tag_facet(PCA, open="", close="", tag_pool = c("(a)"))
PCA

coord <- as.data.frame(cbind(c("14:0", "i15:0", "a15:0", "15:0", "i16:0", "16:1w9", "16:1w5", "16:0", "i17:0",
           "cy17:0", "17:0", "18:2w6", "18:0", "cy19:0"),
           c(-0.05, -0.18, 0.34, -0.47, 0.005, 0.13, 0.29,
             -0.27, 0.14, 0.08, 0.22, -0.32, -0.14, -0.1),
           c(0.27, 0.25, 0.19, 0.095, 0.3, 0.34, 0.26, 0.25, 0.23,
             0.3, 0.26, 0.22, 0.27, 0.28)))
col_headings <- c('label','y','x')
colnames(coord) <- col_headings
coord$y <- as.numeric(coord$y)
coord$x <- as.numeric(coord$x)

PCA2 <- autoplot(For.PCA, shape = F, label = F,
                  frame = F, loadings = T, loadings.label = F, size = 3, scale = 1,
                  loadings.colour = "black") +
  theme_cowplot() +
  ylim(-0.7, 0.5) +
  xlim(-0.5, 0.6) +
  geom_text(data = coord, aes(
    x = x, y = y,
    label = label))
PCA2 <- tag_facet(PCA2, open="", close="", tag_pool = c("(b)"))
PCA2

fig7 <- plot_grid(PCA, PCA2, ncol = 2)
save_plot("fig7.png", fig7, ncol=2, base_height = 5, base_width = 5.5)

# we test the enclosure effect after accounting for layer and environment
# this test terms sequentially
PLFA.test <- adonis2(df.plfa.long[,5:18] ~ Environment*Layer*Treatment, data = df.plfa.long,
                     permutations = 9999, by = "terms", method = "euclidean")
PLFA.test

# another way to test treatment while controlling for environment
PLFA.test <- adonis(df.plfa.long[,5:18] ~ Layer*Treatment, data = df.plfa.long,
                    permutations = 9999, strata= df.plfa.long$Environment ,  method = "euclidean")
PLFA.test

  # average over top and lower layer as layer didnt interact with enclosure effect
df.pfla.long.ave <- aggregate(df.plfa.long[,5:18], df.plfa.long[,c("Treatment","Environment", "Plot_ID")], mean)
PLFA.test <- adonis2(df.pfla.long.ave[,4:17] ~ Environment*Treatment, data = df.pfla.long.ave,
                     permutations = 9999, by = "terms", method = "euclidean")
PLFA.test

  # test individual PFLAs
vars.ind = colnames(df.pfla.long.ave)[4:17]
tests.ind <-list()
for (i in 1:length(vars.ind)) {
  form <- as.formula( paste("log(`",vars.ind[i],"`)"," ~ Environment+Treatment", sep="") )
  tests.ind[[i]] <- lm(form, data = df.pfla.long.ave)
}
res.ind <- lapply(tests.ind, function (x) data.frame(Anova(x))[1:4,] ) %>% bind_rows()
rownames(res.ind) <- paste(rep(vars.ind, each=4),rownames(res.ind), sep=" + ")
res.ind
lapply(tests.ind, plot)

     # Test other explanatory variables on PFLA composition
df.pfla.long.ave$no3 <- df$Nitrate_ug_cm2_day
df.pfla.long.ave$nh4 <- df$Ammonium_ug_cm2_day
df.pfla.long.ave$p <-  rowMeans(cbind(df$LowerPhosphorous_mg_g, df$Top_Phosphorous_mg_g))
df.pfla.long.ave$n <- rowMeans(cbind(df$Lower_Nitrogen_mg_g,df$Top_Nitrogen_mg_g))
df.pfla.long.ave$c <- rowMeans(cbind(df$Lower_Carbon_mg_g, df$Top_Carbon_mg_g))
df.pfla.long.ave$cn <- df.pfla.long.ave$c/df.pfla.long.ave$n 

vars = c("no3", "nh4", "p", "n", "c", "cn")
tests <-list()
for (i in 1:length(vars)) {
  form <- as.formula( paste( "df.pfla.long.ave[,4:17] ~ ", vars[i]) )
tests[[i]] <- adonis2(form, data = df.pfla.long.ave,
                     permutations = 9999, by = "terms", method = "euclidean")
}
res <- lapply(tests, function (x) data.frame(x)[1,] ) %>% bind_rows()

    # would the treatment effect change if we take into account nutrient variables?
tests <-list()
for (i in 1:length(vars)) {
  form <- as.formula( paste( "df.pfla.long.ave[,4:17] ~ ", vars[i], "+ Treatment") )
  tests[[i]] <- adonis2(form, data = df.pfla.long.ave,
                        permutations = 9999, by = "terms", method = "euclidean")
}
res.treat <- lapply(tests, function (x) data.frame(x)[2,] ) %>% bind_rows()
rownames(res.treat) <- paste(vars,rownames(res.treat), sep=" + ")
res.treat

 # Tot PLFA 
df.totpfla.long <- gather(df, key = "Layer", value = "tot.pfla",
                          Top_Total_PLFA, Lower_Total_PLFA)
df.totpfla.long$Layer <- factor(df.totpfla.long$Layer)
df.totpfla.long$Environment <- factor(df.totpfla.long$Environment)
Totpfla.Layer <- lmer(log(tot.pfla) ~ Treatment * Environment*Layer + (1|Plot_number), df.totpfla.long)
summary(Totpfla.Layer)
anova(Totpfla.Layer, type=2)
plot(Totpfla.Layer)
ls_means(Totpfla.Layer)

mean.totplfa <- aggregate(tot.pfla ~ Treatment + Environment + Layer, df.totpfla.long, mean)
mean.totplfa$se <- aggregate(tot.pfla ~ Treatment + Environment + Layer, df.totpfla.long, se)$tot.pfla
mean.totplfa$Layer <- rep(c("Lower", "Top"), each=4)

totplfa <- ggplot(mean.totplfa, aes(x = Layer, y = tot.pfla, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge", color="black") +
  geom_errorbar(aes(ymax = tot.pfla + se, ymin = tot.pfla - se),
                position = position_dodge(0.9), width = 0.25, color = "Gray25") +
  ylab(expression(PLFA~'('~µmol~g^-1~')')) +
  xlab(c("Soil layer")) +
  theme_cowplot() +
  panel_border() +
  theme(legend.position = c(0.05, 0.85),
        legend.title = element_blank()) +
  scale_fill_manual(values = c("white", "gray")) +
  scale_y_continuous(limits = c(0,0.25), expand = c(0, 0)) +
  facet_grid(.~Environment) 
totplfa <- tag_facet(totplfa, open="", close="", tag_pool = c("(a) PINE", "(b) SPRUCE"))
totplfa

  # based on layer average, gives same treatment results
df$PLFA_total <- (df$Top_Total_PLFA + df$Lower_Total_PLFA)/2
PLFA.total <- lm(PLFA_total ~ Treatment*Environment, df)
summary(PLFA.total)
Anova(PLFA.total, type = 2)
plot(PLFA.total)

# Fung:bac 
df.funbac.long <- gather(df, key = "Layer", value = "fun_bac",
                          Top_Fungi_bacteria_ratio, Lower_Fungi_bacteria_ratio)
df.funbac.long$Layer <- factor(df.funbac.long$Layer)
df.funbac.long$Environment <- factor(df.funbac.long$Environment)
funbac.Layer <- lmer(fun_bac ~ Treatment * Environment*Layer + (1|Plot_number), df.funbac.long)
summary(funbac.Layer)
anova(funbac.Layer, type=2)
plot(funbac.Layer)
ls_means(funbac.Layer)

mean.funbac <- aggregate(fun_bac ~ Treatment + Environment + Layer, df.funbac.long, mean)
mean.funbac$se <- aggregate(fun_bac ~ Treatment + Environment + Layer, df.funbac.long, se)$fun_bac
mean.funbac$Layer <- rep(c("Lower", "Top"), each=4)

funbac <- ggplot(mean.funbac, aes(x = Layer, y = fun_bac, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge", color="black") +
  geom_errorbar(aes(ymax = fun_bac + se, ymin = fun_bac - se),
                position = position_dodge(0.9), width = 0.25, color = "Gray25") +
  ylab(c("Fungi:Bacteria")) +
  xlab(c("Soil layer")) +
  theme_cowplot() +
  panel_border() +
  theme(legend.position = c(0.05, 0.85),
        legend.title = element_blank()) +
  scale_fill_manual(values = c("white", "gray")) +
  scale_y_continuous(limits = c(0,0.2), expand = c(0, 0)) +
  facet_grid(.~Environment) 
funbac <- tag_facet(funbac, open="", close="", tag_pool = c("(c) PINE", "(d) SPRUCE"))
funbac

fig8 <- plot_grid(totplfa, funbac, ncol = 1, align = "v")
save_plot("fig8.png", fig8, ncol=2, base_height = 6.71, base_width = 4)

  # based on layer average, gives same treatment results
df$Fungibac <- (df$Top_Fungi_bacteria_ratio + df$Lower_Fungi_bacteria_ratio)/2
Fungibac.test <- lm(Fungibac ~ Treatment*Environment, df)
summary(Fungibac.test)
Anova(Fungibac.test, type = 2)
plot(Fungibac.test)

# Respiration ####
df.resp <- read_excel("Lundgren_etal_soil_resp_data.xlsx")
df.resp.long <- gather(df.resp, key = "Month", value = "Respiration",
       CO2_June, CO2_July, CO2_August, CO2_September, CO2_November)
df.resp.long$Treatment <- factor(df.resp.long$Treatment, levels=c("Reference", "Enclosure"))
df.resp.long$Environment <- factor(df.resp.long$Environment)
df.resp.long$Month <- factor(df.resp.long$Month)
Resp.test <- lmerTest::lmer(log(Respiration) ~ Month*Environment*Treatment + (1|Plot_ID), df.resp.long)
summary(Resp.test)
anova(Resp.test, type = 2)
plot(Resp.test)
ls_means(Resp.test) # used to calculate % diffrences (use exp(enclosure)/exp(reference) because of logged)
difflsmeans(Resp.test)

mean.resp <- aggregate(Respiration~Treatment+Month, df.resp.long, mean)
mean.resp$se <- aggregate(Respiration~Treatment+Month, df.resp.long, se)$Respiration
mean.resp$month <- as.Date(rep(c("2020-08-15","2020-07-15", "2020-06-15", "2020-11-15", "2020-09-15"), each=2)) 

resp <- ggplot(mean.resp, aes(x=month, y=Respiration, shape=Treatment))+
  geom_line(aes(linetype = Treatment)) +
  geom_point()+
  geom_errorbar(aes(ymax = Respiration + se, ymin = Respiration - se), width = 0.1, color = "Gray25") +
  scale_shape_manual(values=c(21,19)) +
  scale_y_continuous(limits = c(0,2000), expand = c(0, 0)) +
  xlab("Month") +
  labs(y = expression(paste("Soil respiration rate"))) +
  theme_cowplot() +
  theme(legend.position = c(0.05, 0.2),
        legend.title = element_blank()) + draw_plot_label("(c)", x= as.Date("2020-06-05"), y=2000 )
resp
fig3 <- plot_grid(tem, moi, resp, ncol = 1, align = "v")
save_plot("fig3.png", fig3, ncol=1, base_height = 6.71, base_width = 4)

    # temp vs resp
        # august
cor.test(df.resp$CO2_August, df.resp$Soil_temperature_August)
cor(df.resp$CO2_November, df.resp$Soil_temperature_November)
cor.test(df.resp$CO2_August, df.resp$Soil_moisture_August)
resp.aug <- aggregate(CO2_August~Plot_ID, df.resp, mean)
temp.aug <- aggregate(Soil_temperature_August~Plot_ID, df.resp, mean)
cor.test(resp.aug[,2], temp.aug[,2])

    # fung:bacteria vs respiration in september
resp.sep <- aggregate(CO2_September~Plot_ID, df.resp, mean)
cor.test(df$Fungibac, resp.sep[,2])
summary(lm(df$Fungibac~resp.sep[,2]))
 


# Tree growth ####
df.tree <- read_excel("Lundgren_etal_tree_core_data.xlsx")

# diff over 2014-2015, ie year 2 and 3 after introduction
df.tree$yr.2.3.post <- rowMeans(df.tree[,12:13])
df.tree$yr.0.1.post <- rowMeans(df.tree[,14:15])
df.tree$Treatment <- factor(df.tree$Treatment, levels=c("Reference", "Enclosure"))
Tree.test2 <- lm(yr.2.3.post ~ scale(yr.0.1.post,scale=F) + Environment*Treatment, 
                 df.tree, )
summary(Tree.test2)
Anova(Tree.test2, type = 2)

# Or use relative growth, which gives the same result
ref.yr = which(colnames(df.tree)=="2006")
end.yr = which(colnames(df.tree)=="2020")
df.tree[,ref.yr:end.yr] <- df.tree[,ref.yr:end.yr]/ df.tree$`2006`
df.tree$RelaTree <- rowMeans(df.tree[,7:15], na.rm=TRUE)

Tree.test <- lm(RelaTree ~ Treatment, df.tree)
summary(Tree.test)
Anova(Tree.test, type = 2)

Tree.Env <- lm(RelaTree ~ Treatment * Environment, df.tree)
Anova(Tree.Env, type = 2)

# Plot relative change over time
df.tree.plot <- gather(df.tree, key = "Year", value = "width",
                       '2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014',
                       '2015', '2016', '2017', '2018', '2019', '2020')

mean.tree <- aggregate(width~Treatment+Year, df.tree.plot, mean)
mean.tree$se <- aggregate(width~Treatment+Year, df.tree.plot, se)$width
mean.tree$Year <- as.Date(rep(c("2006-01-02","2007-01-02", "2008-01-02", "2009-01-02",
                                "2010-01-02", "2011-01-02", "2012-01-02", "2013-01-02",
                                "2014-01-02", "2015-01-02", "2016-01-02", "2017-01-02",
                                "2018-01-02", "2019-01-02", "2020-01-02"), each=2)) 
mean.tree$Treatment <- factor(mean.tree$Treatment, levels=c("Reference", "Enclosure"))

fig9 <- ggplot(mean.tree, aes(x=Year, y=width, shape=Treatment))+
  geom_line(aes(linetype = Treatment)) +
  geom_point()+
  geom_errorbar(aes(ymax = width + se, ymin = width - se), width = 0.1, color = "Gray25") +
  geom_vline(xintercept=as.numeric(mean.tree$Year[14]), linetype=2) +
  scale_shape_manual(values=c(21,19)) +
  scale_y_continuous(limits = c(0.5,2.1)) +
  scale_x_date(date_labels = "%Y", date_breaks = "2 years", limits = c(min(as.Date("2005-09-01")), max = max(as.Date("2020-01-03")))) +
  xlab("Year") +
  ylab(expression(Relative~basal~area~increment~'('~cm^2~')')) +
  theme_cowplot() +
  theme(legend.position = c(0.05, 0.2),
        legend.title = element_blank()) + draw_plot_label(" ", x= as.Date("2006-01-01"), y=0.7 )
fig9

save_plot("fig9.png", fig9, ncol=1, base_height = 3.71, base_width = 5)
