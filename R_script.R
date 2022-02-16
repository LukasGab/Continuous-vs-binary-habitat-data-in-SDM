# Code Authors: Lukas Gabor & Petr Keil, 2021

# ===== Library packages =====
library(PresenceAbsence)
library(xlsx)
library(dplyr)
library(ggplot2)
library(cowplot)


# ===== Load functions =====
# ArcSine function 
ArcSin<-function(x){
  if (!is.vector(x) || !(mode(x)=="numeric")){
    stop("'x' must be a numeric vector")
  }
  if (any(x < 0)){
    stop("Some values of x are < 0")
  }
  if (any(x > 1)){
    stop("Some values of x are > 1")
  }
  else{
    return(asin(sqrt(abs(x))))
  }
}

# GLM function
glm.models <- function(formula, data){
  
  metrics.all <- numeric()
  
  row_number <- nrow(data)
  group <- rep(1:5, length = row_number)
  group <- sample(group)
  D <- numeric(1)
  
  for(i in 1:5) {
    
    m <- glm(formula,family= binomial, data=data[group != i,])
    
    pred <-predict(m, newdata=data[group == i,], type="response")
    
    # D2
    D2 <- (m$null.deviance - m$deviance) / m$null.deviance
    n <- length(m$fitted.values)
    p <- length(m$coefficients)
    D2 <- 1 - ((n - 1) / (n - p)) * (1 - D2)
    
    # Model performance
    D <- as.vector(pred)
    PA <- data[group == i,c("kvadrat","occ")]
    DAT <- cbind(PA,D)
    
    optresh <- optimal.thresholds(DAT, opt.methods="ObsPrev", which.model = 1)
    t <-optresh[,2]
    #cmx <-  cmx(DAT, threshold = t, which.model = 1, na.rm = FALSE);cmx
    #pcc(cmx, st.dev = TRUE)
    performance <- presence.absence.accuracy(DAT, threshold = t, find.auc = TRUE, st.dev = TRUE)
    
    aic <- m$aic
    
    r2<- with(summary(m), 1 - deviance/null.deviance)
    
    tss <- ((performance$sensitivity + performance$specificity)-1)
    
    metrics <- data.frame(performance$sensitivity, performance$specificity, performance$AUC, D2, aic, r2, tss)
    metrics.all <- rbind(metrics.all, metrics)
    
  }
  
  results <- data.frame(sensitivity = mean(metrics.all$performance.sensitivity),
                        specificity = mean(metrics.all$performance.specificity),
                        AUC = mean(metrics.all$performance.AUC),
                        TSS = mean(metrics.all$tss),
                        D2 = mean(metrics.all$D2),
                        aic = mean(metrics.all$aic),
                        R2 = mean(metrics.all$r2),
                        formula = formula,
                        species = spec)
}


# ===== Load data =====
data <- na.omit(read.csv("data.csv", header=T, sep=";"))


# ===== Process data =====
# Summarize occurrences by species and calculate species prevalence
sum_by_species <- data %>% 
  group_by(species) %>% 
  summarise(occ = sum(occ))

sum_by_species <- sum_by_species %>%
  mutate(species_prevalence = sum_by_species$occ / 628)

# Water bodies
data$logwater <- log(data$water + 1) # log

data$sqrtwater <- sqrt(data$water) # sqrt

data$arcsinwater <- ArcSin(data$water / 100) # arcsin

data$bwater[data$water==0] <- 0 # threshold (binary variable)

data$bwater[data$water>0] <- 1

# Forests
data$logforest <- log(data$forest + 1)

data$sqrtforest <- sqrt(data$forest)

data$arcsinforest <- ArcSin(data$forest / 100)

data$forest10 <- ifelse(data$forest<10, 0,1) # 10% threshold (binary variable)

data$forest20 <- ifelse(data$forest<20, 0,1) # 20% threshold 

data$forest30 <- ifelse(data$forest<30, 0,1) # 30% threshold 

data$forest40 <- ifelse(data$forest<40, 0,1) # 40% threshold 

data$forest50 <- ifelse(data$forest<50, 0,1) # 50% threshold

# Filter data for wetlands respectively for forest birds
data.wet <- data[data$habitat == "water", ] # data for wetland birds

data.for <- data[data$habitat == "forest", ] # data for forest birds


# ===== SDM models =====
# Water birds
spec.list <- as.vector(unique(data.wet$species))

results.waterBirds <- numeric()

for(spec in spec.list){
  
  sub.data <- data.wet[data.wet$species == spec, ] # data for one species
  
  sub.results <- rbind(glm.models(formula = "occ~1", data = sub.data),
                       glm.models(formula = "occ~temperature + precipitation + artificalSurfaces + agricultur", data = sub.data),
                       glm.models(formula = "occ~water + temperature + precipitation + artificalSurfaces + agricultur", data = sub.data),
                       glm.models(formula = "occ~bwater + temperature + precipitation + artificalSurfaces + agricultur", data = sub.data),
                       glm.models(formula = "occ~logwater + temperature + precipitation + artificalSurfaces + agricultur", data = sub.data),
                       glm.models(formula = "occ~sqrtwater + temperature + precipitation + artificalSurfaces + agricultur", data = sub.data),
                       glm.models(formula = "occ~arcsinwater + temperature + precipitation + artificalSurfaces + agricultur", data = sub.data))
  
  results.waterBirds <- rbind(results.waterBirds, sub.results)
  
}


# Forest birds
spec.list <- as.vector(unique(data.for$species))

results.forestBirds <- numeric()

for(spec in spec.list){
  
  sub.data <- data.for[data.for$species == spec, ] # data for one species
  
  sub.results <- rbind(glm.models(formula = "occ~1", data = sub.data),
                       glm.models(formula = "occ~temperature + precipitation + artificalSurfaces + agricultur", data = sub.data),
                       glm.models(formula = "occ~forest + temperature + precipitation + artificalSurfaces + agricultur", data = sub.data),
                       glm.models(formula = "occ~logforest + temperature + precipitation + artificalSurfaces + agricultur", data = sub.data),
                       glm.models(formula = "occ~sqrtforest + temperature + precipitation + artificalSurfaces + agricultur", data = sub.data),
                       glm.models(formula = "occ~arcsinforest + temperature + precipitation + artificalSurfaces + agricultur", data = sub.data),
                       glm.models(formula = "occ~forest10 + temperature + precipitation + artificalSurfaces + agricultur", data = sub.data),
                       glm.models(formula = "occ~forest20 + temperature + precipitation + artificalSurfaces + agricultur", data = sub.data),
                       glm.models(formula = "occ~forest30 + temperature + precipitation + artificalSurfaces + agricultur", data = sub.data),
                       glm.models(formula = "occ~forest40 + temperature + precipitation + artificalSurfaces + agricultur", data = sub.data),
                       glm.models(formula = "occ~forest50 + temperature + precipitation + artificalSurfaces + agricultur", data = sub.data))
  
  results.forestBirds <- rbind(results.forestBirds, sub.results)
  
}


# ===== Summarize and save results =====
results.forestBirds$habitat <- "forest"
results.waterBirds$habitat <- "water"

names <- read.table("birds_name.csv", sep = ";", header = T) # upload file with species name

results.forestBirds <- results.forestBirds %>% # add english names
  inner_join(names, by = "species")

results.forestBirds <- results.forestBirds %>%
  inner_join(sum_by_species, by = "species") # add information about prevalence and sample size

results.waterBirds <- results.waterBirds %>% 
  inner_join(names, by = "species")

results.waterBirds <- results.waterBirds %>%
  inner_join(sum_by_species, by = "species") 


results.sum <- rbind(results.forestBirds, results.waterBirds)


write.table(results.sum, "results_sum.csv", sep = ";", row.names = F)


# ===== Plots results =====
# Get difference between models performance
subset_waterBirds <- data.frame(results.waterBirds[results.waterBirds$formula == "occ~bwater + temperature + precipitation + artificalSurfaces + agricultur", ],
                                results.waterBirds[results.waterBirds$formula == "occ~water + temperature + precipitation + artificalSurfaces + agricultur", ])


df_waterBirds <- data.frame(D.sensitivity = subset_waterBirds$sensitivity - subset_waterBirds$sensitivity.1,
                            D.specificity = subset_waterBirds$specificity - subset_waterBirds$specificity.1,
                            D.AUC = subset_waterBirds$AUC - subset_waterBirds$AUC.1,
                            D.TSS = subset_waterBirds$TSS - subset_waterBirds$TSS.1,
                            D.D2 = subset_waterBirds$D2 - subset_waterBirds$D2.1,
                            D.AIC = subset_waterBirds$aic - subset_waterBirds$aic.1,
                            D.R2 = subset_waterBirds$R2 - subset_waterBirds$R2.1,
                            species = subset_waterBirds$species,
                            species2 = subset_waterBirds$species2,
                            habitat = subset_waterBirds$habitat)


subset_forestBirds <- data.frame(results.forestBirds[results.forestBirds$formula == "occ~forest40 + temperature + precipitation + artificalSurfaces + agricultur", ],
                                 results.forestBirds[results.forestBirds$formula == "occ~forest + temperature + precipitation + artificalSurfaces + agricultur", ])


df_forestbirds <- data.frame(D.sensitivity = subset_forestBirds$sensitivity - subset_forestBirds$sensitivity.1,
                             D.specificity = subset_forestBirds$specificity - subset_forestBirds$specificity.1,
                             D.AUC = subset_forestBirds$AUC - subset_forestBirds$AUC.1,
                             D.TSS = subset_forestBirds$TSS - subset_forestBirds$TSS.1,
                             D.D2 = subset_forestBirds$D2 - subset_forestBirds$D2.1,
                             D.AIC = subset_forestBirds$aic - subset_forestBirds$aic.1,
                             D.R2 = subset_forestBirds$R2 - subset_forestBirds$R2.1,
                             species = subset_forestBirds$species,
                             species2 = subset_forestBirds$species2,
                             habitat = subset_forestBirds$habitat)


## Plot results
# AUC
auc1 <- ggplot(df_waterBirds, aes(reorder(species2, D.AUC), D.AUC, fill = species2)) + 
  geom_bar(stat = 'identity', position = 'identity', fill="grey47", width = 0.75) +
  coord_flip() +
  theme_gray() +
  scale_fill_grey() +
  scale_y_continuous(limits = c(-.02, 0.15), breaks = c(0, 0.05, 0.1)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.y = element_text(size = 6)) +
  labs(title="Water Bird species",x="", y = "") +
  guides(fill = 'none') 


auc2 <- ggplot(df_forestbirds, aes(reorder(species2, D.AUC), D.AUC, fill = species2)) + 
  geom_bar(stat = 'identity', position = 'identity', fill="grey47", width = 0.75) +
  coord_flip() +
  theme_gray() +
  scale_fill_grey() +
  scale_y_continuous(limits = c(-.1, 0.05), breaks = c(-.1, -.05, 0)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.y = element_text(size = 6)) +
  labs(title="Forest Bird species",x="", y = "") +
  guides(fill = 'none')


#R2
r1 <- ggplot(df_waterBirds, aes(reorder(species2, D.R2), D.R2, fill = species2)) + 
  geom_bar(stat = 'identity', position = 'identity', fill="grey47", width = 0.75) +
  coord_flip() +
  theme_gray() +
  scale_fill_grey() +
  scale_y_continuous(limits = c(-.008, 0.15), breaks = c(0, .05, 0.1, 0.15)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.y = element_text(size = 6)) +
  labs(x="", y = "") +
  guides(fill = 'none') 


r2 <- ggplot(df_forestbirds, aes(reorder(species2, D.R2), D.R2, fill = species2)) + 
  geom_bar(stat = 'identity', position = 'identity', fill="grey47", width = 0.75) +
  coord_flip() +
  theme_gray() +
  scale_fill_grey() +
  scale_y_continuous(limits = c(-.085, 0.01), breaks = c(-0.06, -.03, 0)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.y = element_text(size = 6)) +
  labs(x="", y = "") +
  guides(fill = 'none')

g <- plot_grid(auc2, auc1, r2, r1, ncol = 2, nrow = 2,
               label_size = 12,
               label_fontfamily = "serif",
               label_fontface = "plain",
               labels = c('AUC', '', 'R2', ''))

# Save plots
ggsave("results.png", g, dpi=300, height = 15, width = 25, units = "cm")
