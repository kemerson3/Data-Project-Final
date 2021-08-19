library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(ggbiplot)
library(moments) #stats
library(psych)
library(pastecs)
library(ggplot2)
library(car)
library(lme4)

baseline <- read_csv("C:\\Users\\kjeme\\OneDrive\\Desktop\\Data Project\\Emerson Microbial Exp1_baseline.csv")
baseline <- select(baseline, -X25)
for (i in 2:ncol(baseline)) {
  if (str_detect(baseline[1,i],"%")) {
    baseline[,i] <-as.numeric(strsplit(c(baseline[[i]]), "%"))
  }
}

base.df <- baseline %>%
  select(`Av. Speed (mm/s)`, `Total Distance (mm)`, `Av. Accel (mm/s^2)`, `Mobility Rate (%)`, `Exploration Rate (%)`, `Total Distance (mm)`, `Time In Center (m:s)`, `Tot. Time Frozen (m:s)`, `Exploration Rate (%)`  )

base.df$`Time In Center (m:s)`<- (as.numeric(base.df$`Time In Center (m:s)`))/3600

base.df$`Tot. Time Frozen (m:s)`<- (as.numeric(base.df$`Tot. Time Frozen (m:s)`))/3600

KMO(base.df)
#Pass assumption

cortest.bartlett(base.df)
#Pass assumption. High chisq

pca1 <- principal(base.df, nfactors = 3, rotate = "varimax")
summary(pca1)

qplot(c(1:7), pca1$values) +
  geom_line() +
  xlab("Principal Component") +
  ylab("Eigenvalue") +
  ylim(0,6)
#scree plot confirms that 3 PCs fit this data set best

Factor <- pca1$scores
df <- cbind(base.df, Factor)
ID <- baseline$ID
MicroTrtmt <- baseline$MicroTrtmt
df <- cbind(df, ID)
df <- cbind(df, MicroTrtmt)
df$MicroTrtmt = factor(df$MicroTrtmt)
#added in some ID tags and treatment to the data frame. may have to remove

# treatment <- c(rep("Colonized", 1),rep("Depleted", 2), rep("Colonized", 1), rep("Depleted", 1),
#                rep("Colonized", 1), rep("Depleted", 1), rep("Colonized", 1), rep("Depleted", 1), rep("Colonized", 1), rep("Depleted", 1),
#                rep("Colonized", 1), rep("Depleted", 1), rep("Colonized", 1), rep("Depleted", 1),
#                rep("Colonized", 1), rep("Depleted", 1), rep("Colonized", 1), rep("Depleted", 1), rep("Colonized", 1), rep("Colonized", 1),
#                rep("Depleted", 1))

pca1$loadings
#Based on our loadings:
#PCA1 = Speed, Distance, Accel, Mobility rate and total time frozen
#PCA2 = Exploration Rate
#PCA3 = Time in Center
#Cumulative 97% of variance explained

PC1.aov <- aov(RC1~MicroTrtmt, data = df)
summary(PC1.aov)

PC2.aov <- aov(RC2~MicroTrtmt, data = df)
summary(PC2.aov)

PC3.aov <- aov(RC3~MicroTrtmt, data = df)
summary(PC3.aov)
#None look significant. Lets graph

#PC1 Boxplot
ggplot(df, aes(x = MicroTrtmt, y = RC1, fill = MicroTrtmt))+
  geom_boxplot() +
  theme_classic() +
  scale_fill_manual(values = c("springgreen4", "dodgerblue1"),
                    name = "Microbial Treatment",
                    labels = c("Colonized", "Depleted")) +
  labs(x = "Microbial Treatment", y = "PC1") +
  theme(panel.background = element_rect(fill = "white", colour                      = "black")) +
  theme(aspect.ratio = 1) +
  theme(axis.text = element_text(face = "bold", size = 12))  +
  theme(axis.title = element_text(face = "bold", size = 14)) +
  theme(legend.position = "bottom") 

#PC2 Box Plot
ggplot(df, aes(x = MicroTrtmt, y = RC2, fill = MicroTrtmt))+
  geom_boxplot() +
  theme_classic() +
  scale_fill_manual(values = c("springgreen4", "dodgerblue1"),
                    name = "Microbial Treatment",
                    labels = c("Colonized", "Depleted")) +
  labs(x = "Microbial Treatment", y = "PC2") +
  theme(panel.background = element_rect(fill = "white", colour                      = "black")) +
  theme(aspect.ratio = 1) +
  theme(axis.text = element_text(face = "bold", size = 12))  +
  theme(axis.title = element_text(face = "bold", size = 14)) +
  theme(legend.position = "bottom")

#PC3 Boxplot
ggplot(df, aes(x = MicroTrtmt, y = RC3, fill = MicroTrtmt))+
  geom_boxplot() +
  theme_classic() +
  scale_fill_manual(values = c("springgreen4", "dodgerblue1"),
                    name = "Microbial Treatment",
                    labels = c("Colonized", "Depleted")) +
  labs(x = "Microbial Treatment", y = "PC3") +
  theme(panel.background = element_rect(fill = "white", colour                      = "black")) +
  theme(aspect.ratio = 1) +
  theme(axis.text = element_text(face = "bold", size = 12))  +
  theme(axis.title = element_text(face = "bold", size = 14)) +
  theme(legend.position = "bottom")

#PC1 + PC2 Biplot
ggplot(df, aes(x= RC1, y = RC2, color = MicroTrtmt)) +
  geom_point() +
  theme_classic() +
  labs(x = "PC1 (67% Variance Explained)", y = "PC2 (15% Variance Explained)") +
  scale_color_manual(values = c("forestgreen", "dodgerblue" ),
                     name = "Microbial Treatment",
                     labels = c("Colonized", "Depleted")) +
  stat_ellipse()+
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(aspect.ratio = 1) +
  theme(axis.text = element_text(face = "bold", size = 12))  +
  theme(axis.title = element_text(face = "bold", size = 14)) +
  theme(legend.position = "bottom")

#PC1 + PC3 Biplot
ggplot(df, aes(x= RC1, y = RC3, color = MicroTrtmt)) +
  geom_point() +
  theme_classic() +
  labs(x = "PC1 (67% Variance Explained)", y = "PC3 (15% Variance Explained)") +
  scale_color_manual(values = c("forestgreen", "dodgerblue" ),
                     name = "Microbial Treatment",
                     labels = c("Colonized", "Depleted")) +
  stat_ellipse()+
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(aspect.ratio = 1) +
  theme(axis.text = element_text(face = "bold", size = 12))  +
  theme(axis.title = element_text(face = "bold", size = 14)) +
  theme(legend.position = "bottom")

###Take home: With the varimax rotation, there were
### no significant impacts of microtrtmt on baseline behavior