library(tidyverse) #package collection for data science
library(readxl) #read excel files
library(GGally) #for correlation matrix
library(ggcorrplot) #for correlation matrix visualization
library(psych) #for tetrachloric function
library(MASS) #stepwise regression
library(ltm) #for point-biserial correlation
library(rcompanion) #Cramer's V
library(broom) #View AIC Variable Selection
library(carData) #for car package
library(car) #VIF 
library(mice) #Multiple Imputation by Chained Equations

#-----------------------------------------------------------------------------------------------------------------------#
#Read in Data Sets
CARD9_Serum <- read_excel("~/Desktop/Computational CARD9/Mycobiome serum_BROAD C. Albican July 19 2022.xlsx")
CARD9_Prism <- read_excel("~/Desktop/Computational CARD9/Mycobiome_sample_tracker_CK.xlsx")
CARD9_Geno1 <- read_excel("~/Desktop/Computational CARD9/Mycobiome serum_BROAD 1 25 22 (3)_with_genotypes.xlsx") 
CARD9_Geno2 <- read_excel("~/Desktop/Computational CARD9/Mycobiome serum_BROAD 3 10 22 (3)_with_genotypes.xlsx")

#Merging into one Data File
CARD9_Serum2 <- CARD9_Serum %>% 
  dplyr::select(`PRISM ID`, DX, `anti-C. Albican IgG (OD relative to control)`, `# days between stool/serum samples`) %>%
  rename(Diagnosis = DX, `anti-C.A_IgG` = `anti-C. Albican IgG (OD relative to control)`, `#days_stool/serum` = `# days between stool/serum samples`) %>%
  dplyr::mutate(Diagnosis = dplyr::recode(Diagnosis,"Crohn's Disease" = "CD", "Ulcerative Colitis" = "UC", "IndeterminateColitis" = "IC", "Indeterminate Colitis" = "IC", "Crohn's disease" = "CD"))
CARD9_Prism2 <- CARD9_Prism %>%
  dplyr::select(`Patient_ID`, `Gender`, `Diagnosis`, `Height`,`Weight`, `BMI`, everything() ) %>%
  rename(`PRISM ID` = `Patient_ID`)
CARD9_Merged <- merge(CARD9_Serum2,CARD9_Prism2, by = "PRISM ID") %>%
  dplyr::select(`PRISM ID`, `Diagnosis.x`, `Diagnosis.y`, everything())

#Converting variables from character to numeric
CARD9_Geno1$CARD9_S12N <-  as.numeric(CARD9_Geno1$CARD9_S12N)
CARD9_Geno2$CARD9_S12N <-  as.numeric(CARD9_Geno2$CARD9_S12N)

#Merging genotype files into one
CARD9_Geno1and2 <- merge(CARD9_Geno1, CARD9_Geno2, by = "PRISM ID", all = TRUE) %>%
  dplyr::select(`PRISM ID`, CARD9_S12N.x, CARD9_S12N.y) %>%
  filter(CARD9_S12N.x==0 | CARD9_S12N.x==1 | CARD9_S12N.x==2 | CARD9_S12N.y==0 | CARD9_S12N.y==1 | CARD9_S12N.y==2) %>%
  mutate(CARD9_S12N.x = replace_na(CARD9_S12N.x, 0)) %>%
  mutate(CARD9_S12N.y = replace_na(CARD9_S12N.y, 0))
CARD9_Geno1and2$CARD9_S12N <- CARD9_Geno1and2$CARD9_S12N.x + CARD9_Geno1and2$CARD9_S12N.y
CARD9_Genotype <- CARD9_Geno1and2 %>%
  dplyr::select(`PRISM ID`, CARD9_S12N)

#Merging into Final Data Set (20 total variables)
CARD9_CAgenotype <- merge(CARD9_Merged,CARD9_Genotype, by = "PRISM ID") %>%
  mutate(Diagnosis =`Diagnosis.y`, `Sample_Collection_Age` = `Age_at_time_of_sample_collection`, `Diagnosis_Age` = `Age_at_diagnosis`,
         `Collection_DiseaseActivity` = `Disease_activity_at_time_of_collection`, `TX_Antidiarrheal` = `TX_antidiarrheal (lomotil, imodium, DTO)`,
         `TX_ASA` = `TX_5_ASA (mesalamine, balasalizide, azulfidine, olsalazine)`, `TX_Biologics`
         =`TX_biologics (infliximab, adalimumab, certalizumab, natalizumab, golimumab, vedolizaumab ustekinumab, tofacitinib)`,
         `TX_Steroid` = `TX_steroid (prednisone, budesonide, medrol, IV steroid, [cortenemas, cortifoam, proctofoam]))`,
         `TX_Immunomodulators` = `Immunomodulators (azathioprine, methotrexate, mercaptopurine)`, 
         `TX_Antibiotic` = `TX_antibiotic (flaygl, cipro, xifaxin, levaquin)`,`Smoking` =`Smoking (Never, Current, Former)` ) %>%
  dplyr::select(`PRISM ID`,`anti-C.A_IgG`, Height, Weight, BMI, `#days_stool/serum`,
                `Sample_Collection_Age`,`Diagnosis_Age`,`Collection_DiseaseActivity`, CARD9_S12N, `Diagnosis`, Gender,`Disease_indices`, `TX_Antidiarrheal`, `TX_ASA`, `TX_Biologics`,
                `TX_Steroid`, `TX_Immunomodulators`,`TX_Antibiotic`,`Smoking` )

#Converting character variables to numeric
CARD9_CAgenotype$Height <- as.numeric(CARD9_CAgenotype$Height)
CARD9_CAgenotype$Weight <- as.numeric(CARD9_CAgenotype$Weight)
CARD9_CAgenotype$BMI <- as.numeric(CARD9_CAgenotype$BMI)
CARD9_CAgenotype$`Sample_Collection_Age`  <- as.numeric(CARD9_CAgenotype$`Sample_Collection_Age`)
CARD9_CAgenotype$`Diagnosis_Age`  <- as.numeric(CARD9_CAgenotype$`Diagnosis_Age`)
CARD9_CAgenotype$`Collection_DiseaseActivity` <- as.numeric(CARD9_CAgenotype$`Collection_DiseaseActivity`)
glimpse(CARD9_CAgenotype)

#-----------------------------------------------------------------------------------------------------------------------#
#look at pattern of missing data
md.pattern(CARD9_CAgenotype)

#complete data and missing data
CARD9_Final <- CARD9_CAgenotype[complete.cases(CARD9_CAgenotype), ]
CARD9_Missing <- CARD9_CAgenotype[!complete.cases(CARD9_CAgenotype), ]

#-----------------------------------------------------------------------------------------------------------------------#
#Overview of data
summary(CARD9_CAgenotype)
str(CARD9_CAgenotype)
#Histograms of Continuous Variables
par(mfrow=c(2,2))
hist(CARD9_CAgenotype$`anti-C.A_IgG`, breaks=25, col=rgb(0,0,1,0.7), xlab="anti-Candida Albicans IgG", 
     ylab="Frequency", main="Distribution of anti-Candida Albicans IgG" )
hist(CARD9_CAgenotype$Height, breaks=25, col=rgb(0,0,1,0.7), xlab="Height", 
     ylab="Frequency", main="Distribution of Height" )
hist(CARD9_CAgenotype$Weight, breaks=25, col=rgb(0,0,1,0.7), xlab="Weight", 
     ylab="Frequency", main="Distribution of Weight" )
hist(CARD9_CAgenotype$BMI, breaks=25, col=rgb(0,0,1,0.7), xlab="BMI", 
     ylab="Frequency", main="Distribution of BMI" )
hist(CARD9_CAgenotype$`Sample_Collection_Age`, breaks=25, col=rgb(0,0,1,0.7), xlab="Age at Sample Collection", 
     ylab="Frequency", main="Distribution of Age at Sample Collection" )
hist(CARD9_CAgenotype$`Diagnosis_Age`, breaks=25, col=rgb(0,0,1,0.7), xlab="Age at Diagnosis", 
     ylab="Frequency", main="Distribution of Age at Diagnosis" )
hist(CARD9_CAgenotype$Collection_DiseaseActivity, breaks=25, col=rgb(0,0,1,0.7), xlab="Disease Activity at Sample Collection", 
     ylab="Frequency", main="Distribution of Disease Activity at Sample Collection")
hist(CARD9_CAgenotype$`#days_stool/serum`, breaks=25, col=rgb(0,0,1,0.7), xlab="# days between stool and serum samples", 
     ylab="Frequency", main="Distribution of # days between stool and serum samples")

#Box-plots of Categorical variables
CARD9_CAgenotype %>% ggplot(aes(Diagnosis, `anti-C.A_IgG`, fill=Diagnosis)) +
  geom_boxplot() + theme_bw() + labs(title="Diagnosis of anti-C. Albican IgG", x="Diagnosis", y="anti-C. Albican_IgG") +
  guides(fill=guide_legend(title="Diagnosis"))
CARD9_CAgenotype %>% ggplot(aes(Gender, `anti-C.A_IgG`, fill=Gender)) +
  geom_boxplot() + theme_bw() + labs(title="Gender Distribution of anti-C. Albican IgG", x="Gender", y="anti-C. Albican_IgG")
CARD9_CAgenotype %>% ggplot(aes(TX_Antidiarrheal, `anti-C.A_IgG`, fill=TX_Antidiarrheal)) +
  geom_boxplot() + theme_bw() + labs(title="Antidiarrheal Treatment Distribution of anti-C. Albican IgG", x="Antidiarrheal Treatment", y="anti-C. Albican_IgG") +
  guides(fill=guide_legend(title="Antidiarrheal Treatment"))

#-----------------------------------------------------------------------------------------------------------------------#
#Scatter-plot Matrix of all Variables, hard to visualize all the plots
ggpairs(CARD9_CAgenotype, upper = list(continuous = wrap("points", alpha = 0.3, size=0.1)),
        lower = list(continuous = wrap('cor', size = 4)))

#Scatter-plot Matrix of continuous Variables
CARD9_ContinuousVariables <- CARD9_CAgenotype %>%
  dplyr::select(`anti-C.A_IgG`,`#days_stool/serum`, Height, Weight, BMI,
                `Sample_Collection_Age`,`Diagnosis_Age`,`Collection_DiseaseActivity`)
ggpairs(CARD9_ContinuousVariables, upper = list(continuous = wrap("points", alpha = 0.3, size=0.1)),
        lower = list(continuous = wrap('cor', size = 4)))

#Correlation Matrix of continuous variables
model.matrix(~0+., data=CARD9_ContinuousVariables) %>% 
  cor(use="pairwise.complete.obs") %>% 
  ggcorrplot(show.diag = F, type="lower", lab=TRUE, lab_size=2)

#VIF of continuous variables, does not work for categorical variables
par(mfrow=c(1,1))
model_cov = lm(data=CARD9_ContinuousVariables, `anti-C.A_IgG` ~ `#days_stool/serum`+ Height + Weight+ BMI+ `Sample_Collection_Age`+`Diagnosis_Age`+`Collection_DiseaseActivity`)
vif_values = vif(model_cov)
barplot(vif_values, main = "VIF Values", col = c("steelblue"))
abline(10,0)

#-----------------------------------------------------------------------------------------------------------------------#
#Tetrachloric Correlation: 2 Binary Variables

CARD9_Tetrachloric <- CARD9_Final[complete.cases(CARD9_Final), ] %>%
  dplyr::select(`Diagnosis`, Gender, `Disease_indices`, `TX_Antidiarrheal`, `TX_ASA`, `TX_Biologics`,
                `TX_Steroid`, `TX_Immunomodulators`,`TX_Antibiotic`)
tetrachoric(table(CARD9_Tetrachloric$Diagnosis,CARD9_Tetrachloric$Gender))
tetrachoric(table(CARD9_Tetrachloric$Diagnosis,CARD9_Tetrachloric$Disease_indices))
tetrachoric(table(CARD9_Tetrachloric$Diagnosis,CARD9_Tetrachloric$TX_Antidiarrheal))
tetrachoric(table(CARD9_Tetrachloric$Diagnosis,CARD9_Tetrachloric$TX_ASA))
tetrachoric(table(CARD9_Tetrachloric$Diagnosis,CARD9_Tetrachloric$TX_Biologics))
tetrachoric(table(CARD9_Tetrachloric$Diagnosis,CARD9_Tetrachloric$TX_Steroid))
tetrachoric(table(CARD9_Tetrachloric$Diagnosis,CARD9_Tetrachloric$TX_Immunomodulators))
tetrachoric(table(CARD9_Tetrachloric$Diagnosis,CARD9_Tetrachloric$TX_Antibiotic))
tetrachoric(table(CARD9_Tetrachloric$Gender,CARD9_Tetrachloric$Disease_indices))
tetrachoric(table(CARD9_Tetrachloric$Gender,CARD9_Tetrachloric$TX_Antidiarrheal))
tetrachoric(table(CARD9_Tetrachloric$Gender,CARD9_Tetrachloric$TX_ASA))
tetrachoric(table(CARD9_Tetrachloric$Gender,CARD9_Tetrachloric$TX_Biologics))
tetrachoric(table(CARD9_Tetrachloric$Gender,CARD9_Tetrachloric$TX_Steroid))
tetrachoric(table(CARD9_Tetrachloric$Gender,CARD9_Tetrachloric$TX_Immunomodulators))
tetrachoric(table(CARD9_Tetrachloric$Gender,CARD9_Tetrachloric$TX_Antibiotic))
tetrachoric(table(CARD9_Tetrachloric$Disease_indices,CARD9_Tetrachloric$TX_Antidiarrheal))
tetrachoric(table(CARD9_Tetrachloric$Disease_indices,CARD9_Tetrachloric$TX_ASA))
tetrachoric(table(CARD9_Tetrachloric$Disease_indices,CARD9_Tetrachloric$TX_Biologics))
tetrachoric(table(CARD9_Tetrachloric$Disease_indices,CARD9_Tetrachloric$TX_Steroid))
tetrachoric(table(CARD9_Tetrachloric$Disease_indices,CARD9_Tetrachloric$TX_Immunomodulators))
tetrachoric(table(CARD9_Tetrachloric$Disease_indices,CARD9_Tetrachloric$TX_Antibiotic))
tetrachoric(table(CARD9_Tetrachloric$TX_Antidiarrheal,CARD9_Tetrachloric$TX_ASA))
tetrachoric(table(CARD9_Tetrachloric$TX_Antidiarrheal,CARD9_Tetrachloric$TX_Biologics))
tetrachoric(table(CARD9_Tetrachloric$TX_Antidiarrheal,CARD9_Tetrachloric$TX_Steroid))
tetrachoric(table(CARD9_Tetrachloric$TX_Antidiarrheal,CARD9_Tetrachloric$TX_Immunomodulators))
tetrachoric(table(CARD9_Tetrachloric$TX_Antidiarrheal,CARD9_Tetrachloric$TX_Antibiotic))
tetrachoric(table(CARD9_Tetrachloric$TX_ASA,CARD9_Tetrachloric$TX_Biologics))
tetrachoric(table(CARD9_Tetrachloric$TX_ASA,CARD9_Tetrachloric$TX_Steroid))
tetrachoric(table(CARD9_Tetrachloric$TX_ASA,CARD9_Tetrachloric$TX_Immunomodulators))
tetrachoric(table(CARD9_Tetrachloric$TX_ASA,CARD9_Tetrachloric$TX_Antibiotic))
tetrachoric(table(CARD9_Tetrachloric$TX_Biologics,CARD9_Tetrachloric$TX_Steroid))
tetrachoric(table(CARD9_Tetrachloric$TX_Biologics,CARD9_Tetrachloric$TX_Immunomodulators))
tetrachoric(table(CARD9_Tetrachloric$TX_Biologics,CARD9_Tetrachloric$TX_Antibiotic))
tetrachoric(table(CARD9_Tetrachloric$TX_Steroid,CARD9_Tetrachloric$TX_Immunomodulators))
tetrachoric(table(CARD9_Tetrachloric$TX_Steroid,CARD9_Tetrachloric$TX_Antibiotic))
tetrachoric(table(CARD9_Tetrachloric$TX_Immunomodulators,CARD9_Tetrachloric$TX_Antibiotic))

#-----------------------------------------------------------------------------------------------------------------------#
#Point Biserial Correlation: Continuous + Binary Variable

biserial.cor(CARD9_Final$`anti-C.A_IgG`, CARD9_Final$Diagnosis)
biserial.cor(CARD9_Final$`anti-C.A_IgG`, CARD9_Final$Gender)
biserial.cor(CARD9_Final$`anti-C.A_IgG`, CARD9_Final$Disease_indices)
biserial.cor(CARD9_Final$`anti-C.A_IgG`, CARD9_Final$TX_Antidiarrheal)
biserial.cor(CARD9_Final$`anti-C.A_IgG`, CARD9_Final$TX_ASA)
biserial.cor(CARD9_Final$`anti-C.A_IgG`, CARD9_Final$TX_Biologics)
biserial.cor(CARD9_Final$`anti-C.A_IgG`, CARD9_Final$TX_Steroid)
biserial.cor(CARD9_Final$`anti-C.A_IgG`, CARD9_Final$TX_Immunomodulators)
biserial.cor(CARD9_Final$`anti-C.A_IgG`, CARD9_Final$TX_Antibiotic)

biserial.cor(CARD9_Final$Height, CARD9_Final$Diagnosis)
biserial.cor(CARD9_Final$Height, CARD9_Final$Gender)
biserial.cor(CARD9_Final$Height, CARD9_Final$Disease_indices)
biserial.cor(CARD9_Final$Height, CARD9_Final$TX_Antidiarrheal)
biserial.cor(CARD9_Final$Height, CARD9_Final$TX_ASA)
biserial.cor(CARD9_Final$Height, CARD9_Final$TX_Biologics)
biserial.cor(CARD9_Final$Height, CARD9_Final$TX_Steroid)
biserial.cor(CARD9_Final$Height, CARD9_Final$TX_Immunomodulators)
biserial.cor(CARD9_Final$Height, CARD9_Final$TX_Antibiotic)

biserial.cor(CARD9_Final$Weight, CARD9_Final$Diagnosis)
biserial.cor(CARD9_Final$Weight, CARD9_Final$Gender)
biserial.cor(CARD9_Final$Weight, CARD9_Final$Disease_indices)
biserial.cor(CARD9_Final$Weight, CARD9_Final$TX_Antidiarrheal)
biserial.cor(CARD9_Final$Weight, CARD9_Final$TX_ASA)
biserial.cor(CARD9_Final$Weight, CARD9_Final$TX_Biologics)
biserial.cor(CARD9_Final$Weight, CARD9_Final$TX_Steroid)
biserial.cor(CARD9_Final$Weight, CARD9_Final$TX_Immunomodulators)
biserial.cor(CARD9_Final$Weight, CARD9_Final$TX_Antibiotic)

biserial.cor(CARD9_Final$BMI, CARD9_Final$Diagnosis)
biserial.cor(CARD9_Final$BMI, CARD9_Final$Gender)
biserial.cor(CARD9_Final$BMI, CARD9_Final$Disease_indices)
biserial.cor(CARD9_Final$BMI, CARD9_Final$TX_Antidiarrheal)
biserial.cor(CARD9_Final$BMI, CARD9_Final$TX_ASA)
biserial.cor(CARD9_Final$BMI, CARD9_Final$TX_Biologics)
biserial.cor(CARD9_Final$BMI, CARD9_Final$TX_Steroid)
biserial.cor(CARD9_Final$BMI, CARD9_Final$TX_Immunomodulators)
biserial.cor(CARD9_Final$BMI, CARD9_Final$TX_Antibiotic)

biserial.cor(CARD9_Final$Sample_Collection_Age, CARD9_Final$Diagnosis)
biserial.cor(CARD9_Final$Sample_Collection_Age, CARD9_Final$Gender)
biserial.cor(CARD9_Final$Sample_Collection_Age, CARD9_Final$Disease_indices)
biserial.cor(CARD9_Final$Sample_Collection_Age, CARD9_Final$TX_Antidiarrheal)
biserial.cor(CARD9_Final$Sample_Collection_Age, CARD9_Final$TX_ASA)
biserial.cor(CARD9_Final$Sample_Collection_Age, CARD9_Final$TX_Biologics)
biserial.cor(CARD9_Final$Sample_Collection_Age, CARD9_Final$TX_Steroid)
biserial.cor(CARD9_Final$Sample_Collection_Age, CARD9_Final$TX_Immunomodulators)
biserial.cor(CARD9_Final$Sample_Collection_Age, CARD9_Final$TX_Antibiotic)

biserial.cor(CARD9_Final$Diagnosis_Age, CARD9_Final$Diagnosis)
biserial.cor(CARD9_Final$Diagnosis_Age, CARD9_Final$Gender)
biserial.cor(CARD9_Final$Diagnosis_Age, CARD9_Final$Disease_indices)
biserial.cor(CARD9_Final$Diagnosis_Age, CARD9_Final$TX_Antidiarrheal)
biserial.cor(CARD9_Final$Diagnosis_Age, CARD9_Final$TX_ASA)
biserial.cor(CARD9_Final$Diagnosis_Age, CARD9_Final$TX_Biologics)
biserial.cor(CARD9_Final$Diagnosis_Age, CARD9_Final$TX_Steroid)
biserial.cor(CARD9_Final$Diagnosis_Age, CARD9_Final$TX_Immunomodulators)
biserial.cor(CARD9_Final$Diagnosis_Age, CARD9_Final$TX_Antibiotic)

biserial.cor(CARD9_Final$Collection_DiseaseActivity, CARD9_Final$Diagnosis)
biserial.cor(CARD9_Final$Collection_DiseaseActivity, CARD9_Final$Gender)
biserial.cor(CARD9_Final$Collection_DiseaseActivity, CARD9_Final$Disease_indices)
biserial.cor(CARD9_Final$Collection_DiseaseActivity, CARD9_Final$TX_Antidiarrheal)
biserial.cor(CARD9_Final$Collection_DiseaseActivity, CARD9_Final$TX_ASA)
biserial.cor(CARD9_Final$Collection_DiseaseActivity, CARD9_Final$TX_Biologics)
biserial.cor(CARD9_Final$Collection_DiseaseActivity, CARD9_Final$TX_Steroid)
biserial.cor(CARD9_Final$Collection_DiseaseActivity, CARD9_Final$TX_Immunomodulators)
biserial.cor(CARD9_Final$Collection_DiseaseActivity, CARD9_Final$TX_Antibiotic)

biserial.cor(CARD9_Final$`#days_stool/serum`, CARD9_Final$Diagnosis)
biserial.cor(CARD9_Final$`#days_stool/serum`, CARD9_Final$Gender)
biserial.cor(CARD9_Final$`#days_stool/serum`, CARD9_Final$Disease_indices)
biserial.cor(CARD9_Final$`#days_stool/serum`, CARD9_Final$TX_Antidiarrheal)
biserial.cor(CARD9_Final$`#days_stool/serum`, CARD9_Final$TX_ASA)
biserial.cor(CARD9_Final$`#days_stool/serum`, CARD9_Final$TX_Biologics)
biserial.cor(CARD9_Final$`#days_stool/serum`, CARD9_Final$TX_Steroid)
biserial.cor(CARD9_Final$`#days_stool/serum`, CARD9_Final$TX_Immunomodulators)
biserial.cor(CARD9_Final$`#days_stool/serum`, CARD9_Final$TX_Antibiotic)

#-----------------------------------------------------------------------------------------------------------------------#
#Polyserial Correlation

polyserial(CARD9_Final$`anti-C.A_IgG`, CARD9_Final$Smoking)
polyserial(CARD9_Final$Height, CARD9_Final$Smoking)
polyserial(CARD9_Final$Weight, CARD9_Final$Smoking)
polyserial(CARD9_Final$BMI, CARD9_Final$Smoking)
polyserial(CARD9_Final$Sample_Collection_Age, CARD9_Final$Smoking)
polyserial(CARD9_Final$Diagnosis_Age, CARD9_Final$Smoking)
polyserial(CARD9_Final$Collection_DiseaseActivity, CARD9_Final$Smoking)
polyserial(CARD9_Final$`#days_stool/serum`, CARD9_Final$Smoking)

polyserial(CARD9_Final$`anti-C.A_IgG`, CARD9_Final$CARD9_S12N)
polyserial(CARD9_Final$Height, CARD9_Final$CARD9_S12N)
polyserial(CARD9_Final$Weight, CARD9_Final$CARD9_S12N)
polyserial(CARD9_Final$BMI, CARD9_Final$CARD9_S12N)
polyserial(CARD9_Final$Sample_Collection_Age, CARD9_Final$CARD9_S12N)
polyserial(CARD9_Final$Diagnosis_Age, CARD9_Final$CARD9_S12N)
polyserial(CARD9_Final$Collection_DiseaseActivity, CARD9_Final$CARD9_S12N)
polyserial(CARD9_Final$`#days_stool/serum`, CARD9_Final$CARD9_S12N)

#-----------------------------------------------------------------------------------------------------------------------#
#Cramer's V 

cramerV(CARD9_Final$Smoking, CARD9_Final$Diagnosis)
cramerV(CARD9_Final$Smoking, CARD9_Final$Gender)
cramerV(CARD9_Final$Smoking, CARD9_Final$Disease_indices)
cramerV(CARD9_Final$Smoking, CARD9_Final$TX_Antidiarrheal)
cramerV(CARD9_Final$Smoking, CARD9_Final$TX_ASA)
cramerV(CARD9_Final$Smoking, CARD9_Final$TX_Biologics)
cramerV(CARD9_Final$Smoking, CARD9_Final$TX_Steroid)
cramerV(CARD9_Final$Smoking, CARD9_Final$TX_Immunomodulators)
cramerV(CARD9_Final$Smoking, CARD9_Final$TX_Antibiotic)

cramerV(CARD9_Final$CARD9_S12N, CARD9_Final$Diagnosis)
cramerV(CARD9_Final$CARD9_S12N, CARD9_Final$Gender)
cramerV(CARD9_Final$CARD9_S12N, CARD9_Final$Disease_indices)
cramerV(CARD9_Final$CARD9_S12N, CARD9_Final$TX_Antidiarrheal)
cramerV(CARD9_Final$CARD9_S12N, CARD9_Final$TX_ASA)
cramerV(CARD9_Final$CARD9_S12N, CARD9_Final$TX_Biologics)
cramerV(CARD9_Final$CARD9_S12N, CARD9_Final$TX_Steroid)
cramerV(CARD9_Final$CARD9_S12N, CARD9_Final$TX_Immunomodulators)
cramerV(CARD9_Final$CARD9_S12N, CARD9_Final$TX_Antibiotic)
cramerV(CARD9_Final$CARD9_S12N, CARD9_Final$Smoking)

#-----------------------------------------------------------------------------------------------------------------------#
#First Linear Model after Variable Reduction, Complete-Case Analysis
lm1 <- lm(data=CARD9_CAgenotype,`anti-C.A_IgG` ~ CARD9_S12N + `Diagnosis` + Gender +
            `Sample_Collection_Age` + `TX_Antidiarrheal` + `TX_Steroid` )
summary(lm1)

attach(CARD9_CAgenotype)
#Standard Residuals vs ISI
StanResMLS <- rstandard(lm1)
dataMLS <- data.frame(`anti-C.A_IgG`,StanResMLS)

# Plotting Standardized Residuals vs ISI
ggplot() + geom_point(data=dataMLS, aes(x=`anti-C.A_IgG`, y=StanResMLS, color = "steelblue"), size = 0.2) +
  geom_hline(yintercept=2,color='blue') + geom_hline(yintercept=-2, color='blue') +
  scale_color_manual(name = element_blank(), labels = c("MLS"), values = c("steelblue")) +
  labs(y = "Standarized Residual") + ggtitle("Standarized Residuals Plot") +
  theme(plot.title = element_text(hjust = 0.5))

# Standarized Residuals vs Fitted
Fitted = fitted(lm1)
dataMLSFitted <- data.frame(Fitted,StanResMLS)
ggplot() + geom_point(data=dataMLSFitted, aes(x=Fitted, y=StanResMLS, color = "steelblue"), size = 0.2) + 
  geom_hline(yintercept=2,color='blue') + geom_hline(yintercept=-2, color='blue') +
  scale_color_manual(name = element_blank(), labels = c("MLS"), values = c("steelblue")) +
  labs(y = "Standarized Residual") + labs(x = "Fitted value") + 
  ggtitle("Standarized Residuals Plot (Fitted) ") +
  theme(plot.title = element_text(hjust = 0.5))

# Q-Q Plot for Normality
par(mfrow=c(1,1))
qqnorm(StanResMLS)
qqline(StanResMLS)

# Outliers
hat=hatvalues(lm1);
res_stud=rstudent(lm1);
dfs=dffits(lm1)
par(mfrow=c(1,3))
plot(hat, main= "hat values")
plot(res_stud, main = "stud. res.")
plot(dfs, main= "dffits")
par(mfrow=c(1,1))

#-----------------------------------------------------------------------------------------------------------------------#
#Step-wise variable selection  
step.model <- stepAIC(lm1, direction = "both", 
                      trace = FALSE)
summary(step.model)
step.model2 <- stepAIC(lm1, direction = "backward", 
                       trace = FALSE)
summary(step.model2)

#-----------------------------------------------------------------------------------------------------------------------#
#Final Linear Model
lm2interaction <- lm(data=CARD9_CAgenotype,`anti-C.A_IgG` ~ `Diagnosis` + `TX_Steroid` + `Diagnosis`*`TX_Steroid`)
lm2 <- lm(data=CARD9_CAgenotype,`anti-C.A_IgG` ~ `Diagnosis` + `TX_Steroid`)
summary(lm2)

#Standard Residuals vs ISI
StanResMLS2 <- rstandard(lm2)
dataMLS2 <- data.frame(`anti-C.A_IgG`,StanResMLS2)

# Plotting Standardized Residuals vs ISI
ggplot() + geom_point(data=dataMLS2, aes(x=`anti-C.A_IgG`, y=StanResMLS2, color = "steelblue"), size = 0.2) +
  geom_hline(yintercept=2,color='blue') + geom_hline(yintercept=-2, color='blue') +
  scale_color_manual(name = element_blank(), labels = c("MLS"), values = c("steelblue")) +
  labs(y = "Standarized Residual") + ggtitle("Standarized Residuals Plot") +
  theme(plot.title = element_text(hjust = 0.5))

# Standarized Residuals vs Fitted
Fitted2 = fitted(lm2)
dataMLSFitted2 <- data.frame(Fitted2,StanResMLS2)
ggplot() + geom_point(data=dataMLSFitted2, aes(x=Fitted2, y=StanResMLS2, color = "steelblue"), size = 0.15) + 
  geom_hline(yintercept=2,color='blue') + geom_hline(yintercept=-2, color='blue') +
  scale_color_manual(name = element_blank(), labels = c("MLS"), values = c("steelblue")) +
  labs(y = "Standarized Residual") + labs(x = "Fitted value") + 
  ggtitle("Standarized Residuals Plot (Fitted) ") +
  theme(plot.title = element_text(hjust = 0.5))

# Q-Q PLot for Normality
par(mfrow=c(1,1))
qqnorm(StanResMLS2)
qqline(StanResMLS2)

# Outliers
hat=hatvalues(lm2);
res_stud=rstudent(lm2);
dfs=dffits(lm2)
par(mfrow=c(1,3))
plot(hat, main= "hat values")
plot(res_stud, main = "stud. res.")
plot(dfs, main= "dffits")
par(mfrow=c(1,1))




