setwd("C:/Temp")

install.packages("lubridate")
library(lubridate)
library(reshape2)
library(ggplot2)

cData = read.csv("CIDP_clinical data_161101.csv", na.strings = c("na"))

cData$gender = factor(cData$gender, labels = c("male", "female"))
cData$date_birth = as.Date(cData$date_birth)
cData$date_onset = as.Date(cData$date_onset)
cData$date_firstVisit = as.Date(cData$date_firstVisit)
cData$date_lastVisit = as.Date(cData$date_lastVisit)
cData$date_EDX = as.Date(cData$date_EDX)
levels(cData$cluster) = list(Typical = "typical", DADS = "DADS", 
                MADSAM = "MADSAM", Pure_sensory = "pure sensory")
cData$date_csf = as.Date(cData$date_csf)
cData$date_disability = as.Date(cData$date_disability)
cData$tx = factor(cData$tx, labels = c("PSL", 
                  "IVIG", "PSL+IVIG", "(PSL or IVIG)+IS", "PSL+IVIG+IS"))

df = tbl_df(cData)
df$age = year(df$date_firstVisit) - year(df$date_birth)
# df$timeOnsetDx = round(difftime(df$date_firstVisit, df$date_onset, unit="")/4)
df$mRS_change = df$mRS_initial - df$mRS_last
df$outcome_improve= ifelse(df$mRS_change > 0, "improved", 
                           "unchanged or worsened")
df$outcome_improve = ifelse(df$mRS_change == 0, "unchanged", df$outcome_improve)
df$outcome_improve = factor(df$outcome_improve, 
                            labels = c("improved", "unchanged", "worsened"))
df$outcome_dichotomy =  ifelse(df$mRS_last > 3, "poor", "favorable")

#################################################################
### statistics: accoridng to phenotype

# age
df %>%
  group_by(cluster) %>%
  summarize(age.mean = mean(age), age.sd = sd(age))

anova(lm(df$age~df$cluster))
boxplot(df$age~df$cluster)
TukeyHSD(aov(df$age~df$cluster)) # DADS vs. pure sensory, p=0.03

# gender 
table(df$cluster, df$gender) # male predominance in CIDP? 특이하네요...
chisq.test(table(df$cluster, df$gender)) # not significant
fisher.test(table(df$cluster, df$gender)) # not significant

# diabetes
table(df$cluster, df$diabetes) 
chisq.test(table(df$cluster, df$diabetes)) # not significant
fisher.test(table(df$cluster, df$diabetes)) # not significant

# ataxia 
table(df$cluster, df$ataxia) # MADSAM 
chisq.test(table(df$cluster, df$ataxia)) # not significant
fisher.test(table(df$cluster, df$ataxia)) # not significant

# cranial nerve involvement 
table(df$cluster, df$cNerve)  
chisq.test(table(df$cluster, df$cNerve)) # not significant
fisher.test(table(df$cluster, df$cNerve)) # not significant

# atrophy 
table(df$cluster, df$atrophy)  
chisq.test(table(df$cluster, df$atrophy)) # not significant
fisher.test(table(df$cluster, df$atrophy)) # not significant

# monoclona gammopathy 
table(df$cluster, df$mgp)  # no monoclonal gammopathy in MADSAM 
# 3/27 in typical (1 uc/uk), 4/14 in DADS (5 uc/uk), 0/10 in madsam (3 uc/uk), 
# 2/7 in pure sensory
table(df$cluster, df$lightch)

# MRC sum score 
df %>%
  group_by(cluster) %>%
  summarize(mrcSum.mean = mean(mrcSum), mrcSum.sd = sd(mrcSum))

anova(lm(df$mrcSum~df$cluster)) # p=0.00017
boxplot(df$mrcSum~df$cluster)
TukeyHSD(aov(df$mrcSum~df$cluster)) # DADS vs Typical, Pure sensory vs. Typical

# modified Rankin disabilty scale 
table(df$cluster, df$mRS_initial)
fisher.test(df$cluster, df$mRS_initial)

table(df$cluster, df$mRS_last)
chisq.test(df$cluster, df$mRS_last)

table(df$cluster, df$mRS_change)
chisq.test(df$cluster, df$mRS_change)

table(df$cluster, df$outcome_dichotomy)
chisq.test(df$cluster, df$outcome_dichotomy)

table(df$cluster, df$outcome_improve)
chisq.test(df$cluster, df$outcome_improve)

# tx 
table(df$cluster, df$tx)
fisher.test(df$cluster, df$tx)

#################################################################
##### statistics according to edx clusters

cluster1 = read.table("cluster1")
cluster1 = as.character(cluster1$V1)

cluster2 = read.table("cluster2")
cluster2 = as.character(cluster2$V1)

cluster3 = read.table("cluster3")
cluster3 = as.character(cluster3$V1)

edx = read.csv("summary_mean_dataset.csv")

temp = strsplit(as.character(edx$hosp_init_id), split="_", fixed=T)
for (i in 1:length(edx[,1])){
  hosp = temp[[i]][1]
  id = temp[[i]][3]
  edx$hosp_id[i] = paste(hosp, id, sep="_")
}

df.m$hosp_id = paste(df.m$hospital, df.m$ID, sep="_")

df.m.m = merge(df.m, edx, by="hosp_id")
df.m.m$id = as.character(df.m.m$id)
df.m.m$cluster = ifelse(df.m.m$id %in% cluster1, "C1", NA)
df.m.m$cluster = ifelse(df.m.m$id %in% cluster2, "C2", df.m.m$cluster)
df.m.m$cluster = ifelse(df.m.m$id %in% cluster3, "C3", df.m.m$cluster)
table(df.m.m$cluster)
df.m$cluster = factor(df.m$cluster)

# age
df.m %>%
  group_by(cluster) %>%
  summarize(age.mean = mean(age), age.sd = sd(age))

anova(lm(df.m$age~df.m$cluster))
boxplot(df.m$age~df.m$cluster)

# gender 
table(df.m$cluster, df.m$gender) # male predominance in CIDP? 특이하네요...
chisq.test(table(df.m$cluster, df.m$gender)) # not significant
fisher.test(table(df.m$cluster, df.m$gender)) # not significant

# subtype
table(df.m$cluster, df.m$subtype) # 
chisq.test(table(df.m$cluster, df.m$subtype)) # p = 0.024
fisher.test(table(df.m$cluster, df.m$subtype)) # p = 0.039

temp = melt(table(df.m$cluster, df.m$subtype))
colnames(temp) = c("cluster", "subtype", "number_of_patients")
p = ggplot(temp, aes(x=subtype, y=number_of_patients, fill=subtype)) +
  geom_bar(stat="identity") + 
  facet_wrap(~cluster) +
  coord_flip()
p

# diabetes
table(df.m$cluster, df.m$diabetes) 
chisq.test(table(df.m$cluster, df.m$diabetes)) # p=0.0052
fisher.test(table(df.m$cluster, df.m$diabetes)) # p=0.0076, more prevalent diabetes in C2 clsuter

# ataxia 
table(df.m$cluster, df.m$ataxia) # MADSAM 
chisq.test(table(df.m$cluster, df.m$ataxia)) # not significant
fisher.test(table(df.m$cluster, df.m$ataxia)) # not significant

# cranial nerve involvement 
table(df.m$cluster, df.m$cNerve)  
chisq.test(table(df.m$cluster, df.m$cNerve)) # not significant
fisher.test(table(df.m$cluster, df.m$cNerve)) # not significant

# atrophy 
table(df.m$cluster, df.m$atrophy)  
chisq.test(table(df.m$cluster, df.m$atrophy)) # not significant
fisher.test(table(df.m$cluster, df.m$atrophy)) # not significant

# monoclona gammopathy 
table(df.m$cluster, df.m$mgp)  # no monoclonal gammopathy in MADSAM 
# no monoclonal gammopaty in C2 clsuter, more prevalent in C1 
table(df.m$cluster, df.m$lightch)

# MRC sum score 
df.m %>%
  group_by(cluster) %>%
  summarize(mrcSum.mean = mean(mrcSum), mrcSum.sd = sd(mrcSum))

anova(lm(df.m$mrcSum~df.m$cluster)) # non significant 
boxplot(df.m$mrcSum~df.m$cluster)

# modified Rankin disabilty scale 
table(df.m$cluster, df.m$mRS_initial)
fisher.test(df.m$cluster, df.m$mRS_initial)

table(df.m$cluster, df.m$mRS_last)
chisq.test(df.m$cluster, df.m$mRS_last)

table(df.m$cluster, df.m$mRS_change)
chisq.test(df.m$cluster, df.m$mRS_change)

table(df.m$cluster, df.m$outcome_dichotomy)
chisq.test(df.m$cluster, df.m$outcome_dichotomy) # non significant, but C1 & C3 more favorable than C2

temp = melt(table(df.m$cluster, df.m$outcome_dichotomy))
colnames(temp) = c("cluster", "disability_latest_visit", "number_of_patients")
p = ggplot(temp, aes(x=disability_latest_visit, y=number_of_patients, fill=disability_latest_visit)) +
  geom_bar(stat="identity") + 
  facet_wrap(~cluster) + 
  theme(axis.text.x=element_blank())
p



table(df.m$cluster, df.m$outcome_improve)
chisq.test(df.m$cluster, df.m$outcome_improve) # p = 0.019

temp = melt(table(df.m$cluster, df.m$outcome_improve))
colnames(temp) = c("cluster", "treatment_outcome", "number_of_patients")
p = ggplot(temp, aes(x=treatment_outcome, y=number_of_patients, fill=treatment_outcome)) +
  geom_bar(stat="identity") + 
  facet_wrap(~cluster) + 
  theme(axis.text.x=element_blank())
p


# tx 
table(df.m$cluster, df.m$tx)
fisher.test(df.m$cluster, df.m$tx)

temp1 = df.m[,c(47,40:46)]
temp2 = melt(temp1)

png("boxplot_cluster.png",
    width=15*300,
    height=5*300,
    res=500,
    pointsize = 12)

p = ggplot(temp2, aes(cluster, value)) +
  geom_boxplot(aes(fill=cluster)) +
  facet_grid(~variable)
p

dev.off()

kruskal.test(value~cluster, data=temp2[temp2$variable=="DML",]) # p = 8.229e-11
kruskal.test(value~cluster, data=temp2[temp2$variable=="DUR",]) # p-value = 0.000106
kruskal.test(value~cluster, data=temp2[temp2$variable=="NCV",]) # p-value = 0.0001185
kruskal.test(value~cluster, data=temp2[temp2$variable=="CB",]) # p-value = 0.8754
kruskal.test(value~cluster, data=temp2[temp2$variable=="TD",]) # p-value = 0.1773
kruskal.test(value~cluster, data=temp2[temp2$variable=="FL",]) # p-value = 7.901e-05
kruskal.test(value~cluster, data=temp2[temp2$variable=="CMAP",]) # p-value = 4.572e-10

C3 = df.m[df.m$cluster == "C3",]
plot(C3$DML, C3$CMAP)
plot(df.m$DML, df.m$CMAP)
plot(df.m$DUR, df.m$CMAP)
plot(df.m$DML, df.m$DUR)
# favorable outcome, enrichment analysis... 
