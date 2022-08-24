###########################################
# Meta-analysis of outcomes following CBT #
###########################################

module add general/R;R

# Script 1. Phenotype processing

rm /users/k1620287/brc_scratch/meta.analysis/pheno/scripts/make.phenotype.datasets.sh

####cat <<'EOT'>> /users/k1620287/brc_scratch/meta.analysis/pheno/scripts/make.phenotype.datasets.sh

#!/bin/sh
#$-S /bin/sh
#$ -V
#$ -pe smp 8
#$ -l h_vmem=9G
#$ -l h_rt=2:00:00
#$ -m be
#$ -o /users/k1620287/brc_scratch/meta.analysis/output
#$ -e /users/k1620287/brc_scratch/meta.analysis/output

awk '{print $1}' /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/ki_mdd.postimputation.maf05.hwe.keeps.udsex.pheno.locf.fam > /users/k1620287/brc_scratch/meta.analysis/pheno/pheno.data/ki_mdd.gid

###echo '

#### PACKAGES ####
library(dplyr)
library(lme4)

#### ACRONYMS ####

# gxte = bochum and braunschweig sample (exposure-cbt genes for treatment study)
# gxtc = genes for treatment child study
# ki = karolinska panic disorder sample (icbt study)
# wb = wurzburg sample (panic.net consortium)
# bl = baseline
# pt = pot-treatment
# fu = follow-up
# z = z-score
# res. = residualised
# st. = standardised
# sd = standard deviation
# csr = clinical severity rating
# cgi = clinical global improvement
# madrs = montgomery-asberg depression rating scale
# pdssr = panic disorder severity scale self-report
# pas = panic agaoraphobia scale

#### DATA ####

### child sample ###
# primary outcome measure = csr

gxtc <- read.table("/users/k1620287/brc_scratch/meta.analysis/pheno/original.data/GxT_phenotypes_long_format_labels.csv", ",", head=T)

names(gxtc)[1] <- "IID"

gxtcfu <- gxtc[which(gxtc$time==5),]
gxtcfu <- gxtcfu[,c("IID","prisevlocf")]

gxtc <- gxtc[which(gxtc$time==1),]
gxtc <- gxtc[,c("IID","sex","age","diag","siter","treatment","trial","aprisev","prisev")]
colnames(gxtc) <- c("IID","sex","age","primary.diagnosis","site","treatment","trial","csr.bl","csr.pt")

gxtc <- merge(gxtc, gxtcfu, by = "IID", type="left")
colnames(gxtc) <- c("IID","sex","age","primary.diagnosis","site","treatment","trial","csr.bl","csr.pt","csr.fu.locf")

gxtc[gxtc==999] <- NA
gxtc[gxtc==998] <- NA
gxtc[gxtc==" "] <- NA

gxtc <- gxtc %>%
  mutate(
    IID = as.character(IID),
    age = as.integer(age),
    sex = as.factor(sex),
    site = as.factor(site),
    treatment = as.factor(treatment),
    primary.diagnosis = as.factor(primary.diagnosis),
    trial = as.factor(trial)
  )

gxtc1 <- read.table("/users/k1620287/brc_scratch/meta.analysis/pheno/original.data/old.gxtc/gxtc_cr_250917.csv", head=T, sep=",")
gxtc1 <- gxtc1[,c("IID","IID","Sex","Age","Primary_diagnosis","Site","Treatment_type","N_Sessions","Med_pre","N_Comorbidities","CSR_pre","CSR_post")]
colnames(gxtc1) <- c("IID","FID","sex","age","primary.dsm","site","treatment","n.sessions","psych.med.bl","n.comorbidities","csr.bl","csr.pt")

gxtc1$psych.med.bl <- ifelse(gxtc1$psych.med.bl==0, 0, ifelse(gxtc1$psych.med.bl>=1, 1, NA))

gxtc1 <- gxtc1[,c("IID","n.sessions","psych.med.bl","n.comorbidities")]

gxtc <- gxtc %>%
  mutate(
    IID = as.character(IID),
    psych.med.bl = as.factor(psych.med.bl),
    n.comorbidities = as.integer(n.comorbidities),
    n.sessions = as.integer(n.sessions)
  )

gxtc <- merge(gxtc, gxtc1, by = "IID", all.x=TRUE)

# covariates used in Jonis paper: baseline severity, age, gender, treatment type, diagnosis and trial
# lots of missing data for number of sessions (half), 20 missing for n.comorbidities, psych.med
# n.sessions does not have an significant effect in the model when assessed in individuals who have data, or when imputing using mean imputation (mean taken from each treatment group)
# therefore n.sessions excluded from the model - because of excessive missing data!

# identify individuals in data with QC'd genotype data

gid <- read.table("/users/k1620287/brc_scratch/meta.analysis/analysis/gxtc/binary/gxtc.prs.id",head=F, sep=" ")

names(gid)[1] <- "IID"

gid$V2 <- NULL

gid$qc.genotpyes <- 1

gxtc <- merge(gxtc, gid, by = "IID", all.x=TRUE)

gxtc$qc.genotpyes <- ifelse(is.na(gxtc$qc.genotpyes), 0, gxtc$qc.genotpyes)

gxtc_g <- gxtc[which(gxtc$qc.genotpyes==1 & !is.na(gxtc$csr.pt)),]

# gxtc_g <- gxtc[which(gxtc$qc.genotpyes==1 & !is.na(gxtc$csr.fu.locf)),]

gxtc_g$csr.bl.z <- (gxtc_g$csr.bl/(sd(gxtc_g$csr.bl, na.rm=TRUE)))
gxtc_g$csr.pt.z <- (gxtc_g$csr.pt/(sd(gxtc_g$csr.bl, na.rm=TRUE)))
gxtc_g$csr.fu.z <- (gxtc_g$csr.fu.locf/(sd(gxtc_g$csr.bl, na.rm=TRUE)))

gxtc_g$primary.cbt.outcomes.res <- resid(lm(csr.pt.z ~ csr.bl.z + age + sex + site + treatment + primary.diagnosis + trial, data = gxtc_g, na.action = na.exclude))

gxtc_g$FID <- gxtc_g$IID

gxtc.outcomes <- gxtc_g[,c("FID","IID", "primary.cbt.outcomes.res")]

write.table(gxtc.outcomes, "/users/k1620287/brc_scratch/meta.analysis/pheno/pheno.data/gxtc.primary.z.outcomes.res.txt", col.names=T, row.names=F, quote=F, sep=" ")

g.bl <- lmer(csr.bl.z ~ (1 | trial) + (1 | site) + (1 | treatment) + age + sex + primary.diagnosis + n.comorbidities + psych.med.bl, data = gxtc_g, na.action = na.exclude)
g.pt <- lmer(csr.pt.z ~ (1 | trial) + (1 | site) + (1 | treatment) + csr.bl.z + age + sex + primary.diagnosis + n.comorbidities + psych.med.bl, data = gxtc_g, na.action = na.exclude)

summary(g.bl)
summary(g.pt)


g.bl <- lmer(csr.bl.z ~ (1 | trial) + (1 | site) + (1 | treatment) + (1 | primary.diagnosis) + age + sex + n.comorbidities + psych.med.bl, data = gxtc_g, na.action = na.exclude)
g.pt <- lmer(csr.pt.z ~ (1 | trial) + (1 | site) + (1 | treatment) + (1 | primary.diagnosis) + csr.bl.z + age + sex + n.comorbidities + psych.med.bl, data = gxtc_g, na.action = na.exclude)

summary(g.bl)
summary(g.pt)


write.table(gxtc_g, "/users/k1620287/brc_scratch/meta.analysis/pheno/pheno.data/gxtc.genotyped.phenotype.data.txt", col.names=T, row.names=F, quote=F, sep=" ")

# FULL DATA (INCLUDING THOSE WITHOUT GENETIC DATA)

gxtc$csr.bl.z <- (gxtc$csr.bl/(sd(gxtc$csr.bl, na.rm=TRUE)))
gxtc$csr.pt.z <- (gxtc$csr.pt/(sd(gxtc$csr.bl, na.rm=TRUE)))
gxtc$csr.fu.z <- (gxtc$csr.fu.locf/(sd(gxtc$csr.bl, na.rm=TRUE)))

full.bl <- lmer(csr.bl.z ~ (1 | trial) + (1 | site) + (1 | treatment) + age + sex + primary.diagnosis + n.comorbidities + psych.med.bl, data = gxtc, na.action = na.exclude)
full.pt <- lmer(csr.pt.z ~ (1 | trial) + (1 | site) + (1 | treatment) + csr.bl + age + sex + primary.diagnosis + n.comorbidities + psych.med.bl, data = gxtc, na.action = na.exclude)
summary(full.bl)
summary(full.pt)

write.table(gxtc_g, "/users/k1620287/brc_scratch/meta.analysis/pheno/pheno.data/gxtc.phenotype.data.txt", col.names=T, row.names=F, quote=F, sep=" ")

###########################################
### adult samples ###

# karolinska mdd iCBT sample
# primary outcome measure = madrs

kimd <- read.table("/users/k1620287/brc_scratch/meta.analysis/pheno/original.data/SWE_iCBT_phenotypes.csv", head=T, sep=",")

kimd <- kimd[,c(1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)]

colnames(kimd) <- c("IID","age","site","primary.dsm","n.sessions","treatment","sex","n.comorbidities","psych.med.bl","primary.bl","primary.pt",
"primary.w1","primary.w2","primary.w3","primary.w4","primary.w5","primary.w6","primary.w7","primary.w8","primary.w9","primary.w10")

kimd$primary.locf <- ifelse(is.na(kimd$primary.pt) & !is.na(kimd$primary.w10), kimd$primary.w10,
 ifelse(is.na(kimd$primary.pt) & !is.na(kimd$primary.w9), kimd$primary.w9,
 ifelse(is.na(kimd$primary.pt) & !is.na(kimd$primary.w8), kimd$primary.w8,
 ifelse(is.na(kimd$primary.pt) & !is.na(kimd$primary.w7), kimd$primary.w7,
 ifelse(is.na(kimd$primary.pt) & !is.na(kimd$primary.w6), kimd$primary.w6, kimd$primary.pt)))))

kimd$sex <- ifelse(kimd$sex==1, "female", "male")

kimd$FID <- kimd$IID

kimd <- kimd[,c("IID","FID","age","site","primary.dsm","n.sessions","treatment","sex","n.comorbidities","psych.med.bl","primary.bl","primary.locf","primary.pt")]

kimd <- kimd %>%
	mutate(
		FID = as.character(FID), IID = as.character(IID), age = as.integer(age), sex = as.factor(sex),
		site = as.factor(site), treatment = as.factor(treatment), primary.dsm = as.factor(primary.dsm),
		psych.med.bl = as.factor(psych.med.bl), n.comorbidities = as.integer(n.comorbidities), n.sessions = as.integer(n.sessions))

gid <- read.table("/users/k1620287/brc_scratch/meta.analysis/pheno/pheno.data/ki_mdd.gid", head=F)

names(gid)[1] <- "IID"

gid$qc.genotpyes <- 1

kimd <- merge(kimd, gid, by = "IID", all.x=TRUE)

kimd$qc.genotpyes <- ifelse(is.na(kimd$qc.genotpyes), 0, kimd$qc.genotpyes)

kimd_g <- kimd[which(kimd$qc.genotpyes==1 & !is.na(kimd$primary.locf) & !is.na(kimd$primary.bl) & !is.na(kimd$age)),]

kimd_g$primary.bl.z <- (kimd_g$primary.bl/(sd(kimd_g$primary.bl, na.rm=TRUE)))
kimd_g$primary.pt.z <- (kimd_g$primary.locf/(sd(kimd_g$primary.bl, na.rm=TRUE)))

g.bl <- lm(primary.bl.z ~ age + sex  + n.comorbidities + n.sessions + psych.med.bl, data = kimd_g, na.action = na.exclude)
g.pt <- lm(primary.pt.z ~ primary.bl.z + age + sex + n.comorbidities + n.sessions + psych.med.bl, data = kimd_g, na.action = na.exclude)

summary(g.bl)
summary(g.pt)

write.table(kimd_g, "/users/k1620287/brc_scratch/meta.analysis/pheno/pheno.data/kimd.genotype.phenotype.data.txt", col.names=T, row.names=F, quote=F, sep=" ")

kimd_g$primary.z.outcomes.res <- resid(lm(primary.pt.z ~ primary.bl.z + age + sex + n.comorbidities + n.sessions + psych.med.bl, data = kimd_g, na.action = na.exclude))

kimd.outcomes <- kimd[,c("FID","IID","primary.z.outcomes.res")]

write.table(kimd.outcomes, "/users/k1620287/brc_scratch/meta.analysis/pheno/pheno.data/kimd.primary.z.pt.res.txt", col.names=T, row.names=F, quote=F, sep=" ")

# FULL DATA (INCLUDING THOSE WITHOUT GENETIC DATA)

kimd$primary.bl.z <- (kimd$primary.bl/(sd(kimd$primary.bl, na.rm=TRUE)))
kimd$primary.pt.z <- (kimd$primary.pt/(sd(kimd$primary.bl, na.rm=TRUE)))

full.bl <- lmer(primary.bl.z ~ age + sex  + n.comorbidities + n.sessions + psych.meds.bl,, data = kimd, na.action = na.exclude)
full.pt <- lmer(primary.pt.z ~ primary.bl.z + age + sex  + n.comorbidities + n.sessions + psych.meds.bl, data = kimd, na.action = na.exclude)
summary(full.bl)
summary(full.pt)

write.table(kimd_g, "/users/k1620287/brc_scratch/meta.analysis/pheno/pheno.data/kimd.phenotype.data.txt", col.names=T, row.names=F, quote=F, sep=" ")

##########################################################

gxte <- read.table("/users/k1620287/brc_scratch/meta.analysis/pheno/new.data/gxte.analysis.variables.030418.txt", head=T)

wb <- read.table("/users/k1620287/brc_scratch/meta.analysis/pheno/new.data/panicnet.analysis.variables.030418.txt", head=T)

ki <- read.table("/users/k1620287/brc_scratch/meta.analysis/pheno/new.data/karolinska.merged.analysis.variables.030418.txt", head=T)

gid <- read.table("/users/k1620287/brc_scratch/meta.analysis/analysis/adult_merged/binary/adult_merged.keeps", head=F)

names(gid)[1] <- "iid"

gid$V2 <- NULL

library(dplyr)

ki <- merge(gid, ki, by = "iid", type="left")

gxte <- merge(gid, gxte, by = "iid", type="left")

wb <- merge(gid, wb, by = "iid", type="left")

# MICE Variables imputed in above files:
# GxTE: gxte$psych.med.bl.imputed ("logreg"); gxte$primary.dsm.imputed ("polyreg"); gxte$n.sessions.imputed ("polr")
# KI: ki$psych.med.bl.imputed ("logreg"); ki$n.sessions.imputed ("polr")

gxte$primary.dsm <- gxte$primary.dsm.imputed

gxte$psych.med.bl <- gxte$psych.med.bl.imputed

gxte$n.sessions <- gxte$n.sessions.imputed

ki$psych.med.bl <- ki$psych.med.bl.imputed

ki$n.sessions <- ki$n.sessions.imputed

gxte$gender <- ifelse(gxte$sex==2,"female","male")

ki$gender <- ifelse(ki$sex==2,"female","male")

wb$gender <- ifelse(wb$sex==2,"female","male")

# recode cgi from 1 - 7 to 0 - 6
# cgi 0 = missing data recode by -1 then recode -1=NA
gxte$cgi.bl <- gxte$cgi.bl - 1

gxte$cgi.pt <- gxte$cgi.pt - 1

gxte$cgi.bl <- ifelse(gxte$cgi.bl==-1,NA,gxte$cgi.bl)

gxte$cgi.pt <- ifelse(gxte$cgi.pt==-1,NA,gxte$cgi.pt)

#### ADD MISSING VARIBLES - WINSORIZE OUTLIERS ####
gxte$n.sessions <- ifelse(gxte$n.sessions >30, 30, gxte$n.sessions)

gxte$treatment <- ifelse(gxte$site=="bochum.dental","dcbt",ifelse(gxte$n.sessions>7 & !gxte$site=="bochum.dental", "ecbt", "exp"))

wb$treatment <- "ecbt"

ki$treatment <- "icbt"

# LOCF imputation

ki$pdsssr.pt <- ifelse(is.na(ki$pdsssr.pt),ki$pdss.locf,ki$pdsssr.pt)

wb$pas.pt <- ifelse(is.na(wb$pas.pt),wb$pas.ia,wb$pas.pt)

#### GXTE ####
#### primary = cgi ####

gxte$primary.bl.z <- (gxte$cgi.bl/(sd(gxte$cgi.bl, na.rm=TRUE)))

gxte$primary.pt.z <- (gxte$cgi.pt/(sd(gxte$cgi.bl, na.rm=TRUE)))

gxte$primary.z.outcomes.res <- resid(lm(primary.pt.z ~ primary.bl.z + site + primary.dsm + age + sex + n.comorbidities + n.sessions + psych.med.bl, data = gxte, na.action = na.exclude))

#### WURZBURG #####
#### primary = pas ####

wb$primary.bl.z <- (wb$pas.bl/(sd(wb$pas.bl, na.rm=TRUE)))

wb$primary.pt.z <- (wb$pas.pt/(sd(wb$pas.bl, na.rm=TRUE)))

wb$primary.z.outcomes.res <- resid(lm(primary.pt.z ~ primary.bl.z + site + primary.dsm + age + sex + n.comorbidities + n.sessions, data = wb, na.action = na.exclude))

#### KAROLINSKA #####
#### primary = pdss ####

ki$primary.bl.z <- (ki$pdsssr.bl/(sd(ki$pdsssr.bl, na.rm=TRUE)))

ki$primary.pt.z <- (ki$pdsssr.pt/(sd(ki$pdsssr.bl, na.rm=TRUE)))

ki$primary.z.outcomes.res <- resid(lm(primary.pt.z ~ primary.bl.z + site + primary.dsm + age + sex + n.comorbidities + n.sessions + psych.med.bl, data = ki, na.action = na.exclude))

# Merge datasets

# z score data

gxte <- gxte[,c("fid","iid","age","sex","gender","cohort","site","treatment","primary.dsm","psych.med.bl","n.sessions","n.comorbidities","cgi.bl","cgi.pt","primary.bl.z","primary.pt.z","primary.z.outcomes.res")]

colnames(gxte) <- c("FID","IID","age","sex","gender","cohort","site","treatment","primary.dsm","psych.med.bl","n.sessions","n.comorbidities","primary.bl","primary.pt","primary.bl.z","primary.pt.z","primary.z.outcomes.res")

wb <- wb[,c("fid","iid","age","sex","gender","cohort","site","treatment","primary.dsm","psych.med.bl","n.sessions","n.comorbidities","pas.bl","pas.pt","primary.bl.z","primary.pt.z","primary.z.outcomes.res")]

colnames(wb) <- c("FID","IID","age","sex","gender","cohort","site","treatment","primary.dsm","psych.med.bl","n.sessions","n.comorbidities","primary.bl","primary.pt","primary.bl.z","primary.pt.z","primary.z.outcomes.res")

ki <- ki[,c("fid","iid","age","sex","gender","cohort","site","treatment","primary.dsm","psych.med.bl","n.sessions","n.comorbidities","pdsssr.bl","pdsssr.pt","primary.bl.z","primary.pt.z","primary.z.outcomes.res")]

colnames(ki) <- c("FID","IID","age","sex","gender","cohort","site","treatment","primary.dsm","psych.med.bl","n.sessions","n.comorbidities","primary.bl","primary.pt","primary.bl.z","primary.pt.z","primary.z.outcomes.res")

wb <- wb %>%
	mutate(
		FID = as.character(FID),IID = as.character(IID), age = as.integer(age),sex = as.factor(sex),
		cohort = as.factor(cohort),site = as.factor(site), treatment = as.factor(treatment), primary.dsm = as.factor(primary.dsm),
		psych.med.bl = as.factor(psych.med.bl),  n.comorbidities = as.integer(n.comorbidities), n.sessions = as.integer(n.sessions))

ki <- ki %>%
	mutate(
		FID = as.character(FID), IID = as.character(IID), age = as.integer(age), sex = as.factor(sex),
		cohort = as.factor(cohort), site = as.factor(site), treatment = as.factor(treatment), primary.dsm = as.factor(primary.dsm),
		psych.med.bl = as.factor(psych.med.bl), n.comorbidities = as.integer(n.comorbidities), n.sessions = as.integer(n.sessions))

gxte <- gxte %>%
	mutate(
		FID = as.character(FID), IID = as.character(IID), age = as.integer(age), sex = as.factor(sex),
		cohort = as.factor(cohort), site = as.factor(site), treatment = as.factor(treatment), primary.dsm = as.factor(primary.dsm),
		psych.med.bl = as.factor(psych.med.bl), n.comorbidities = as.integer(n.comorbidities), n.sessions = as.integer(n.sessions))

wc <- dplyr::bind_rows(ki, wb, gxte)

wc$primary.z.outcomes.res.whole <-  resid(lm(primary.pt.z ~ primary.bl.z + cohort + site + primary.dsm + age + sex + n.comorbidities + n.sessions + psych.med.bl, data = wc, na.action = na.exclude))

wc <- wc %>%
	mutate(
		FID = as.character(FID), IID = as.character(IID), age = as.integer(age), sex = as.factor(sex), gender = as.factor(gender),
		cohort = as.factor(cohort), site = as.factor(site), treatment = as.factor(treatment),primary.dsm = as.factor(primary.dsm),
		psych.med.bl = as.factor(psych.med.bl), n.comorbidities = as.integer(n.comorbidities), n.sessions = as.integer(n.sessions))

# This is the phenotype I used:

wc$primary.z.outcomes.res <-  resid(lm(primary.z.outcomes.res ~ cohort, data = wc, na.action = na.exclude))

write.table(wc, "/users/k1620287/brc_scratch/meta.analysis/pheno/pheno.data/adult.cohorts.for.pheno.analysis", col.names=T, row.names=F, quote=F, sep=" ")

gxte.outcomes <- gxte[,c("FID","IID","primary.z.outcomes.res")]

wb.outcomes <- wb[,c("FID","IID","primary.z.outcomes.res")]

ki.outcomes <- ki[,c("FID","IID","primary.z.outcomes.res")]

wc.outcomes <- wc[,c("FID","IID","primary.z.outcomes.res")]

write.table(gxte.outcomes, "/users/k1620287/brc_scratch/meta.analysis/pheno/pheno.data/gxte.primary.z.outcomes.res.txt", col.names=T, row.names=F, quote=F, sep=" ")

write.table(ki.outcomes, "/users/k1620287/brc_scratch/meta.analysis/pheno/pheno.data/ki.primary.z.outcomes.res.txt", col.names=T, row.names=F, quote=F, sep=" ")

write.table(wb.outcomes, "/users/k1620287/brc_scratch/meta.analysis/pheno/pheno.data/wb.primary.z.outcomes.res.txt", col.names=T, row.names=F, quote=F, sep=" ")

write.table(wc.outcomes, "/users/k1620287/brc_scratch/meta.analysis/pheno/pheno.data/wc.primary.z.outcomes.res.txt", col.names=T, row.names=F, quote=F, sep=" ")

#' > /users/k1620287/brc_scratch/meta.analysis/pheno/scripts/make.phenotype.datasets.R

module load general/R

R --file=/users/k1620287/brc_scratch/meta.analysis/pheno/scripts/make.phenotype.datasets.R

EOT

qsub /users/k1620287/brc_scratch/meta.analysis/pheno/scripts/make.phenotype.datasets.sh


######
# MERGED ANALYSIS
library(data.table)
library(dplyr)

ki <- fread("/users/k1620287/brc_scratch/meta.analysis/pheno/pheno.data/kimd.primary.z.pt.res.txt", head=T)
ad <- fread("/users/k1620287/brc_scratch/meta.analysis/pheno/pheno.data/wc.primary.z.outcomes.res.txt", head=T)
ch <- fread("/users/k1620287/brc_scratch/meta.analysis/pheno/pheno.data/gxtc.primary.z.outcomes.res.txt", head=T)

ki$cohort <- "KI"
ad$cohort <- "AD"
ch$cohort <- "CH"

wc <- bind_rows(ki, ad, ch)

wc <- wc %>% mutate(FID = as.character(FID), IID = as.character(IID), cohort = as.factor(cohort))

wc$primary.z.outcomes.res <-  resid(lm(primary.z.outcomes.res ~ cohort, data = wc, na.action = na.exclude))

wc$cohort <- NULL

write.table(wc, "/users/k1620287/brc_scratch/meta.analysis/pheno/pheno.data/merged.primary.z.outcomes.res.txt", col.names=T, row.names=F, quote=F, sep=" ")

adult <- bind_rows(ki, ad)

adult <- adult %>% mutate(FID = as.character(FID), IID = as.character(IID), cohort = as.factor(cohort))

adult$primary.z.outcomes.res <-  resid(lm(primary.z.outcomes.res ~ cohort, data = adult, na.action = na.exclude))

adult$cohort <- NULL

write.table(adult, "/users/k1620287/brc_scratch/meta.analysis/pheno/pheno.data/adult_merged.primary.z.outcomes.res.txt", col.names=T, row.names=F, quote=F, sep=" ")

anx <- bind_rows(ad, ch)

anx <- anx %>% mutate(FID = as.character(FID), IID = as.character(IID), cohort = as.factor(cohort))

anx$primary.z.outcomes.res <-  resid(lm(primary.z.outcomes.res ~ cohort, data = anx, na.action = na.exclude))

anx$cohort <- NULL

write.table(anx, "/users/k1620287/brc_scratch/meta.analysis/pheno/pheno.data/anx_merged.primary.z.outcomes.res.txt", col.names=T, row.names=F, quote=F, sep=" ")



##############################################################################################
# DATASETS FOR TABLES:
##############################################################################################

gid <- read.table("/users/k1620287/brc_scratch/meta.analysis/analysis/adult_merged/binary/adult_merged.keeps", head=F)

names(gid)[1] <- "iid"

gid$V2 <- NULL

gxte <- read.table("/users/k1620287/brc_scratch/meta.analysis/pheno/new.data/gxte.analysis.variables.030418.txt", head=T)

wb <- read.table("/users/k1620287/brc_scratch/meta.analysis/pheno/new.data/panicnet.analysis.variables.030418.txt", head=T)

ki <- read.table("/users/k1620287/brc_scratch/meta.analysis/pheno/original.data/Clinical variables_PD_TOCHRIS_171009_CR_EDIT.csv", ",", head=T)

colnames(ki) <- c("id","gp","fid","iid","sex","age","primary.diagnosis","co_1","co_2","co_3","co_4","n.comorbidities","site.incorrect","site","pdsssr.bl","pdsssr.pt","pdsssr.locf","pdssr.fu",
"madrs.bl","madrs.pt","madrs.fu","n.sessions","psych.med.bl","qc.gw.data")

ki <- ki[,c("iid","sex","age","primary.diagnosis","n.comorbidities","site","pdsssr.bl","pdsssr.pt","pdsssr.locf","n.sessions","psych.med.bl","qc.gw.data")]

ki$primary.dsm <- ifelse(ki$primary.diagnosis=="AGORAPHOBIA",30022,ifelse(ki$primary.diagnosis=="GAD",30002,ifelse(ki$primary.diagnosis=="MDD",311,ifelse(ki$primary.diagnosis=="OCD",3003,
ifelse(ki$primary.diagnosis=="PD",30001,ifelse(ki$primary.diagnosis=="SAD",30023,NA))))))

ki <- ki[which(ki$qc.gw.data==1),]
ki$fid <- ki$iid
ki$cohort <- "karolinska"
ki$sex <- ifelse(ki$sex==0, 2, ki$sex)

library(dplyr)

ki <- merge(gid, ki, by = "iid", type="left")

gxte <- merge(gid, gxte, by = "iid", type="left")

wb <- merge(gid, wb, by = "iid", type="left")

gxte$gender <- ifelse(gxte$sex==2,"female","male")
ki$gender <- ifelse(ki$sex==2,"female","male")
wb$gender <- ifelse(wb$sex==2,"female","male")

gxte$cgi.bl <- gxte$cgi.bl - 1
gxte$cgi.pt <- gxte$cgi.pt - 1
gxte$cgi.bl <- ifelse(gxte$cgi.bl==-1,NA,gxte$cgi.bl)
gxte$cgi.pt <- ifelse(gxte$cgi.pt==-1,NA,gxte$cgi.pt)

#### ADD MISSING VARIBLES - WINSORIZE OUTLIERS ####
gxte$n.sessions <- ifelse(gxte$n.sessions >30, 30, gxte$n.sessions)
gxte$treatment <- ifelse(gxte$site=="bochum.dental","dcbt",ifelse(gxte$n.sessions>7 & !gxte$site=="bochum.dental", "ecbt", "exp"))
wb$treatment <- "ecbt"
ki$treatment <- "icbt"

# LOCF imputation

ki$pdsssr.pt <- ifelse(is.na(ki$pdsssr.pt),ki$pdsssr.locf,ki$pdsssr.pt)

wb$pas.pt <- ifelse(is.na(wb$pas.pt),wb$pas.ia,wb$pas.pt)

#### GXTE ####
#### primary = cgi ####

gxte$primary.bl.z <- (gxte$cgi.bl/(sd(gxte$cgi.bl, na.rm=TRUE)))
gxte$primary.pt.z <- (gxte$cgi.pt/(sd(gxte$cgi.bl, na.rm=TRUE)))

#### WURZBURG #####
#### primary = pas ####

wb$primary.bl.z <- (wb$pas.bl/(sd(wb$pas.bl, na.rm=TRUE)))
wb$primary.pt.z <- (wb$pas.pt/(sd(wb$pas.bl, na.rm=TRUE)))

#### KAROLINSKA #####
#### primary = pdss ####

ki$primary.bl.z <- (ki$pdsssr.bl/(sd(ki$pdsssr.bl, na.rm=TRUE)))
ki$primary.pt.z <- (ki$pdsssr.pt/(sd(ki$pdsssr.bl, na.rm=TRUE)))

# Merge datasets
# z score data

gxte <- gxte[,c("fid","iid","age","sex","gender","cohort","site","treatment","primary.dsm","psych.med.bl","n.sessions","n.comorbidities","cgi.bl","cgi.pt","primary.bl.z","primary.pt.z")]
colnames(gxte) <- c("fid","iid","age","sex","gender","cohort","site","treatment","primary.dsm","psych.med.bl","n.sessions","n.comorbidities","primary.bl","primary.pt","primary.bl.z","primary.pt.z")

wb <- wb[,c("fid","iid","age","sex","gender","cohort","site","treatment","primary.dsm","psych.med.bl","n.sessions","n.comorbidities","pas.bl","pas.pt","primary.bl.z","primary.pt.z")]
colnames(wb) <- c("fid","iid","age","sex","gender","cohort","site","treatment","primary.dsm","psych.med.bl","n.sessions","n.comorbidities","primary.bl","primary.pt","primary.bl.z","primary.pt.z")

ki <- ki[,c("fid","iid","age","sex","gender","cohort","site","treatment","primary.dsm","psych.med.bl","n.sessions","n.comorbidities","pdsssr.bl","pdsssr.pt","primary.bl.z","primary.pt.z")]
colnames(ki) <- c("fid","iid","age","sex","gender","cohort","site","treatment","primary.dsm","psych.med.bl","n.sessions","n.comorbidities","primary.bl","primary.pt","primary.bl.z","primary.pt.z")

wb <- wb %>%
  mutate(
    fid = as.character(fid), iid = as.character(iid),  age = as.integer(age),  sex = as.factor(sex),
    cohort = as.factor(cohort), site = as.factor(site), treatment = as.factor(treatment), primary.dsm = as.factor(primary.dsm),
	psych.med.bl = as.factor(psych.med.bl),  n.comorbidities = as.integer(n.comorbidities), n.sessions = as.integer(n.sessions)
  )

ki <- ki %>%
  mutate(
    fid = as.character(fid), iid = as.character(iid), age = as.integer(age), sex = as.factor(sex),
    cohort = as.factor(cohort), site = as.factor(site), treatment = as.factor(treatment), primary.dsm = as.factor(primary.dsm),
    psych.med.bl = as.factor(psych.med.bl), n.comorbidities = as.integer(n.comorbidities), n.sessions = as.integer(n.sessions)
  )

gxte <- gxte %>%
  mutate(
    fid = as.character(fid),iid = as.character(iid), age = as.integer(age), sex = as.factor(sex),
    cohort = as.factor(cohort), site = as.factor(site), treatment = as.factor(treatment), primary.dsm = as.factor(primary.dsm),
    psych.med.bl = as.factor(psych.med.bl), n.comorbidities = as.integer(n.comorbidities), n.sessions = as.integer(n.sessions)
  )

wc <- dplyr::bind_rows(ki, wb, gxte)

wc$primary.z.outcomes.res.whole <-  resid(lm(primary.pt.z ~ primary.bl.z + cohort + site + primary.dsm + age + sex + n.comorbidities + n.sessions + psych.med.bl, data = wc, na.action = na.exclude))

wc <- wc %>%
  mutate(
    fid = as.character(fid), iid = as.character(iid), age = as.integer(age), sex = as.factor(sex),
    gender = as.factor(gender), cohort = as.factor(cohort),  site = as.factor(site), treatment = as.factor(treatment),
    primary.dsm = as.factor(primary.dsm), psych.med.bl = as.factor(psych.med.bl), n.comorbidities = as.integer(n.comorbidities),
    n.sessions = as.integer(n.sessions))

wc$site <- factor(wc$site, levels=c("bochum.clinic","bochum.dental","braunschweig","CLINIC","RCT","panicnet.1","panicnet.2"), labels=c("BC","BD","BS","KI T","KI C","PNC I", "PNC II"))

wc$cohort <- factor(wc$cohort, levels=c("gxte","karolinska","panicnet"), labels=c("Bochum & Braunschweig","Karolinska","Panicnet Consortium"))

wc$Female <- ifelse(wc$sex==2,1,0)

write.table(wc, "/users/k1620287/brc_scratch/meta.analysis/pheno/pheno.data/adult.cohorts.for.table", col.names=T, row.names=F, quote=F, sep=",")


# Summarise missing categorical data as 0 for % prevalence of 1

wc$primary.dsm <- with(wc, ifelse(is.na(primary.dsm),0,primary.dsm))

wc$psych.med.bl <- with(wc, ifelse(is.na(psych.med.bl),0,psych.med.bl))


# Create tables
sink("/users/k1620287/brc_scratch/meta.analysis/pheno/pheno.data/analysis.sample.descriptives")

cat("Female (N/[%]) by Site")
f <- with(wc, table(Female, site))
f.p <- prop.table(f, 2)
print(with(wc, table(Female, site)))
print(f.p*100, digits = 3)

cat("Female (N/[%])")
f <- with(wc, table(Female))
f.p <- prop.table(f)
print(with(wc, table(Female)))
print(f.p*100, digits = 3)

cat("Female Chi-square")
print(chisq.test(wc$sex,wc$site))

cat("Main diagnosis (N/[%]) by Site")
d <- with(wc, table(primary.dsm, site))
d.p <- prop.table(d, 2)
print(with(wc, table(primary.dsm, site)))
print(d.p*100, digits = 3)

cat("Main diagnosis (N/[%])")
d <- with(wc, table(primary.dsm))
d.p <- prop.table(d)
print(with(wc, table(primary.dsm)))
print(d.p*100, digits = 3)

cat("Main diagnosis Chi-square")
print(chisq.test(wc$primary.dsm,wc$site))

cat("Taking psych. meds (N/[%]) by Site")
m <- with(wc, table(psych.med.bl, site))
m.p <- prop.table(m, 2)
print(with(wc, table(psych.med.bl, site)))
print(m.p*100, digits = 3)

cat("Taking psych. meds (N/[%])")
m <- with(wc, table(psych.med.bl))
m.p <- prop.table(m)
print(with(wc, table(psych.med.bl)))
print(m.p*100, digits = 3)

cat("Taking psych. meds Chi-square")
print(chisq.test(wc$psych.med.bl,wc$site))

cat("Means by site")
m.age.site <- aggregate(age ~ site, wc, function(x) c(mean = round(mean(x),1), sd = round(sd(x), 1)))
m.co.site <- aggregate(n.comorbidities ~ site, wc, function(x) c(mean = round(mean(x),1), sd = round(sd(x), 1)))
m.sess.site <- aggregate(n.sessions ~ site, wc, function(x) c(mean = round(mean(x),1), sd = round(sd(x), 1)))
m.bl.z.site <- aggregate(primary.bl.z ~ site, wc, function(x) c(mean = round(mean(x),2), sd = round(sd(x), 2)))
m.pt.z.site <- aggregate(primary.pt.z ~ site, wc, function(x) c(mean = round(mean(x),2), sd = round(sd(x), 2)))
mean.site <- plyr::join_all(list(m.age.site,m.co.site,m.sess.site,m.bl.z.site,m.pt.z.site ), by = "site", type="left")
print(mean.site)

cat("Means")
print(with(wc, c(mean = round(mean(age, na.rm=T),1), sd = round(sd(age, na.rm=T), 1))))
print(with(wc, c(mean = round(mean(n.comorbidities, na.rm=T),1), sd = round(sd(n.comorbidities, na.rm=T), 1))))
print(with(wc, c(mean = round(mean(n.sessions, na.rm=T),1), sd = round(sd(n.sessions, na.rm=T), 1))))
print(with(wc, c(mean = round(mean(primary.bl.z, na.rm=T),2), sd = round(sd(primary.bl.z, na.rm=T), 2))))
print(with(wc, c(mean = round(mean(primary.pt.z, na.rm=T),2), sd = round(sd(primary.pt.z, na.rm=T), 2))))

cat("Anovas")
print(summary(aov(age ~ site, data=wc)))
print(summary(aov(n.sessions ~ site, data=wc)))
print(summary(aov(n.comorbidities ~ site, data=wc)))
print(summary(aov(primary.bl.z ~ site, data=wc)))
print(summary(aov(primary.pt.z ~ site, data=wc)))

sink()



dat <- read.table("/users/k1620287/brc_scratch/meta.analysis/pheno/pheno.data/adult.cohorts.for.pheno.analysis", head=T)

library(lme4)
library(dplyr)

dat$primary.dsm <- ifelse(dat$primary.dsm==29620, 50000, dat$primary.dsm)

dat <- dat %>%
  mutate(
    age = as.integer(age), sex = as.factor(sex),
    gender = as.factor(gender), cohort = as.factor(cohort),  site = as.factor(site), treatment = as.factor(treatment),
    primary.dsm = as.factor(primary.dsm), psych.med.bl = as.factor(psych.med.bl), n.comorbidities = as.integer(n.comorbidities),
    n.sessions = as.integer(n.sessions))

full.bl <- lmer(primary.bl.z ~ (1 | cohort) + (1 | site) + age + gender + n.comorbidities + n.sessions + psych.med.bl + primary.dsm, data = dat, na.action = na.exclude)
full.pt <- lmer(primary.pt.z ~ (1 | cohort) + (1 | site) + primary.bl.z + age + gender + n.comorbidities + n.sessions + psych.med.bl + primary.dsm, data = dat, na.action = na.exclude)
summary(full.bl)
summary(full.pt)

full.bl <- lmer(primary.bl.z ~ (1 | cohort) + (1 | site) + (1 | primary.dsm) + age + gender + n.comorbidities + n.sessions + psych.med.bl, data = dat, na.action = na.exclude)
full.pt <- lmer(primary.pt.z ~ (1 | cohort) + (1 | site) + (1 | primary.dsm) + primary.bl.z + age + gender + n.comorbidities + n.sessions + psych.med.bl, data = dat, na.action = na.exclude)
summary(full.bl)
summary(full.pt)

full.bl1 <- lmer(primary.bl.z ~ (1 | cohort) + (1 | site) + age + sex + n.comorbidities + n.sessions + psych.med.bl + primary.dsm, data = dat1, na.action = na.exclude)
full.pt1 <- lmer(primary.pt.z ~ (1 | cohort) + (1 | site) + primary.bl.z + age + sex + n.comorbidities + n.sessions + psych.med.bl + primary.dsm, data = dat1, na.action = na.exclude)
summary(full.bl1)
summary(full.pt1)
