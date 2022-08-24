####################################
# Meta-analysis of response to CBT #
####################################

# Final attempt:

cd /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files

for i in adult_anx.loco.mlma.results  child_anx.loco.mlma.results  ki_mdd.loco.mlma.results
do
awk '{print $1,$2,$3}' $i | tail -n +2 > $i.snplist
done

cat *anx.loco.mlma.results.snplist > anx.meta.snplist

awk '{print $1,$2,$3}' anx.meta.snplist | sort | uniq -d > anx.snplist.keeps

cat anx.snplist.keeps ki_mdd.loco.mlma.results.snplist > full.meta.snplist

awk '{print $1,$2,$3}' full.meta.snplist | sort | uniq -d > full.meta.snplist.keeps

module add general/R;R

ad <- read.table("adult_anx.loco.mlma.results",head=T)
ch <- read.table("child_anx.loco.mlma.results",head=T)
ki <- read.table("ki_mdd.loco.mlma.results",head=T)

snps <- read.table("full.meta.snplist.keeps",head=F)
snps$V1 <- NULL
snps$V3 <- NULL
names(snps)[1] <- "SNP"

ad <- merge(snps, ad, by = "SNP", type="left")
ch <- merge(snps, ch, by = "SNP", type="left")
ki <- merge(snps, ki, by = "SNP", type="left")

write.table(ad, "adult_anx.loco.mlma.results.common_snps", col.names=T, row.names=F, quote=F, sep=" ")
write.table(ch, "child_anx.loco.mlma.results.common_snps", col.names=T, row.names=F, quote=F, sep=" ")
write.table(ki, "ki_mdd.loco.mlma.results.common_snps", col.names=T, row.names=F, quote=F, sep=" ")


cd  /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/sign_tests/
rm sign.test.sh
cat <<'EOT'>> sign.test.sh
#!/bin/sh
#$ -V
#$ -m ea
#$ -pe smp 1
#$ -l h_vmem=4G
#$ -S /bin/sh
#$ -m be
#$ -o /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/sign_tests/
#$ -e /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/sign_tests/

export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

module add general/R

cd  /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/sign_tests/

path=/mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files
base=$(basename $1 | sed -e 's/.loco.mlma.results.common_snps//g')
target=$(basename $2 | sed -e 's/.loco.mlma.results.common_snps//g')
mkdir /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/sign_tests/$base.$target
out=sign_tests/$base.$target

sh /mnt/lustre/groups/ukbiobank/usr/KLP/GAD7/UKBB/scripts/SignTest.sh \
-b $path/$1 \
-t $path/$2 \
-l $path/$out/Clumped.P1/$1.WG.clump \
-p "0.00005 0.0005 0.005 0.05 0.5" \
-o $path/$out/$base.$target \
-v

#-k 1000 -r 0.01 -s \

EOT

echo 'ki_mdd.loco.mlma.results.common_snps adult_anx.loco.mlma.results.common_snps
ki_mdd.loco.mlma.results.common_snps child_anx.loco.mlma.results.common_snps
adult_anx.loco.mlma.results.common_snps ki_mdd.loco.mlma.results.common_snps
adult_anx.loco.mlma.results.common_snps child_anx.loco.mlma.results.common_snps
child_anx.loco.mlma.results.common_snps ki_mdd.loco.mlma.results.common_snps
child_anx.loco.mlma.results.common_snps adult_anx.loco.mlma.results.common_snps' > sign_test_list

while read sign_test_list; do qsub sign.test.sh $sign_test_list;  done < sign_test_list

# qsub sign.test.sh ki_mdd.loco.mlma.results.common_snps adult_anx.loco.mlma.results.common_snps

qsub sign.test.sh ki_mdd.loco.mlma.results.common_snps child_anx.loco.mlma.results.common_snps
qsub sign.test.sh adult_anx.loco.mlma.results.common_snps ki_mdd.loco.mlma.results.common_snps

#qsub sign.test.sh adult_anx.loco.mlma.results.common_snps child_anx.loco.mlma.results.common_snps
#qsub sign.test.sh child_anx.loco.mlma.results.common_snps ki_mdd.loco.mlma.results.common_snps
#qsub sign.test.sh child_anx.loco.mlma.results.common_snps adult_anx.loco.mlma.results.common_snps



########################################################################

cd /users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/assoc/

awk '{print $1,$2,$3}' adult_zscore_gwas_results | tail -n +2 > adult.snplist

awk '{print $1,$2,$3}' child_zscore_gwas_results | tail -n +2 > child.snplist

cat adult.snplist child.snplist > all.meta.snplist

awk '{print $1,$2,$3}' all.meta.snplist | sort | uniq -d > meta.snplist.keeps


adult <- read.table("/users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/assoc/adult_zscore_gwas_results",head=T)
child <- read.table("/users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/assoc/child_zscore_gwas_results",head=T)

snps <- read.table("/users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/assoc/meta.snplist.keeps",head=F)
snps$V1 <- NULL
snps$V3 <- NULL
names(snps)[1] <- "SNP"

adult <- merge(snps, adult, by = "SNP", type="left")
child <- merge(snps, child, by = "SNP", type="left")

write.table(adult, "/users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/assoc/adult_zscore_gwas_results_common_snps", col.names=T, row.names=F, quote=F, sep=" ")
write.table(child, "/users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/assoc/child_zscore_gwas_results_common_snps", col.names=T, row.names=F, quote=F, sep=" ")



wb <- read.table("/users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/assoc/wb_zscore_gwas_results",head=T)

ki <- read.table("/users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/assoc/ki_zscore_gwas_results",head=T)

gxte <- read.table("/users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/assoc/gxte_zscore_gwas_results",head=T)

snps <- read.table("/users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/assoc/adult.snplist",head=F)

snps$V1 <- NULL
snps$V3 <- NULL

names(snps)[1] <- "SNP"

wb <- merge(snps, wb, by = "SNP", type="left")
ki <- merge(snps, ki, by = "SNP", type="left")
gxte <- merge(snps, gxte, by = "SNP", type="left")

write.table(wb, "/users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/assoc/wb_zscore_gwas_results_common_snps", col.names=T, row.names=F, quote=F, sep=" ")
write.table(ki, "/users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/assoc/ki_zscore_gwas_results_common_snps", col.names=T, row.names=F, quote=F, sep=" ")
write.table(gxte, "/users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/assoc/gxte_zscore_gwas_results_common_snps", col.names=T, row.names=F, quote=F, sep=" ")


# Script 12. SIGN TESTS #

rm /users/k1620287/brc_scratch/meta.analysis/scripts/adult.child.sign.test.sh
cat <<'EOT'>> /users/k1620287/brc_scratch/meta.analysis/scripts/adult.child.sign.test.sh
#!/bin/sh
#$ -V
#$ -m ea
#$ -pe smp 8
#$ -l h_vmem=9G
#$ -S /bin/sh
#$ -m be
#$ -o /users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/sign/adult/
#$ -e /users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/sign/adult/

export MKL_NUM_THREADS=8
export NUMEXPR_NUM_THREADS=8
export OMP_NUM_THREADS=8

module add general/R

mkdir /users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/sign/adult

cd /users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/sign/adult

sh /mnt/lustre/groups/ukbiobank/usr/KLP/GAD7/UKBB/scripts/SignTest.sh \
-b /users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/assoc/adult_zscore_gwas_results_common_snps \
-t "/users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/assoc/child_zscore_gwas_results_common_snps" \
-k 1000 -r 0.01 -s \
-p "0.00005 0.0005 0.005 0.05 0.5" \
-o /users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/sign/adult/adult.child.sign.test \
-v
EOT

qsub /users/k1620287/brc_scratch/meta.analysis/scripts/adult.child.sign.test.sh

rm /users/k1620287/brc_scratch/meta.analysis/scripts/child.adult.sign.test.sh
cat <<'EOT'>> /users/k1620287/brc_scratch/meta.analysis/scripts/child.adult.sign.test.sh
#!/bin/sh
#$ -V
#$ -m ea
#$ -pe smp 8
#$ -l h_vmem=9G
#$ -S /bin/sh
#$ -m be
#$ -o /users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/sign/child/
#$ -e /users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/sign/child/

export MKL_NUM_THREADS=8
export NUMEXPR_NUM_THREADS=8
export OMP_NUM_THREADS=8

module add general/R

mkdir /users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/sign/child
cd /users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/sign/child

sh /mnt/lustre/groups/ukbiobank/usr/KLP/GAD7/UKBB/scripts/SignTest.sh \
-b /users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/assoc/child_zscore_gwas_results_common_snps \
-t "/users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/assoc/adult_zscore_gwas_results_common_snps" \
-k 1000 -r 0.01 -s \
-p "0.00005 0.0005 0.005 0.05 0.5" \
-o /users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/sign/child/child.adult.sign.test \
-v

EOT
qsub /users/k1620287/brc_scratch/meta.analysis/scripts/child.adult.sign.test.sh

rm /users/k1620287/brc_scratch/meta.analysis/scripts/ki.sign.test.sh
cat <<'EOT'>> /users/k1620287/brc_scratch/meta.analysis/scripts/ki.sign.test.sh
#!/bin/sh
#$ -V
#$ -m ea
#$ -pe smp 8
#$ -l h_vmem=9G
#$ -S /bin/sh
#$ -m be
#$ -o /users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/sign/ki/
#$ -e /users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/sign/ki/

export MKL_NUM_THREADS=8
export NUMEXPR_NUM_THREADS=8
export OMP_NUM_THREADS=8

module add general/R

mkdir /users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/sign/ki
cd /users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/sign/ki

sh /mnt/lustre/groups/ukbiobank/usr/KLP/GAD7/UKBB/scripts/SignTest.sh \
-b /users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/assoc/ki_zscore_gwas_results_common_snps \
-t "/users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/assoc/wb_zscore_gwas_results_common_snps \
/users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/assoc/gxte_zscore_gwas_results_common_snps" \
-k 1000 -r 0.01 -s \
-p "0.0005 0.005 0.05 0.5" \
-o /users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/sign/ki/ki.sign.test \
-v

EOT
qsub /users/k1620287/brc_scratch/meta.analysis/scripts/ki.sign.test.sh

rm /users/k1620287/brc_scratch/meta.analysis/scripts/wb.sign.test.sh
cat <<'EOT'>> /users/k1620287/brc_scratch/meta.analysis/scripts/wb.sign.test.sh
#!/bin/sh
#$ -V
#$ -m ea
#$ -pe smp 8
#$ -l h_vmem=9G
#$ -S /bin/sh
#$ -m be
#$ -o /users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/sign/wb/
#$ -e /users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/sign/wb/

export MKL_NUM_THREADS=8
export NUMEXPR_NUM_THREADS=8
export OMP_NUM_THREADS=8

module add general/R

mkdir /users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/sign/wb
cd /users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/sign/wb

sh /mnt/lustre/groups/ukbiobank/usr/KLP/GAD7/UKBB/scripts/SignTest.sh \
-b /users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/assoc/wb_zscore_gwas_results_common_snps \
-t "/users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/assoc/ki_zscore_gwas_results_common_snps \
/users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/assoc/gxte_zscore_gwas_results_common_snps" \
-k 1000 -r 0.01 -s \
-p "0.0005 0.005 0.05 0.5" \
-o /users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/sign/wb/wb.sign.test \
-v

EOT
qsub /users/k1620287/brc_scratch/meta.analysis/scripts/wb.sign.test.sh


rm /users/k1620287/brc_scratch/meta.analysis/scripts/gxte.sign.test.sh
cat <<'EOT'>> /users/k1620287/brc_scratch/meta.analysis/scripts/gxte.sign.test.sh
#!/bin/sh
#$ -V
#$ -m ea
#$ -pe smp 8
#$ -l h_vmem=9G
#$ -S /bin/sh
#$ -m be
#$ -o /users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/sign/gxte/
#$ -e /users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/sign/gxte/

export MKL_NUM_THREADS=8
export NUMEXPR_NUM_THREADS=8
export OMP_NUM_THREADS=8

module add general/R

mkdir /users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/sign/gxte
cd /users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/sign/gxte

sh /mnt/lustre/groups/ukbiobank/usr/KLP/GAD7/UKBB/scripts/SignTest.sh \
-b /users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/assoc/gxte_zscore_gwas_results_common_snps \
-t "/users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/assoc/ki_zscore_gwas_results_common_snps \
/users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/assoc/wb_zscore_gwas_results_common_snps" \
-k 1000 -r 0.01 -s \
-p "0.0005 0.005 0.05 0.5" \
-o /users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/sign/gxte/gxte.sign.test \
-v

EOT
qsub /users/k1620287/brc_scratch/meta.analysis/scripts/gxte.sign.test.sh






























# Munge summary statistics for consistency in snps across all traits

rm /users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/munge.for.sign.test.sh
cat <<'EOT'>> /users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/munge.for.sign.test.sh
#!/bin/sh
#$-S /bin/sh
#$ -V
#$ -pe smp 8
#$ -l h_vmem=9G
#$ -l h_rt=72:00:00
#$ -m be
#$ -o /mnt/lustre/groups/ukbiobank/usr/crayner/ukb18177/assoc/output
#$ -e /mnt/lustre/groups/ukbiobank/usr/crayner/ukb18177/assoc/output

export MKL_NUM_THREADS=8
export NUMEXPR_NUM_THREADS=8
export OMP_NUM_THREADS=8

cd /mnt/lustre/groups/ukbiobank/usr/crayner/ukb18177/ldsc/
module add utilities/anaconda/2.5.0

wb=/users/k1620287/brc_scratch/meta.analysis/analysis/wb/results/wb.primary.z.response.res.assoc.linear.a2.space
ki=/users/k1620287/brc_scratch/meta.analysis/analysis/ki/results/ki.primary.z.response.res.assoc.linear.a2.space
gxte=/users/k1620287/brc_scratch/meta.analysis/analysis/gxte/results/gxte.primary.z.response.res.assoc.linear.a2.space
gxtc=/users/k1620287/brc_scratch/meta.analysis/analysis/gxtc/results/gxtc.primary.z.response.res.assoc.linear.a2.space
adult_merged=/users/k1620287/brc_scratch/meta.analysis/analysis/adult_merged/results/adult.primary.z.response.res.assoc.linear.a2.space
adhd=/users/k1620287/brc_scratch/meta.analysis/analysis/pgs/final_scores/bin/adhd_demontis_2017_cc_pgs_edit
asd=/users/k1620287/brc_scratch/meta.analysis/analysis/pgs/final_scores/bin/asd_pgc_2017_cc_pgs_edit
mdd2=/users/k1620287/brc_scratch/meta.analysis/analysis/pgs/final_scores/bin/mdd2_pgc_2016_cc_pgs_edit
mdd3=/users/k1620287/brc_scratch/meta.analysis/analysis/pgs/final_scores/bin/mdd_pgc23_2018_cc_pgs_edit
scz=/users/k1620287/brc_scratch/meta.analysis/analysis/pgs/final_scores/bin/scz_ripke_2015_cc_pgs_edit
bmi=/users/k1620287/brc_scratch/meta.analysis/analysis/pgs/final_scores/quant/bmi_ukb_2017_q_pgs_edit
dnt=/users/k1620287/brc_scratch/meta.analysis/analysis/pgs/final_scores/quant/dep.not.trauma_coleman_2018_q_pgs_edit
ea2=/users/k1620287/brc_scratch/meta.analysis/analysis/pgs/final_scores/quant/edu.ea2_okbay_2016_pgs_edit
ea3=/users/k1620287/brc_scratch/meta.analysis/analysis/pgs/final_scores/quant/edu.ea3_lee_2018_q_pgs_edit
env=/users/k1620287/brc_scratch/meta.analysis/analysis/pgs/final_scores/quant/env_keers_2016_q_pgs_edit
gad=/users/k1620287/brc_scratch/meta.analysis/analysis/pgs/final_scores/quant/gad_purves_2018_pgs_edit
iq=/users/k1620287/brc_scratch/meta.analysis/analysis/pgs/final_scores/quant/iq_coleman_2017_q_pgs_edit
neu=/users/k1620287/brc_scratch/meta.analysis/analysis/pgs/final_scores/quant/neur_okbay_2016_q_pgs_edit
swb=/users/k1620287/brc_scratch/meta.analysis/analysis/pgs/final_scores/quant/swb_okbay_2016_q_pgs_edit
trsk=/users/k1620287/brc_scratch/meta.analysis/analysis/pgs/final_scores/quant/tr.seek_rayner_2018_pgs_edit
output_path=/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats

for i in BODY08 ALCO01 CONS01 DEPR02 EXTR01 INCO03 INSO01 NEUR02B OPEN01 SUBJ01 TIRE01 EDUC03 INTE01 INTE03 ADHD05 ALCD03 ANOR02 ANXI03 AUTI05 BIPO01 BORD01 CANU01 CROS01 DEPR05 SCHI02 MIGR01
do
scp /mnt/lustre/groups/ukbiobank/sumstats/munged_noMHC/${i}_noMHC.sumstats.gz ./

for i in ${wb} ${ki} ${gxte} ${gxtc} ${adult_merged} ${adhd} ${asd} ${mdd2} ${mdd3} ${scz} ${bmi} ${dnt} ${ea2} ${ea3} ${env} ${gad} ${iq} ${neu} ${swb} ${trsk}
do
out=$(basename $i)

/opt/apps/general/anaconda/2.5.0/anaconda2/bin/python /mnt/lustre/groups/ukbiobank/usr/crayner/ldsc/munge_sumstats.py \
--sumstats ${i} \
--out ${output_path}/${out} \
--merge-alleles /mnt/lustre/groups/ukbiobank/usr/crayner/ldsc/w_hm3.snplist \
--N

done
EOT
qsub /users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/munge.for.sign.test.sh


# Reformat sumstats for use in sign tests

ls /users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/munged > /users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/scripts/sumstats_list

while read list
do echo 'dat <- read.table(gzfile("/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/munged/'$list'.gz"), header=T, na.strings=c(""," ","NA"), fill = TRUE)
dat$P <- 2*pnorm(-abs(dat$Z))
dat$BETA <- dat$Z
dat <- na.omit(dat)
write.table(dat, "/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/'$list'.p", quote=FALSE, row.names=FALSE, col.names=TRUE, sep=" ")
' > /users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/scripts/$list.R
done  < /users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/scripts/sumstats_list

cat /users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/scripts/*.sumstats.gz.R > /users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/scripts/sumstats_list.R

R --file=/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/scripts/sumstats_list.R

rm /users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/scripts/*.sumstat.gz.R

echo '/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/CBTRESPZ_ACMETA.munged.sumstat.gz.p
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/CBTRESPZ_ADULTMETA.munged.sumstat.gz.p
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/CBTRESPZ_GXTC.munged.sumstat.gz.p
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/CBTRESPZ_GXTE.munged.sumstat.gz.p
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/CBTRESPZ_KI.munged.sumstat.gz.p
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/CBTRESPZ_WB.munged.sumstat.gz.p' \
> /users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/scripts/base_sign_test_list

cd /users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/
WB=CBTRESPZ_WB.munged.sumstat.gz.p
GXTE=CBTRESPZ_GXTE.munged.sumstat.gz.p
KI=CBTRESPZ_KI.munged.sumstat.gz.p

for i in $WB $GXTE $KI
do
awk '{print $1,$2,$3,$4,$5,$6,$7*-1}' $i > $i.tmp && mv $i.tmp $i
done
# CHANGE HEADER TO BETA


rm /users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/scripts/internal.sign.test.sh
cat <<'EOT'>> /users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/scripts/internal.sign.test.sh
## Summary: Run sign test script
## Author: KLP
## Date: 08/18
#!/bin/sh
#$-S /bin/sh
#$ -V
#$ -o /users/k1620287/brc_scratch/meta.analysis/analysis/sign_test
#$ -pe smp 8
#$ -l h_vmem=4G
#$ -l h_rt=72:00:00
#$ -m be

export MKL_NUM_THREADS=8
export NUMEXPR_NUM_THREADS=8
export OMP_NUM_THREADS=8

module add general/R

output_path=/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test
out=$(basename ${1} | sed -e 's/.munged.sumstat.gz.p//g' -e 's/CBTRESPZ_//g')
mkdir ${output_path}/internal/
mkdir ${output_path}/internal/${out}

sh /mnt/lustre/groups/ukbiobank/usr/KLP/GAD7/UKBB/scripts/SignTest.sh \
-b ${1} \
-t "/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/CBTRESPZ_ACMETA.munged.sumstat.gz.p \
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/CBTRESPZ_ADULTMETA.munged.sumstat.gz.p \
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/CBTRESPZ_GXTC.munged.sumstat.gz.p \
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/CBTRESPZ_GXTE.munged.sumstat.gz.p \
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/CBTRESPZ_KI.munged.sumstat.gz.p \
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/CBTRESPZ_WB.munged.sumstat.gz.p" \
-p "0.0005 0.05 0.5" \
-k 1000 \
-r 0.5 \
-o ${output_path}/internal/${out}/${out}.sign.test \
-v

EOT

while read base_sign_test_list
do qsub /users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/scripts/internal.sign.test.sh $base_sign_test_list
done < /users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/scripts/base_sign_test_list


for i in CBTRESPZ_ACMETA CBTRESPZ_ADULTMETA CBTRESPZ_GXTC CBTRESPZ_GXTE CBTRESPZ_KI CBTRESPZ_WB
do
sed -i  -e 's/.munged.sumstats.gz.p//g' -e 's/pthreshold,Base_sample,Target_sample,NA,NA,NA,NA//g' -e 's/pthreshold,Base_sample,Target_sample,Total_SNPs,Total_Shared,Proportion,Binomial_test//g' -e 's/.munged.sumstat.gz.p//g' -e 's/_noMHC.sumstats.gz.p.effect//g' -e 's/CBTRESPZ_//g' -e 's/.effect//g' $i/$i.sign.test_results.csv
done

rm NEWRESULTS
cat *.csv > NEWRESULTS

sed -i '/NA/d' NEWRESULTS
sed -i '/,,/d' NEWRESULTS
sed -i '/^\s*$/d' NEWRESULTS
sed -i '/1,GXTE,-/d' NEWRESULTS


# EXTERNAL

/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/ADHD05_noMHC.sumstats.gz.p \
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/ALCD03_noMHC.sumstats.gz.p \
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/ALCO01_noMHC.sumstats.gz.p \
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/ANOR02_noMHC.sumstats.gz.p \
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/ANXI03_noMHC.sumstats.gz.p \
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/AUTI05_noMHC.sumstats.gz.p \
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/BIPO01_noMHC.sumstats.gz.p \
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/BODY08_noMHC.sumstats.gz.p \
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/BORD01_noMHC.sumstats.gz.p \
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/CANU01_noMHC.sumstats.gz.p \
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/CONS01_noMHC.sumstats.gz.p \
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/CROS01_noMHC.sumstats.gz.p \
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/DEPR02_noMHC.sumstats.gz.p \
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/DEPR05_noMHC.sumstats.gz.p \
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/EDUC03_noMHC.sumstats.gz.p \
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/EXTR01_noMHC.sumstats.gz.p \
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/INCO03_noMHC.sumstats.gz.p \
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/INSO01_noMHC.sumstats.gz.p \
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/INTE01_noMHC.sumstats.gz.p \
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/INTE03_noMHC.sumstats.gz.p \
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/MIGR01_noMHC.sumstats.gz.p \
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/NEUR02B_noMHC.sumstats.gz.p \
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/OPEN01_noMHC.sumstats.gz.p \
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/SCHI02_noMHC.sumstats.gz.p \
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/SUBJ01_noMHC.sumstats.gz.p \
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/TIRE01_noMHC.sumstats.gz.p \
/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/sumstats/for_test/TRSEEK.munged.sumstats.gz.p" \











dat <- read.table("NEWRESULTS",head=F, sep=",", na.strings=c(""," ","NA"))

colnames(dat) <- c('pthreshold','Base_sample','Target_sample','Total_SNPS','Total_Shared','Proportion','Bionomial_test')

library(dplyr)
dat <- dat %>%
  mutate(
    pthreshold = as.factor(pthreshold),
    Base_sample = as.factor(Base_sample),
    Target_sample = as.factor(Target_sample),
    Total_SNPS = as.integer(Total_SNPS),
    Total_Shared = as.integer(Total_Shared),
    Proportion = as.numeric(Proportion),
    Bionomial_test = as.numeric(Bionomial_test)
  )

# dat <- dat[which(!is.na(dat$Target_sample)),]

dat$Total_Unshared <- dat$Total_SNPS - dat$Total_Shared

datp <- subset(dat,select=c('Target_sample','Base_sample','Proportion','pthreshold','Total_Shared','Total_Unshared','Bionomial_test'))

datp$sig <- ifelse(datp$Bionomial_test < 0.05, "*", " ")

dats <- datp[which(datp$Bionomial_test < 0.05),]

S.EXT <- dats[which(dats$Target_sample!="KI" & dats$Target_sample!="WB" & dats$Target_sample!="GXTE" & dats$Target_sample!="GXTC" & dats$Target_sample!="ADULTMETA" & dats$Target_sample!="datsMETA"),]
S.INT <- dats[which(dats$Target_sample=="KI" | dats$Target_sample=="WB" | dats$Target_sample=="GXTE" | dats$Target_sample=="GXTC" | dats$Target_sample=="ADULTMETA" | dats$Target_sample=="datsMETA"),]

write.table(S.INT, "/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/singificant.sign.tests.internal.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep=" ")
write.table(S.EXT, "/users/k1620287/brc_scratch/meta.analysis/analysis/sign_test/singificant.sign.tests.external.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep=" ")

AC <- datp[which(datp$Target_sample!="ACMETA" & datp$Base_sample=="ACMETA"),]
ACEXT <- AC[which(AC$Target_sample!="KI" & AC$Target_sample!="WB" & AC$Target_sample!="GXTE" & AC$Target_sample!="GXTC" & AC$Target_sample!="ADULTMETA" & AC$Target_sample!="ACMETA"),]
ACINT <- AC[which(AC$Target_sample=="KI" | AC$Target_sample=="WB" | AC$Target_sample=="GXTE" | AC$Target_sample=="GXTC" | AC$Target_sample=="ADULTMETA" | AC$Target_sample=="ACMETA"),]

AD <- datp[which(datp$Target_sample!="ADULTMETA" & datp$Base_sample=="ADULTMETA"),]
ADEXT <- AD[which(AD$Target_sample!="KI" & AD$Target_sample!="WB" & AD$Target_sample!="GXTE" & AD$Target_sample!="GXTC" & AD$Target_sample!="ADULTMETA" & AD$Target_sample!="ACMETA"),]
ADINT <- AD[which(AD$Target_sample=="KI" | AD$Target_sample=="WB" | AD$Target_sample=="GXTE" | AD$Target_sample=="GXTC" | AD$Target_sample=="ADULTMETA" | AD$Target_sample=="ACMETA"),]

GXTE <- datp[which(datp$Target_sample!="GXTE" & datp$Base_sample=="GXTE"),]
GXTEEXT <- GXTE[which(GXTE$Target_sample!="KI" & GXTE$Target_sample!="WB" & GXTE$Target_sample!="GXTE" & GXTE$Target_sample!="GXTC" & GXTE$Target_sample!="ADULTMETA" & GXTE$Target_sample!="ACMETA"),]
GXTEINT <- GXTE[which(GXTE$Target_sample=="KI" | GXTE$Target_sample=="WB" | GXTE$Target_sample=="GXTE" | GXTE$Target_sample=="GXTC" | GXTE$Target_sample=="ADULTMETA" | GXTE$Target_sample=="ACMETA"),]

GXTC <- datp[which(datp$Target_sample!="GXTC" & datp$Base_sample=="GXTC"),]
GXTCEXT <- GXTC[which(GXTC$Target_sample!="KI" & GXTC$Target_sample!="WB" & GXTC$Target_sample!="GXTC" & GXTC$Target_sample!="GXTE" & GXTC$Target_sample!="ADULTMETA" & GXTC$Target_sample!="ACMETA"),]
GXTCINT <- GXTC[which(GXTC$Target_sample=="KI" | GXTC$Target_sample=="WB" | GXTC$Target_sample=="GXTC" | GXTC$Target_sample=="GXTE" | GXTC$Target_sample=="ADULTMETA" | GXTC$Target_sample=="ACMETA"),]

KI <- datp[which(datp$Target_sample!="KI" & datp$Base_sample=="KI"),]
KIEXT <- KI[which(KI$Target_sample!="KI" & KI$Target_sample!="WB" & KI$Target_sample!="GXTE" & KI$Target_sample!="GXTC" & KI$Target_sample!="ADULTMETA" & KI$Target_sample!="ACMETA"),]
KIINT <- KI[which(KI$Target_sample=="KI" | KI$Target_sample=="WB" | KI$Target_sample=="GXTE" | KI$Target_sample=="GXTC" | KI$Target_sample=="ADULTMETA" | KI$Target_sample=="ACMETA"),]

WB <- datp[which(datp$Target_sample!="WB" & datp$Base_sample=="WB"),]
WBEXT <- WB[which(WB$Target_sample!="KI" & WB$Target_sample!="WB" & WB$Target_sample!="GXTE" & WB$Target_sample!="GXTC" & WB$Target_sample!="ADULTMETA" & WB$Target_sample!="ACMETA"),]
WBINT <- WB[which(WB$Target_sample=="KI" | WB$Target_sample=="WB" | WB$Target_sample=="GXTE" | WB$Target_sample=="GXTC" | WB$Target_sample=="ADULTMETA" | WB$Target_sample=="ACMETA"),]

dat1	<-	ACEXT
dat2	<-	ADEXT
dat3	<-	GXTCEXT
dat4	<-	GXTEEXT
dat5	<-	KIEXT
dat6	<-	WBEXT
dat7	<-	ACINT
dat8	<-	ADINT
dat9	<-	GXTEINT
dat10	<-	GXTCINT
dat11	<-	KIINT
dat12	<-	WBINT

dat1 <- data.table::melt(dat1,id.vars=c('Target_sample','Base_sample','pthreshold','Proportion', 'sig'))
dat2 <- data.table::melt(dat2,id.vars=c('Target_sample','Base_sample','pthreshold','Proportion', 'sig'))
dat3 <- data.table::melt(dat3,id.vars=c('Target_sample','Base_sample','pthreshold','Proportion', 'sig'))
dat4 <- data.table::melt(dat4,id.vars=c('Target_sample','Base_sample','pthreshold','Proportion', 'sig'))
dat5 <- data.table::melt(dat5,id.vars=c('Target_sample','Base_sample','pthreshold','Proportion', 'sig'))
dat6 <- data.table::melt(dat6,id.vars=c('Target_sample','Base_sample','pthreshold','Proportion', 'sig'))
dat7 <- data.table::melt(dat7,id.vars=c('Target_sample','Base_sample','pthreshold','Proportion', 'sig'))
dat8 <- data.table::melt(dat8,id.vars=c('Target_sample','Base_sample','pthreshold','Proportion', 'sig'))
dat9 <- data.table::melt(dat9,id.vars=c('Target_sample','Base_sample','pthreshold','Proportion', 'sig'))
dat10 <- data.table::melt(dat10,id.vars=c('Target_sample','Base_sample','pthreshold','Proportion', 'sig'))
dat11 <- data.table::melt(dat11,id.vars=c('Target_sample','Base_sample','pthreshold','Proportion', 'sig'))
dat12 <- data.table::melt(dat12,id.vars=c('Target_sample','Base_sample','pthreshold','Proportion', 'sig'))

dat1$variable <- factor(dat1$variable,levels=c('Total_Unshared','Total_Shared'))
dat2$variable <- factor(dat2$variable,levels=c('Total_Unshared','Total_Shared'))
dat3$variable <- factor(dat3$variable,levels=c('Total_Unshared','Total_Shared'))
dat4$variable <- factor(dat4$variable,levels=c('Total_Unshared','Total_Shared'))
dat5$variable <- factor(dat5$variable,levels=c('Total_Unshared','Total_Shared'))
dat6$variable <- factor(dat6$variable,levels=c('Total_Unshared','Total_Shared'))
dat7$variable <- factor(dat7$variable,levels=c('Total_Unshared','Total_Shared'))
dat8$variable <- factor(dat8$variable,levels=c('Total_Unshared','Total_Shared'))
dat9$variable <- factor(dat9$variable,levels=c('Total_Unshared','Total_Shared'))
dat10$variable <- factor(dat10$variable,levels=c('Total_Unshared','Total_Shared'))
dat11$variable <- factor(dat11$variable,levels=c('Total_Unshared','Total_Shared'))
dat12$variable <- factor(dat12$variable,levels=c('Total_Unshared','Total_Shared'))


library(ggplot2)

xlabels <- c('5.10-6','5.10-4','5.10-2','5.10-1','1')

plot <- ggplot(dat8)  +
		geom_bar(aes(x=as.factor(pthreshold),y=value,fill=variable),
		stat='identity',position='fill') +
		facet_grid(Base_sample ~ Target_sample) +
		scale_x_discrete(labels=xlabels) +
		theme(axis.text.x = element_text(size = 8, angle=45, hjust=1)) +
		scale_y_continuous(expand = c(0, 0)) +
		xlab("p-threshold") +
		ylab("Proportion of total SNPs shared") +
		theme(legend.position='none',
		panel.border = element_blank(),
		panel.grid = element_blank(),
		strip.text = element_text(face='bold',size=8)) +
		scale_fill_manual(values = c("azure3","#636363"))
show(plot)

pdf(‘path/and/filename.pdf’,width=10)
plot
dev.off()
