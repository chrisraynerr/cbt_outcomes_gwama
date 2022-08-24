####################################
# Meta-analysis of response to CBT #
####################################

# Script 8. Post Imputation Quality Control #

# Remove genotyped: snps info (<1); Info threshold (<.75); all indels
# Include only biallelic variants

cd /mnt/lustre/users/k1620287/meta.analysis/imputation
rm hrc.info.indels.sh
cat <<'EOT'>> hrc.info.indels.sh
#!/bin/sh
#$ -V
#$ -m ea
#$ -pe smp 10
#$ -l h_vmem=19G
#$ -l s_vmem=19G
#$ -l h_rt=04:00:00
#$ -S /bin/sh
#$ -m be
#$ -t 1-23
#$ -tc 10
#$ -o /users/k1620287/brc_scratch/meta.analysis/imputation/logfiles
#$ -e /users/k1620287/brc_scratch/meta.analysis/imputation/logfiles

module load bioinformatics/bcftools/1.3.1
module load bioinformatics/htslib/1.3.1

for site in $*
do
# cd ~/brc_scratch/meta.analysis/imputation/$site
cd /mnt/lustre/groups/ukbiobank/usr/crayner/
bcftools view -e 'INFO/INFO<0.75' -m2 -M2 -v snps ${SGE_TASK_ID}.vcf.gz | bgzip -c > ${SGE_TASK_ID}.info75.indels.multi.gz
done
EOT

qsub hrc.info.indels.sh ki_mdd


lfs quota -h -u k1620287 /mnt/lustre/
lfs quota -gh k1620287 /users/k1620287/brc_scratch/

##########################################
# NEED TO REPLACE '.' WITH: chr:bp_A1_A2 #
##########################################

cd /mnt/lustre/users/k1620287/meta.analysis/imputation
rm replace.dot.sh
cat <<'EOT'>> replace.dot.sh
#!/bin/sh
#$ -V
#$ -m ea
#$ -pe smp 10
#$ -l h_vmem=4G
#$ -l s_vmem=4G
#$ -l h_rt=04:00:00
#$ -S /bin/sh
#$ -m be
#$ -t 1-23
#$ -tc 10
#$ -o /users/k1620287/brc_scratch/meta.analysis/imputation/logfiles/
#$ -e /users/k1620287/brc_scratch/meta.analysis/imputation/logfiles/
module load bioinformatics/htslib/1.3.1

# for site in $*
# do
# cd ~/brc_scratch/meta.analysis/imputation/$site
cd /mnt/lustre/groups/ukbiobank/usr/crayner/
zcat ${SGE_TASK_ID}.info75.indels.multi.gz | awk '$3=="." {$3= $1 ":" $2"_"$4"_"$5} {print}' | tr ' ' '\t' | bgzip -c > ${SGE_TASK_ID}.info75.indels.multi.dot.gz
done
EOT

qsub replace.dot.sh
# qsub replace_dot.sh 'wb1' 'wb2' 'ki' 'gxte' 'gxtc'

##########################
# QUALITY CONTROL CHECKS #
##########################
# Generate info file
cd /mnt/lustre/users/k1620287/meta.analysis/imputation
rm hrc.post.qc.sh
cat <<'EOT'>> hrc.post.qc.sh
#!/bin/sh
#$ -V
#$ -m ea
#$ -pe smp 10
#$ -l h_vmem=3G
#$ -l s_vmem=3G
#$ -l h_rt=08:00:00
#$ -S /bin/sh
#$ -m be
#$ -t 1-23
#$ -tc 10
#$ -o /users/k1620287/brc_scratch/meta.analysis/imputation/logfiles
#$ -e /users/k1620287/brc_scratch/meta.analysis/imputation/logfiles

module load bioinformatics/bcftools/1.3.1
module load bioinformatics/htslib/1.3.1

# for site in $*
# do
# cd ~/brc_scratch/meta.analysis/imputation/$site
cd /mnt/lustre/groups/ukbiobank/usr/crayner/

	bcftools query -f '%CHROM %POS %ID %REF %ALT %RefPanelAF %AN %AC %INFO/INFO\n' ${SGE_TASK_ID}.info75.indels.multi.dot.gz \
	| awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$8/$7}' | bgzip -c > ${SGE_TASK_ID}.info75.indels.multi.dot.INFO.gz;
done
EOT
qsub hrc.post.qc.sh

-hold_jid replace.dot.sh hrc.post.qc.sh
qsub hrc.post.qc.sh wb1 wb2 ki gxte gxtc

############################################
# Write script to add header to INFO files #
############################################

cd /mnt/lustre/users/k1620287/meta.analysis/imputation
for site in 'wb1' 'wb2' 'ki' 'gxte' 'gxtc'
do
for chr in {1..23};
do
echo 'library(data.table)
setwd("/users/k1620287/brc_scratch/meta.analysis/imputation/'$site'")
file <- fread("zcat '$chr'.info75.indels.multi.dot.INFO.gz", data.table=F)
colnames(file) <- c("CHR","BP","SNP","A0","A1","HRC_AF_A1","AN","AC","Info","AF_A1")
genot <- read.table("/users/k1620287/brc_scratch/meta.analysis/imputation/'$site'/'$site'.Hardy_Imiss",h=T)[,c(1,2,14)]
m <- merge(file,genot,by=c("CHR","SNP"),all.x=T)
m$genotyped[is.na(m$genotyped)] <- 0
write.table(m,gzfile("/users/k1620287/brc_scratch/meta.analysis/imputation/'$site'/'$chr'.info75.indels.multi.dot.INFOh.gz"),col.names=T,row.names=F,quote=F,sep=" ")
print("Done")
q(save="no")' > /mnt/lustre/users/k1620287/meta.analysis/imputation/$site/$chr.$site.header.R; done
done

cd /mnt/lustre/users/k1620287/meta.analysis/imputation
for site in 'wb1' 'wb2' 'ki' 'gxte' 'gxtc'
do
cd ~/brc_scratch/meta.analysis/imputation/$site
for chr in {1..23};
do
echo $chr.$site.header.R; done > /mnt/lustre/users/k1620287/meta.analysis/imputation/$site/$site.all.header.R
done


####################################

for chr in {1..23};
do
echo 'library(data.table)
setwd("/mnt/lustre/groups/ukbiobank/usr/crayner")
file <- fread("zcat '$chr'.info75.indels.multi.dot.INFO.gz", data.table=F)
colnames(file) <- c("CHR","BP","SNP","A0","A1","HRC_AF_A1","AN","AC","Info","AF_A1")
genot <- read.table("ki_mdd.Hardy_Imiss",h=T)[,c(1,2,14)]
m <- merge(file,genot,by=c("CHR","SNP"),all.x=T)
m$genotyped[is.na(m$genotyped)] <- 0
write.table(m,gzfile("/mnt/lustre/groups/ukbiobank/usr/crayner/'$chr'.info75.indels.multi.dot.INFOh.gz"),col.names=T,row.names=F,quote=F,sep=" ")
print("Done")
q(save="no")' > /mnt/lustre/groups/ukbiobank/usr/crayner/$chr.header.R; done
done

for chr in {1..23};
do
echo $chr.header.R; done > ./all.header.R
done


##########################################
# Run script to add header to INFO files #
##########################################

cd /mnt/lustre/users/k1620287/meta.analysis/imputation
rm add.header.hrc.sh
cat <<'EOT'>> add.header.hrc.sh
#!/bin/sh
#$ -V
#$ -m ea
#$ -pe smp 10
#$ -l h_vmem=3G
#$ -l s_vmem=3G
#$ -l h_rt=05:00:00
#$ -S /bin/sh
#$ -m be
#$ -t 1-23
#$ -tc 10
#$ -o /users/k1620287/brc_scratch/meta.analysis/imputation/logfiles
#$ -e /users/k1620287/brc_scratch/meta.analysis/imputation/logfiles

module load general/R

# for site in $*
# do
# cd ~/brc_scratch/meta.analysis/imputation/$site
cd /mnt/lustre/groups/ukbiobank/usr/crayner/

#script=$(awk "NR==$SGE_TASK_ID" /mnt/lustre/users/k1620287/meta.analysis/imputation/$site/$site.all.header.R)
script=$(awk "NR==$SGE_TASK_ID" ./all.header.R)

R CMD BATCH ${script};

done
EOT
qsub -N add.header.hrc.sh -hold_jid hrc.post.qc.sh add.header.hrc.sh


# Remove old files #
for site in 'wb1' 'wb2' 'ki' 'gxte' 'gxtc'
do
	rm $site/*info.gz
	rm $site/*.R $site/*.Rout
done

####################################
# Check bad snps have been removed #
####################################
cd /mnt/lustre/users/k1620287/meta.analysis/imputation
for site in 'wb1' 'wb2' 'ki' 'gxte' 'gxtc'
do
for chr in {1..23};
do
echo '
library(data.table)
library(ggplot2)
setwd("/users/k1620287/brc_scratch/meta.analysis/imputation/'$site'")
file <- fread("zcat '$chr'.info75.indels.multi.dot.INFOh.gz",data.table=F)

png("chr'$chr'.'$site'.HRC_AF.AF_2.png")
plot(file$HRC_AF_A1,file$AF_A1)
dev.off()

png("chr'$chr'.'$site'.InfoDistr_2.png")
hist(file$Info)
dev.off()

scat <- ggplot(file,aes(x=Info,y=AF_A1)) +
geom_point(aes(colour=factor(genotyped)))+
facet_wrap("genotyped")
ggsave("chr.'$chr'.'$site'.Info.AF_2.png",scat)

print("Done")
q(save="no")' > /mnt/lustre/users/k1620287/meta.analysis/imputation/$site/$site.$chr.post.plots.R; done
done

# WRITE SCRIPT TO INCLUDE ALL CHR PLOTS FOR EACH SITE #

cd /mnt/lustre/users/k1620287/meta.analysis/imputation
for site in 'wb1' 'wb2' 'ki' 'gxte' 'gxtc'
do
cd /mnt/lustre/users/k1620287/meta.analysis/imputation/$site
for chr in {1..23};
do
echo ${site}.${chr}.post.plots.R; done > /mnt/lustre/users/k1620287/meta.analysis/imputation/$site/$site.imputation.post.plots.all.chr
done


# Generate Plots  #

cd /mnt/lustre/users/k1620287/meta.analysis/imputation
rm info.plots.hrc.sh
cat <<'EOT'>> info.plots.hrc.sh
#!/bin/sh
#$ -V
#$ -m ea
#$ -pe smp 10
#$ -l h_vmem=3G
#$ -l s_vmem=3G
#$ -l h_rt=04:00:00
#$ -S /bin/sh
#$ -m be
#$ -t 1-23
#$ -tc 10
#$ -o /users/k1620287/brc_scratch/meta.analysis/imputation/logfiles/output
#$ -e /users/k1620287/brc_scratch/meta.analysis/imputation/logfiles/error
module load bioinformatics/R/3.3.3
for site in gxtc
do
cd ~/brc_scratch/meta.analysis/imputation/$site
script=$(awk "NR==$SGE_TASK_ID" /mnt/lustre/users/k1620287/meta.analysis/imputation/$site/$site.imputation.post.plots.all.chr)
R CMD BATCH ${script};
done
EOT
qsub -N info.plots.hrc.sh -hold_jid add.header.hrc.sh info.plots.hrc.sh
qsub info.plots.hrc.sh wb1 wb2 ki gxte gxtc

# Remove old files and move figures into folder #
for site in wb1 wb2 ki gxte gxtc
do
	cd ~/brc_scratch/meta.analysis/imputation/
	rm $site/*.R $site/*.Rout
	mkdir $site/figures
	mv $site/*.png $site/figures
done


**********************************************************************************
####################################
# Create INFO file for all gw snps #
####################################
for site in wb1 wb2 ki gxte gxtc
do
	cd ~/brc_scratch/meta.analysis/imputation/$site
# copy header
zcat 1.info75.indels.multi.dot.INFOh.gz > 1.info75
cat 1.info75 | head -1 > info.header
# remove header for
for chr in {1..23};
do
zcat $chr.info75.indels.multi.dot.INFOh.gz | tail -n +2 > $chr.info75
done
# combine chr snps
cat info.header *.info75 > $site.info75.snps
cat *.keepSNPs > qc.snps
awk -v FS="[ =]" 'NR==FNR{rows[$1]++;next}(substr($NF,1,length($NF)-1) in rows)' gxte_merged.qc.snps gxte_merged.info75.snps
grep -F -f gxte_merged.qc.snps gxte_merged.info75.snps > gxte_merged.info75.qc.snps
grep -F -f file1 file2

awk '{print $2,$9}' gxte_merged.info75.snps | column -t > gxte_merged.info75.snps.less


rm info.header
rm *.info75
rm *.info.AF.gz
rm *.info75.indels.multi.dot.INFOh.gz


#############################
#  Convert to Plink format  #
#############################

# FOR SOME REASON WORKING ON A FEW FILES AND GENERATING TEMPORARY FILES
# TRIED RUNNING ONE CHROMOSOME AT A TIME BUT STILL GENERATES TEMP FILES
# USE ECHO VERSION BELOW - SLOWER BUT SEEMS TO WORK

cd /mnt/lustre/users/k1620287/meta.analysis/imputation
rm post.qc.convert.plink.sh
cat <<'EOT'>> post.qc.convert.plink.sh
#!/bin/sh
#$ -V
#$ -m ea
#$ -pe smp 10
#$ -l h_vmem=19G
#$ -l s_vmem=19G
#$ -l h_rt=04:00:00
#$ -S /bin/sh
#$ -m be
#$ -t 1-23
#$ -tc 10
#$ -o /mnt/lustre/groups/ukbiobank/usr/crayner
#$ -e /mnt/lustre/groups/ukbiobank/usr/crayner

################ #$ -o /users/k1620287/brc_scratch/meta.analysis/imputation/logfiles/output
################ #$ -e /users/k1620287/brc_scratch/meta.analysis/imputation/logfiles/error

module load bioinformatics/plink/1.90b3.31

# for site in gxtc
# do
# cd ~/brc_scratch/meta.analysis/imputation/$site

cd /mnt/lustre/groups/ukbiobank/usr/crayner

plink \
--vcf ${SGE_TASK_ID}.info75.indels.multi.dot.gz \
--geno 0.02 \
--maf 0.01 \
--make-bed \
--id-delim '|' \
--out ${SGE_TASK_ID}.callr02.maf01

# done

EOT

#qsub -N post.qc.convert.plink.sh -hold_jid replace.dot.sh post.qc.convert.plink.sh

qsub post.qc.convert.plink.sh


##########################
#  Drop duplicate SNPs   #
##########################

# Cat version doesn't work :s
cd /mnt/lustre/users/k1620287/meta.analysis/imputation
rm find.duplicates.sh
cat <<'EOT'>> find.duplicates.sh
#!/bin/sh
#$ -V
#$ -m ea
#$ -pe smp 10
#$ -l h_vmem=4G
#$ -l s_vmem=4G
#$ -l h_rt=04:00:00
#$ -S /bin/sh
#$ -m be
#$ -t 1-23
#$ -tc 10
#$ -o /mnt/lustre/groups/ukbiobank/usr/crayner
#$ -e /mnt/lustre/groups/ukbiobank/usr/crayner

################ #$ -o /users/k1620287/brc_scratch/meta.analysis/imputation/logfiles/output
################ #$ -e /users/k1620287/brc_scratch/meta.analysis/imputation/logfiles/error

module load bioinformatics/plink/1.90b3.31

# for site in gxtc
# do
# cd ~/brc_scratch/meta.analysis/imputation/$site

cd /mnt/lustre/groups/ukbiobank/usr/crayner

awk '{print $2}' ${SGE_TASK_ID}.callr02.maf01.bim | sort | uniq -d > ${SGE_TASK_ID}.duplicates.to.remove

#done

EOT

qsub find.duplicates.sh

rm drop.duplicates.sh
cat <<'EOT'>> drop.duplicates.sh
#!/bin/sh
#$ -V
#$ -m ea
#$ -pe smp 10
#$ -l h_vmem=9G
#$ -l s_vmem=9G
#$ -l h_rt=04:00:00
#$ -S /bin/sh
#$ -m be
#$ -t 1-23
#$ -tc 10
#$ -o /mnt/lustre/groups/ukbiobank/usr/crayner
#$ -e /mnt/lustre/groups/ukbiobank/usr/crayner

################ #$ -o /users/k1620287/brc_scratch/meta.analysis/imputation/logfiles/output
################ #$ -e /users/k1620287/brc_scratch/meta.analysis/imputation/logfiles/error

module load bioinformatics/plink/1.90b3.31

# for site in gxtc
# do
# cd ~/brc_scratch/meta.analysis/imputation/$site

cd /mnt/lustre/groups/ukbiobank/usr/crayner

# awk '{print $2}' ${SGE_TASK_ID}.callr02.maf01.bim | sort | uniq -d > ${SGE_TASK_ID}.duplicates.to.remove

plink \
--bfile ${SGE_TASK_ID}.callr02.maf01 \
--exclude ${SGE_TASK_ID}.duplicates.to.remove \
--make-bed \
--out ${SGE_TASK_ID}.callr02.maf01.1

#done

EOT

qsub drop.duplicates.sh

#qsub -N drop_duplicates.sh -hold_jid post.qc.convert.plink.echo.sh drop.duplicates.sh

###############################
#  Drop duplicate positions   #
###############################

cat <<'EOT'>> find.pos.duplicates.sh
#!/bin/sh
#$ -V
#$ -m ea
#$ -pe smp 10
#$ -l h_vmem=4G
#$ -l s_vmem=4G
#$ -l h_rt=04:00:00
#$ -S /bin/sh
#$ -m be
#$ -t 1-23
#$ -tc 10
#$ -o /mnt/lustre/groups/ukbiobank/usr/crayner
#$ -e /mnt/lustre/groups/ukbiobank/usr/crayner

################ #$ -o /users/k1620287/brc_scratch/meta.analysis/imputation/logfiles/output
################ #$ -e /users/k1620287/brc_scratch/meta.analysis/imputation/logfiles/error

module load bioinformatics/plink/1.90b3.31

# for site in gxtc
# do
# cd ~/brc_scratch/meta.analysis/imputation/$site

cd /mnt/lustre/groups/ukbiobank/usr/crayner

sort -uk4,4 ${SGE_TASK_ID}.callr02.maf01.1.bim | awk '{print $2}' > ${SGE_TASK_ID}.keepSNPs

EOT

qsub find.pos.duplicates.sh

cd /mnt/lustre/users/k1620287/meta.analysis/imputation
rm pos.duplicates.sh
cat <<'EOT'>> pos.duplicates.sh
#!/bin/sh
#$ -V
#$ -m ea
#$ -pe smp 10
#$ -l h_vmem=9G
#$ -l s_vmem=9G
#$ -l h_rt=04:00:00
#$ -S /bin/sh
#$ -m be
#$ -t 1-23
#$ -tc 10
#$ -o /mnt/lustre/groups/ukbiobank/usr/crayner
#$ -e /mnt/lustre/groups/ukbiobank/usr/crayner

################ #$ -o /users/k1620287/brc_scratch/meta.analysis/imputation/logfiles/output
################ #$ -e /users/k1620287/brc_scratch/meta.analysis/imputation/logfiles/error

module load bioinformatics/plink/1.90b3.31

# for site in gxtc
# do
# cd ~/brc_scratch/meta.analysis/imputation/$site

cd /mnt/lustre/groups/ukbiobank/usr/crayner

# sort -uk4,4 ${SGE_TASK_ID}.callr02.maf01.1.bim | awk '{print $2}' > ${SGE_TASK_ID}.keepSNPs

plink \
--bfile ${SGE_TASK_ID}.callr02.maf01.1 \
--extract ${SGE_TASK_ID}.keepSNPs \
--make-bed \
--out ${SGE_TASK_ID}.callr02.maf01.2

# done

EOT

qsub pos.duplicates.sh

#qsub -N pos.duplicates.sh -hold_jid drop.duplicates.sh pos.duplicates.sh


####################################
#  MERGE CHR INTO WHOLE DATA SET   #
####################################
# merge to  chr1 genotypes #

# cd ~/brc_scratch/meta.analysis/imputation
# for site in gxtc
# do
# cd ~/brc_scratch/meta.analysis/imputation/$site
for chr in {2..23};
do
echo ''$chr'.callr02.maf01.2.bed '$chr'.callr02.maf01.2.bim '$chr'.callr02.maf01.2.fam'; done > merge.list
#done

cd ~/brc_scratch/meta.analysis/imputation
rm merge.chrs.sh
cat <<'EOT'>> merge.chrs.sh
#!/bin/sh
#$ -V
#$ -m ea
#$ -l h_vmem=9G
#$ -l s_vmem=9G
#$ -l h_rt=04:00:00
#$ -S /bin/sh
#$ -m be
#$ -o /mnt/lustre/groups/ukbiobank/usr/crayner
#$ -e /mnt/lustre/groups/ukbiobank/usr/crayner

################ #$ -o /users/k1620287/brc_scratch/meta.analysis/imputation/logfiles/output
################ #$ -e /users/k1620287/brc_scratch/meta.analysis/imputation/logfiles/error

module load bioinformatics/plink/1.90b3.31

for site in $*
do
# cd ~/brc_scratch/meta.analysis/imputation/$site

cd /mnt/lustre/groups/ukbiobank/usr/crayner

plink \
--bfile 1.callr02.maf01.2 \
--merge-list merge.list \
--make-bed \
--out $site.postimputation

done

EOT

qsub merge.chrs.sh ki_mdd


**********************************************************************
rm *callr02.maf01.*
rm *indels.multi.gz
rm *indels.multi.dot.INFO.gz
rm *duplicates.to.remove
rm *keepSNPs
rm *.sh.*

**********************************************************************

###########################################
# Quality control - MAF 0.05 - HWE 0.0005 #
###########################################

cd /mnt/lustre/users/k1620287/meta.analysis/imputation
rm post.imputation.plink.qc.sh
cat <<'EOT'>> post.imputation.plink.qc.sh
#!/bin/sh
#$ -V
#$ -m ea
#$ -l h_vmem=9G
#$ -l s_vmem=9G
#$ -l h_rt=04:00:00
#$ -S /bin/sh
#$ -m be
#$ -o /users/k1620287/brc_scratch/meta.analysis/imputation/logfiles/
#$ -e /users/k1620287/brc_scratch/meta.analysis/imputation/logfiles/

module load bioinformatics/plink2/1.90b3.38

for site in $*
do
#cd /mnt/lustre/users/k1620287/meta.analysis/imputation/$site/
cd /mnt/lustre/groups/ukbiobank/usr/crayner
(
plink \
--bfile $site.postimputation \
--maf 0.05 \
--make-bed \
--out $site.postimputation.maf05
)&
done
wait

for site in $*
do
#cd /mnt/lustre/users/k1620287/meta.analysis/imputation/$site/
cd /mnt/lustre/groups/ukbiobank/usr/crayner
(
plink \
--bfile $site.postimputation.maf05 \
--hwe 0.00001 \
--make-bed \
--out $site.postimputation.maf05.hwe
)&
done
wait

for site in $*
do
#cd /mnt/lustre/users/k1620287/meta.analysis/imputation/$site/
cd /mnt/lustre/groups/ukbiobank/usr/crayner
(
plink \
--bfile $site.postimputation.maf05.hwe \
--hardy \
--freq \
--out $site.postimputation.maf05.hwe
)&
done
done

EOT
qsub post.imputation.plink.qc.sh ki_mdd



sed -i -e 's/cas_mdd_icbt_eur_rk_Il5M//g' ki_mdd.postimputation.maf01.hwe.fam
sed -i -e 's/*//g' ki_mdd.postimputation.maf01.hwe.fam

# ******************************************************************************

STEP 1

cat <<'EOT'>> post.imputation.PCA.1.sh
#!/bin/sh
#$ -V
#$ -m ea
#$ -l h_vmem=9G
#$ -l s_vmem=9G
#$ -l h_rt=04:00:00
#$ -S /bin/sh
#$ -m be
#$ -o /users/k1620287/brc_scratch/meta.analysis/imputation/logfiles/
#$ -e /users/k1620287/brc_scratch/meta.analysis/imputation/logfiles/

module load bioinformatics/plink2/1.90b3.38

cd /mnt/lustre/groups/ukbiobank/usr/crayner

#note: high-LD list for build h37 retreived from https://github.com/hou/GWAS/blob/master/GWAS_QC.sh
for site in $*
do
plink \
--bfile $site.postimputation.maf05.hwe \
--exclude range highLD.list \
--indep-pairwise 50 5 0.2 \
--out $site.postimputation.maf05.hwe.rmhLD.pruned
done
EOT

qsub post.imputation.PCA.1.sh ki_mdd

STEP 2
rm post.imputation.PCA.2.sh
cat <<'EOT'>> post.imputation.PCA.2.sh
#!/bin/sh
#$ -V
#$ -m ea
#$ -l h_vmem=9G
#$ -l s_vmem=9G
#$ -l h_rt=04:00:00
#$ -S /bin/sh
#$ -m be
#$ -o /users/k1620287/brc_scratch/meta.analysis/imputation/logfiles/
#$ -e /users/k1620287/brc_scratch/meta.analysis/imputation/logfiles/

cd /mnt/lustre/groups/ukbiobank/usr/crayner

for site in $*
do
/mnt/lustre/groups/ukbiobank/Edinburgh_Data/Software/tools/plink2 \
--bfile $site.postimputation.maf05 \
--extract $site.postimputation.maf05.hwe.rmhLD.pruned.prune.in \
--freq \
--pca var-wts \
--out $site.postimputation.maf05.hwe.rmhLD.pruned.pca_wts
done
EOT

qsub post.imputation.PCA.2.sh ki_mdd






########################
# Merging imputed data #
########################

site1=wb1
site2=wb2
site3=gxte
site4=ki

{
module add bioinformatics/plink2/1.90b3.38

plink \
--bfile /mnt/lustre/users/k1620287/meta.analysis/imputation/$site1/$site1.postimputation.maf05.hwe \
--bmerge /mnt/lustre/users/k1620287/meta.analysis/imputation/$site2/$site2.postimputation.maf05.hwe.bed \
	/mnt/lustre/users/k1620287/meta.analysis/imputation/$site2/$site2.postimputation.maf05.hwe.bim \
	/mnt/lustre/users/k1620287/meta.analysis/imputation/$site2/$site2.postimputation.maf05.hwe.fam \
--make-bed \
--out /mnt/lustre/users/k1620287/meta.analysis/adult_merged/$site1.$site2.postimputation.maf05.hwe.merged

# TO DO:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
^
REPEAT THIS FOR:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#plink --bfile ./$site2/$site2.hrc.imputed.qc \
#--flip ./merge/$site1.$site2.hrc.imputed.qc.merged-merge.missnp \
#--make-bed --out ./$site2/$site2.hrc.imputed.qc.flipped

#plink --bfile $site1/$site1.hrc.imputed.qc \
#--bmerge $site2/$site2.hrc.imputed.qc.flipped.bed $site2/$site2.hrc.imputed.qc.flipped.bim $site2/$site2.hrc.imputed.qc.flipped.fam \
#--make-bed --out ./merge/$site1.$site2.hrc.imputed.qc.merged

plink --bfile ./merge/$site1.$site2.hrc.imputed.qc.merged \
--bmerge $site3/$site3.hrc.imputed.qc.bed $site3/$site3.hrc.imputed.qc.bim $site3/$site3.hrc.imputed.qc.fam \
--make-bed --out ./merge/$site1.$site2.$site3.hrc.imputed.qc.merged

#plink --bfile ./$site3/$site3.hrc.imputed.qc \
#--flip ./merge/$site1.$site2.$site3.hrc.imputed.qc.merged-merge.missnp \
#--make-bed --out ./$site3/$site3.hrc.imputed.qc.flipped

plink \
--bfile ./$site3/$site3.hrc.imputed.qc \
--exclude ./merge/wb1.wb2.gxte.hrc.imputed.qc.merged-merge.missnp \
--make-bed \
--out ./$site3/$site3.hrc.imputed.qc.mis.exc

plink --bfile ./merge/$site1.$site2.hrc.imputed.qc.merged \
--bmerge $site3/$site3.hrc.imputed.qc.mis.exc.bed $site3/$site3.hrc.imputed.qc.mis.exc.bim $site3/$site3.hrc.imputed.qc.mis.exc.fam \
--make-bed --out ./merge/$site1.$site2.$site3.hrc.imputed.qc.merged

#plink --bfile ./merge/$site1.$site2.hrc.imputed.qc.merged \
#--bmerge $site3/$site3.hrc.imputed.qc.flipped.bed $site3/$site3.hrc.imputed.qc.flipped.bim $site3/$site3.hrc.imputed.qc.flipped.fam \
#--make-bed --out ./merge/$site1.$site2.$site3.hrc.imputed.qc.merged

plink --bfile ./merge/$site1.$site2.$site3.hrc.imputed.qc.merged \
--bmerge $site4/$site4.hrc.imputed.qc.bed $site4/$site4.hrc.imputed.qc.bim $site4/$site4.hrc.imputed.qc.fam \
--make-bed --out ./merge/$site1.$site2.$site3.$site4.hrc.imputed.qc.merged

plink \
--bfile ./$site4/$site4.hrc.imputed.qc \
--exclude ./merge/wb1.wb2.gxte.ki.hrc.imputed.qc.merged-merge.missnp \
--make-bed \
--out ./$site4/$site4.hrc.imputed.qc.mis.exc

plink --bfile ./merge/$site1.$site2.$site3.hrc.imputed.qc.merged \
--bmerge $site4/$site4.hrc.imputed.qc.mis.exc.bed $site4/$site4.hrc.imputed.qc.mis.exc.bim $site4/$site4.hrc.imputed.qc.mis.exc.fam \
--make-bed --out ./merge/$site1.$site2.$site3.$site4.hrc.imputed.qc.merged


awk '{print $2}' ./merge/$site1.$site2.$site3.$site4.hrc.imputed.qc.merged.bim | sort | uniq -d > ./merge/imputed.merged.duplicates.to.remove
wc -l ./merge/imputed.merged.duplicates.to.remove
#0 ./merge/merged.duplicates.to.remove

sort -uk4,4 ./merge/$site1.$site2.$site3.$site4.hrc.imputed.qc.merged.bim | awk '{print $2}' > ./merge/imputed.merged.keep.SNPs
wc -l ./merge/imputed.merged.keep.SNPs
#5686961 ./merge/imputed.merged.keep.SNPs

plink \
--bfile ./merge/$site1.$site2.$site3.$site4.hrc.imputed.qc.merged \
--extract ./merge/imputed.merged.keep.SNPs \
--make-bed \
--out ./merge/$site1.$site2.$site3.$site4.hrc.imputed.qc.merged.no.dup

plink \
--bfile ./merge/$site1.$site2.$site3.$site4.hrc.imputed.qc.merged.no.dup \
--maf 0.05 \
--make-bed \
--out ./merge/adult_merged.hrc.imputed.maf.05

sh Iterative_Missingness.sh 90 99 1

}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##########################################################################
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

################################
#### Adult_merged analysis ##### !!!!!!!!!!!!!!!! NEEDS FIXING: TO DO!!!
################################

cd /users/k1620287/brc_scratch/meta.analysis/analysis/adult_merged
plink --bfile ../../imputation/adult_merged/adult_merged.postimputation.maf05.hwe --update-ids update.files/adult_merged.update.ids --make-bed --out binary/adult_merged.postimputation.maf05.hwe.udid
plink --bfile binary/adult_merged.postimputation.maf05.hwe.udid --keep binary/adult_merged.keeps --make-bed --out binary/adult_merged.postimputation.maf05.hwe.udid.keeps
plink --bfile binary/adult_merged.postimputation.maf05.hwe.udid.keeps --update-sex update.files/adult_merged.genders.txt --make-bed --out binary/adult_merged.postimputation.maf05.hwe.udid.keeps.udsex

plink --bfile binary/adult_merged.postimputation.maf05.hwe.udid.keeps.udsex --pheno pheno/adult_merged.pc.post --make-bed --out binary/adult_merged.postimputation.maf05.hwe.updated.2002.pc.post
plink --bfile binary/adult_merged.postimputation.maf05.hwe.udid.keeps.udsex --pheno pheno/adult_merged.pc.fu --make-bed --out binary/adult_merged.postimputation.maf05.hwe.updated.2002.pc.fu
plink --bfile binary/adult_merged.postimputation.maf05.hwe.udid.keeps.udsex --pheno pheno/adult_merged.zscore.post --make-bed --out binary/adult_merged.postimputation.maf05.hwe.updated.2002.zscore.post
plink --bfile binary/adult_merged.postimputation.maf05.hwe.udid.keeps.udsex --pheno pheno/adult_merged.zscore.fu --make-bed --out binary/adult_merged.postimputation.maf05.hwe.updated.2002.zscore.fu

adult_merged = 12 significant PCs
gcta --bfile binary/adult_merged.postimputation.maf05.hwe.updated.2002.pc.fu --autosome --maf 0.01 --thread-num 2 --make-grm --out grm/adult_merged.postimputation.maf05.hwe.updated.2002
gcta --grm grm/adult_merged.postimputation.maf05.hwe.updated.2002 --pca --out cov/pcs/adult_merged.pca
awk '{print $1,$2,$3}' cov/pcs/adult_merged.pca.eigenvec > cov/pcs/adult_merged.one.pc
awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' cov/pcs/adult_merged.pca.eigenvec > cov/adult_merged.ten.pc
cd /users/k1620287/brc_scratch/meta.analysis/analysis/adult_merged/cov
R
ten <- read.table("pcs/adult_merged.ten.pc", head =F)
one <- read.table("pcs/adult_merged.one.pc", head =F)
conc <- read.table("batch/adult_merged.conc.txt", head=F)
names(ten)[1] <- "FID"
tenconc <- merge(ten, conc, by = "FID")
tenconc$IID <- NULL
names(one)[1] <- "FID"
oneconc <- merge(one, conc, by = "FID")
oneconc$IID <- NULL
write.table(oneconc,"final/one.pc.conc.qcov", col.names=F, row.names=F, quote=F, sep=" ")
write.table(tenconc,"final/ten.pc.conc.qcov", col.names=F, row.names=F, quote=F, sep=" ")
q()
n

tail -n +2 batch/adult_merged.type.filt.txt > final/adult_merged.type.filt.cov

covariate files:
/users/k1620287/brc_scratch/meta.analysis/analysis/adult_merged/cov/final
adult_merged.one.pc.conc.qcov
adult_merged.ten.pc.conc.qcov
adult_merged.type.filt.cov


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

###############
#### GXTC #####
###############

cd /users/k1620287/brc_scratch/meta.analysis/analysis/gxtc
plink --bfile ../../imputation/gxtc/gxtc.postimputation.maf05.hwe --update-ids update.files/gxtc.update.ids --make-bed --out binary/gxtc.postimputation.maf05.hwe.udid
plink --bfile binary/gxtc.postimputation.maf05.hwe.udid --keep update.files/gxtc.keeps.txt --make-bed --out binary/gxtc.postimputation.maf05.hwe.udid.keeps
plink --bfile binary/gxtc.postimputation.maf05.hwe.udid.keeps --update-sex update.files/gxtc.genders.txt --make-bed --out binary/gxtc.postimputation.maf05.hwe.udid.keeps.udsex

plink --bfile binary/gxtc.postimputation.maf05.hwe.udid.keeps.udsex --pheno pheno/gxtc.pc.post --make-bed --out binary/gxtc.postimputation.maf05.hwe.updated.2002.pc.post
plink --bfile binary/gxtc.postimputation.maf05.hwe.udid.keeps.udsex --pheno pheno/gxtc.pc.fu --make-bed --out binary/gxtc.postimputation.maf05.hwe.updated.2002.pc.fu
plink --bfile binary/gxtc.postimputation.maf05.hwe.udid.keeps.udsex --pheno pheno/gxtc.zscore.post --make-bed --out binary/gxtc.postimputation.maf05.hwe.updated.2002.zscore.post
plink --bfile binary/gxtc.postimputation.maf05.hwe.udid.keeps.udsex --pheno pheno/gxtc.zscore.fu --make-bed --out binary/gxtc.postimputation.maf05.hwe.updated.2002.zscore.fu

gxtc = 1 significant PCs
gcta --bfile binary/gxtc.postimputation.maf05.hwe.updated.2002.pc.fu --autosome --maf 0.01 --thread-num 2 --make-grm --out grm/gxtc.postimputation.maf05.hwe.updated.2002
gcta --grm grm/gxtc.postimputation.maf05.hwe.updated.2002 --pca --out cov/pcs/gxtc.pca
awk '{print $1,$2,$3}' cov/pcs/gxtc.pca.eigenvec > cov/pcs/gxtc.one.pc
awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' cov/pcs/gxtc.pca.eigenvec > cov/gxtc.ten.pc
cd /users/k1620287/brc_scratch/meta.analysis/analysis/gxtc/cov
R
ten <- read.table("pcs/gxtc.ten.pc", head =F)
one <- read.table("pcs/gxtc.one.pc", head =F)
conc <- read.table("batch/gxtc.conc.txt", head=F)
names(ten)[1] <- "FID"
tenconc <- merge(ten, conc, by = "FID")
tenconc$IID <- NULL
names(one)[1] <- "FID"
oneconc <- merge(one, conc, by = "FID")
oneconc$IID <- NULL
write.table(oneconc,"final/one.pc.conc.qcov", col.names=F, row.names=F, quote=F, sep=" ")
write.table(tenconc,"final/ten.pc.conc.qcov", col.names=F, row.names=F, quote=F, sep=" ")
q()
n

tail -n +2 batch/gxtc.type.filt.txt > final/gxtc.type.filt.cov

covariate files:
/users/k1620287/brc_scratch/meta.analysis/analysis/gxtc/cov/final
gxtc.one.pc.conc.qcov
gxtc.ten.pc.conc.qcov
gxtc.type.filt.cov


###########################
# MERGING THREE COHORTS
###########################

snplist /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt_v1/analysis_files/snps_shared_between_all_cohorts
adult anx /users/k1620287/brc_scratch/meta.analysis/analysis/adult_merged/binary/adult_merged.postimputation.maf05.hwe.udid.keeps.udsex
adult md /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/ki_mdd.postimputation.maf05.hwe.keeps.udsex.pheno.locf
child anx /users/k1620287/brc_scratch/meta.analysis/analysis/gxtc/binary/gxtc.postimputation.maf05.hwe.udid.keeps.udsex

module add bioinformatics/plink2/1.90b3.38

plink \
--bfile /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/ki_mdd.postimputation.maf05.hwe.keeps.udsex.pheno.locf \
--extract /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt_v1/analysis_files/snps_shared_between_all_cohorts \
--make-bed \
--out /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/ki_mdd.postimputation.maf05.hwe.keeps.udsex.pheno.locf.shared.snps

plink \
--bfile /users/k1620287/brc_scratch/meta.analysis/analysis/adult_merged/binary/adult_merged.postimputation.maf05.hwe.udid.keeps.udsex \
--extract /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt_v1/analysis_files/snps_shared_between_all_cohorts \
--make-bed \
--out /users/k1620287/brc_scratch/meta.analysis/analysis/adult_merged/binary/adult_merged.postimputation.maf05.hwe.udid.keeps.udsex.shared.snps

plink \
--bfile /users/k1620287/brc_scratch/meta.analysis/analysis/gxtc/binary/gxtc.postimputation.maf05.hwe.udid.keeps.udsex \
--extract /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt_v1/analysis_files/snps_shared_between_all_cohorts \
--make-bed \
--out /users/k1620287/brc_scratch/meta.analysis/analysis/gxtc/binary/gxtc.postimputation.maf05.hwe.udid.keeps.udsex.shared.snps

# MERGE

plink \
--bfile /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/ki_mdd.postimputation.maf05.hwe.keeps.udsex.pheno.locf.shared.snps \
--bmerge /users/k1620287/brc_scratch/meta.analysis/analysis/adult_merged/binary/adult_merged.postimputation.maf05.hwe.udid.keeps.udsex.shared.snps.bed \
/users/k1620287/brc_scratch/meta.analysis/analysis/adult_merged/binary/adult_merged.postimputation.maf05.hwe.udid.keeps.udsex.shared.snps.bim \
/users/k1620287/brc_scratch/meta.analysis/analysis/adult_merged/binary/adult_merged.postimputation.maf05.hwe.udid.keeps.udsex.shared.snps.fam \
--make-bed \
--out /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/adult_anx.adult_mdd.shared_snps


plink \
--bfile /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/adult_anx.adult_mdd.shared_snps \
--bmerge /users/k1620287/brc_scratch/meta.analysis/analysis/gxtc/binary/gxtc.postimputation.maf05.hwe.udid.keeps.udsex.shared.snps.bed \
/users/k1620287/brc_scratch/meta.analysis/analysis/gxtc/binary/gxtc.postimputation.maf05.hwe.udid.keeps.udsex.shared.snps.bim \
/users/k1620287/brc_scratch/meta.analysis/analysis/gxtc/binary/gxtc.postimputation.maf05.hwe.udid.keeps.udsex.shared.snps.fam \
--make-bed \
--out /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/adult_anx.adult_mdd.child_anx.shared_snps
