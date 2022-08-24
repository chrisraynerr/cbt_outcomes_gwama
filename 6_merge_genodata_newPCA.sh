####################################
# Meta-analysis of response to CBT #
####################################

# Script 5. Merging data & PCA of merged data

mkdir /mnt/lustre/users/k1620287/meta.analysis/qc/adult.merged
cd /mnt/lustre/users/k1620287/meta.analysis/qc/adult.merged
module add bioinformatics/plink2/1.90b3.38

# TO DO: might be worth swaping these around: moving from most SNPs to least SNPS

site1=wb1
site2=wb2
site3=gxte
site4=ki

plink \
--bfile /mnt/lustre/users/k1620287/meta.analysis/qc/$site1/$site1.full.qc.build.updated \
--bmerge /mnt/lustre/users/k1620287/meta.analysis/qc/$site2/$site2.full.qc.build.updated.bed \
		/mnt/lustre/users/k1620287/meta.analysis/qc/$site2/$site2.full.qc.build.updated.bim \
		/mnt/lustre/users/k1620287/meta.analysis/qc/$site2/$site2.full.qc.build.updated.fam \
--make-bed \
--out /mnt/lustre/users/k1620287/meta.analysis/qc/adult.merged/$site1.$site2.qc.merged

plink \
--bfile /mnt/lustre/users/k1620287/meta.analysis/qc/$site2/$site2.full.qc.build.updated \
--flip /mnt/lustre/users/k1620287/meta.analysis/adult.merged/$site1.$site2.qc.merged-merge.missnp \
--make-bed \
--out /mnt/lustre/users/k1620287/meta.analysis/qc/$site2/$site2.full.qc.build.updated.flipped

plink \
--bfile /mnt/lustre/users/k1620287/meta.analysis/qc/$site1/$site1.qc \
--bmerge /mnt/lustre/users/k1620287/meta.analysis/qc/$site2/$site2.full.qc.build.updated.flipped.bed \
		/mnt/lustre/users/k1620287/meta.analysis/qc/$site2/$site2.full.qc.build.updated.flipped.bim \
		/mnt/lustre/users/k1620287/meta.analysis/qc/$site2/$site2.full.qc.build.updated.flipped.fam \
--out /mnt/lustre/users/k1620287/meta.analysis/qc/adult.merged/$site1.$site2.qc.merged

plink \
--bfile /mnt/lustre/users/k1620287/meta.analysis/qc/adult.merged/$site1.$site2.qc.merged \
--bmerge /mnt/lustre/users/k1620287/meta.analysis/qc/$site3/$site3.full.qc.build.updated.bed \
		/mnt/lustre/users/k1620287/meta.analysis/qc/$site3/$site3.full.qc.build.updated.bim \
		/mnt/lustre/users/k1620287/meta.analysis/qc/$site3/$site3.full.qc.build.updated.fam \
--make-bed \
--out /mnt/lustre/users/k1620287/meta.analysis/qc/adult.merged/$site1.$site2.$site3.qc.merged

plink \
--bfile /mnt/lustre/users/k1620287/meta.analysis/qc/$site3/$site3.full.qc.build.updated \
--flip /mnt/lustre/users/k1620287/meta.analysis/adult.merged/$site1.$site2.$site3.qc.merged-merge.missnp \
--make-bed \
--out /mnt/lustre/users/k1620287/meta.analysis/qc/$site3/$site3.full.qc.build.updated.flipped

plink \
--bfile /mnt/lustre/users/k1620287/meta.analysis/qc/adult.merged/$site1.$site2.qc.merged \
--bmerge /mnt/lustre/users/k1620287/meta.analysis/qc/$site3/$site3.full.qc.build.updated.flipped.bed \
		/mnt/lustre/users/k1620287/meta.analysis/qc/$site3/$site3.full.qc.build.updated.flipped.bim \
		/mnt/lustre/users/k1620287/meta.analysis/qc/$site3/$site3.full.qc.build.updated.flipped.fam \
--make-bed \
--out /mnt/lustre/users/k1620287/meta.analysis/qc/adult.merged/$site1.$site2.$site3.qc.merged

plink \
--bfile /mnt/lustre/users/k1620287/meta.analysis/qc/adult.merged/$site1.$site2.$site3.qc.merged \
--bmerge /mnt/lustre/users/k1620287/meta.analysis/qc/$site4/$site4.full.qc.build.updated.bed \
		/mnt/lustre/users/k1620287/meta.analysis/qc/$site4/$site4.full.qc.build.updated.bim \
		/mnt/lustre/users/k1620287/meta.analysis/qc/$site4/$site4.full.qc.build.updated.fam \
--make-bed \
--out /mnt/lustre/users/k1620287/meta.analysis/qc/adult.merged/$site1.$site2.$site3.$site4.qc.merged

plink \
--bfile /mnt/lustre/users/k1620287/meta.analysis/qc/$site4/$site4.full.qc.build.updated \
--flip /mnt/lustre/users/k1620287/meta.analysis/adult.merged/$site1.$site2.$site3.$site4.qc.merged-merge.missnp \
--make-bed \
--out /mnt/lustre/users/k1620287/meta.analysis/qc/$site4/$site4.full.qc.build.updated.flipped

plink \
--bfile /mnt/lustre/users/k1620287/meta.analysis/qc/adult.merged/$site1.$site2.$site3.qc.merged \
--bmerge /mnt/lustre/users/k1620287/meta.analysis/qc/$site4/$site4.full.qc.build.updated.flipped.bed \
		/mnt/lustre/users/k1620287/meta.analysis/qc/$site4/$site4.full.qc.build.updated.flipped.bim \
		/mnt/lustre/users/k1620287/meta.analysis/qc/$site4/$site4.full.qc.build.updated.flipped.fam \
--make-bed \
--out /mnt/lustre/users/k1620287/meta.analysis/qc/adult.merged/$site1.$site2.$site3.$site4.qc.merged

# If flipping the snps doesnt work and there is a low number of missnps, you can just remove them:
plink \
--bfile /mnt/lustre/users/k1620287/meta.analysis/qc/adult.merged/$site1.$site2.$site3.qc.merged \
--exclude /mnt/lustre/users/k1620287/meta.analysis/adult.merged/$site1.$site2.$site3.$site4.qc.merged-merge.missnp \
--make-bed \
--out /mnt/lustre/users/k1620287/meta.analysis/qc/adult.merged/$site1.$site2.$site3.qc.merged.missnp.exc

plink \
--bfile /mnt/lustre/users/k1620287/meta.analysis/qc/adult.merged/$site1.$site2.$site3.qc.merged.missnp.exc \
--bmerge /mnt/lustre/users/k1620287/meta.analysis/qc/$site4/$site4.full.qc.build.updated.flipped.bed \
		/mnt/lustre/users/k1620287/meta.analysis/qc/$site4/$site4.full.qc.build.updated.flipped.bim \
		/mnt/lustre/users/k1620287/meta.analysis/qc/$site4/$site4.full.qc.build.updated.flipped.fam \
--make-bed \
--out /mnt/lustre/users/k1620287/meta.analysis/qc/adult.merged/$site1.$site2.$site3.$site4.qc.merged

awk '{print $2}' /mnt/lustre/users/k1620287/meta.analysis/qc/adult.merged/$site1.$site2.$site3.$site4.qc.merged \
| sort | uniq -d > /mnt/lustre/users/k1620287/meta.analysis/qc/adult.merged/all.sites.merged.duplicates.to.remove

# wc -l /mnt/lustre/users/k1620287/meta.analysis/qc/adult.merged/merge/all.sites.merged.duplicates.to.remove
# 0 /mnt/lustre/users/k1620287/meta.analysis/qc/adult.merged/all.sites.merged.duplicates.to.remove

#plink \
#--bfile /mnt/lustre/users/k1620287/meta.analysis/qc/adult.merged/$site1.$site2.$site3.$site4.qc.merged \
#--exclude /mnt/lustre/users/k1620287/meta.analysis/qc/adult.merged/all.sites.merged.duplicates.to.remove \
#--make-bed \
#--out /mnt/lustre/users/k1620287/meta.analysis/qc/adult.merged/$site1.$site2.$site3.$site4.qc.merged.snps

sort -uk4,4 /mnt/lustre/users/k1620287/meta.analysis/qc/adult.merged/$site1.$site2.$site3.$site4.qc.merged.bim \
| awk '{print $2}' > /mnt/lustre/users/k1620287/meta.analysis/qc/adult.merged/merged.keep.snps

# wc -l /mnt/lustre/users/k1620287/meta.analysis/qc/adult.merged/merged.keep.SNPs
#664456 /mnt/lustre/users/k1620287/meta.analysis/qc/adult.merged/merged.keep.SNPs

rm */*flipped*
rm */*missnp*
rm */*mis.exc.*

plink \
--bfile /mnt/lustre/users/k1620287/meta.analysis/qc/adult.merged/$site1.$site2.$site3.$site4.qc.merged \
--extract /mnt/lustre/users/k1620287/meta.analysis/qc/adult.merged/merged.keep.snps \
--make-bed \
--out /mnt/lustre/users/k1620287/meta.analysis/qc/adult.merged/$site1.$site2.$site3.$site4.qc.merged.dup.rm

plink \
--bfile /mnt/lustre/users/k1620287/meta.analysis/qc/adult.merged/$site1.$site2.$site3.$site4.qc.merged.dup.rm \
--maf 0.05 \
--make-bed \
--out /mnt/lustre/users/k1620287/meta.analysis/qc/adult.merged/$site1.$site2.$site3.$site4.qc.merged.dup.rm.maf.05


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# TO DO: Check this:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
site=./adult.merged/$site1.$site2.$site3.$site4.qc.merged.dup.rm
site=./merge/$site1.$site2.$site3.$site4
sh Iterative_Missingness.sh 90 99 1

mv .............. adult_merged.maf05.mind01.geno1

update phenotype
NEED TO EDIT:
for site in 'ki' 'wb1' 'wb2' 'gxte' 'adult_merged'
do
cd /users/k1620287/brc_scratch/meta.analysis/analysis/WAVE2/$site
plink --bfile $site.hrc.imputed.qc --keep $keeps --make-bed --out $site.hrc.qc.keeps
plink --bfile $site.hrc.qc.keeps --update-sex $genders --make-bed --out $site.hrc.qc.updated_gen
plink --bfile $site.hrc.qc.updated_gen --pheno ../$percent_no_outliers --make-bed --out $site.hrc.qc.updated
done


###############################################
# Principal Component Analysis in Merged Data #
###############################################

cat <<'EOT'>> pca.for.all.sh
#!/bin/sh
#$ -V
#$ -m ea
#$ -pe smp 10
#$ -l h_vmem=4G
#$ -l s_vmem=4G
#$ -l h_rt=04:00:00
#$ -S /bin/sh
#$ -m be
#$ -o /users/k1620287/brc_scratch/meta.analysis/adult.merged/
#$ -e /users/k1620287/brc_scratch/meta.analysis/adult.merged/

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
UPDATED UPTO HERE -ALSO DONE ALL IMUTATION STUFF!

module add bioinformatics/R/3.3.3
module add bioinformatics/EIGENSOFT/6.1.4

for site in 'ki' 'wb1' 'wb2' 'gxte' 'adult'
do
cd /users/k1620287/brc_scratch/meta.analysis/analysis/WAVE2/$site

convertf -p <(printf "genotypename: $site.qc.updated_gen_phen.bed
snpname: $site.qc.updated_gen_phen.bim
indivname: $site.qc.updated_gen_phen.fam
outputformat: EIGENSTRAT
genotypeoutname: $site.PCS_for_covariates.eigenstratgeno
snpoutname: $site.PCS_for_covariates.snp
indivoutname: $site.PCS_for_covariates.ind")
done

for site in 'ki' 'wb1' 'wb2' 'gxte' 'adult'
do
cd /users/k1620287/brc_scratch/meta.analysis/analysis/WAVE2/$site
module add bioinformatics/EIGENSOFT/6.1.4

smartpca.perl \
-i $site.PCS_for_covariates.eigenstratgeno \
-a $site.PCS_for_covariates.snp \
-b $site.PCS_for_covariates.ind \
-o $site.PCS_for_covariates.pca \
-p $site.PCS_for_covariates.plot \
-e $site.PCS_for_covariates.eval \
-l $site.PCS_for_covariates_smartpca.log \
-m 0 \
-t 100 \
-k 100 \
-s 6
done

for site in 'ki' 'wb1' 'wb2' 'gxte' 'adult'

do

cd /users/k1620287/brc_scratch/meta.analysis/analysis/$site/cov

sed -i -e 's/^[ \t]*//' -e 's/:/ /g' $site.PCS_for_covariates.pca.evec

R --file=/mnt/lustre/users/k1620287/meta.analysis/gwas.scripts/PC-VS-OUTCOME_IN_R_SHORT.R \
	--args $site.PCS_for_covariates $pheno

awk '$3 < 0.05 {print $1,$2,$3,$4}' $site.PCS_for_covariates.PC_Output_Associations_SHORT.txt > $site.significant.pcs.txt

sig_assoc_pcs=`wc -l $site.significant.pcs.txt | awk '{print $1}'`

sed -i '1d' $site.PCS_for_covariates.pca.evec

echo 'IID FID PC1  PC2  PC3  PC4  PC5  PC6  PC7  PC8  PC9  PC10  PC11  PC12  PC13  PC14  PC15  PC16  PC17  PC18  PC19  PC20  PC21  PC22  PC23  PC24  PC25  PC26  PC27  PC28  PC29  PC30  PC31  PC32  PC33  PC34  PC35  PC36  PC37  PC38  PC39  PC40  PC41  PC42  PC43  PC44  PC45  PC46  PC47  PC48  PC49  PC50  PC51  PC52  PC53  PC54  PC55  PC56  PC57  PC58  PC59  PC60  PC61  PC62  PC63  PC64  PC65  PC66  PC67  PC68  PC69  PC70  PC71  PC72  PC73  PC74  PC75  PC76  PC77  PC78  PC79  PC80  PC81  PC82  PC83  PC84  PC85  PC86  PC87  PC88  PC89  PC90  PC91  PC92  PC93  PC94  PC95  PC96  PC97  PC98  PC99  PC100' > pc.header

cat pc.header $site.PCS_for_covariates.pca.evec | column -t > /users/k1620287/brc_scratch/meta.analysis/analysis/$site/cov/$site.PCS_for_covariates

n=$sig_assoc_pcs+2

awk '{for (i=1; i>=NF-'$sig_assoc_pcs'; i++) $i = $(i+'$sig_assoc_pcs'); NF='$sig_assoc_pcs'; print}' \
/users/k1620287/brc_scratch/meta.analysis/analysis/$site/cov/$site.PCS_for_covariates > /users/k1620287/brc_scratch/meta.analysis/analysis/$site/cov/$site.pc.covariates.txt

awk '{for (i=1; i>=NF-'$n'; i++) $i = $(i+'$n'); NF='$n'; print}' \
/users/k1620287/brc_scratch/meta.analysis/analysis/$site/cov/$site.PCS_for_covariates > /users/k1620287/brc_scratch/meta.analysis/analysis/$site/cov/$site.pc.covariates.txt



EOT
qsub pca.for.all.sh
