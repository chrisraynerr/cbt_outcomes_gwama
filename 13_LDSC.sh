####################################
# Meta-analysis of response to CBT #
####################################

# Script 13. LD SCORE REGRESSION #

#######################
# LD score regression #
#######################
# git clone https://github.com/bulik/ldsc.git
# cd ldsc
# wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
# wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
# tar -jxvf eur_w_ld_chr.tar.bz2
# bunzip2 w_hm3.snplist.bz2


cd /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files

scp /users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/assoc/child_anx.loco.mlma.results ./
scp /users/k1620287/brc_scratch/meta.analysis/analysis/final_analyses/assoc/adult_anx.loco.mlma.results ./

/mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/ki_mdd.loco.mlma.results

/mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/adult_anx.child_anx.ki_mdd.se.gwama.results1.tbl
/mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/adult_anx.ki_mdd.se.gwama.results1.tbl
/mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/adult_anx.child_anx.se.gwama.results1.tbl

python=/opt/apps/general/anaconda/2.5.0/anaconda2/bin/python
ldsc=/mnt/lustre/groups/ukbiobank/usr/crayner/ldsc

for i in \
child_anx.loco.mlma.results \
adult_anx.loco.mlma.results \
ki_mdd.loco.mlma.results \
adult_anx.child_anx.ki_mdd.se.gwama.results1.tbl \
adult_anx.ki_mdd.se.gwama.results1.tbl \
adult_anx.child_anx.se.gwama.results1.tbl
do
$python $ldsc/munge_sumstats.py \
--sumstats $i \
--out $i \
--merge-alleles $ldsc/w_hm3.snplist
done

for i in \
child_anx.loco.mlma.results \
adult_anx.loco.mlma.results \
ki_mdd.loco.mlma.results \
adult_anx.child_anx.ki_mdd.se.gwama.results1.tbl \
adult_anx.ki_mdd.se.gwama.results1.tbl \
adult_anx.child_anx.se.gwama.results1.tbl
do
$python $ldsc/ldsc.py \
--h2 $i.sumstats.gz  \
--ref-ld-chr $ldsc/eur_w_ld_chr/ \
--w-ld-chr $ldsc/eur_w_ld_chr/ \
--out $i.ldsc.h2
done



# After ancestry outliers removed
python=/opt/apps/general/anaconda/2.5.0/anaconda2/bin/python
ldsc=/mnt/lustre/groups/ukbiobank/usr/crayner/ldsc
for i in \
ki_mdd.postimputation.maf05.hwe.keeps.udsex.pheno.locf.3.29S_POP.loco.mlma
do
$python $ldsc/munge_sumstats.py \
--sumstats $i \
--out $i \
--merge-alleles $ldsc/w_hm3.snplist \
--N 778
done

for i in \
ki_mdd.postimputation.maf05.hwe.keeps.udsex.pheno.locf.3.29S_POP.loco.mlma
do
$python $ldsc/ldsc.py \
--h2 $i.sumstats.gz  \
--ref-ld-chr $ldsc/eur_w_ld_chr/ \
--w-ld-chr $ldsc/eur_w_ld_chr/ \
--out $i.ldsc.h2
done


# REMOVING MHC


for i in \
ki_mdd.loco.mlma.results \
adult_anx.child_anx.ki_mdd.se.gwama.results1.tbl \
adult_anx.ki_mdd.se.gwama.results1.tbl
do
if [ -f $i.sumstats.gz ]
then zcat $i.sumstats.gz > $i.tmp
gawk -F"\t" '{if(f==1){r[$1]}else if(!($1 in r)){print}}' f=1 /mnt/lustre/groups/ukbiobank/sumstats/scripts/hapmap_hmc_uniq.txt f=2 ${i}.tmp > ${i}.noMHC.sumstats; gzip ${i}.noMHC.sumstats; fi;
done

for i in \
ki_mdd.loco.mlma.results \
adult_anx.child_anx.ki_mdd.se.gwama.results1.tbl \
adult_anx.ki_mdd.se.gwama.results1.tbl
do
zcat $i.noMHC.sumstats.gz > $i.tmp
awk '$2!=""' $i.tmp > $i.noMHC.sumstats.rmblank; gzip $i.noMHC.sumstats.rmblank
done

for i in \
ki_mdd.loco.mlma.results \
adult_anx.child_anx.ki_mdd.se.gwama.results1.tbl \
adult_anx.ki_mdd.se.gwama.results1.tbl
do
$python $ldsc/ldsc.py \
--h2 $i.noMHC.sumstats.rmblank.gz  \
--ref-ld-chr $ldsc/eur_w_ld_chr/ \
--w-ld-chr $ldsc/eur_w_ld_chr/ \
--out $i.noMHC.ldsc.h2
done






# Therapy x AD treatment response

python=/opt/apps/general/anaconda/2.5.0/anaconda2/bin/python
munge=/mnt/lustre/groups/ukbiobank/Edinburgh_Data/Software/ldscore/ldsc/munge_sumstats.py
weights=/mnt/lustre/groups/ukbiobank/Edinburgh_Data/Software/ldscore/ldsc/weights_hm3_no_hla/

$python $munge \
--sumstats meta_improv_SSRIs_new_plink.meta \
--out munged/SSRIRESP.gwas \
--merge-alleles $weights \
--N 1828

$python $munge \
--sumstats meta_rem_SSRIs_new_plink.meta \
--out munged/SSRIREM.gwas \
--merge-alleles $weights \
--N 1828

$python /mnt/lustre/groups/ukbiobank/usr/crayner/ldsc/munge_sumstats.py \
--sumstats meta_rem_SSRIs_new_plink.meta \
--out munged/SSRIREM.gwas \
--merge-alleles /mnt/lustre/groups/ukbiobank/usr/crayner/ldsc/w_hm3.snplist \
--N 1828

$python /mnt/lustre/groups/ukbiobank/usr/crayner/ldsc/munge_sumstats.py \
--sumstats /mnt/lustre/groups/ukbiobank/usr/crayner/ukb18177/assoc/results_5/meta_improv_new_plink.meta \
--out /mnt/lustre/groups/ukbiobank/usr/crayner/ukb18177/assoc/results_5/munged/ANTIDEPRESP.gwas \
--merge-alleles /mnt/lustre/groups/ukbiobank/usr/crayner/ldsc/w_hm3.snplist \
--N 2145

$python /mnt/lustre/groups/ukbiobank/usr/crayner/ldsc/munge_sumstats.py \
--sumstats /mnt/lustre/groups/ukbiobank/usr/crayner/ukb18177/assoc/results_5/meta_rem_new_plink.meta \
--out /mnt/lustre/groups/ukbiobank/usr/crayner/ukb18177/assoc/results_5/munged/ANTIDEPREM.gwas \
--merge-alleles /mnt/lustre/groups/ukbiobank/usr/crayner/ldsc/w_hm3.snplist \
--N 2145
