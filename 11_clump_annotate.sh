####################################
# Meta-analysis of response to CBT #
####################################

# Script 11. Clump and annotate all results #

bfile_1=/users/k1620287/brc_scratch/meta.analysis/analysis/gxtc/binary/gxtc.postimputation.maf05.hwe.udid.keeps.udsex
bfile_2=/users/k1620287/brc_scratch/meta.analysis/analysis/adult_merged/binary/adult_merged.postimputation.maf05.hwe.udid.keeps.udsex
bfile_3=/mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/ki_mdd.postimputation.maf05.hwe.keeps.udsex.pheno.locf
bfile_4=/mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/adult_anx.adult_mdd.child_anx.shared_snps
bfile_5=/mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/adult_anx.adult_mdd.child_anx.shared_snps
bfile_6=/mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/adult_anx.adult_mdd.child_anx.shared_snps

results_1=/mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/child_anx.loco.mlma.results
results_2=/mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/adult_anx.loco.mlma.results
results_3=/mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/ki_mdd.loco.mlma.results
results_4=/mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/adult_anx.child_anx.ki_mdd.se.gwama.results1.tbl
results_5=/mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/adult_anx.ki_mdd.se.gwama.results1.tbl
results_6=/mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/adult_anx.child_anx.se.gwama.results1.tbl

rm /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/scripts/clump.sh
cat <<'EOT'>> /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/scripts/clump.sh
#!/bin/sh
#$ -V
#$ -pe smp 8
#$ -l h_vmem=20G
#$ -l h_rt=04:00:00
#$ -S /bin/sh
#$ -m ea
#$ -m be
#$ -o /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files
#$ -e /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files

module load bioinformatics/plink/1.90b3.31

cd /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files

plink \
--bfile $1 \
--clump $2 \
--clump-snp-field SNP MarkerName \
--clump-field P P-value \
--clump-p1 0.0001 \
--clump-p2 0.0001 \
--clump-r2 0.1 \
--clump-kb 3000 \
--clump-range /mnt/lustre/groups/ukbiobank/Edinburgh_Data/Resources/glist-hg19 \
--out annot/$2.clumped.for.RA

# plink \
# --bfile $1 \
# --clump $2 \
# --clump-snp-field SNP MarkerName \
# --clump-field P P-value \
# --clump-p1 1 \
# --clump-p2 1 \
# --clump-r2 0.25 \
# --clump-kb 250 \
# --clump-range /mnt/lustre/groups/ukbiobank/Edinburgh_Data/Resources/glist-hg19 \
# --out $2.clump
EOT

while read clump_file_list; do qsub /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/scripts/clump.sh $clump_file_list; done < clump_file_list



rm /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/scripts/anot.sh
cat <<'EOT'>> /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/scripts/anot.sh
#!/bin/sh
#$ -V
#$ -pe smp 8
#$ -l h_vmem=20G
#$ -l h_rt=04:00:00
#$ -S /bin/sh
#$ -m ea
#$ -m be
#$ -o /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/
#$ -e /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/

module load bioinformatics/plink/1.90b3.31

cd /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files

for i in \
child_anx.loco.mlma.results \
adult_anx.loco.mlma.results \
ki_mdd.loco.mlma.results \
adult_anx.child_anx.ki_mdd.se.gwama.results1.tbl \
adult_anx.ki_mdd.se.gwama.results1.tbl \
adult_anx.child_anx.se.gwama.results1.tbl
do

plink \
--annotate annot/$i.clumped.for.RA.clumped \
ranges=/mnt/lustre/groups/ukbiobank/Edinburgh_Data/Resources/glist-hg19 \
--border 250 \
--out $i.clump.annot.2

done
EOT
qsub /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/scripts/anot.sh


for i in \
child_anx.loco.mlma.results \
adult_anx.loco.mlma.results \
ki_mdd.loco.mlma.results \
adult_anx.child_anx.ki_mdd.se.gwama.results1.tbl \
adult_anx.ki_mdd.se.gwama.results1.tbl \
adult_anx.child_anx.se.gwama.results1.tbl
do
awk '{print $1,$3,$4,$5,$13}' $i.clump.annot.2.annot | column -t > $i.clump.annot.2.annot.less
done

for i in \
child_anx.loco.mlma.results \
adult_anx.loco.mlma.results \
ki_mdd.loco.mlma.results \
adult_anx.child_anx.ki_mdd.se.gwama.results1.tbl \
adult_anx.ki_mdd.se.gwama.results1.tbl \
adult_anx.child_anx.se.gwama.results1.tbl
do
awk '{print $3}' $i.clump.annot.2.annot | column -t > $i.clumped.snp.list
done



for i in \
child_anx.loco.mlma.results \
adult_anx.loco.mlma.results \
ki_mdd.loco.mlma.results
do
echo '

library(data.table)

dat <- fread("'$i'", head=T)
snp <- fread("'$i'.clumped.snp.list", head=T)
annot <- fread("'$i'.clump.annot.2.annot.less", head=T)
dat2 <- merge(snp, dat, by="SNP")
dat2 <- merge(dat2, annot, by="SNP")
write.table(dat2, "'$i'.clumped.significant.snps", col.names=T, row.names=F, quote=F, sep=" ")

'
done

for i in \
adult_anx.child_anx.ki_mdd.se.gwama.results1.tbl \
adult_anx.ki_mdd.se.gwama.results1.tbl \
adult_anx.child_anx.se.gwama.results1.tbl
do
echo '
library(data.table)

dat <- fread("'$i'", head=T)
snp <- fread("'$i'.clumped.snp.list", head=T)
names(dat)[1] <- "SNP"
annot <- fread("'$i'.clump.annot.2.annot.less", head=T)
dat2 <- merge(snp, dat, by="SNP")
dat2 <- merge(dat2, annot, by="SNP")
write.table(dat2, "'$i'.clumped.significant.snps", col.names=T, row.names=F, quote=F, sep=" ")

'
done



while read $i.clumped.snp.list
do
grep $i.clumped.snp.list $i > $i.clumped.snp.list.2
done < $i.clumped.snp.list
done










for i in \
child_anx.loco.mlma.results \
adult_anx.loco.mlma.results \
ki_mdd.loco.mlma.results \
adult_anx.child_anx.ki_mdd.se.gwama.results1.tbl \
adult_anx.ki_mdd.se.gwama.results1.tbl \
adult_anx.child_anx.se.gwama.results1.tbl
do
awk ' {if(NR > 1){print $5,$3}}' ${i}.clumped.for.RA.clumped.ranges > ${i}.clumped.for.RA.clumped.ranges.1
done


for i in \
child_anx.loco.mlma.results \
adult_anx.loco.mlma.results \
ki_mdd.loco.mlma.results \
adult_anx.child_anx.ki_mdd.se.gwama.results1.tbl \
adult_anx.ki_mdd.se.gwama.results1.tbl \
adult_anx.child_anx.se.gwama.results1.tbl
do
sed -e 's/\:/ /g' -e 's/\.\./ /g' ${i}.clumped.for.RA.clumped.ranges.1 > ${i}.clumped.for.RA.clumped.ranges.2
done

for i in \
child_anx.loco.mlma.results \
adult_anx.loco.mlma.results \
ki_mdd.loco.mlma.results \
adult_anx.child_anx.ki_mdd.se.gwama.results1.tbl \
adult_anx.ki_mdd.se.gwama.results1.tbl \
adult_anx.child_anx.se.gwama.results1.tbl
do
sort -k1,1 -k2,2n ${i}.clumped.for.RA.clumped.ranges.2 | \
awk -F $'\t' '{print $1, $2, $3, $4}' | \
tr ' ' '\t' >> ${i}.clumped.bed
done

for i in \
child_anx.loco.mlma.results \
adult_anx.loco.mlma.results \
ki_mdd.loco.mlma.results \
adult_anx.child_anx.ki_mdd.se.gwama.results1.tbl \
adult_anx.ki_mdd.se.gwama.results1.tbl \
adult_anx.child_anx.se.gwama.results1.tbl
do
/opt/apps/bioinformatics/bedtools2/2.25.0/bin/bedtools \
merge -d 50000 -c 4 -o min \
-i ${i}.clumped.bed > ${i}.clumped.bedtooled
done

echo 'CHR BP1 BP2 P' > bedtooled.header

for i in \
child_anx.loco.mlma.results \
adult_anx.loco.mlma.results \
ki_mdd.loco.mlma.results \
adult_anx.child_anx.ki_mdd.se.gwama.results1.tbl \
adult_anx.ki_mdd.se.gwama.results1.tbl \
adult_anx.child_anx.se.gwama.results1.tbl
do
cat bedtooled.header ${i}.clumped.bedtooled | column -t > ${i}.clumped.bedtooled.headered
done

rm /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/annot/annotate.results.sh
cat <<'EOT'>> /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/annot/annotate.results.sh
#!/bin/sh
#$-S /bin/sh
#$ -V
#$ -pe smp 4
#$ -l h_vmem=20G
#$ -l h_rt=24:00:00
#$ -m be
#$ -o /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/annot/

export MKL_NUM_THREADS=4
export NUMEXPR_NUM_THREADS=4
export OMP_NUM_THREADS=4

cd /mnt/lustre/groups/ukbiobank/Edinburgh_Data/Software/RegionAnnotator/

#Load the database with gene data. One file. If .txt we need to specify the input format.

sh /mnt/lustre/groups/ukbiobank/Edinburgh_Data/Software/RegionAnnotator/RegionAnnotator.sh \
-input /mnt/lustre/groups/ukbiobank/Edinburgh_Data/Software/RegionAnnotator/inputGene/gencode.genes.txt \
-gene \
-iformat TSV

#Load the database with reference data. Since it's many files, we can point to a directory. If .txt-files we need to specify the input format.

sh /mnt/lustre/groups/ukbiobank/Edinburgh_Data/Software/RegionAnnotator/RegionAnnotator.sh \
-input /mnt/lustre/groups/ukbiobank/Edinburgh_Data/Software/RegionAnnotator/inputReference \
-reference \
-iformat TSV

#Run the program on an input file with regions/genes.
#The program operations will be run, using the previously loaded gene and reference data, and an output-file will be created (standardized EXCEL-file as default).
#If the input file is a .txt we need to specify the input format. If risk of the database being in use - you can use a timeout-flag
#- the program will then wait the specified number of seconds for the database to become available.

printf running phenotype: $i \n

sh /mnt/lustre/groups/ukbiobank/Edinburgh_Data/Software/RegionAnnotator/RegionAnnotator.sh \
-input /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/annot/$1.clumped.bedtooled.headered \
-iformat TSV \
-timeout 60000000000 \
-output /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/annot/$1.RA.xlsx

EOT
for i in \
child_anx.loco.mlma.results \
adult_anx.loco.mlma.results \
ki_mdd.loco.mlma.results \
adult_anx.child_anx.ki_mdd.se.gwama.results1.tbl \
adult_anx.ki_mdd.se.gwama.results1.tbl \
adult_anx.child_anx.se.gwama.results1.tbl
do
qsub /mnt/lustre/groups/ukbiobank/usr/crayner/ki_mdd_cbt/analysis_files/annot/annotate.results.sh $i
done








rm /users/k1620287/brc_scratch/meta.analysis/analysis/meta/scripts/clump.sh
cat <<'EOT'>> /users/k1620287/brc_scratch/meta.analysis/analysis/meta/scripts/clump.sh
#!/bin/sh
#$ -V
#$ -pe smp 8
#$ -l h_vmem=20G
#$ -l h_rt=04:00:00
#$ -S /bin/sh
#$ -m ea
#$ -m be
#$ -o /users/k1620287/brc_scratch/meta.analysis/analysis/meta/output/
#$ -e /users/k1620287/brc_scratch/meta.analysis/analysis/meta/output/

module add bioinformatics/plink2/1.90b3.38
for data in gxtc adult_merged
{
	for pheno in percent zscore
	{
	plink2 \
	--bfile /users/k1620287/brc_scratch/meta.analysis/analysis/$data/binary/$data.postimputation.maf05.hwe.udid.keeps.udsex \
	--clump /users/k1620287/brc_scratch/meta.analysis/analysis/meta/results/$data.$pheno.post.mlma.for.meta \
	--clump-p1 1 \
	--clump-p2 1 \
	--clump-r2 0.25 \
	--clump-kb 250 \
	--out /users/k1620287/brc_scratch/meta.analysis/analysis/meta/results/$data.$pheno
	}
}
EOT
qsub /users/k1620287/brc_scratch/meta.analysis/analysis/meta/scripts/clump.sh

rm /users/k1620287/brc_scratch/meta.analysis/analysis/meta/scripts/clump.meta.sh
cat <<'EOT'>> /users/k1620287/brc_scratch/meta.analysis/analysis/meta/scripts/clump.meta.sh
#!/bin/sh
#$ -V
#$ -pe smp 8
#$ -l h_vmem=20G
#$ -l h_rt=04:00:00
#$ -S /bin/sh
#$ -m ea
#$ -m be
#$ -o /users/k1620287/brc_scratch/meta.analysis/analysis/meta/output/
#$ -e /users/k1620287/brc_scratch/meta.analysis/analysis/meta/output/

module add bioinformatics/plink2/1.90b3.38
for data in gxtc adult_merged
{
	for pheno in percent zscore
	{
	plink \
	--bfile /users/k1620287/brc_scratch/meta.analysis/analysis/$data/binary/$data.postimputation.maf05.hwe.udid.keeps.udsex \
	--clump /users/k1620287/brc_scratch/meta.analysis/analysis/meta/results/adult.x.child.se.gwama.$pheno.post.snps \
	--clump-kb 250 \
	--clump-p1 1 \
	--clump-p2 1 \
	--clump-r2 0.25 \
	--out /users/k1620287/brc_scratch/meta.analysis/analysis/meta/results/adult.x.child.se.gwama.$pheno.post.snps.$data
	}
}
EOT
qsub /users/k1620287/brc_scratch/meta.analysis/analysis/meta/scripts/clump.meta.sh

rm /users/k1620287/brc_scratch/meta.analysis/analysis/meta/scripts/anot.sh
cat <<'EOT'>> /users/k1620287/brc_scratch/meta.analysis/analysis/meta/scripts/anot.sh
#!/bin/sh
#$ -V
#$ -pe smp 8
#$ -l h_vmem=20G
#$ -l h_rt=04:00:00
#$ -S /bin/sh
#$ -m ea
#$ -m be
#$ -o /users/k1620287/brc_scratch/meta.analysis/analysis/meta/output/
#$ -e /users/k1620287/brc_scratch/meta.analysis/analysis/meta/output/

module add bioinformatics/plink2/1.90b3.38
for data in gxtc adult_merged
{
	for pheno in percent zscore
	{
	plink \
	--annotate /users/k1620287/brc_scratch/meta.analysis/analysis/meta/results/$data.$pheno.clumped \
	ranges=/users/k1620287/brc_scratch/meta.analysis/analysis/meta/scripts/glist_hg19 \
	--border 250 \
	--out /users/k1620287/brc_scratch/meta.analysis/analysis/meta/results/$data.$pheno.clumped.clumped
	}
}
EOT
qsub /users/k1620287/brc_scratch/meta.analysis/analysis/meta/scripts/anot.sh

rm /users/k1620287/brc_scratch/meta.analysis/analysis/meta/scripts/anot.meta.sh
cat <<'EOT'>> /users/k1620287/brc_scratch/meta.analysis/analysis/meta/scripts/anot.meta.sh
#!/bin/sh
#$ -V
#$ -pe smp 8
#$ -l h_vmem=20G
#$ -l h_rt=04:00:00
#$ -S /bin/sh
#$ -m ea
#$ -m be
#$ -o /users/k1620287/brc_scratch/meta.analysis/analysis/meta/output/
#$ -e /users/k1620287/brc_scratch/meta.analysis/analysis/meta/output/

module add bioinformatics/plink2/1.90b3.38
for data in gxtc adult_merged
{
	for pheno in percent zscore
	{
	plink \
	--annotate /users/k1620287/brc_scratch/meta.analysis/analysis/meta/results/adult.x.child.se.gwama.$pheno.post.snps.$data.clumped \
	ranges=/users/k1620287/brc_scratch/meta.analysis/analysis/meta/scripts/glist_hg19 \
	--border 250 \
	--out /users/k1620287/brc_scratch/meta.analysis/analysis/meta/results/adult.x.child.se.gwama.$pheno.post.snps.$data.clumped
	}
}
EOT
qsub /users/k1620287/brc_scratch/meta.analysis/analysis/meta/scripts/anot.meta.sh
