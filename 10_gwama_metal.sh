
# Collection of scripts used for meta-analyses

#########################################################

#!/bin/sh
#$ -V
#$ -l h_vmem=20G
#$ -l h_rt=04:00:00
#$ -S /bin/sh
#$ -m ea
#$ -m be
#$ -o /mnt/lustre/groups/ukbiobank/usr/crayner/output/
#$ -e /mnt/lustre/groups/ukbiobank/usr/crayner/err/

source $1

cd $working_dir

module add bioinformatics/metal/2011-03-25

metal << EOT
SCHEME STDERR
AVERAGEFREQ ON

CUSTOMVARIABLE N
LABEL N as N

MARKER SNP
ALLELE A1 A2
EFFECT BETA
PVALUE P
STDERRLABEL SE
FREQLABEL MAF

PROCESS $sumstats_1 # adult_anx_ki
PROCESS $sumstats_2 # adult_anx_gxte
PROCESS $sumstats_3 # adult_anx_wb

OUTFILE $output_1   # adult_anx_meta
ANALYZE HETEROGENEITY

QUIT

EOT

#########################################################

#!/bin/sh
#$ -V
#$ -l h_vmem=20G
#$ -l h_rt=04:00:00
#$ -S /bin/sh
#$ -m ea
#$ -m be
#$ -o /mnt/lustre/groups/ukbiobank/usr/crayner/output/
#$ -e /mnt/lustre/groups/ukbiobank/usr/crayner/err/

source $1

cd $working_dir

module add bioinformatics/metal/2011-03-25

metal << EOT
SCHEME STDERR
AVERAGEFREQ ON

CUSTOMVARIABLE N
LABEL N as N

MARKER MarkerName
ALLELE Allele1 Allele2
EFFECT Effect
PVALUE P-value
STDERRLABEL StdErr
PROCESS $output_1 # adult_anx_meta

MARKER SNP
ALLELE A1 A2
EFFECT BETA
PVALUE P
STDERRLABEL SE
FREQLABEL MAF
PROCESS $sumstats_3 # adult_dep_ki

OUTFILE $output_2 # adult_anxdep_meta
ANALYZE HETEROGENEITY
QUIT
EOT

#########################################################

#!/bin/sh
#$ -V
#$ -l h_vmem=20G
#$ -l h_rt=04:00:00
#$ -S /bin/sh
#$ -m ea
#$ -m be
#$ -o /mnt/lustre/groups/ukbiobank/usr/crayner/output/
#$ -e /mnt/lustre/groups/ukbiobank/usr/crayner/err/

source $1

cd $working_dir

module add bioinformatics/metal/2011-03-25

metal << EOT
SCHEME STDERR
AVERAGEFREQ ON

CUSTOMVARIABLE N
LABEL N as N

MARKER SNP
ALLELE A1 A2
EFFECT BETA
PVALUE P
STDERRLABEL SE
FREQLABEL MAF
PROCESS $sumstats_4 # child_anx_meta

MARKER MarkerName
ALLELE Allele1 Allele2
EFFECT Effect
PVALUE P-value
STDERRLABEL StdErr
PROCESS $output_1 # adult_anx_meta

OUTFILE $output_3
ANALYZE HETEROGENEITY
QUIT
EOT 
