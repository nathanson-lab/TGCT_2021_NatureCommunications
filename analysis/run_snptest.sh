#!/bin/bash
# pluta 4/20/17
# script to run association testing
# additive model, caucasian samples

NPARAM=$#

if [ $NPARAM -lt 1 ]
then
	echo ""
	echo "USAGE :: "
	echo "./run_snptest.sh CHR"
	echo "CHR is the chr number"
	echo ""
	exit
fi

CHR=$1


GEN=/project/knathanslab/TECAC/chr${CHR}/HRC/Imputed/chr${CHR}.qc.gen
OUTDIR=/project/knathanslab/TECAC/chr${CHR}/Association/
SAMPLE=/project/knathanslab/TECAC/case_ctl_cauc_3ev.sample
LOG=${OUTDIR}/snptest_chr${CHR}_3ev.log
OUTFILE=${OUTDIR}/snptest_chr${CHR}_3ev.out


# STOP if any of the input files are not found
if [ ! -e $GEN ]
then
	echo "$GEN :: file not found."
	exit 1
fi

if [ ! -e $SAMPLE ]
then
	echo "$SAMPLE :: file not found."
	exit 1
fi

# frequentist
# 1 = Additive
# 2 = Dominant
# 3 = Recessive
# 4 = General
# 5 = Heterozygote

# EVs come from runPCA.sh

# snptest v2.5.2 has some kind of error with the file format that i couldnt resolve
# the exact same file works with v2.5
# 2.5.4 has a library issue
CMD="/home/jpluta/snptest_v2.5/snptest_v2.5     -data $GEN $SAMPLE 
                                                -log $LOG 
                                                -o $OUTFILE
						                                    -cov_names Site EV1 EV2 EV3	
                                                -frequentist 1
                                                -method ml
                                                -pheno phenotype"
eval $CMD



if [ $? -ne 0 ]
then
	echo "$CMD failed. aborting"
	exit 1
fi
