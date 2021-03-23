#!/bin/sh
# pluta 3/23/18

# for checking independence against multiple snps
# modified version of conditionalAnalysis.sh and plotCond.r
# pluta 3/21/18
NPARAM=$#

if [ $NPARAM -lt 2 ]
then
	echo ""
	echo "USAGE :: "
	echo "conditionalAnalysis.sh CONDSNPLST WINDOW"
	echo "SNPLST: list of snsp to condition on SIMULTANEOUSLY, in chr:bp for
mat"
	echo "eg., if the list is snpA and snpB, script will check for independe
nce from both in a joint model."
	echo "this assume snpA and snpB are independent."
	echo ""
	echo "WINDOW: flanking length in kbp"
	echo ""
	exit 
fi


# in chr:pos format
CONDSNPLST=$1
WINDOW=$2

### CONSTANTS ###
TECACDIR="/project/knathanslab/TECAC"
METATBL="${TECACDIR}/Meta/TECAC_META_8site1.tbl"
COJOFILE="${TECACDIR}/Meta/GCTA-input3.ma"
#################

echo "1. Set up variables..."

# the default file to store snp data
# remove the output from previous runs
if [ -e snps-pvals ]
then
	rm snps-pvals
fi

# read snplist into an array variable
IFS=$'\n' read -d '' -r -a snparr < $CONDSNPLST

# name files based on the first snp in the list
CHR=$(echo ${snparr[0]} | tr ":" " " | awk '{print $1}')
BP=$(echo ${snparr[0]} | tr ":" " " | awk '{print $2}')


BFILE=${TECACDIR}/chr${CHR}/HRC/Imputed/chr${CHR} 

# get number of snps in snplist
nsnp=$(echo ${snparr[@]} | awk '{print NF}')

# get the p-value of the reference snp in the meta analysis
# but first have to make sure it exists in the data
# vector version of the same process in conditionalAnalysis.sh
for ((i=0; i < ${nsnp}; i++))
do

	tmp=$(LC_ALL=C grep -w "${snparr[$i]}" $METATBL)
	
	if [ $? != 0 ]
	then
		echo "${snparr[$i]} not found in $METATBL."
		echo "exiting"
		exit 1
	fi

	PVAL=$(echo $tmp | awk '{print $8}')
	echo "${snparr[$i]} $PVAL" >> snps-pvals
done

echo "done!"
echo ""


# the window should remain on the reference snp
# TO DO: simplify this to reflect edit
BPlo=$(echo ${snparr[0]} | tr ":" " " | awk '{print $2}')
BPhi=$(echo ${snparr[0]} | tr ":" " " | awk '{print $2}')



echo "2. Extracting snp list..."
# extract the list of snps with WINDOW flanking the reference snp
awk -v bphi="$BPhi" -v bplo="$BPlo" -v win="$WINDOW" '$4 > bplo - win * 1000 && 
$4 < bphi + win * 1000 {print $2}' ${BFILE}.bim > SNPLST

echo "done!"
echo ""

OUTNAME=${CHR}_${BP}_${nsnp}_joint

if [ $CHR == "X" ]
then
	CHR=23
fi


echo "3. Conditional analysis..."
CMD="/home/jpluta/gcta64 --bfile $BFILE
			 --chr $CHR 
			 --cojo-file $COJOFILE 
			 --cojo-cond $CONDSNPLST 
			 --cojo-collinear 0.8
			 --extract SNPLST
			 --out $OUTNAME"

eval $CMD > ${OUTNAME}_cojo.log

if [ $? != 0 ]
then
	echo "$CMD failed!"
	echo "exiting"
	exit 1
fi

echo "done!"
echo ""





echo "5. Plotting data..."
CMD="Rscript plotCondJoint.R snps-pvals ${OUTNAME}.cma.cojo"

eval $CMD > ${OUTNAME}_jt_Rplot.log
if [ $? != 0 ]
then
	echo "$CMD failed!"
	echo "exiting"
	exit 1
fi

echo "done!"

echo ""
echo "conditional analysis - successfully complete."

rm SNPLST
