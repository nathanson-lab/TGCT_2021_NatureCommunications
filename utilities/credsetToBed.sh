#pluta 5/8/19

#!/bin/bash

# convert a credible set to bed file format
CREDFILE=$1

SNP=$( basename $CREDFILE "_credSet.txt")
BP=$( echo $SNP | tr "_" " " | awk '{print $2}')
tail -n +2 $CREDFILE | awk -v snp="$SNP" -v bp="$BP" '{print $1, $3, $3+=1, snp,
 bp, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15}' | tr " " "\t" > tmp
echo "Chr Start End RefSNP RefBP refA freq b se p n freq_geno bC bC_se pC r2 ind
p" | tr " " "\t" > header
cat header tmp > ${SNP}.bed
rm tmp
rm header
