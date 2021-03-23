Annotation Pipeline

1. addFunctionalAnnotation.R: determines which features overlap with CRVs
2. paintor-extractLDMatrix.sh: extracts snps in LD with the ref snp and retains LD values. removes any duplicate SNPs or SNPs with LD = 1 (a singular LD matrix will cause PAINTOR to crash). invokes removeRedundantSnps.R
3. PAINTOR-modelSelect.sh: runs univariate association of each annotated feature. Check these for independence, then create PAINTOR model. invokes LRT.BF.R.
