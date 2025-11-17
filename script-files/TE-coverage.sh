#!/bin/bash

#SBATCH --cpus-per-task=20
#SBATCH --mem=40g
#SBATCH --partition=c
#SBATCH -e "%x.e.%j.txt"
#SBATCH -o "%x.o.%j.txt"
#SBATCH --qos=short
#SBATCH --time=8:00:00


hostname
set -ux

###################################################################################################
#extract variables
VARI=$(echo "$1" | sed 's/,/\t/g;s/"//g')
eval "$VARI"

TIME=$(date "+%s")

###################################################################################################
#setup-phase
CORES=$(( $SLURM_CPUS_PER_TASK * 2 ))
MEM=$(scontrol show job $SLURM_JOBID | grep TRES | awk '{ split($NF, X, /,|=|G/);{print X[5]-5}}' | head -n 1)

APPTAINERdir=$APPTAINERdir
locTMP=${TMPdir}/TEcoverage/
mkdir -p ${locTMP}

source ${SCRIPTdir}tools

#***********************************no commenting above************************************************************
#build TE bowtie index
if [[ ! -s ${locTMP}TEindex.1.ebwt ]]; then 
  bowtieBuild ${TEconsensus} ${locTMP}TEindex
fi


nOSC=$(grep -n ${assemblyNAME} ${assemblyFILE} | tr ':' '\t' | cut -f 1)

currASSEMBLYline=$(sed -n ${nOSC}p ${assemblyFILE})
ASSEMBLYseq=$(echo $currASSEMBLYline |tr ' ' '\t' | cut -f 1)  #remove file extension
ASSEMBLYname=$(echo $currASSEMBLYline |tr ' ' '\t' | cut -f 2)  #remove file extension
ASSEMBLYflam=$(echo $currASSEMBLYline |tr ' ' '\t' | cut -f 3)  #remove file extension




#determine sRNA coverage on the TEs
seqkit fx2tab ${wt_sRNA} | 
  mawk -v OFS="\t" '{if(length($2)>23) print $1,substr($2,1,25) }' |
  seqkit tab2fx --line-width  0|
  bowtie -fS -p $CORES -a -M 1 --best --strata -x ${locTMP}TEindex -  | 
  samtools view -bS - | 
  samtools sort -T ${locTMP} -@ $CORES - | 
  bedtools bamtobed -i - | 
  mawk -v OFS="\t" '
  {
    split($4,splitNAME,/:|=/)
    for(i=1; i<=splitNAME[4]; i++) print 
  }'> ${locTMP}TEsRNA.bed

bedtools genomecov -strand + -d -i ${locTMP}TEsRNA.bed -g ${TEconsensus}.fai | LC_COLLATE=C sort -k1,1 -k2,2n |  mawk -v OFS="\t" '{print $1":"$2,$3}' > ${locTMP}TEcoverage.sRNA.sense.bedgraph
bedtools genomecov -strand - -d -i ${locTMP}TEsRNA.bed -g ${TEconsensus}.fai | LC_COLLATE=C sort -k1,1 -k2,2n | mawk -v OFS="\t" '{print $1":"$2,-$3}' > ${locTMP}TEcoverage.sRNA.antisense.bedgraph

BGname=TEcoverage.sRNA.sense.bedgraph
BGname="${BGname} TEcoverage.sRNA.antisense.bedgraph"

#determine cluster coverage on the TEs
for GENOME in $ASSEMBLYseq $refFASTA ; do
  NAME=$(basename $GENOME .fa)
  echo $NAME

  #generate 25mers from the cluster regions defined in the cluster-bed
  bedtools getfasta -fi $GENOME -bed ${UTILITYdir}cluster-coordinates.${NAME}.bed -fo - | 
    seqkit sliding -s 1 -W 25 - > ${locTMP}${NAME}.25mers.fa

  #align the 25mers to the TE consensus and generate bedgraph
  bowtie -fS -p $CORES -v 1 -a -x ${locTMP}TEindex ${locTMP}${NAME}.25mers.fa  | 
    samtools view -bS - | 
    samtools sort -T ${locTMP} -@ $CORES - | 
    bedtools bamtobed -i -  > ${locTMP}TEcoverage.${NAME}.bedgraph

  bedtools genomecov -strand + -d -i ${locTMP}TEcoverage.${NAME}.bedgraph -g ${TEconsensus}.fai | LC_COLLATE=C sort -k1,1 -k2,2n | mawk -v OFS="\t" '{print $1":"$2,$3}'> ${locTMP}TEcoverage.${NAME}.sense.bedgraph
  bedtools genomecov -strand - -d -i ${locTMP}TEcoverage.${NAME}.bedgraph -g ${TEconsensus}.fai | LC_COLLATE=C sort -k1,1 -k2,2n | mawk -v OFS="\t" '{print $1":"$2,-$3}' > ${locTMP}TEcoverage.${NAME}.antisense.bedgraph
  
  BGname="${BGname} TEcoverage.${NAME}.sense.bedgraph"
  BGname="${BGname} TEcoverage.${NAME}.antisense.bedgraph"  
done

cut -f 1 ${locTMP}TEcoverage.sRNA.sense.bedgraph > ${locTMP}TEcoverage.merged.txt
for BG in ${BGname}; do
  echo $BG
  join ${locTMP}TEcoverage.merged.txt ${locTMP}${BG} > ${locTMP}TEcoverage.merged.tmp
  mv ${locTMP}TEcoverage.merged.tmp ${locTMP}TEcoverage.merged.txt
done

echo CHR:POS $BGname | tr ' ' '\t' > ${OPENdir}TEcoverage.merged.txt
cat ${locTMP}TEcoverage.merged.txt | tr ' ' '\t' >> ${OPENdir}TEcoverage.merged.txt
#plot individual TE coverage


#quantify accuracy of TE coverage