#!/bin/bash

#SBATCH --cpus-per-task=20
#SBATCH --mem=50g
#SBATCH --partition=c
#SBATCH -e "%x.e.%A-%a.txt"
#SBATCH -o "%x.o.%A-%a.txt"
#SBATCH --qos=medium


hostname
set -x

###################################################################################################
#extract variables
VARI=$(echo "$1" | sed 's/,/\t/g;s/"//g')
echo $VARI
eval "$VARI"

TIME=$(date "+%s")

###################################################################################################
#setup-phase

#create path variables
locTMP=${TMPdir}TEanalysis-genome/

#create directories
mkdir $locTMP
mkdir -p ${OPENdir}TEanalysis-genome/

#load tools
source ${SCRIPTdir}tools

THREADS=$(( $SLURM_CPUS_PER_TASK * 2 ))
MEM=$(scontrol show job $SLURM_JOBID | grep TRES | awk '{ split($NF, X, /,|=|G/); {print X[5]-5}}' | head -n 1)

###################################################################################################

if [[ $SLURM_ARRAY_TASK_ID -gt 0 ]]; then
  clonalChIPinput=$(echo $clonalChIPinput | tr '~' '\t')
  currFILE=$(echo $clonalChIPinput | tr ' ' '\t' | cut -f $SLURM_ARRAY_TASK_ID)
  currNAME=$(basename $currFILE | cut -d '.' -f 1)
else
  #process genome uniqueness
  nOSC=$(grep -n ${assemblyNAME} ${assemblyFILE} | tr ':' '\t' | cut -f 1)

  currASSEMBLYline=$(sed -n ${nOSC}p ${assemblyFILE})
  ASSEMBLYseq=$(echo $currASSEMBLYline |tr ' ' '\t' | cut -f 1)  #remove file extension
  ASSEMBLYname=$(echo $currASSEMBLYline |tr ' ' '\t' | cut -f 2)  #remove file extension

  #create genome sequences for bowtie
  if [[ ! -s ${locTMP}uniquely_aligned.fa ]]; then
    seqkit sliding -s 1 -W 50 -g ${ASSEMBLYseq} | seqkit seq --min-len 50 | seqkit fx2tab | 
    sort -k2,2 --parallel=${THREADS} -S${MEM}g |
    uniq -c -f 1 | 
    mawk -v OFS="\t" '
    {
      if($1==1){
        print $2, $3
      }
    }' | 
    seqkit tab2fx --line-width 0 |
    gzip > ${locTMP}genomic_reads.fa.gz

    bowtieBuild --threads ${THREADS} --noref ${ASSEMBLYseq} ${locTMP}${ASSEMBLYname}
    bowtie --threads ${THREADS} -fS -m 1 -v 2 -x ${locTMP}${ASSEMBLYname} ${locTMP}genomic_reads.fa.gz | 
    mawk 'BEGIN{OFS="\n"} $1 !~ /^@/ && $3 != "*" {print ">"$1, $10}' > "${locTMP}uniquely_aligned.fa"
  fi
 
  currFILE=${locTMP}uniquely_aligned.fa
  currNAME=genomeUNIQUENESS
fi

#pre-filter reads for uniqueness
bowtie --threads ${THREADS} -fS -m 1 -v 2 -x ${locTMP}${ASSEMBLYname} ${currFILE} | 
  mawk 'BEGIN{OFS="\n"} $1 !~ /^@/ && $3 != "*" {print ">"$1, $10}' > "${locTMP}${currNAME}_uniquely_aligned.fa"

currFILE=${locTMP}${currNAME}_uniquely_aligned.fa

#start analysis
bowtie --threads ${THREADS} -fS -m 1 -v 2 -x ${locTMP}deletion_flanks_combined ${currFILE} | 
samtools view -bS | 
bedtools bamtobed -i - | 
sort -k1,1 -k2,2n --parallel=${THREADS} -S${MEM}g |
mawk -v OFS="\t" -v SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID '
{
  if(SLURM_ARRAY_TASK_ID==0 || $4!~"count="){
    print $0
  } else {
    split($4, splitNAME, /:|=/)
    if($4!~"mapping=m"){
      for(i=1; i<=splitNAME[4]; i++) {
        print $0 
      }
    }
  }
}' > ${locTMP}${currNAME}.expanded.bed

#filter away reads that are only contained in the TE and not anchored in the flanks
mawk -v OFS="\t" -v SEQextension=$SEQextension '
{
  if($1~"~us"){
    if($2<SEQextension-20 ){
      print $0
    }
  }
  if($1~"~ds"){
    if($3> SEQextension+20){
      print $0
    }
  }
  if($1~"~dm6"){
    if( ( $2<SEQextension-15 && $3> SEQextension+15 ) || $1~"~dm6" ){
      print $0
    }
  }
}' ${locTMP}${currNAME}.expanded.bed > ${locTMP}${currNAME}.expanded.filtered.bed


bedtools genomecov -d -i ${locTMP}${currNAME}.expanded.filtered.bed -g ${locTMP}deletion_flanks_combined.size  | 
  mawk -v OFS="\t" '
  {
    if($1~"~us"){
      sub("~us","",$1)
      print $1,"us",$2,$3
    }
    if($1~"~ds"){
      sub("~ds","",$1)
      print $1,"ds",$2,$3
    }
    if($1~"~dm6"){
      sub("~dm6","",$1)
      print $1,"dm6",$2,$3
    }

  }' > ${locTMP}${currNAME}_coverage.bedgraph


mawk -v OFS="\t" -v SVconversion=${locTMP}variants.vcf -v SVannotation=${OPENdir}TEanalysis-genome/TEsummary.${ASSEMBLYname}.txt '
BEGIN{
  while((getline < SVconversion) > 0) {
    if($1 !~ "^#"){
      CONVERSION[$1"_"$2]=$3
    }
  }
  while((getline < SVannotation) > 0) {
    if($1 !~ "^ID") {
      split($1, splitID, /!:!/)
      TEstatus[CONVERSION[splitID[1]]]=$3
      ZYGOSITYstatus[CONVERSION[splitID[1]]]=$6
    }
  }
  print "SV VERSION POS COUNT TEstatus ZYGOSITYstatus"
}
{
  print $0, TEstatus[$1], ZYGOSITYstatus[$1]

}' ${locTMP}${currNAME}_coverage.bedgraph | tr ' ' '\t' > ${OPENdir}TEanalysis-genome/${currNAME}_coverage.bedgraph

################################################################################
#hetTEs
bowtie --threads ${THREADS} -fS -m 1 -v 2 -x ${locTMP}deletion_hetTE_flanks_combined ${currFILE} | 
samtools view -bS | 
bedtools bamtobed -i - | 
sort -k1,1 -k2,2n --parallel=${THREADS} -S${MEM}g |
mawk -v OFS="\t" -v SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID '
{
  if(SLURM_ARRAY_TASK_ID==0 || $4!~"count="){
    print $0
  } else {
    split($4, splitNAME, /:|=/)
    if($4!~"mapping=m"){
      for(i=1; i<=splitNAME[4]; i++) {
        print $0 
      }
    }
  }
}' > ${locTMP}${currNAME}.expanded.hetTE.bed

#filter away reads that are only contained in the TE and not anchored in the flanks
mawk -v OFS="\t" -v SEQextension=$SEQextension '
{
  if($1~"~us"){
    if($2<SEQextension-20 ){
      print $0
    }
  }
  if($1~"~ds"){
    if($3> SEQextension+20){
      print $0
    }
  }
  if($1~"~dm6"){
    if( ( $2<SEQextension-15 && $3> SEQextension+15 ) || $1~"~dm6" ){
      print $0
    }
  }
}' ${locTMP}${currNAME}.expanded.hetTE.bed > ${locTMP}${currNAME}.expanded.hetTE.filtered.bed


bedtools genomecov -d -i ${locTMP}${currNAME}.expanded.hetTE.filtered.bed -g ${locTMP}deletion_hetTE_flanks_combined.size  | 
  mawk -v OFS="\t" '
  {
    if($1~"~us"){
      sub("~us","",$1)
      print $1,"us",$2,$3
    }
    if($1~"~ds"){
      sub("~ds","",$1)
      print $1,"ds",$2,$3
    }
    if($1~"~dm6"){
      sub("~dm6","",$1)
      print $1,"dm6",$2,$3
    }

  }' > ${OPENdir}TEanalysis-genome/${currNAME}_coverage.hetTE.bedgraph


#? mawk -v OFS="\t" -v SVconversion=${locTMP}variants.vcf -v SVannotation=${OPENdir}TEanalysis-genome/TEsummary.${ASSEMBLYname}.txt '
#? BEGIN{
#?   while((getline < SVconversion) > 0) {
#?     if($1 !~ "^#"){
#?       CONVERSION[$1"_"$2]=$3
#?     }
#?   }
#?   while((getline < SVannotation) > 0) {
#?     if($1 !~ "^ID") {
#?       split($1, splitID, /!:!/)
#?       TEstatus[CONVERSION[splitID[1]]]=$3
#?       ZYGOSITYstatus[CONVERSION[splitID[1]]]=$6
#?     }
#?   }
#?   print "SV VERSION POS COUNT TEstatus ZYGOSITYstatus"
#? }
#? {
#?   print $0, TEstatus[$1], ZYGOSITYstatus[$1]
#? 
#? }' ${locTMP}${currNAME}_coverage.hetTE.bedgraph | tr ' ' '\t' > ${OPENdir}TEanalysis-genome/${currNAME}_coverage.hetTE.bedgraph

################################################################################
#TE insertions

  bowtie --threads ${THREADS} -fS -m 1 -v 2 -x ${locTMP}TE_flanks_combined ${currFILE} | 
  samtools view -bS | 
  bedtools bamtobed -i - | 
  sort -k1,1 -k2,2n --parallel=${THREADS} -S${MEM}g |
  mawk -v OFS="\t" -v SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID '
  {
    if(SLURM_ARRAY_TASK_ID==0 || $4!~"count="){
      print $0
    } else {
      split($4, splitNAME, /:|=/)
      if($4!~"mapping=m"){
        for(i=1; i<=splitNAME[4]; i++) {
          print $0 
        }
      }
    }
  }' > ${locTMP}${currNAME}.TE_flanks_expanded.bed


#filter away reads that are only contained in the TE and not anchored in the flanks
mawk -v OFS="\t" -v SEQextension=$SEQextension '
{
  if($1~"~us"){
    if($2<SEQextension-20 ){
      print $0
    }
  }
  if($1~"~ds"){
    if($3> SEQextension+20){
      print $0
    }
  }
  if($1~"~dm6"){
      print $0
  }

}' ${locTMP}${currNAME}.TE_flanks_expanded.bed > ${locTMP}${currNAME}.TE_flanks_expanded.filtered.bed

bedtools genomecov -d -i ${locTMP}${currNAME}.TE_flanks_expanded.filtered.bed -g ${locTMP}TE_flanks_combined.size  | 
  mawk -v OFS="\t" '
  {
    if($1~"~us"){
      sub("~us","",$1)
      print $1,"us",$2,$3
    }
    if($1~"~ds"){
      sub("~ds","",$1)
      print $1,"ds",$2,$3
    }
    if($1~"~dm6"){
      sub("~dm6","",$1)
      print $1,"dm6",$2,$3
    }
  }' > ${locTMP}${currNAME}_TE_flanks_coverage.bedgraph

mawk -v OFS="\t"  -v H3K9classes=${locTMP}TEs_H3K9me3_matrix_region_clusters.txt '
BEGIN{
  while((getline < H3K9classes) > 0) {
    if($1 !~ "^ID") {
      H3K9status[$4]=$6
    }
  }
  print "SV VERSION POS COUNT H3K9cluster"
}
{
  print $0, H3K9status[$1]

}' ${locTMP}${currNAME}_TE_flanks_coverage.bedgraph | tr ' ' '\t' > ${OPENdir}TEanalysis-genome/${currNAME}_TE_flank_coverage.bedgraph
