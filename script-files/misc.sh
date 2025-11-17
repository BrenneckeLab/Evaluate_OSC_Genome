#!/bin/bash

#SBATCH --cpus-per-task=20
#SBATCH --mem=40g
#SBATCH --partition=c
#SBATCH -e "%x.e.%j.txt"
#SBATCH -o "%x.o.%j.txt"
#SBATCH --qos=short
#SBATCH --time=6:00:00


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
locTMP=${TMPdir}/misc/
mkdir -p ${locTMP}

source ${SCRIPTdir}tools

OPENdir=${OPENdir}/misc/
mkdir -p ${OPENdir}
#***********************************no commenting above************************************************************

###################################################################################################
#determine 3' ends for region
#CG12535
echo X_RagTag 3301558 3355000 CG12535 0 - | tr ' ' '\t' > ${locTMP}CG12535.bed


for GENO in siGFP siPiwi ; do
  if [[ $GENO == siPiwi ]]; then 
    #path to the siPIWI direct RNA seq data
    BAMfile=
  elif [[ $GENO == siGFP ]]; then
    #path to the siGFP direct RNA seq data
    BAMfile=
  fi

  #determine total number of reads
  totalREADS=$(samtools view -c -F 0x904 $BAMfile)

  bedtools bamtobed -i $BAMfile | 
    mawk -v OFS="\t" '{
      if ($6 == "+") {print $1,$3-1,$3,$4,$5,$6}
      else if ($6 == "-") {print $1,$2,$2+1,$4,$5,$6}
    }' |
    bedtools intersect -wao -s -a ${locTMP}CG12535.bed -b  -  | 
    sort -k1,1 -k2,2n > ${locTMP}CG12535.ONT_3prime.${GENO}.bed 

  echo endline |
  cat ${locTMP}CG12535.ONT_3prime.${GENO}.bed - | 
  awk -v OFS="\t" -v TMP=$locTMP -v CHRfile=$CHRsizes -v PEAKextension=25 -v TOTALreads=$totalREADS '
  BEGIN{
    PROCINFO["sorted_in"] = "@val_num_desc"
    while (getline LINE < CHRfile) {
      split(LINE,splitLINE,/\t| /)
      CHRsize[splitLINE[1]]=splitLINE[2]
    }
  }
  {
    ID=$4":!:"$1":"$2":"$3
    CHR=$1
    STRAND=$6
    START=$2
    STOP=$3

    if(oldID==ID){
      # Count ends
      TOTAL+=1
      X[$9]+=1
    }else{
      if(NR>1 && SKIP=="N"){
        for(POS in X){
          if(POS in ASSIGNEDalready){
            a=b
          }else{
            if(POS+0 > PEAKextension){XY=PEAKextension }else{XY=POS }
            if(oldSTRAND=="+"){ YZ=PEAKextension; if(POS-(oldSTART+1) > PEAKextension) {XY=PEAKextension;}else{XY=POS-(oldSTART+1)}}
            if(oldSTRAND=="-"){ XY=PEAKextension; if(oldSTOP-POS > PEAKextension) {YZ=PEAKextension;}else{YZ=oldSTOP-POS}}

            for(i=POS-XY; i<=POS+YZ; i++){
              if(i in ASSIGNEDalready){
                a=b
              }else{
                if(i+0 < CHRsize[oldCHR]){
                  TILEcount+=X[i]
                  ASSIGNEDalready[i]=1
                  if(X[i]>0){
                    COVERcount+=1
                  }
                  n=n+1
                  currREGION[n]=i
                }
              }
            }

            if(oldSTRAND == "sense"){
              COLOR="50,149,124"
              lightCOLOR="198,236,226"
            }else{
              COLOR="149,124,50"
              lightCOLOR="230,217,179"
            }

            print oldCHR,currREGION[1]-1,currREGION[n],TILEcount*1000000/TOTALreads

            TILEcount=0
            COVERcount=0
            delete currREGION
            n=0
          }
        }
      }
      # Reset to new region
      SKIP="N"
      if($9 != "-1"){
        oldID=ID
        oldSTART=START 
        oldSTOP=STOP
        oldCHR=CHR
        oldSTRAND=STRAND
        delete X
        TOTAL=1
        X[$9]=1
        delete ASSIGNEDalready
      }else{
        SKIP="Y"
      }
    }
  }' | LC_COLLATE=C sort -k1,1 -k2,2n > ${locTMP}CG12535.ONT_3prime.${GENO}.minus.depth.bg
  
  bedGraphToBigWig ${locTMP}CG12535.ONT_3prime.${GENO}.minus.depth.bg ${CHRsizes} ${OPENdir}CG12535.ONT_3prime.${GENO}.minus.depth.norm.bw
done


###################################################################################################
#SNP number in SoYb

for GENOME in dm6 OSC_r1.01; do
  if [[ $GENOME == dm6 ]]; then
    SNPfile=${TMPdir}/SNV/SNV_Illumina/dm6/output.vcf.gz
    echo 2L 9996932 9999043 SOYB 0 - | tr ' ' '\t' > ${locTMP}SOYB.bed
    echo 2L 9999110 10001428 SOYB2 0 - | tr ' ' '\t' >> ${locTMP}SOYB.bed
  elif [[ $GENOME == OSC_r1.01 ]]; then
    SNPfile=${TMPdir}/SNV/SNV_Illumina/OSC_r1.01/output.vcf.gz
    echo 2L_RagTag 10508893 10511005 SOYB 0 - | tr ' ' '\t' > ${locTMP}SOYB.bed
    echo 2L_RagTag 10511071 10513391 SOYB2 0 - | tr ' ' '\t' >> ${locTMP}SOYB.bed
  fi

  bedtools intersect -wao -a ${locTMP}SOYB.bed -b ${SNPfile} | sort -k1,1 -k2,2n > ${locTMP}SOYB.SNPs.bed
  
  for i in 1/1 0/1 0/0; do
    X=$(grep $i ${locTMP}SOYB.SNPs.bed | wc -l)
    echo $i $X
  done > ${OPENdir}SOYB.SNPnumber.${GENOME}.txt
done