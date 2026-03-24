#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=20g
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



#create path variables
topOPENdir=$OPENdir
locTMPtop=${TMPdir}SNV/

#create directories
mkdir -p $locTMPtop
mkdir -p ${OPENdir}/SNV_ONT/

#load tools
source ${SCRIPTdir}tools


if [[ $SLURM_ARRAY_TASK_ID -gt $nASSEMBLY ]]; then
  EXT=.Siomi
  SLURM_ARRAY_TASK_ID=$(( $SLURM_ARRAY_TASK_ID - $nASSEMBLY ))
  READS=${PacBio_SIOMI}
else
  EXT=
  READS=${ONT_DNA}
fi


if [[ $SLURM_ARRAY_TASK_ID -eq 0 ]]; then
  #if no assembly is provided, use the dm6 genome
  currASSEMBLYseq=${refFASTA}
  currASSEMBLYname="dm6"
  currASSEMBLYflam=""
  assemblyFILE=${assemblyFILE:-${SCRIPTdir}dm6_assembly.txt}
else
  currASSEMBLYline=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${assemblyFILE})
  currASSEMBLYseq=$(echo $currASSEMBLYline |tr ' ' '\t' | cut -f 1)  #remove file extension
  currASSEMBLYname=$(echo $currASSEMBLYline |tr ' ' '\t' | cut -f 2)  #remove file extension
  currASSEMBLYflam=$(echo $currASSEMBLYline |tr ' ' '\t' | cut -f 3)  #remove file extension
fi
locTMP=${locTMPtop}SNV_ONT/${currASSEMBLYname}/
mkdir -p $locTMP


###################################################################################################
THREADS=$(( $SLURM_CPUS_PER_TASK * 2 ))

#map transcripts to contigs
cd ${LOG}SNV/

samtools faidx $currASSEMBLYseq

###################################################################################################
#---------------------------------------------------------------------------------------------------------
#map reads
if [[ ! -s ${locTMP}reads.mapped${EXT}.bam ]]; then
  if [[ $COMPUTING == C ]]; then
    
    sbatch --parsable --wait --job-name=minimap_variations -o "%x.o.%A-%a.txt" -e "%x.e.%A-%a.txt" --cpus-per-task=20 --mem=60g --wrap="
        set -ux
        THREADS=\$(( \$SLURM_CPUS_PER_TASK * 2 ))
        APPTAINERdir=$APPTAINERdir
        TMPdir=$TMPdir
        source ${SCRIPTdir}tools
        minimap2 -Lax map-ont -t \$THREADS --secondary=no $currASSEMBLYseq $READS | samtools sort -T $locTMP -@ \$THREADS -m 750M -O BAM - > ${locTMP}reads.mapped${EXT}.bam

        samtools index -@ \$THREADS ${locTMP}reads.mapped${EXT}.bam

      " 
  else
    minimap2 -Lax map-ont -t $THREADS --secondary=no $currASSEMBLYseq $READS | samtools sort -T $locTMP -@ $THREADS -m 750M -O BAM - > ${locTMP}reads.mapped${EXT}.bam 

    samtools index -@ $THREADS ${locTMP}reads.mapped${EXT}.bam
  fi

fi

#sniffles
if [[ ! -s ${locTMP}SV.sniffles${EXT}.vcf ]]; then
  if [[ $COMPUTING == C ]]; then

    sbatch --wait  --job-name=SV_sniffles -o "%x.o.%A-%a.txt" -e "%x.e.%A-%a.txt" --cpus-per-task=20 --mem=50g --qos=short --time=2:00:00 --wrap="
      set -ux
      THREADS=\$(( \$SLURM_CPUS_PER_TASK * 2 ))
      APPTAINERdir=$APPTAINERdir
      TMPdir=$TMPdir
      source ${SCRIPTdir}tools

      sniffles --input ${locTMP}reads.mapped${EXT}.bam --vcf ${locTMP}SV.sniffles.${EXT}vcf --reference $currASSEMBLYseq --threads \$THREADS --allow-overwrite
      rtg vcfstats ${locTMP}SV.sniffles.${EXT}vcf > ${locTMP}${currASSEMBLYname}.vcfstats${EXT}.txt

      " 
  else
    sniffles --input ${locTMP}reads.mapped${EXT}.bam --vcf ${locTMP}SV.sniffles.vcf --reference $currASSEMBLYseq --threads $THREADS --allow-overwrite
    rtg vcfstats ${locTMP}SV.sniffles.vcf > ${locTMP}${currASSEMBLYname}.vcfstats${EXT}.txt
  fi
fi

if [[ $currASSEMBLYname != OSC_r1.01 ]]; then
  exit
fi


###################################################################################################
#determine if SVs are TE or not
mawk -v OFS="\t"  '
{ 
  if($8 ~"SVTYPE=INS" && length($5)>50) {
    print ">"$3"\n"$5;
  }else{
    if($8 ~"SVTYPE=DEL" && length($4)>50) {
      print ">"$3"\n"$4
    }
  }
}' ${locTMP}SV.sniffles.vcf  | seqkit seq -w 0 -i - > ${locTMP}SVs.fasta

minimap2 --paf-no-hit ${TEconsensus} ${locTMP}SVs.fasta > ${locTMP}variants_aligned_to_TE.paf


#! inverting deletion/insertion as I want to state their status in the OSC genome but for later stages I had to use OSC as the reference genome
echo TYPE noTE TE OSCzygo | tr ' ' '\t' > ${OPENdir}/SNV_ONT/TE_summary.${currASSEMBLYname}.txt


mawk -v OFS="\t" '
{
  print $1,$3,$4,$6,$2,$5
}' ${locTMP}variants_aligned_to_TE.paf | sort -k1,1 -k2,2n  > ${locTMP}variants_aligned_to_TE.bed

bedtools merge -i ${locTMP}variants_aligned_to_TE.bed -c 4,5,6 -o distinct,distinct,distinct > ${locTMP}variants_aligned_to_TE.merge.bed

awk -v OFS="\t" '
{
    X[$1]+=$3-$2
    Y[$1]=$5
    Z[$1][$4]+=1
  
}
END{
  for(i in X) {
    if( X[i] > Y[i]*0.8 ) {
      n=0
      for(j in Z[i]) {
        n++
        if(j ~",") {
          n+=10
        }
      }
      if(n>1) {
        currTE="multipleTE"
      }else{
        currTE=j
      }
    }else{
      currTE="noTE"
    }
    print i,currTE,X[i],Y[i]
  }
}' ${locTMP}variants_aligned_to_TE.merge.bed > ${locTMP}SV_TE.txt


#determine TEs in LOH if file is available
if [[ -s ${OPENdir}/SNV_Illumina/${currASSEMBLYname}_LOH.bed ]]; then

  mawk -v OFS="\t" '
  {
    if($2 !="noTE" && $4 > 5000){
      print $1,$2 
    }
  }' ${locTMP}SV_TE.txt  > ${locTMP}onlyTE-SV.txt

  mawk -v OFS="\t" -v TE_SVs=${locTMP}onlyTE-SV.txt '
    BEGIN{
      while( getline < TE_SVs ) {
        TE[$1]=$2
      }
    }
    {
      if(TE[$3] != "") {
        print $1,$2,$2+1, TE[$3] 
      }
    }' ${locTMP}SV.sniffles.vcf >  ${locTMP}onlyTE-SV.bed

  bedtools intersect -a ${OPENdir}/SNV_Illumina/${currASSEMBLYname}_LOH.bed -b ${locTMP}onlyTE-SV.bed -wao | 
    mawk -v OFS="\t" '
    BEGIN{
      print "TE","COUNT"
    }
    {
      if($5 != ".") {
        X[$8]++
      }
    }
    END{
      for(i in X) {
        print i,X[i]
      }
    }' >  ${OPENdir}/SNV_ONT/${currASSEMBLYname}_LOH_TE.txt
fi

###################################################################################################
#generate data for the chromosome-plot

bigBedToBed $OSCuniqueness_100nt ${locTMP}uniqueness.bed


cat ${locTMP}SV.sniffles.vcf | 
  grep -v '#' |
  grep PASS |  
  mawk -v OFS="\t" -v CHRsize=${CHRsizes} -v TEfile=${locTMP}SV_TE.txt '
  BEGIN{
    print "CHR","POS","VAL","SVlength","TYPE","TE"

    while( getline < CHRsize ) {
      print $1,$2,0,100,"INS","NA"
      print $1,$2,0,100,"DEL","NA"
      print $1,$2,0,10000,"INS","NA"
      print $1,$2,0,10000,"DEL","NA"
    }
    while( getline < TEfile ) {
      TE[$1]=$2
    }
  }
  {
    split($10,X,/:/)
    split($8,Y,/;|=/)

    if(TE[$3] == "") {
      currTE="noTE"
    }else{
      currTE=TE[$3]
    }
    

    if( $8 ~ "SVTYPE=INS"){
      POS=$2+(Y[5]/2)
      TYPE="INS"
    }else{
      POS=$2
      TYPE="DEL"
    }
    print $1,POS,X[4]/(X[3]+X[4]),Y[5],TYPE,currTE
  }'  > ${OPENdir}/SNV_ONT/${currASSEMBLYname}_SVs_along_chromosome.txt


exit


#unblock file
rm -rf ${TMPdir}wait.txt
exit
###################################################################################################
#finish script

#clean up
if [[ $DEBUG == N ]]; then
  rm -rf $locTMP
fi

#report processing time
PROCESSED_TIME=$(echo -e $(date "+%s") $TIME | awk '{ print ($1-$2)/60 }')
echo "map transcripts using BLAT=" ${PROCESSED_TIME} >>${LOG}time-log.txt





exit

