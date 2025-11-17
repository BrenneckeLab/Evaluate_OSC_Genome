#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=20g
#SBATCH --partition=c
#SBATCH -e "%x.e.%A-%a.txt"
#SBATCH -o "%x.o.%A-%a.txt"
#SBATCH --qos=medium
#SBATCH --time=2-00:00:00


hostname
set -ux

###################################################################################################
#extract variables
VARI=$(echo "$1" | sed 's/,/\t/g;s/"//g')
eval "$VARI"

TIME=$(date "+%s")

###################################################################################################
#setup-phase
locTMP=$TMPdir

#create path variables
topOPENdir=$OPENdir
locTMP=${TMPdir}HiC/

#create directories
mkdir $locTMP

#load tools
source ${SCRIPTdir}tools

###################################################################################################
###################################################################################################
THREADS=$(( $SLURM_CPUS_PER_TASK * 2 ))

###################################################################################################
#ONT Clair3

currASSEMBLYline=$(grep OSC_r1.01  ${assemblyFILE})
currASSEMBLYseq=$(echo $currASSEMBLYline |tr ' ' '\t' | cut -f 1)  #remove file extension
currASSEMBLYname=$(echo $currASSEMBLYline |tr ' ' '\t' | cut -f 2)  #remove file extension
currASSEMBLYflam=$(echo $currASSEMBLYline |tr ' ' '\t' | cut -f 3)  #remove file extension

#---------------------------------------------------------------------------------------------------------
#map reads


cp $currASSEMBLYseq ${locTMP}input.fa
assemblyFASTA=${locTMP}input.fa
VARI=${VARI},assemblyFASTA=${locTMP}input.fa

if [[ ! -s ${locTMP}output.sam ]]; then

    if [[ ! -s ${locTMP}input.fa.sa ]]; then
      bwa index ${assemblyFASTA}
    fi

  if [[ $COMPUTING == C ]]; then
    cd $LOG
    sbatch --wait --job-name=HiC_bwa -o "%x.o.%j.txt" -e "%x.e.%j.txt" --cpus-per-task=38 --mem=0g --qos=medium --wrap="
      set -ux
      APPTAINERdir=${APPTAINERdir}
      TMPdir=$TMPdir
      locTMP=${locTMP}
      source ${SCRIPTdir}tools
      THREADS=\$(( \$SLURM_CPUS_PER_TASK * 2  ))
      mamba activate hic_tools
      
      bwa mem -P -v 2 -t \${THREADS} ${assemblyFASTA}  -p $HiCreads  > ${locTMP}output.sam
  " 
  else
    echo not coded
    exit 1
  fi
fi


if [[ ! -s ${locTMP}output.nodups.pairs.gz ]]; then
  if [[ $COMPUTING == C ]]; then
    cd $LOG
    sbatch --wait --job-name=HiC_pairtools -o "%x.o.%j.txt" -e "%x.e.%j.txt" --cpus-per-task=5 --mem=30g --qos=medium --wrap="
      set -ux
      APPTAINERdir=${APPTAINERdir}
      TMPdir=$TMPdir
      locTMP=${locTMP}
      source ${SCRIPTdir}tools
      THREADS=\$(( \$SLURM_CPUS_PER_TASK * 2 -2 ))

      pairtools parse -c ${CHRsizes} --assembly OSC_r1.01 ${locTMP}output.sam | \
      pairtools sort --tmpdir $locTMP --memory 20G | \
      pairtools dedup \
          --output ${locTMP}output.nodups.pairs.gz \
          --output-dups ${locTMP}output.dups.pairs.gz \
          --output-unmapped ${locTMP}output.unmapped.pairs.gz \
          --output-stats ${locTMP}output.dedup.stats
    "
  fi
fi


cat ${CHRsizes} | grep -E '^(2L|2R|3L|3R|4|X)_RagTag' > ${locTMP}CHRsizes.filtered.txt

cooler makebins ${locTMP}CHRsizes.filtered.txt 1000 > ${locTMP}bins.1k.bed

cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 --assembly OSC_r1.01 ${locTMP}bins.1k.bed ${locTMP}output.nodups.pairs.gz ${locTMP}output.1k.cool

cooler balance --mad-max 5 ${locTMP}output.1k.cool

cooler zoomify --resolutions 1000,2000,5000,10000,25000,50000,100000 --balance --balance-args "--mad-max 5" ${locTMP}output.1k.cool

cooler info ${locTMP}output.1k.mcool

cooler show ${locTMP}output.1k.mcool::resolutions/10000
