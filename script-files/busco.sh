#!/bin/bash

#SBATCH --cpus-per-task=15
#SBATCH --mem=40g
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
locTMP=${TMPdir}busco/


#create directories
mkdir -p $locTMP

#load tools
source ${SCRIPTdir}tools

THREADS=$(( SLURM_CPUS_PER_TASK * 2 ))
###################################################################################################

#run busco
rm -rf ${OPENdir}busco
#mkdir -p  ${OPENdir}busco
#cd ${OPENdir}busco


mkdir -p ${locTMP}/augustus/config

printf "
  cp -r /augustus/config/ ${locTMP}/augustus/
  AUGUSTUS_CONFIG_PATH="${locTMP}/augustus/config/"
  export AUGUSTUS_CONFIG_PATH="${locTMP}/augustus/config/"
  busco --in $assemblyFASTA -c $THREADS -o busco --out_path ${OPENdir} --mode geno --lineage_dataset drosophila_odb12 --force  " > ${locTMP}busco_run.sh
singularity exec --cleanenv --contain -B /groups -B /scratch -B /scratch-cbe ${APPTAINERdir}busco.app ${locTMP}busco_run.sh


nASSEMBLY=$(wc -l ${assemblyFILE} | cut -f 1 -d ' ')
mkdir $locTMP/input/
rm -rf $locTMP/input/*

seqkit seq  ${refFASTA} > ${locTMP}input/dm6.fa
for currASSEMBLY in $(seq 1 $nASSEMBLY); do
  currASSEMBLYline=$(sed -n ${currASSEMBLY}p ${assemblyFILE})
  currASSEMBLYseq=$(echo $currASSEMBLYline |tr ' ' '\t' | cut -f 1)  #remove file extension
  currASSEMBLYname=$(echo $currASSEMBLYline |tr ' ' '\t' | cut -f 2)  #remove file extension

  seqkit seq  ${currASSEMBLYseq} > ${locTMP}input/${currASSEMBLYname}.fa
done

  if [[  ! -s ${locTMP}${currASSEMBLYname}.fa ]]; then
    seqkit seq --line-width 0 ${currASSEMBLYseq}  > ${locTMP}${currASSEMBLYname}.fa
   fi


singularity exec --cleanenv --contain -B /groups -B /scratch -B /scratch-cbe --env BBTOOLS_JAVA_OPTIONS="-Xmx35g" ${APPTAINERdir}busco.app busco --in ${locTMP}input/ -c $THREADS -o busco --out_path ${OPENdir} --download_path ${locTMP} --mode geno --lineage_dataset drosophila_odb12 --force 
singularity exec --cleanenv --contain -B /groups -B /scratch -B /scratch-cbe --env BBTOOLS_JAVA_OPTIONS="-Xmx35g" ${APPTAINERdir}busco.app busco --plot ${OPENdir}/busco

###################################################################################################
#finish script

#clean up
if [[ $DEBUG == N ]]; then
  rm -rf $locTMP
fi

#report processing time
PROCESSED_TIME=$(echo -e $(date "+%s") $TIME | awk '{ print ($1-$2)/60 }')
echo "map ONT reads=" ${PROCESSED_TIME} >>${LOG}time-log.txt

