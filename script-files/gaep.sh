#!/bin/bash

#SBATCH --cpus-per-task=20
#SBATCH --mem=40g
#SBATCH --partition=c
#SBATCH -e "%x.e.%j.txt"
#SBATCH -o "%x.o.%j.txt"
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
CORES=$(( $SLURM_CPUS_PER_TASK * 2 ))
MEM=$(scontrol show job $SLURM_JOBID | grep TRES | awk '{ split($NF, X, /,|=|G/);{print X[5]-5}}' | head -n 1)


if [[ $SLURM_ARRAY_TASK_ID == 0 ]]; then
  GENOME=${refFASTA}
  GENOMEname="dm6"

else
  currASSEMBLYline=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${assemblyFILE})
  GENOME=$(echo $currASSEMBLYline |tr ' ' '\t' | cut -f 1)  #remove file extension
  GENOMEname=$(echo $currASSEMBLYline |tr ' ' '\t' | cut -f 2)  #remove file extension
fi

APPTAINERdir=$APPTAINERdir
locTMP=${TMPdir}/gaep/${GENOMEname}/
mkdir -p ${locTMP}

source ${SCRIPTdir}tools

#***********************************no commenting above************************************************************

#download lineage if not already present
wget -O ${locTMP}diptera_odb10.2024-01-08.tar.gz https://busco-data.ezlab.org/v5/data/lineages/diptera_odb10.2024-01-08.tar.gz
tar -xvzf  ${locTMP}diptera_odb10.2024-01-08.tar.gz -C ${locTMP}

GAEP pipe -o ${GENOMEname} -d ${locTMP} -r $GENOME -t $CORES  \
  -l ${locTMP}diptera_odb10   \
  --sr1 ${TMPdir}illumina_1.fa --sr2 ${TMPdir}illumina_2.fa \
  --lr ${ONT_DNA} -x ont 


exit
#gaep
GAEP busco -o ${GENOMEname} -d ${locTMP} -r $GENOME -t $CORES  \
  -l ${locTMP}diptera_odb10   \
  --lr ${ONT_DNA} -x ont \
 --sr1 ${TMPdir}illumina_1.fa --sr2 ${TMPdir}illumina_2.fa

cp ${locTMP}*html ${OPENdir}

