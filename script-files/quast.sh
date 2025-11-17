#!/bin/bash

#SBATCH --cpus-per-task=20
#SBATCH --mem=60g
#SBATCH --partition=c
#SBATCH -e "%x.e.%j.txt"
#SBATCH -o "%x.o.%j.txt"
#SBATCH --qos=medium
#SBATCH --time=1-00:00:00


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
locTMP=${TMPdir}/quast/
mkdir -p ${locTMP}

source ${SCRIPTdir}tools

#***********************************no commenting above************************************************************

if [[ ! -s ${locTMP}flybase.gff ]]; then
  wget -O ${locTMP}flybase.gff.gz https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/dmel_r6.64_FB2025_03/gff/dmel-all-r6.64.gff.gz
  if [[ ! -s ${locTMP}flybase.gff.gz ]]; then
    printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    printf "       Download failed  !\n"
    printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n"
    exit 1
  fi

  gunzip ${locTMP}flybase.gff.gz

fi

#assemble all assemblies
allASSEMBLIES=$(cat ${assemblyFILE} | cut -f 1 | tr '\n' ' ')
allASSEMBLIES="${refFASTA} ${allASSEMBLIES}"  #add reference genome to assemblies
#minimal quast

quast --threads $CORES -o ${OPENdir}quast_quick_module/ --circos --large --features ${locTMP}flybase.gff --split-scaffolds --conserved-genes-finding --eukaryote $allASSEMBLIES  


#assemble all assemblies
allASSEMBLIES=$(cat ${assemblyFILE} | cut -f 1 | tr '\n' ' ')

quast --pe1 ${TMPdir}illumina_1.fa --pe2 ${TMPdir}illumina_2.fa --nanopore ${ONT_DNA} --pacbio ${PacBio} --threads $CORES -o ${OPENdir}quast/ --circos --large --features ${locTMP}flybase.gff --split-scaffolds -r $refFASTA --conserved-genes-finding --eukaryote $allASSEMBLIES  
quast --pe12 ${ILLUMINA_DNAseq} --threads $CORES -o ${OPENdir}quast_onlyILL/ --circos --large --features ${locTMP}flybase.gff --split-scaffolds -r $refFASTA --conserved-genes-finding --eukaryote $allASSEMBLIES --debug

#--nanopore ${ONT_DNA} --pacbio ${PacBio} 
#--pe12 ${ILLUMINA_DNAseq} 