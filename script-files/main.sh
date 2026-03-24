#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=20g
#SBATCH -e "%x.e.%j.txt"
#SBATCH -o "%x.o.%j.txt"
#SBATCH --qos=short
#SBATCH --time=1:00:00

hostname
set -ux

###################################################################################################
#extract variables
VARI=$(echo "$1" | sed 's/,/\t/g;s/"//g')
eval "$VARI"
VARI=$1

#report all variables to log
echo $1 | tr ',' '\n'

locTMP=$TMPdir

#initiate time-variable
TIME=$(date "+%s")

#determine CORES and MEM
CORES=$(( $SLURM_CPUS_PER_TASK * 2 ))
MEM=$(scontrol show job $SLURM_JOBID | grep TRES | awk '{ split($NF, X, /,|=|G/);{print X[5]-5}}' | head -n 1)

#load tools 
source ${SCRIPTdir}tools
ID=""

###################################################################################################
#download dm6 and remove Y-chromosome if required
if [[ $refFASTA == dm6 ]]; then
  if [[ ! -s ${TMPdir}dm6.fa ]]; then
    #download genome from flybase
    #@ wget  -O ${TMPdir}dm6.fa.gz https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/current/fasta/dmel-all-chromosome-r6.63.fasta.gz
    cat $refFASTAseq | gzip >  ${TMPdir}dm6.fa.gz

    if [[ ! -s ${TMPdir}dm6.fa.gz ]]; then
      printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
      printf "       Download failed  !\n"
      printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n"
      exit 1
    fi

    #remove Y-chromosome using seqkit if required 
    if [[ $Ychrom == N ]]; then
      seqkit seq -i ${TMPdir}dm6.fa.gz  | seqkit grep -v -r -p "Y"   > ${TMPdir}dm6.fa
      rm -rf ${TMPdir}dm6.fa.gz
    else
      gunzip -c ${TMPdir}dm6.fa.gz | seqkit seq -i > ${TMPdir}dm6.fa
    fi
  fi
  refFASTA=${TMPdir}dm6.fa
  VARI="${VARI},refFASTA=${refFASTA}"
fi

###################################################################################################
#generate synteny plots
if [[ ( -z $STAGE || $STAGE == *synteny* ) ]]; then
  
  #run
  COMMAND="${SCRIPTdir}create_synteny.sh"

  if [[ $COMPUTING == C ]]; then
    newID=$(sbatch --parsable $COMMAND ${VARI})
    ID="${ID}:${newID}"
  else
    $COMMAND ${VARI}
  fi
fi

###################################################################################################
#analyse busco
if [[ ( -z $STAGE || $STAGE == *busco* ) ]]; then
  
  #run
  COMMAND="${SCRIPTdir}busco.sh"

  if [[ $COMPUTING == C ]]; then
    newID=$(sbatch --parsable $COMMAND ${VARI})
    ID="${ID}:${newID}"
  else
    $COMMAND ${VARI}
  fi
fi


###################################################################################################
#prepare paired-end reads from interleaved
if [[ ! -s ${TMPdir}illumina_2.fa && ! -z ${ILLUMINA_DNAseq} ]]; then 
  seqkit fx2tab  ${ILLUMINA_DNAseq}  |  
    mawk -v OFS="\t" '
    {
      
      sub("^>", "", $1)
      #split paired end reads into two files using seqkit to convert with tab2fx
      if (NR%2==1) {
        FW=$0
        sub("/1$", "", FW)
        nFW=NF
      }else{
        RV=$0
        sub("/2$", "", RV)
        nRV=NF
      
        if(nRV==3 && nFW==3){
          print FW > "'${TMPdir}'illumina_1.tab"
          print RV > "'${TMPdir}'illumina_2.tab"
        }
      }
    }'

  seqkit tab2fx ${TMPdir}illumina_1.tab --line-width 0 | sed 's/\/1//g' > ${TMPdir}illumina_1.fa &
  seqkit tab2fx ${TMPdir}illumina_2.tab --line-width 0 | sed 's/\/2//g' > ${TMPdir}illumina_2.fa &
  wait
fi




###################################################################################################
#run quast to evaluate genome

if [[ ( -z $STAGE || $STAGE == *quast* ) ]]; then

  #run gaep
  COMMAND="${SCRIPTdir}quast.sh"

  if [[ $COMPUTING == C ]]; then
    newID=$(sbatch --parsable $COMMAND ${VARI})
    #@ ID="${ID}:${newID}"
  else
    $COMMAND ${VARI}
  fi
fi




###################################################################################################
#run gaep to evaluate genome

if [[ ( -z $STAGE || $STAGE == *gaep* ) ]]; then

  #run gaep
  COMMAND="${SCRIPTdir}gaep.sh"

  nASSEMBLY=$(wc -l ${assemblyFILE} | cut -f 1 -d ' ')
  if [[ $COMPUTING == C ]]; then
    newID=$(sbatch --array=0-$nASSEMBLY --parsable $COMMAND ${VARI})
    #@ ID="${ID}:${newID}"
  else
    $COMMAND ${VARI},SLURM_ARRAY_TASK_ID=1
  fi
fi


#---------------------------------------------------------------------------------------------------------
#deterime TE coverage in the clusters and compare to piRNA covaerage

if [[ ( -z $STAGE || $STAGE == *TEcoverage* ) ]]; then

  #run TE-coverage
  COMMAND="${SCRIPTdir}TE-coverage.sh"

  if [[ $COMPUTING == C ]]; then
    newID=$(sbatch  --parsable $COMMAND ${VARI})
    ID="${ID}:${newID}"
  else
    $COMMAND ${VARI}
  fi
fi



###################################################################################################
#SNP prediction
if [[ ( -z $STAGE || $STAGE == *VARIANTS_ONT* || $STAGE == *VARIANTS_ILLUMINA ) ]]; then
  #prepare log directory
  mkdir $LOG/SNV/
  #@ rm -rf $LOG/SNV/*
  cd $LOG/SNV/
  
  nASSEMBLY=$(wc -l ${assemblyFILE} | cut -f 1 -d ' ')
  nASSEMBLYdouble=$(( $nASSEMBLY * 2 )) 

  #ONT
  if [[ ( -z $STAGE || $STAGE == *VARIANTS_ONT* ) ]]; then
    COMMAND="${SCRIPTdir}predict-SNPs.ONT.sh"

    if [[ $COMPUTING == C ]]; then
      DEPEND=""
#!      DEPEND=$(sbatch --array=0-$nASSEMBLYdouble --parsable $COMMAND ${VARI},nASSEMBLY=$nASSEMBLY )
      DEPEND=$(sbatch --array=0-$nASSEMBLY --parsable $COMMAND ${VARI},nASSEMBLY=$nASSEMBLY )

      #collect all VCFs
      sbatch --dependency=afterany:$DEPEND --parsable --wrap="
        set -ux
        APPTAINERdir=$APPTAINERdir
        TMPdir=$TMPdir
        source ${SCRIPTdir}tools
        export TMPDIR=$locTMP



      "
    else
      $COMMAND ${VARI},SLURM_ARRAY_TASK_ID=1,nASSEMBLY=$nASSEMBLY
      exit
      #@ #collect all VCFs
      source ${SCRIPTdir}functions
      combinVCF ${assemblyFILE} ${TMPdir}SNV/SNV_ONT/
       
    fi


  fi

  #Illumina
  if [[ ( -z $STAGE || $STAGE == *VARIANTS_ILLUMINA* ) ]]; then

    COMMAND="${SCRIPTdir}predict-SNPs.illumina.sh"

    if [[ $COMPUTING == C ]]; then
      nASSEMBLY=1
      sbatch --array=0-$nASSEMBLY --parsable $COMMAND ${VARI}
      #@ sbatch --parsable $COMMAND ${VARI},SLURM_ARRAY_TASK_ID=0
    else
      $COMMAND ${VARI},nASSEMBLY=$nASSEMBLY,SLURM_ARRAY_TASK_ID=1
    #@   source ${SCRIPTdir}functions
    #@  combinVCF ${assemblyFILE} ${TMPdir}SNV/SNV_Illumina/
    fi
  fi
fi



###################################################################################################
#TE analysis genome_vs_genome
if [[ ( -z $STAGE || $STAGE == *TEanalysisGENOME* ) ]]; then
  #run TE analysis
  COMMAND="${SCRIPTdir}TE-analysis.genome.sh"
  if [[ $COMPUTING == C ]]; then
    newID=$(sbatch --parsable $COMMAND ${VARI})
    ID="${ID}:${newID}"
  else
    $COMMAND ${VARI}
  fi
fi  

###################################################################################################
#HiC analysis 
if [[ ( -z $STAGE || $STAGE == *HiC* ) ]]; then
  #run TE analysis
  COMMAND="${SCRIPTdir}HiC.sh"
  if [[ $COMPUTING == C ]]; then
    newID=$(sbatch --parsable $COMMAND ${VARI})
    ID="${ID}:${newID}"
  else
    $COMMAND ${VARI}
  fi
fi  


###################################################################################################
#cluster analysis
if [[ ( -z $STAGE || $STAGE == *clusterAnalysis* ) ]]; then
  #run cluster analysis
  COMMAND="${SCRIPTdir}cluster-analysis.sh"
  if [[ $COMPUTING == C ]]; then
    newID=$(sbatch --parsable $COMMAND ${VARI})
    ID="${ID}:${newID}"
  else
    $COMMAND ${VARI}
  fi
fi  

###################################################################################################
#misc analysis
if [[ ( -z $STAGE || $STAGE == *misc* ) ]]; then
  #run cluster analysis
  COMMAND="${SCRIPTdir}misc.sh"
  if [[ $COMPUTING == C ]]; then
    newID=$(sbatch --parsable $COMMAND ${VARI})
    ID="${ID}:${newID}"
  else
    $COMMAND ${VARI}
  fi
fi  
