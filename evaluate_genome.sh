#!/bin/bash


#!#############################################################################################
#!#############################################################################################
#hard-coded things 

STAGEavail="synteny gaep SNP TEcoverage quast VARIANTS_ONT VARIANTS_ILLUMINA busco TEanalysisGENOME clusterAnalysis HiC misc"

#add path to the OSC assembly sequence file here
OSCassembly=

#path to tmp-storage
TMP=""

#!#############################################################################################
#!#############################################################################################
###############################################################################################
set -u

# Argument = -i input -c chunksize -D blast-database -v
usage() {
  cat <<EOF
  usage: $0 options
  
  ###############################################################################
  
  usage: [PATH]/selectUTRs [options] -F 
  
  OPTIONS:
      -h  Show this message
      -A  assembly file
        path to the file that contains all assemblies to be analysed
        [can be left empty to analyse default OSC genome]
      -R reference fasta file to be used for annotation
        [can be left empty to analyse dm6]
      -N  Name for analysis
        [can be left empty = default]
      -O  Results directory
      -Y  include Y chromosome in analysis

      -r  ONT raw DNA read file for aligning to the genome

      -T  TE consensus file - if unset then the AP TE file will be used
      -s  wt-sRNA file for TE coverage evaluation

      -i  Illumina DNA reads for SNP calling
      -S  run only particular stages
            $STAGEavail

      -C  set flag for local processing (use only if multiple cores available)
      -D  sed debug mode - does not trigger git commit and tmp-files not deleted
      -F  force re-generation of bowtie-indexes and other fixed files
      -W  wipe all data and start fresh
EOF
}

#define setup
refFASTA=dm6
assemblyNAME=OSC_r1.01
analysisNAME=default
outPATH=
Ychrom=N

#define paths to input files

#add path to the chromosome sizes file here
CHRsizes=

#path to the nanopore reads used for assembly
ONT_DNA=

#path to the  HiC reads used for scaffolding
HiCreads=

#path to the illumina DNA reads used for polishing and SNP calling
ILLUMINA_DNAseq=

#path to the annotations directory in the OSC genome hub
OSC_annotations=

#path to the RNAseq PE data used for flamenco splicing analysis (normalization factor is for the GEO deposited data)
OSC_RNAseq_PE=
OSC_RNAseq_PE_NORM=2.5774

#path to the 100nt uniqueness track for OSC genome (3MM)
OSCuniqueness_100nt=

#path to the INPUT data from the clonal cell line
clonalChIPinput=

#path to the dm6 reference genome sequence
refFASTAseq=

wt_sRNA=
CHIPdataH3K9=
CHIPdataH3K9_SIENSKI=
CHIPdataH3K9_Saito=
PacBio_SIOMI=
TEconsensus=

#define computing parameters
STAGE=
COMPUTING=C
DEBUG=N
FORCE=N
WIPE=N

###############################################################################################
while getopts ÒhA:R:N:O:Yr:i:T:s:S:CDFW,Ó OPTION; do
  case $OPTION in
  h)
    usage
    exit 1
    ;;
  A)
    assemblyFILE=$OPTARG
    ;;
  R)
    refFASTA=$OPTARG
    ;;
  N)
    analysisNAME=$OPTARG
    ;;
  O)
    outPATH=$OPTARG
    ;;
  Y)
    Ychrom=Y
    ;;
  r)
    ONT_DNA=$OPTARG
    ;;
  i)
    ILLUMINA_DNAseq=$OPTARG
    ;;
  T)
    TEconsensus=$OPTARG
    ;;
  s)
    wt_sRNA=$OPTARG
    ;;
  S)
    STAGE=$OPTARG
    ;;
  C)
    COMPUTING=L
    ;;
  D)
    DEBUG=Y
    ;;
  F)
    FORCE=Y
    ;;
  W)
    WIPE=Y
    ;;
  ?)
    usage
    exit
    ;;
  esac
done

###################################################################################################
#fixed variables


#@!@#
###################################################################################################
#variable-setup



#dates
DATE_OF_DAY=$(date +%F)
FULL_DATE=$(date)
#DATE_OF_DAY=2017-02-14


if [[ -z $outPATH ]]; then
  #usage
  printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
  printf "       Please provide the path to the output directory in option O!\n"
  printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n"
  exit
fi

if [[ $analysisNAME == *"+"* ]]; then
  #usage
  printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
  printf "       no special characters like + allowed in assembly name!\n"
  printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n"
  exit
fi

if [[ ! -z $STAGE ]]; then
  STAGEavail=$(echo $STAGEavail | tr ' ' ',')

  STAGE=$(echo $STAGE | tr ',' '\t')
  for i in $STAGE; do
    if [[ $STAGEavail != *$i* ]];  then
      #usage
      printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
      printf "       Stage $i does not exist!\n"
      printf "$STAGEavail \n"
      printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n"
      exit
    fi

  done
fi

TMPdir="${TMP}${analysisNAME}/"
OPENdir="${outPATH}/${analysisNAME}/"

###################################################################################################
#test if all files required for the executed stages are present

if [[ $STAGE == *synteny* || -z $STAGE ]]; then
  if [[ -z $refFASTA ]]; then
    #usage
    printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    printf "     Please provide refFASTA in option R to run synteny-stage!\n"
    printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n"
    exit
  fi
fi
if [[ $STAGE == *gaep* || -z $STAGE || $STAGE == *VARIANTS_ILLUMINA* ]]; then
  if [[ -z $ONT_DNA ]]; then
    #usage
    printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    printf "     Please provide ONT_DNA in option r to run gaep-stage!\n"
    printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n"
    exit
  fi
  if [[ -z $ILLUMINA_DNAseq ]]; then
    #usage
    printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    printf "     Please provide ILLUMINA_DNAseq in option i to run gaep-stage!\n"
    printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n"
    exit
  fi
fi
if [[ STAGE == "TEcoverage" || -z $STAGE ]]; then
  if [[ -z $TEconsensus ]]; then
    #usage
    printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    printf "     Please provide TEconsensus in option T to run TEcoverage-stage!\n"
    printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n"
    exit
  fi
  if [[ -z $wt_sRNA ]]; then
    #usage
    printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    printf "     Please provide wt_sRNA in option s to run TEcoverage-stage!\n"
    printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n"
    exit
  fi
fi

###################################################################################################
#setup directories
if [[ $WIPE == Y ]]; then
  rm -rf $OPENdir
  rm -rf $TMPdir
fi

if [[ $FORCE == Y ]]; then
  rm -rf ${TMPdir}indeces/
fi

mkdir -p $OPENdir
mkdir -p $TMPdir

if [[ -z $assemblyFILE ]]; then
  #default assembly
  echo $OSCassembly OSC_r1.01 | tr ' ' '\t' > ${TMP}assemblyFILE.txt
  assemblyFILE=${TMP}assemblyFILE.txt
else
  assemblyFASTA=$assemblyFILE
fi

###################################################################################################
#preset scripts

#determine script-location
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$(cd -P "$(dirname "$SOURCE")" && pwd)"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done

#generate final SCRIPT_DIR variables
SCRIPTdir_raw="$(cd -P "$(dirname "$SOURCE")" && pwd)"
SCRIPTdir_raw="${SCRIPTdir_raw}/"
SCRIPTdir="${SCRIPTdir_raw}script-files/"
UTILITYdir="${SCRIPTdir_raw}utility-files/"

#move scripts to TMP-directory
cp -r ${SCRIPTdir} ${TMPdir}
SCRIPTdir=${TMPdir}script-files/
mv ${SCRIPTdir}main.sh ${SCRIPTdir}genomeEVAL_${analysisNAME}.sh

###################################################################################################
#push git version and commit

cd ${SCRIPTdir_raw}
echo ${SCRIPTdir_raw}

if [[ $DEBUG != Y ]]; then

  #ask for commit-message
  while true; do
    read -r -p "Plese specify a commit-message: " msg
    case $msg in
    [Nn]) break ;;
    *)
      commitMESSAGE=$msg
      break
      ;;
    esac
  done

  #commit all changes
  git add .
  if [[ -z $commitMESSAGE ]]; then
    git commit -m "automatic commit on submission"
  else
    git commit -m "$commitMESSAGE"
  fi
  #git push --all
fi

#@ commitID=$(git log -1 --pretty=format:"%h")

###################################################################################################
#download containers

APPTAINERdir=${TMP}apptainer/

mkdir -p ${APPTAINERdir}
cd ${APPTAINERdir}

wget -O basicTools.app https://brenneckelab.imba.oeaw.ac.at/Publication_Data/2026_Handler_OSC-genome/Apptainer/basicTools.app
wget -O genome_evaluation.app https://brenneckelab.imba.oeaw.ac.at/Publication_Data/2026_Handler_OSC-genome/Apptainer/genome_evaluation.app
wget -O ragtag.app https://brenneckelab.imba.oeaw.ac.at/Publication_Data/2026_Handler_OSC-genome/Apptainer/ragtag.app
wget -O quast.app https://brenneckelab.imba.oeaw.ac.at/Publication_Data/2026_Handler_OSC-genome/Apptainer/quast.app
wget -O AP_R.app https://brenneckelab.imba.oeaw.ac.at/Publication_Data/2026_Handler_OSC-genome/Apptainer/AP_R.app
wget -O deepvariant.app https://brenneckelab.imba.oeaw.ac.at/Publication_Data/2026_Handler_OSC-genome/Apptainer/deepvariant.app
wget -O rtg-tools.app https://brenneckelab.imba.oeaw.ac.at/Publication_Data/2026_Handler_OSC-genome/Apptainer/rtg-tools.app
wget -O sniffles.app https://brenneckelab.imba.oeaw.ac.at/Publication_Data/2026_Handler_OSC-genome/Apptainer/sniffles.app
wget -O busco.app https://brenneckelab.imba.oeaw.ac.at/Publication_Data/2026_Handler_OSC-genome/Apptainer/busco.app
wget -O HiC_tools.app https://brenneckelab.imba.oeaw.ac.at/Publication_Data/2026_Handler_OSC-genome/Apptainer/HiC_tools.app
wget -O repeatmasker.app https://brenneckelab.imba.oeaw.ac.at/Publication_Data/2026_Handler_OSC-genome/Apptainer/repeatmasker.app

###################################################################################################
#submit main-run script
LOG=${OPENdir}/LOGs/
#@ rm -rf ${LOG}
mkdir -p ${LOG}
#@ rm -rf ${LOG}FINann_*.txt

#add settings to LOGs
mkdir -p ${LOG}/settings/
cat ${SCRIPTdir_raw}*.sh |
  awk -v RS="#@!@#" '{if (NR==1) print }' >${LOG}/settings/SETTINGS_${commitID}.log

#clear left-over stop-commands
rm -rf ${TMPdir}wait.txt

#cd to log directory for correct log deposition
cd $LOG

#convert variables for submission
ILLUMINA_DNAseq=$(echo $ILLUMINA_DNAseq | tr ',' '~')
ONT_DNA=$(echo $ONT_DNA | tr ',' '~')
STAGE=$(echo $STAGE | tr '\t' '~' | tr ' ' '~')


COMMAND=${SCRIPTdir}genomeEVAL_${analysisNAME}.sh
VARI="OPENdir=${OPENdir},TMPdir=${TMPdir},LOG=${LOG},COMPUTING=${COMPUTING},DEBUG=${DEBUG},FORCE=${FORCE},SCRIPTdir=${SCRIPTdir},UTILITYdir=${UTILITYdir},APPTAINERdir=${APPTAINERdir},refFASTA=${refFASTA},assemblyNAME=${assemblyNAME},assemblyFASTA=${assemblyFASTA},assemblyFILE=${assemblyFILE},analysisNAME=${analysisNAME},CHRsizes=${CHRsizes},ONT_DNA=${ONT_DNA},PacBio_SIOMI=${PacBio_SIOMI},ILLUMINA_DNAseq=${ILLUMINA_DNAseq},STAGE=${STAGE},Ychrom=${Ychrom},TEconsensus=${TEconsensus},wt_sRNA=${wt_sRNA},OSC_RNAseq_PE=${OSC_RNAseq_PE},HiCreads=${HiCreads},OSCuniqueness_100nt=${OSCuniqueness_100nt},CHIPdataH3K9=${CHIPdataH3K9},OSC_annotations=${OSC_annotations},OSC_RNAseq_PE_NORM=${OSC_RNAseq_PE_NORM},clonalChIPinput=${clonalChIPinput},CHIRdataH3K9_SIENSKI=${CHIPdataH3K9_SIENSKI},CHIPdataH3K9_Saito=${CHIPdataH3K9_Saito}"


if [[ $COMPUTING == C ]]; then
  sbatch $COMMAND ${VARI}
else
  if [[ -z ${SLURM_CPUS_PER_TASK+x} ]]; then
    srun --cpus-per-task=10 --mem=20g --time=1:00:00 --qos=short $COMMAND ${VARI}
  else
    $COMMAND ${VARI}
  fi
fi

exit
