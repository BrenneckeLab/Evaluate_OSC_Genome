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
locTMPtop=${TMPdir}SNV/

#create directories
mkdir $locTMPtop
mkdir -p ${OPENdir}/SNV_Illumina/

#load tools
source ${SCRIPTdir}tools

###################################################################################################
###################################################################################################
THREADS=$(( $SLURM_CPUS_PER_TASK * 2 ))

#map transcripts to contigs
cd ${LOG}SNV/



###################################################################################################
#ONT Clair3

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
locTMP=${locTMPtop}SNV_Illumina/${currASSEMBLYname}/
mkdir -p $locTMP

#---------------------------------------------------------------------------------------------------------
#map reads
TEST=""
TEST=$(seqkit fx2tab $ILLUMINA_DNAseq | head -n 1 | grep "/1" )
if [[ -z $TEST ]]; then 
  PAIRED=N
else
  PAIRED=Y
fi

TEST=""
if [[ -s ${locTMP}reads.mapped.bam ]]; then
  TEST=$(samtools view ${locTMP}reads.mapped.bam | head -n 100 | wc -l)
else
  TEST=0
fi

cp $currASSEMBLYseq ${locTMP}input.fa
assemblyFASTA=${locTMP}input.fa
VARI=${VARI},assemblyFASTA=${locTMP}input.fa


if [[ $TEST -lt 10 ]]; then

    if [[ ! -s ${locTMP}input.fa.sa ]]; then
      bwa index ${assemblyFASTA}
    fi

    if [[ $PAIRED == Y ]]; then
      COMMAND="-p $ILLUMINA_DNAseq | samtools view --threads 3 -F  0x904 -b - | samtools fixmate -m --threads 3  - - | samtools sort -m 3g -T $locTMP --threads 20 - | samtools markdup -r  --threads 5 - ${locTMP}reads.mapped.bam "
      COMMAND="-p $ILLUMINA_DNAseq | samtools sort -m 3G -T $locTMP --threads 20  -O BAM -o ${locTMP}reads.mapped.bam - "
    else
      COMMAND="$ILLUMINA_DNAseq | samtools view --threads 3 -F  0x904 -b -| samtools sort -m 5g -T $locTMP -o ${locTMP}reads.mapped.bam  --threads 5 - "
    fi

  if [[ $COMPUTING == C ]]; then
    
    sbatch --wait --job-name=SNV_ILL_bwa -o "%x.o.%j.txt" -e "%x.e.%j.txt" --cpus-per-task=38 --mem=0g --qos=medium --wrap="
      set -ux
      APPTAINERdir=${APPTAINERdir}
      TMPdir=$TMPdir
      locTMP=${locTMP}
      source ${SCRIPTdir}tools
      THREADS=\$(( \$SLURM_CPUS_PER_TASK * 2 ))
      bwa mem -v 2 -t \${THREADS} ${assemblyFASTA} $COMMAND  
  " 
  else
    eval bwa mem -v 2 -t ${THREADS} ${assemblyFASTA} $COMMAND 
  fi

  samtools index -@ $THREADS ${locTMP}reads.mapped.bam
fi

TEST=""
if [[ -s ${locTMP}reads.mapped.unique.bam ]]; then
  TEST=$(samtools view ${locTMP}reads.mapped.unique.bam | head -n 100 | wc -l)
else
  TEST=0
fi

#filter for unique mappings
if [[ $TEST -lt 10 ]]; then
  cd $locTMP
  samtools view -h -q 1 ${locTMP}reads.mapped.bam | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -b - >${locTMP}reads.mapped.unique.bam
  samtools index -@ $THREADS ${locTMP}reads.mapped.unique.bam
fi


#predict SNPs using deepvariant
if [[ ! -s ${locTMP}/output.vcf.gz ]]; then
  if [[ $COMPUTING == C ]]; then
    sbatch --wait --job-name=SNV_ILL_deepvariant --partition=c  -o "%x.o.%A-%a.txt" -e "%x.e.%A-%a.txt" --cpus-per-task=15 --mem=30g --wrap="
      set -ux
      THREADS=\$(( \$SLURM_CPUS_PER_TASK * 2 - 5 ))
      APPTAINERdir=$APPTAINERdir
      TMPdir=$TMPdir
      source ${SCRIPTdir}tools
      cd $locTMP
      export TMPDIR=$locTMP

      #index fasta
      samtools faidx $assemblyFASTA

      #run deepvariant
      deepvariant /opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref=${assemblyFASTA} --reads=${locTMP}reads.mapped.bam --output_vcf=${locTMP}/output.vcf.gz --output_gvcf=${locTMP}/output.g.vcf.gz --num_shards=\$THREADS
      
      #?@ deepvariant /opt/deepvariant/bin/postprocess_variants --ref ${assemblyFASTA} --infile ${locTMP}/intermediate_results_dir/call_variants_output@\${THREADS}.tfrecord.gz --outfile ${locTMP}/output.vcf.gz --cpus \$THREADS --small_model_cvo_records ${locTMP}/intermediate_results_dir/make_examples_call_variant_outputs.tfrecord@\${THREADS}.gz --gvcf_outfile=${locTMP}/output.g.vcf.gz --nonvariant_site_tfrecord_path ${locTMP}/intermediate_results_dir/gvcf.tfrecord@\${THREADS}.gz

      rtg vcfstats ${locTMP}/output.vcf.gz > ${locTMP}${currASSEMBLYname}.deepvariant.vcfstats.txt

    "
  else
    #deleted due to divergence from the cluster variant - not tested
    echo not coded
    exit 1
  fi 
fi
wait

if [[ $currASSEMBLYname != OSC_r1.01 ]] && [[ $currASSEMBLYname != dm6 ]]; then
  exit
fi

###################################################################################################
#calculate average read coverage

#make 1kb windows
bedtools makewindows -g ${assemblyFASTA}.fai -w 1000 >${locTMP}genome.1kb.bed
bedtools coverage  -sorted  -a ${locTMP}genome.1kb.bed -b ${locTMP}reads.mapped.unique.bam >${locTMP}genome.1kb.cov.txt

mawk -v OFS="\t" '
{
  if($NF > 0.5) print $0
}' ${locTMP}genome.1kb.cov.txt | cut -f 1-3 > ${locTMP}genome.1kb.cov.filtered.txt

bcftools view -i 'GT="1/1"' ${locTMP}output.vcf.gz -Oz -o ${locTMP}hom.vcf.gz
tabix -p vcf ${locTMP}hom.vcf.gz
bcftools view -i 'GT="0/1"' ${locTMP}output.vcf.gz -Oz -o ${locTMP}het.vcf.gz
tabix -p vcf ${locTMP}het.vcf.gz

bedtools intersect -sorted -c \
  -a ${locTMP}genome.1kb.cov.filtered.txt \
  -b ${locTMP}hom.vcf.gz \
  -g ${assemblyFASTA}.fai \
  > ${OPENdir}/SNV_Illumina/${currASSEMBLYname}_SNPs_per_1kb.hom.txt

# Heterozygous SNP counts
bedtools intersect -sorted -c \
  -a ${locTMP}genome.1kb.cov.filtered.txt \
  -b ${locTMP}het.vcf.gz \
  -g ${assemblyFASTA}.fai \
  > ${OPENdir}/SNV_Illumina/${currASSEMBLYname}_SNPs_per_1kb.het.txt

paste ${OPENdir}/SNV_Illumina/${currASSEMBLYname}_SNPs_per_1kb.hom.txt ${OPENdir}/SNV_Illumina/${currASSEMBLYname}_SNPs_per_1kb.het.txt \
 | awk 'BEGIN{OFS="\t"; print "chrom","start","end","hom_count","het_count"} \
        {print $1,$2,$3,$4,$8}' \
 > ${OPENdir}/SNV_Illumina/${currASSEMBLYname}_SNPs_per_1kb.txt

if [[ $currASSEMBLYname != OSC_r1.01 ]]; then
  exit
fi

###################################################################################################
#generate data for the chromosome-plot

bigBedToBed $OSCuniqueness_100nt ${locTMP}uniqueness.bed

gunzip -c ${locTMP}output.vcf.gz | 
  grep -v '#' |
  grep PASS |  
  mawk -v OFS="\t" '{if ( length($4)<5 && length($5)<5) {split($10,X,/:/); print $1,$2-1,$2,X[5] }}' | 
  bedtools intersect -c -b ${locTMP}uniqueness.bed -a - | 
  mawk -v OFS="\t" '
    BEGIN{
      print "CHR","POS","POSend","VAL","unique100nt"
    }
    { 
      if($5>=0) print $1,$2,$3,$4,$5
    }' > ${OPENdir}/SNV_Illumina/${currASSEMBLYname}_SNPs_along_chromosome.txt




###################################################################################################
###################################################################################################
#add to trackDb.txt

#do not continue until file is unblocked by other process
while [[ -f ${TMPdir}wait.txt ]]; do
  sleep 10s
done

#block trackDb from other processes
touch ${TMPdir}wait.txt

#remove old lines if present
awk -v FS="\n" -v RS="\n\n" -v OFS="\t" -v ORS="\n\n" -v NAME="SNV_Illumina_" '
  {
    if( $0 !~ NAME ) print
  }
' ${HUBdir}/trackDb.txt >${TMPdir}trackDB.tmp

mv ${TMPdir}trackDB.tmp ${HUBdir}/trackDb.txt

printf "
  track SNV_Illumina_deepvariant_hom
  shortLabel homSNV_deepvariant_ILL
  longLabel homozygous SNVs called using deepvariant and Illumina reads
  priority 30
  visibility squish
  maxWindowToDraw 10000000
  type vcfTabix
  parent VARIATIONS
  bigDataUrl annotations/VARIATIONS/ILL.deepvariant.hom.vcf.gz

  track SNV_Illumina_deepvariant_het
  shortLabel _deepvariant_ILL
  longLabel heterozygous SNVs called using deepvariant and Illumina reads
  priority 33
  visibility squish
  maxWindowToDraw 10000000
  type vcfTabix
  parent VARIATIONS
  bigDataUrl annotations/VARIATIONS/ILL.deepvariant.het.vcf.gz

" >>${HUBdir}/trackDb.txt


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



###################################################################################################
