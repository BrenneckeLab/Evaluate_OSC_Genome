#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --mem=20g
#SBATCH --partition=c
#SBATCH -e "%x.e.%A-%a.txt"
#SBATCH -o "%x.o.%A-%a.txt"
#SBATCH --qos=short


hostname
set -x

###################################################################################################
#extract variables
VARI=$(echo "$1" | sed 's/,/\t/g;s/"//g')
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
###################################################################################################
#compare assemblies to reference genome and determine TE content in SVs

#define variables
cd ${LOG}


refNAME="dm6"


nOSC=$(grep -n ${assemblyNAME} ${assemblyFILE} | tr ':' '\t' | cut -f 1)

currASSEMBLYline=$(sed -n ${nOSC}p ${assemblyFILE})
ASSEMBLYseq=$(echo $currASSEMBLYline |tr ' ' '\t' | cut -f 1)  #remove file extension
ASSEMBLYname=$(echo $currASSEMBLYline |tr ' ' '\t' | cut -f 2)  #remove file extension



###################################################################################################
#determine SV differences between OSC and dm6

#align assembly to reference genome
minimap2 -a -x asm5 --cs -r2k -t ${THREADS} $ASSEMBLYseq $refFASTA |
  samtools sort -o ${locTMP}alignments.sorted.bam
samtools index ${locTMP}alignments.sorted.bam

#determine SVs
source ~/.bashrc
mamba activate svimasm_py39
svim-asm haploid ${locTMP} ${locTMP}alignments.sorted.bam $ASSEMBLYseq

cp ${locTMP}sv-lengths.png ${OPENdir}TEanalysis-genome/sv-lengths.${ASSEMBLYname}.png
rtg vcfstats ${locTMP}/variants.vcf > ${OPENdir}TEanalysis-genome/sv-stats.${ASSEMBLYname}.txt


###################################################################################################
#determine if the SVs are homozygous or heterozygous in the OSC genome


sniffles --input ${TMPdir}/SNV/SNV_ONT/OSC_r1.01/reads.mapped.bam --genotype-vcf ${locTMP}/variants.vcf --vcf ${locTMP}SV.sniffles.OSC.vcf --reference $ASSEMBLYseq --threads $THREADS --allow-overwrite 
rtg vcfstats ${locTMP}SV.sniffles.OSC.vcf > ${OPENdir}TEanalysis-genome/sv-stats.zygosity_in_OSC.txt


###################################################################################################
#determine TE fraction in SVs and plot the results

#determine if SVs are TE or not
mawk -v OFS="\t" -v INFILE=${locTMP}SV.sniffles.OSC.vcf '
BEGIN{
  while((getline < INFILE) > 0) {
    if($1 !~ "^#") {
      split($NF, splitTAG, ":")
      ZYGO[$3]=splitTAG[1]
    }

  }
}
{ 
  if($8 ~"SVTYPE=INS" && length($5)>50) {
    print ">"$1"_"$2"!:!LENGTH="length($5)"!:!insertion:!:"ZYGO[$3]"\n"$5;
  }else{
    if($8 ~"SVTYPE=DEL" && length($4)>50){
      print ">"$1"_"$2"!:!LENGTH="length($4)"!:!deletion:!:"ZYGO[$3]"\n"$4
    }
  }
}' ${locTMP}variants.vcf |
seqkit seq -w 0 -i - > ${locTMP}SVs.fasta

minimap2 --paf-no-hit ${TEconsensus} ${locTMP}SVs.fasta > ${locTMP}variants_aligned_to_TE.paf

#! inverting deletion/insertion as I want to state their status in the OSC genome but for later stages I had to use OSC as the reference genome

mawk -v OFS="\t" '
{
  print $1,$3,$4,$6,$2,$5
}' ${locTMP}variants_aligned_to_TE.paf | sort -k1,1 -k2,2n  > ${locTMP}variants_aligned_to_TE.bed

bedtools merge -i ${locTMP}variants_aligned_to_TE.bed -c 4,5,6 -o distinct,distinct,distinct > ${locTMP}variants_aligned_to_TE.merge.bed


awk -v OFS="\t" '
BEGIN{
  print "ID TYPE TE TElength SVlength OSCzygo"
}
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

    if(i~"0/0"){
      ZYGO="hom"
    }else{
      ZYGO="het"
    }

    if(i~"deletion"){
      TYPE="insertion"
    }else{
      TYPE="deletion"
    }

    print i,TYPE,currTE,X[i],Y[i],ZYGO
  }
}' ${locTMP}variants_aligned_to_TE.merge.bed  | tr -s ' ' | tr ' ' '\t' > ${OPENdir}TEanalysis-genome/TEsummary.${ASSEMBLYname}.txt



#plot in R
ml build-env/f2022
ml  r/4.5.1-gfbf-2023b
  R -e "
    library(tidyverse)
    library(cowplot)
    library(ggalluvial)
    theme_set(theme_cowplot())

  RAW = read_tsv('${OPENdir}TEanalysis-genome/TEsummary.${ASSEMBLYname}.txt', col_names = TRUE)%>%
    mutate(
      TEclass = case_when(
        TE == 'noTE' ~ 'noTE',
        TRUE ~ 'TE'
      )
    )

  p = RAW %>%
    select(TEclass, ID, TYPE)%>%
    group_by(TEclass, TYPE) %>%
    summarise(COUNT = n()) %>%
    ggplot(aes(axis1=TYPE, axis2=TEclass, y = COUNT, fill = TYPE)) +
    geom_alluvium() +
    geom_stratum() +
    scale_y_continuous(breaks = seq(0, 5500, by = 1000)) +
    scale_x_discrete(expand = c(.1, .1))+
    geom_text(stat = 'stratum', aes(label = after_stat(stratum)), size = 3, na.rm = TRUE) +
    scale_fill_manual(
      values = c('insertion' = '#081626', 'deletion' = '#e69f00'),
      name = 'TE type'
    ) +
    labs(
      title = 'TE analysis of SVs in ${ASSEMBLYname}',
      x = 'TE class',
      y = 'Number of SVs'
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = 'bottom',
      plot.title = element_text(hjust = 0.5, size = 16, face = 'bold')
    ) 

    ggsave(p, filename = '${OPENdir}TEanalysis-genome/TEsummary.${ASSEMBLYname}.pdf', width = 4, height = 8, dpi = 300)


  p = RAW %>%
    mutate(
      SVlength = abs(SVlength),
      SVlength = ifelse(SVlength > 12000, 12000, SVlength),
      ALPHA = case_when(
        TEclass == 'noTE' ~ 0.3,
        TRUE ~ 1
      )
    ) %>%
    ggplot(aes(x = SVlength, fill = TYPE, alpha=TEclass)) +
    # Insertions (positive y)
    geom_histogram(
      data = ~filter(., TYPE == 'insertion'),
      bins = 50, color = 'black', linewidth = 0.1, position = 'stack'
    ) +
    # Deletions (negative y)
    geom_histogram(
      data = ~filter(., TYPE == 'deletion'),
      aes(y = -after_stat(count)),
      bins = 50, color = 'black', linewidth = 0.1, position = 'stack'
    ) +
    scale_fill_manual(
      values = c('insertion' = '#081626', 'deletion' = '#e69f00'),
      name = 'TE type'
    ) +
    scale_alpha_manual(values=c('noTE'=0.3, 'TE'=1)) +
    scale_x_log10() +
    theme_cowplot(14) +
    theme(
      panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.3),
    ) +
    labs(
      x = 'SV length (capped at 12 kb)',
      y = 'Count (insertions positive, deletions negative)'
    )

    ggsave(p, filename = '${OPENdir}TEanalysis-genome/SV-histogram.${ASSEMBLYname}.pdf', width = 16, height = 4, dpi = 300)

  p = RAW %>% 
    filter(TE != 'noTE', SVlength > 1000) %>%
    group_by(TE,TYPE) %>%
    summarise(COUNT = n()) %>%
    mutate(
      COUNT = case_when(
        TYPE == 'deletion' ~ -COUNT,
        TRUE ~ COUNT
      )
    )%>%
    group_by(TE) %>%
    mutate(ORDERVAL = max(ifelse(TYPE == 'insertion', COUNT, 0))) %>%
    ungroup() %>%
    mutate(TE = fct_reorder(TE, ORDERVAL, .desc = FALSE)) %>%
    ggplot(aes(x = TE, y = COUNT, fill=TYPE)) +
      geom_bar(stat = 'identity', width = 0.5, position = 'stack') +
      coord_flip() +
      scale_fill_manual(
        values = c('insertion' = '#081626', 'deletion' = '#e69f00'),
        name = 'TE type'
      ) +
    theme_cowplot(14) +
    theme(
      panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.3),
      axis.title.y = element_blank(),
      legend.position.inside = c(0.8, 0.8)
    ) +
    labs(
      x = 'number of SVs containing indicated TE'
    )

    ggsave(p, filename = '${OPENdir}TEanalysis-genome/TE-splitup.${ASSEMBLYname}.pdf', width = 9, height = 15, dpi = 300)

"

###################################################################################################
#analyze TE insertions for H3K9me3

#determine full-length TE insertions
minimap2 -x map-ont -N 1000 -p 0.5 -k 15 -w 2 --end-bonus 100 --secondary=yes -P -t ${THREADS}  $ASSEMBLYseq $TEconsensus > ${locTMP}TEs_in_genome.paf

mawk -v OFS="\t"  '{
  if($3<100 && $4>$2-100 && $10 > $2*0.8){
    print 
  }
}' ${locTMP}TEs_in_genome.paf | sort -k1,1 -k2,2n > ${locTMP}TEs_in_genome.filtered.paf

wc -l ${locTMP}TEs_in_genome.filtered.paf 


blat $ASSEMBLYseq $TEconsensus ${locTMP}TEs_in_genome.psl \
    -tileSize=11 -minIdentity=80 -minScore=100 -out=psl



mawk -v OFS="\t" -v CHRsizes=${CHRsizes} '
BEGIN{
  while(getline < CHRsizes){
    chrSizes[$1]=$2
  }
  WINDOWsize=5000
}{
  # PSL fields:

  qStart=$12
  qEnd=$13
  qSize=$11
  tStart=$16
  tEnd=$17
  matches=$1

  if(qStart < 50 && qEnd > qSize - 50 && matches > qSize * 0.8 && tEnd-tStart < qSize*1.2){
    if(tStart>WINDOWsize && tEnd < chrSizes[$14]-WINDOWsize){
      print $14,tStart-WINDOWsize,tStart,$10":!:"$14"_"tStart"_"tEnd"_us",0,$9
      print $14,tEnd,tEnd+WINDOWsize,$10":!:"$14"_"tStart"_"tEnd"_ds",0,$9
    }
  }
}' ${locTMP}TEs_in_genome.psl | sort -k1,1 -k2,2n > ${locTMP}TEs_in_genome.filtered.bed

bigWigToBedGraph ${CHIPdataH3K9}siLuc_H3K9me3_uniq.bw ${locTMP}siLuc_H3K9me3_uniq.bedgraph

bedtools coverage \
  -a ${locTMP}TEs_in_genome.filtered.bed \
  -b ${locTMP}siLuc_H3K9me3_uniq.bedgraph \
  -d > ${locTMP}interval_signal_perbase.tsv

exit
#!##################################################################################################
#deeptools aproach does not really work for ggplot plotting
mawk -v OFS="\t" '{
  # PSL fields:

  qStart=$12
  qEnd=$13
  qSize=$11
  tStart=$16
  tEnd=$17
  matches=$1

  # Check for near end-to-end alignment on query
  if(qStart < 50 && qEnd > qSize - 50 && matches > qSize * 0.8 && tEnd-tStart < qSize*1.2){
    print $14,tStart,tEnd,$10":!:"$14"_"tStart"_"tEnd,0,$9
  }
}' ${locTMP}TEs_in_genome.psl | sort -k1,1 -k2,2n > ${locTMP}TEs_in_genome.filtered.bed


printf "
2L_RagTag 0 23171000
2R_RagTag 5972000 26061764
3L_RagTag 0 24176000
3R_RagTag 5501500 34930624
X_RagTag 0 22542000
" | tr ' ' '\t' > ${locTMP}euchromatin.bed

bedtools intersect -a ${locTMP}TEs_in_genome.filtered.bed -b ${locTMP}euchromatin.bed -u > ${locTMP}TEs_in_genome.filtered.euchromatin.bed


set +x 
ml build-env/f2022
ml deeptools/3.5.4-foss-2022a
set -x 

WINDOWsize=5000
BINsize=10


#path to the bb-file containing expressed transcript isoforms determined by stringtie using Illumina and direct RNA seq data
expressedTRANSCRIPTS=

bigBedToBed $expressedTRANSCRIPTS ${locTMP}HQmergedAnnotations.genepred
cut -f 1-12 ${locTMP}HQmergedAnnotations.genepred | sort -k1,1 -k2,2n > ${locTMP}HQmergedAnnotations.bed12

bedparse introns ${locTMP}HQmergedAnnotations.bed12 > ${locTMP}HQmergedAnnotations.annotations.bed
bedparse 3pUTR ${locTMP}HQmergedAnnotations.bed12 >> ${locTMP}HQmergedAnnotations.annotations.bed
bedparse 5pUTR ${locTMP}HQmergedAnnotations.bed12 >> ${locTMP}HQmergedAnnotations.annotations.bed
bedparse cds ${locTMP}HQmergedAnnotations.bed12 >> ${locTMP}HQmergedAnnotations.annotations.bed

sort -k1,1 -k2,2n ${locTMP}HQmergedAnnotations.annotations.bed > ${locTMP}HQmergedAnnotations.annotations.sorted.bed
bedtools intersect -wao -a ${locTMP}TEs_in_genome.filtered.euchromatin.bed -b ${locTMP}HQmergedAnnotations.annotations.sorted.bed |
  mawk -v OFS="\t" '
  {
    if($7 == ".") {
      $7="NONE"
    }else{
      if($6==$12){
        $7="SENSE"
      }else{
        $7="ANTISENSE"
      }
    }
    print
  }' | cut -f 1-7 | sort -k4,4 |
   bedtools groupby -i - -g 1,2,3,4,5,6  -c 7 -o distinct -i - | sort -k1,1 -k2,2n  |
  mawk -v OFS="\t" '{gsub(",","~",$7); $4=$4":!:"$7; print}' | cut -f 1-6  | sort -k1,1 -k2,2n > ${locTMP}TEs_in_genome.filtered.euchromatin.strandOverlap.bed 


  computeMatrix scale-regions \
    -S 2025-09-26-ChIP_Sienski-2012/siGFP_H3K9me3_Sienski-2012_uniq.bw 2025-09-26-ChIP_Sienski-2012/siPiwi_H3K9me3_Sienski-2012_uniq.bw \
    -R ${locTMP}TEs_in_genome.filtered.euchromatin.strandOverlap.bed \
    --beforeRegionStartLength 5000 \
    --afterRegionStartLength 5000 \
    --regionBodyLength 10 \
    --binSize 10 \
    -o ${locTMP}matrix_TE.gz \
    --outFileNameMatrix ${OPENdir}TEanalysis-genome/TEs_H3K9me3_matrix.txt \
    --skipZeros \
    --numberOfProcessors ${THREADS} \
    --missingDataAsZero 


  plotHeatmap -m ${locTMP}matrix_TE.gz \
    --heatmapHeight 20 \
    --colorMap Greens \
    --kmeans 3 \
    --heatmapWidth 5 \
    --outFileName ${OPENdir}TEanalysis-genome/TEs_H3K9me3.${ASSEMBLYname}.Sienski-2012.pdf \
    --outFileNameMatrix ${OPENdir}TEanalysis-genome/TEs_H3K9me3_matrix.plotHeatmap.txt.gz  \
    --samplesLabel siLuc_H3K9me3 siPiwi_H3K9me3

exit
#determine copia and intron overlap
bigBedToBed ${OSC_annotations}HQmergedAnnotations.bb ${locTMP}HQmergedAnnotations.genepred
cut -f 1-12 ${locTMP}HQmergedAnnotations.genepred | sort -k1,1 -k2,2n > ${locTMP}HQmergedAnnotations.bed12


bedparse introns ${locTMP}HQmergedAnnotations.bed12 > ${locTMP}HQmergedAnnotations.introns.bed
bedparse exons ${locTMP}HQmergedAnnotations.bed12 > ${locTMP}HQmergedAnnotations.exons.bed

grep copia ${locTMP}TEs_in_genome.filtered.euchromatin.bed > ${locTMP}copia_TEs.bed

echo TE STRAND REGION COUNT > ${OPENdir}TEanalysis-genome/copia_counts.txt
for i in sense antisense; do
  if [[ $i == sense ]]; then
    STRANDoption="-s"
  else
    STRANDoption="-S"
  fi

  for j in introns exons; do
    bedtools intersect $STRANDoption -a ${locTMP}HQmergedAnnotations.${j}.bed -b ${locTMP}copia_TEs.bed -wa > ${locTMP}copia_in_${j}.${i}.bed
    COUNT=$(cat ${locTMP}copia_in_${j}.${i}.bed | wc -l) 
    echo copia ${i} ${j} ${COUNT} >> ${OPENdir}TEanalysis-genome/copia_counts.txt
  done
done

exit

#!##################################################################################################
#!##################################################################################################
#!##################################################################################################
#! #old stuff
###################################################################################################
#determine SV length distribution

#! inverting deletion/insertion as I want to state their status in the OSC genome but for later stages I had to use OSC as the reference genome
mawk -v OFS="\t" '
{
  if($8 ~"SVTYPE=INS" ) {
    print $3,"DEL",length($5);
  }else{
    if($8 ~"SVTYPE=DEL" ){
      print $3,"INS",length($4);
    }
  }
}' ${locTMP}variants.vcf > ${locTMP}SVlengths.txt

set +x
ml build-env/f2022
ml  r/4.5.1-gfbf-2023b
set -x
#plot in R
R -e "
  library(ggplot2)
  library(dplyr)
  library(cowplot)
  library(colorblindr)

SVlengths <- read.table('${locTMP}SVlengths.txt', header = FALSE, col.names = c('ID', 'TYPE', 'LENGTH'))%>%
  filter(LENGTH < 20000) 

  # Create separate datasets for positive and negative
  pos_data <- SVlengths %>% filter(TYPE != 'DEL')
  neg_data <- SVlengths %>% filter(TYPE == 'DEL')


  ggplot() +
    geom_histogram(data = pos_data, aes(x = LENGTH, fill = TYPE), 
                   bins = 50,  position = 'identity') +
    geom_histogram(data = neg_data, aes(x = LENGTH, y = -after_stat(count), fill = TYPE), 
                   bins = 50, position = 'identity') +
    scale_x_log10()+
    geom_hline(yintercept = 0, color = 'black', linewidth = 0.5) +
    labs( x = 'Length', y = 'Count') +
    scale_fill_manual(values = c('INS' = '#081626', 'DEL' = '#28728e')) +
    scale_fill_OkabeIto() +
    theme_cowplot(16)+
    theme(
      legend.position = 'inside',
      legend.position.inside = c(0.2, 0.8),      
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank()
    )


  ggsave('${OPENdir}TEanalysis-genome/SVlengths.${ASSEMBLYname}.pdf', width = 5, height = 4)
"

###################################################################################################
#determine which SVs are present in W1118
#download ONT data

ONT_DNA_W1118=${locTMP}ONT_W1118.fastq.gz

#align the reads to the OSC genome
if [[ ! -s ${locTMP}W1118.mapped.bam ]]; then
  if [[ ! -s ${locTMP}ONT_W1118.fastq.gz ]]; then
    curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR136/009/ERR13697009/ERR13697009.fastq.gz > ${locTMP}ONT_W1118.fastq.gz
  fi

  minimap2 -Lax map-ont -t $THREADS --secondary=no $ASSEMBLYseq $ONT_DNA_W1118  | samtools sort -@ $THREADS -m 1G -O BAM - > ${locTMP}W1118.mapped.bam
  samtools index -@ $THREADS ${locTMP}W1118.mapped.bam
fi

#sniffles
sniffles --input ${locTMP}W1118.mapped.bam --genotype-vcf ${locTMP}/variants.vcf --vcf ${locTMP}SV.sniffles.W1118.vcf --reference $ASSEMBLYseq --threads $THREADS --allow-overwrite 
sniffles --input ${TMPdir}/SNV/SNV_ONT/OSC_r1.01/reads.mapped.bam --genotype-vcf ${locTMP}/variants.vcf --vcf ${locTMP}SV.sniffles.OSC.vcf --reference $ASSEMBLYseq --threads $THREADS --allow-overwrite 

bgzip -f ${locTMP}SV.sniffles.W1118.vcf
bgzip -f ${locTMP}SV.sniffles.OSC.vcf

bcftools index -f ${locTMP}SV.sniffles.W1118.vcf.gz
bcftools index -f ${locTMP}SV.sniffles.OSC.vcf.gz
#compare the two VCF files
bcftools merge -O z -o ${locTMP}SV.sniffles.merged.vcf.gz  ${locTMP}SV.sniffles.OSC.vcf.gz ${locTMP}SV.sniffles.W1118.vcf.gz --force-samples
bcftools index -f ${locTMP}SV.sniffles.merged.vcf.gz

bgzip -d -c  ${locTMP}SV.sniffles.merged.vcf.gz > ${locTMP}SV.sniffles.merged.vcf

mawk -v OFS="\t" '
{
    if($8 ~"SVTYPE=INS" && length($5)>100) {
      print ">"$3"\n"$5;
    }else{
      if($8 ~"SVTYPE=DEL" && length($4)>100){
        print ">"$3"\n"$4
      }
    }
  
}' ${locTMP}SV.sniffles.merged.vcf > ${locTMP}SV.sniffles.merged.fasta

minimap2 --paf-no-hit ${TEconsensus} ${locTMP}SV.sniffles.merged.fasta > ${locTMP}variants_aligned_to_TE.classified.paf


mawk -v OFS="\t" -v TEs=${locTMP}variants_aligned_to_TE.classified.paf '
BEGIN{
  while ((getline < TEs) > 0) {
    if($7==0){
      TE[$1] = "noTE"
    }else{
      if($4-$3 > $2*0){}
        TE[$1] = $6
    }
  }
}
{
  print $0, TE[$3]
} ' ${locTMP}SV.sniffles.W1118.vcf > ${locTMP}variants_aligned_to_TE.classified.annotated.paf


mawk -v OFS="\t" -v  INFILE=${locTMP}SV.sniffles.W1118.vcf -v TEs=${locTMP}variants_aligned_to_TE.classified.paf  '
BEGIN{
  while ((getline < INFILE) > 0) {
    split($NF,splitNAME,":")
    W1118[$3]= splitNAME[1]

    split($10,splitNAME,":")
    OSC[$3]= splitNAME[1]
  }
  while ((getline < TEs) > 0) {
    if($7==0){
      TE[$1] = "noTE"
    }else{
      if($4-$3 > $2*0){}
        TE[$1] = $6
    }
  }

}
{
  print $0, W1118[$3], OSC[$3], TE[$3]
}' ${locTMP}variants.vcf > ${locTMP}variants.inclW1118.vcf

awk -v OFS="\t" '
{
  if($1 ~"^#") {
    next
  }

  if($(NF-2) ~"1/1") {
    W1118 = "noW1118"
  }else{
    if($(NF-2) ~"0/1") {
      W1118 = "het"
    }else{
      W1118 = "hom"
    }
  }
  if($(NF-1) ~"0/0" || $(NF-1) ~"0/1") {
    if($8 ~"SVTYPE=INS" && length($5)>100) {
      X["INS"][$NF][W1118] +=-1
    }
    if($8 ~"SVTYPE=DEL" && length($4)>100){
      X["DEL"][$NF][W1118] +=1
    }
  }
}
END{
  print "CLASS", "TE", "W1118", "COUNT"
  for (CLASS in X) {
    for (TE in X[CLASS]) {
      for (W1118 in X[CLASS][TE]) {
        print CLASS, TE, W1118, X[CLASS][TE][W1118]
      }
    }
  }
}' ${locTMP}variants.inclW1118.vcf > ${OPENdir}TEanalysis-genome/TEquantification.${ASSEMBLYname}inclW1118.txt


#plot in R
ml build-env/f2022
ml  r/4.5.1-gfbf-2023b
R -e "
  library(tidyverse)
  library(dplyr)

  RAW <- read_tsv('TEquantification.OSC_r1.01inclW1118.txt', col_names = TRUE) %>%
    pivot_wider(names_from = c(CLASS), values_from = 'COUNT', values_fill = 0)

  # Calculate total INS, DEL, CLASS, and sortSUM per TE (across all W1118)
  TE_summary <- RAW %>%
    group_by(TE) %>%
    summarise(
      total_INS = sum(INS, na.rm = TRUE),
      total_DEL = sum(DEL, na.rm = TRUE)
    ) %>%
    mutate(
      CLASS = case_when(
        abs(total_INS / total_DEL) < 0.5 ~ 'INSERTED',
        abs(total_INS / total_DEL) > 2 ~ 'DELETED',
        TRUE ~ 'BALANCED'
      ),
      TOTAL_SUM = abs(total_INS) + abs(total_DEL)
    )

  # Join back to main data
  RAW <- RAW %>%
    left_join(TE_summary %>% select(TE, CLASS, TOTAL_SUM, total_DEL), by = 'TE') %>%
    filter(!str_detect(TE, 'random'), !str_detect(TE, 'noTE'), TOTAL_SUM > 15)

  # Pivot longer and reorder TE for plotting
  RAW_long <- RAW %>%
    pivot_longer(cols = c(INS, DEL), names_to = 'SVclass', values_to = 'COUNT') %>%
    mutate(
      COUNT = as.numeric(COUNT),
      TE = fct_reorder(TE, total_DEL, .desc = FALSE)
    )

  # Now you can plot, facet by CLASS if you want
  p <- RAW_long %>%
    ggplot(aes(x = TE, y = COUNT, fill = W1118)) +
    geom_bar(stat = 'identity', width = 0.5, position = 'stack') +
    coord_flip() +
    theme_minimal() +
    # ggforce::facet_col(~ CLASS, scales = 'free_y', space = 'free')+
  labs(y = 'Count', x = 'TE')

  
  ggsave(p, filename = '${OPENdir}TEanalysis-genome/TEquantification.${ASSEMBLYname}.incl-W1118.pdf', width = 6, height = 4)
"

###################################################################################################
#test non W1118 SVs for loss of heterozygosity

#create extended bed-ranges from vcf
mawk -v OFS="\t"  -v TEs=${locTMP}variants_aligned_to_TE.classified.paf -v CHRsizes=${CHRsizes} '
BEGIN{
  while ((getline < TEs) > 0) {
    if($7==0){
      TE[$1] = "noTE"
    }else{
      if($4-$3 > $2*0){}
        TE[$1] = "TE"
    }
  }
  while ((getline < CHRsizes) > 0) {
    if($1 !~ "^#") {
      CHR[$1]=$2
    }
  }
  EXTNSION=20000
}
{
  if($1 !~ "^#") {
    if($8 ~ "SVTYPE=DEL" && TE[$3] == "TE") {
      split($8, splitTAG, /=|;/)
      split($NF, splitTAG2, ":")
      if($2-EXTNSION >= 0 ) {
        print $1, $2-EXTNSION,$2, $3":"splitTAG2[1]":us", 0, "+"
      }
      if($2+splitTAG[6]+EXTNSION <= CHR[$1]) {
        print $1, $2+splitTAG[6], $2+splitTAG[6]+EXTNSION, $3":"splitTAG2[1]":ds", 0, "+"
      }
    }
  }
}' ${locTMP}SV.sniffles.OSC.vcf | LC_COLLATE=C sort -k1,1 -k2,2n > ${locTMP}TE_insertions.extended.bed


#determine analyzable regions in the window
#filter bam
samtools view -F 0x904 -q 1 -b ${TMPdir}SNV/SNV_Illumina/OSC_r1.01/reads.mapped.bam > ${TMPdir}SNV/SNV_Illumina/OSC_r1.01/reads.mapped.filtered.bam
samtools index ${TMPdir}SNV/SNV_Illumina/OSC_r1.01/reads.mapped.filtered.bam
bedtools genomecov -ibam ${TMPdir}SNV/SNV_Illumina/OSC_r1.01/reads.mapped.filtered.bam -bg > ${locTMP}coverage.bed

LC_COLLATE=C sort -k1,1 -k2,2n $CHRsizes > ${locTMP}chr.sizes.bed 
LC_COLLATE=C sort -k1,1 -k2,2n --parallel=$THREADS -S${MEM}g ${locTMP}coverage.bed > ${locTMP}coverage.sorted.bed

bedtools intersect -sorted -wao -g ${locTMP}chr.sizes.bed -a ${locTMP}TE_insertions.extended.bed -b ${locTMP}coverage.sorted.bed | 
  mawk -v OFS="\t" '
  {
    if($4 == currID) {
      if($10>20){
        SUM+=$NF
      }
    }else{
      if(NR>1){
        print CHR,START,STOP,currID":!:"SUM, 0, "+"
      }
      currID=$4
      CHR=$1
      START=$2
      STOP=$3
  
      SUM=0
      if($10>20){
        SUM=$NF
      }
    }

  }' | sort -k1,1 -k2,2n > ${locTMP}TE_insertions.extended.uniq.bed


bgzip -dc ${TMPdir}SNV/SNV_Illumina/OSC_r1.01/output.vcf.gz > ${locTMP}OSC_SNPs.vcf 
bedtools intersect -wao -a ${locTMP}TE_insertions.extended.uniq.bed -b ${locTMP}OSC_SNPs.vcf | 
  awk -v OFS="\t" -v W1118vcf=${locTMP}SV.sniffles.W1118.vcf '
  BEGIN{
    while ((getline < W1118vcf) > 0) {
      if($1 !~ "^#") {
        split($NF, splitTAG, ":")
        W1118[$3] = splitTAG[1]
      }
    }
  }
  {
    if($7 != "." && $0!~ "RefCall" ){
      split($(NF-1), splitTAG, ":")
      X[$4][splitTAG[1]] += 1
    }
  }
  END{
    print "RegionID","RegionZygosity","HET", "HOM", "NO_CALL","UNIQUE_POS","W1118"
    for (i in X) {
      if(!("0/1" in X[i])) { X[i]["0/1"]=0 }
      if(!("1/1" in X[i])) { X[i]["1/1"]=0 }
      if(!("." in X[i])) { X[i]["."]=0 }
      split(i, splitID,/:!:/)
      split(splitID[1],splitZYGO , ":")
      print i, splitZYGO[2],X[i]["0/1"]*1000/splitID[2], X[i]["1/1"]*1000/splitID[2], X[i]["."]*1000/splitID[2],splitID[2],W1118[splitZYGO[1]]
    }
  }' > ${OPENdir}TEanalysis-genome/SNPs_in_TE_insertions.txt
exit
#!##################################################################################################
#!####################    not functional due to missint W1118 data.    #############################   
#!##################################################################################################
exit
#compare OSC to W1118 SVs determined on the OSC genome

#align the reads to the OSC genome
#map reads
if [[ ! -s ${locTMP}reads.mapped.bam ]]; then
  if [[ $COMPUTING == C ]]; then
    
    sbatch --parsable --array=0-1 --wait --job-name=minimap_variations -o "%x.o.%A-%a.txt" -e "%x.e.%A-%a.txt" --cpus-per-task=20 --mem=50g --wrap="
        THREADS=\$(( \$SLURM_CPUS_PER_TASK * 2 ))
        APPTAINERdir=$APPTAINERdir
        TMPdir=$TMPdir
        locTMP=${locTMP}
        source ${SCRIPTdir}tools

        if [[ \$SLURM_ARRAY_TASK_ID -eq 0 ]]; then
          inputDNA=$ONT_DNA
          NAMEext="OSC"
        else
          inputDNA=${ONT_DNA_W1118}
          NAMEext="W1118"
        fi

        minimap2 -Lax map-ont -t \$THREADS --secondary=no $ASSEMBLYseq \$inputDNA | samtools sort -@ \$THREADS -m 1G -O BAM - > ${locTMP}${NAMEext}.mapped.bam
        samtools index -@ \$THREADS ${locTMP}${NAMEext}.mapped.bam

      " 
  else
    minimap2 -Lax map-ont -t $THREADS --secondary=no $ASSEMBLYseq $ONT_DNA  | samtools sort -@ $THREADS -m 1G -O BAM - > ${locTMP}OSC.mapped.bam 
    samtools index -@ $THREADS ${locTMP}OSC.mapped.bam
    minimap2 -Lax map-ont -t $THREADS --secondary=no $ASSEMBLYseq $ONT_DNA_W1118  | samtools sort -@ $THREADS -m 1G -O BAM - > ${locTMP}W1118.mapped.bam
    samtools index -@ $THREADS ${locTMP}W1118.mapped.bam
  fi

fi
exit
#sniffles
  if [[ $COMPUTING == C ]]; then

    sbatch --wait  --job-name=SV_sniffles -o "%x.o.%A-%a.txt" -e "%x.e.%A-%a.txt" --cpus-per-task=20 --mem=50g --qos=short --time=2:00:00 --wrap="
      set -ux
      THREADS=\$(( \$SLURM_CPUS_PER_TASK * 2 ))
      APPTAINERdir=$APPTAINERdir
      TMPdir=$TMPdir
      source ${SCRIPTdir}tools

      sniffles --input ${locTMP}reads.mapped.bam --vcf ${locTMP}SV.sniffles.vcf --reference $currASSEMBLYseq --threads \$THREADS --allow-overwrite
      rtg vcfstats ${locTMP}SV.sniffles.vcf > ${locTMP}${currASSEMBLYname}.vcfstats.txt

      " 
  else
    sniffles --input ${locTMP}reads.mapped.bam --vcf ${locTMP}SV.sniffles.vcf --reference $currASSEMBLYseq --threads $THREADS --allow-overwrite
    rtg vcfstats ${locTMP}SV.sniffles.vcf > ${locTMP}${currASSEMBLYname}.vcfstats.txt
  
  fi
