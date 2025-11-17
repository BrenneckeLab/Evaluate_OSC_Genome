#!/bin/bash

#SBATCH --cpus-per-task=30
#SBATCH --mem=80g
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
locTMP=${TMPdir}cluster-analysis/

#create directories
mkdir $locTMP
mkdir -p ${OPENdir}cluster-analysis/

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
ASSEMBLYflam=$(echo $currASSEMBLYline |tr ' ' '\t' | cut -f 3)  #remove file extension


###################################################################################################
#determine cluster synteny



###################################################################################################
#determine cluster length
TILEsize=100

#path to the 25mer 1MM uniqueness track for OSC genome hub
UNIQUENESSfile=

#path to chromosome sizes file
CHRsizes=

#directory containing the sRNA-seq data from the Annotation Pipeline run
## structure within is /individual-libraries/LIBRARYNAME/LIBRARYNAME_annotated.bed.gz and /individual-libraries/LIBRARYNAME/normalization.txt
## The reads in the bed files are collapsed and a name structure of NR_484128_1_UMIcount_:GACTCGCGGCAGTTGGAGACGGATTT:count=1 is expected with the count-field indicating the number of identical reads
APdir=

###################################################################################################

#convert cluster coordinates to bed
echo $ASSEMBLYflam | mawk -v OFS="\t" '{
  split($1,splitCOORD,":|-")
  print splitCOORD[1],splitCOORD[2],splitCOORD[3]+800000,"flam",0,"+"
}' > ${locTMP}flam.bed

echo X_RagTag 22553405 22591523 20A 0 + | tr ' ' '\t' >> ${locTMP}flam.bed
echo 3L_RagTag	21354864	21379509	77B	0	+ | tr ' ' '\t' >> ${locTMP}flam.bed
echo 2L_RagTag 20471755 20473266 tj 0 + | tr ' ' '\t' >> ${locTMP}flam.bed
echo X_RagTag 3296526 3300341 myc 0 + | tr ' ' '\t' >> ${locTMP}flam.bed

#create cluster 1kb tiles
mawk -v OFS="\t" -v TILEsize=$TILEsize '
  BEGIN {
    print "chr","start","end","name","score","strand"
  }
{
  N=1
  for(i=$2; i<=$3-TILEsize; i+=TILEsize) {
    print $1,i,i+TILEsize,$4"_"N,0,$6
    N++
  }
}' ${locTMP}flam.bed | sort -k1,1 -k2,2n > ${locTMP}flam_tiles.bed

#---------------------------------------------------------------------------------------------------------
#determine uniqueness of the tiles

bigBedToBed $UNIQUENESSfile  ${locTMP}uniqueness.sRNA.25mer.1MM.bed

bedtools intersect -wao -sorted -nobuf -g ${CHRsizes} -s -a ${locTMP}flam_tiles.bed -b ${locTMP}uniqueness.sRNA.25mer.1MM.bed | 
  mawk -v OFS="\t" -v TILEsize=$TILEsize '
  {
    if($4 == currID){
      SUM+=$NF
    }else{
      if(currID != "") {
        print currCHROM, currSTART, currEND, currID":!:"SUM":!:"SUM*100/TILEsize,0,currSTRAND
      }
      currID=$4
      SUM=$NF
      currCHROM=$1
      currSTART=$2
      currEND=$3  
      currSTRAND=$6
    }
  }
  END{
    if(currID != "") {
      print currCHROM, currSTART, currEND, currID":!:"SUM":!:"SUM*100/TILEsize,0,currSTRAND
    }
  }' | sort -k1,1 -k2,2n > ${locTMP}flam_tiles.uniqueness.bed


#---------------------------------------------------------------------------------------------------------
#quantify sRNAs
rm -rf ${locTMP}flam_tiles.sRNA.counts.txt
for currLIB in $(ls ${APdir}individual-libraries/); do
  echo "Processing library: ${currLIB}"
  NORM=$(tail -n 1 ${APdir}individual-libraries/${currLIB}/normalization.txt | tr ' ' '\t' | cut -f 1 )
  gunzip -c ${APdir}individual-libraries/${currLIB}/${currLIB}_annotated.bed.gz | 
    grep "mapping=u" |
    bedtools intersect -wao -sorted -nobuf -g ${CHRsizes} -s -a ${locTMP}flam_tiles.uniqueness.bed -b stdin | 
    mawk -v NORM=$NORM -v currLIB=$currLIB -v OFS="\t" -v TILEsize=$TILEsize '
    {
      split($10, splitREAD, /:|=/)

      if($4 == currID){

        SUM+=splitREAD[4]
      }else{
        if(currID != "") {
          split(currID, splitID,/:!:/)
          print currID, currLIB, SUM/NORM, (SUM/NORM)/splitID[2]*TILEsize
        }
        currID=$4
        SUM=splitREAD[4]
      }
    }
    END{
      if(currID != "") {
          split(currID, splitID,/:!:/)
          print currID, currLIB, SUM/NORM, (SUM/NORM)/splitID[2]*TILEsize
      }
    }' >> ${locTMP}flam_tiles.sRNA.counts.txt
done 

mawk -v OFS="\t" '
BEGIN{
  print "ID","uniqPOS","uniqPERC","LIBRARY","sRNA.counts","sRNA.counts.perPOS"
}
{
  split($1,splitNAME,/:!:/)
  print splitNAME[1],splitNAME[2],splitNAME[3],$2,$3,$4
}' ${locTMP}flam_tiles.sRNA.counts.txt > ${OPENdir}cluster-analysis/flam_tiles.sRNA.counts.txt

#---------------------------------------------------------------------------------------------------------

###################################################################################################
###################################################################################################
#determine cluster length full chromosome
TILEsize=1000


#convert cluster coordinates to bed
echo $ASSEMBLYflam | mawk -v OFS="\t" '{
  split($1,splitCOORD,":|-")
  print splitCOORD[1],splitCOORD[2],splitCOORD[3]+800000,"flam",0,"+"
}' > ${locTMP}cluster.bed

echo X_RagTag 22553405 22591523 20A 0 + | tr ' ' '\t' >> ${locTMP}cluster.bed
echo 2L_RagTag 20471755 20473266 tj 0 + | tr ' ' '\t' >> ${locTMP}flam.bed
echo X_RagTag 3296526 3300341 myc 0 + | tr ' ' '\t' >> ${locTMP}flam.bed

#create X-chromosome 1kb tiles
grep X_RagTag $CHRsizes | 
  mawk -v OFS="\t" -v TILEsize=$TILEsize '
  {
    N=1
    for(i=0; i<=$2-TILEsize; i+=TILEsize) {
      print $1,i,i+TILEsize,"X_"i"_+",0,"+"
      print $1,i,i+TILEsize,"X_"i"_-",0,"-"
      N++
    }
  }'  | LC_COLLATE=C sort -k1,1 -k2,2n --parallel=$THREADS -S${MEM}g > ${locTMP}chrX.tiles.bed

#---------------------------------------------------------------------------------------------------------
#determine uniqueness of the tiles

bigBedToBed $UNIQUENESSfile  ${locTMP}uniqueness.sRNA.25mer.1MM.orig.bed
LC_COLLATE=C sort -k1,1 -k2,2n --parallel=$THREADS -S${MEM}g ${locTMP}uniqueness.sRNA.25mer.1MM.orig.bed > ${locTMP}uniqueness.sRNA.25mer.1MM.bed
LC_COLLATE=C sort -k1,1 -k2,2n --parallel=$THREADS -S${MEM}g $CHRsizes > ${locTMP}chr.sizes
bedtools intersect -wao -sorted -nobuf -g ${locTMP}chr.sizes -s -a ${locTMP}chrX.tiles.bed -b ${locTMP}uniqueness.sRNA.25mer.1MM.bed | 
  mawk -v OFS="\t" -v TILEsize=$TILEsize '
  {
    if($4 == currID){
      SUM+=$NF
    }else{
      if(currID != "") {
        print currCHROM, currSTART, currEND, currID":!:"SUM":!:"SUM*100/TILEsize,0,currSTRAND
      }
      currID=$4
      SUM=$NF
      currCHROM=$1
      currSTART=$2
      currEND=$3  
      currSTRAND=$6
    }
  }
  END{
    if(currID != "") {
      if(SUM*100/TILEsize > 50){
        print currCHROM, currSTART, currEND, currID":!:"SUM":!:"SUM*100/TILEsize,0,currSTRAND
      } 
    }
  }' | sort -k1,1 -k2,2n > ${locTMP}chrX.uniqueness.bed


#---------------------------------------------------------------------------------------------------------
#quantify sRNAs
rm -rf ${locTMP}chrX.sRNA.counts${EXT}.txt
for currLIB in $(ls ${APdir}individual-libraries/); do
  echo "Processing library: ${currLIB}"
  NORM=$(tail -n 1 ${APdir}individual-libraries/${currLIB}/normalization.txt | tr ' ' '\t' | cut -f 1 )
  gunzip -c ${APdir}individual-libraries/${currLIB}/${currLIB}_annotated.bed.gz | 
    grep "mapping=u" |
    bedtools intersect -wao -sorted -nobuf -g ${CHRsizes} -s -a ${locTMP}chrX.uniqueness.bed -b stdin | 
    mawk -v NORM=$NORM -v currLIB=$currLIB -v OFS="\t" -v TILEsize=$TILEsize '
    {
      split($10, splitREAD, /:|=/)

      if($4 == currID){

        SUM+=splitREAD[4]
      }else{
        if(currID != "") {
          split(currID, splitID,/:!:/)
          print currID, currLIB, SUM/NORM, (SUM/NORM)/splitID[2]*TILEsize
        }
        currID=$4
        SUM=splitREAD[4]
      }
    }
    END{
      if(currID != "") {
          split(currID, splitID,/:!:/)
          print currID, currLIB, SUM/NORM, (SUM/NORM)/splitID[2]*TILEsize
      }
    }' >> ${locTMP}chrX.sRNA.counts${EXT}.txt
done

mawk -v OFS="\t" '
BEGIN{
  print "ID","uniqPOS","uniqPERC","LIBRARY","sRNA.counts","sRNA.counts.perPOS"
}
{
  split($1,splitNAME,/:!:/)
  print splitNAME[1],splitNAME[2],splitNAME[3],$2,$3,$4
}' ${locTMP}chrX.sRNA.counts${EXT}.txt > ${OPENdir}cluster-analysis/chrX.sRNA.counts${EXT}.txt


###################################################################################################
###################################################################################################
#flamenco splicing

cp ${ASSEMBLYseq} ${locTMP}OSC.genome.fa
#cutoff for read counts to consider a splice junction
CUTOFF=3
#minimal US-region cutoff
CUTOFFus=10
#upstream region used to quantify the donor site
US=1



#index assembly using STAR
if [[ ! -s ${locTMP}STARindex/SAindex ]]; then
  STAR --runThreadN 5 --runMode genomeGenerate --genomeDir ${locTMP}STARindex \
    --genomeFastaFiles ${locTMP}OSC.genome.fa \
    --genomeSAindexNbases 12 \
    --limitGenomeGenerateRAM 10000000000 
fi
  

#split interleaved paired-end reads
if [[ ! -s ${locTMP}R1.fastq && ! -s ${locTMP}R2.fastq ]]; then
  seqkit fx2tab ${OSC_RNAseq_PE} | 
  mawk -v OFS="\t" -v locTMP=${locTMP} '
  BEGIN{
    COMMAND1="seqkit tab2fx -w 0 > "locTMP"R1.fastq"
    COMMAND2="seqkit tab2fx -w 0 > "locTMP"R2.fastq"
  }
  {
    sub("/1$", "", $1)
    sub("/2$", "", $1)
    if(NR%2==1) {X=$0; Y=length($2);}
    else if(NR%2==0) {
      if(length($2)>90 && Y>90) {
        print X | COMMAND1
        print $0 | COMMAND2
      }
    }
  }' 

fi

cp ${locTMP}R2.fastq ${locTMP}unpaired.fastq
seqkit seq --reverse --complement ${locTMP}R1.fastq >> ${locTMP}unpaired.fastq

#map reads to assembly using STAR 
STAR --runThreadN $THREADS \
    --genomeDir ${locTMP}STARindex \
    --readFilesIn ${locTMP}unpaired.fastq \
    --outFileNamePrefix ${locTMP}STAR_unpaired/  \
    --outFilterMultimapNmax 1 --winAnchorMultimapNmax 100 \
    --alignEndsType EndToEnd --twopassMode Basic \
    --outSAMtype BAM Unsorted --outFilterIntronMotifs RemoveNoncanonical


#sort all mappers
currTHREADS=$(( $THREADS / 2 ))
MEMperThread=$(( ( $MEM - 10 ) / $currTHREADS))
samtools sort --output-fmt BAM -m ${MEMperThread}G -@ $currTHREADS -o ${locTMP}all_mappers_unpaired.bam ${locTMP}STAR_unpaired/Aligned.out.bam 


samtools view -h ${locTMP}all_mappers_unpaired.bam | 
    awk '$0 ~ /^@/ || $0 ~ /NH:i:1(\s|$)/' |
    samtools view -b -o ${locTMP}unique_mappers_unpaired.bam

samtools index -@ $THREADS ${locTMP}unique_mappers_unpaired.bam


samtools view -h ${locTMP}unique_mappers_unpaired.bam | awk '
BEGIN {OFS="\t"}
/^@/ {print; next}  # Print headers
$6 ~ /[0-9]+N/ {print}  # Print reads with N in CIGAR
' | samtools view -b -o ${locTMP}spliced_reads_unpaired.bam

samtools index -@ $THREADS ${locTMP}spliced_reads_unpaired.bam
  
cp ${locTMP}spliced_reads_unpaired.bam* ${OPENdir}cluster-analysis/

#create bw file

bedtools genomecov -ibam ${locTMP}unique_mappers_unpaired.bam -split -bg -g ${CHRsizes} -strand + | 
  mawk -v OFS="\t" -v NORM=$OSC_RNAseq_PE_NORM '{$4=$4/NORM; print}'  | sort -k1,1 -k2,2n > ${locTMP}unique_mappers.sense_unpaired.bedGraph

#create bw file
bedtools genomecov -ibam ${locTMP}unique_mappers_unpaired.bam -split -bg -g ${CHRsizes} -strand - |
  mawk -v OFS="\t" -v NORM=$OSC_RNAseq_PE_NORM '{$4=-$4/NORM; print}'  | sort -k1,1 -k2,2n  > ${locTMP}unique_mappers.antisense_unpaired.bedGraph

bedGraphToBigWig ${locTMP}unique_mappers.sense_unpaired.bedGraph ${CHRsizes} ${OPENdir}cluster-analysis/unique_mappers.sense_unpaired.bw
bedGraphToBigWig ${locTMP}unique_mappers.antisense_unpaired.bedGraph ${CHRsizes} ${OPENdir}cluster-analysis/unique_mappers.antisense_unpaired.bw


mawk -v OFS="\t" -v CUTOFF=$CUTOFF '
{
  if($4>CUTOFF){
    print $1,$2,$3,NR,0,"+"
  }
}'  ${locTMP}unique_mappers.sense_unpaired.bedGraph|
  sort -k1,1 -k2,2n  > ${locTMP}analyzable-regions_for_spliceAnalysis.sense.bed
bedToBigBed ${locTMP}analyzable-regions_for_spliceAnalysis.sense.bed ${CHRsizes} ${OPENdir}cluster-analysis/analyzable-regions_for_spliceAnalysis.sense.bb


mawk -v OFS="\t" -v CUTOFF=$CUTOFFus '
{
  if(-$4>CUTOFF){
    print $1,$2,$3,NR,0,"-"
  }
}'  ${locTMP}unique_mappers.antisense_unpaired.bedGraph| 
  sort -k1,1 -k2,2n > ${locTMP}unique_mappers.antisense_unpaired.filtered.antisense.bedGraph
  bedToBigBed ${locTMP}unique_mappers.antisense_unpaired.filtered.antisense.bedGraph ${CHRsizes}  ${OPENdir}cluster-analysis/analyzable-regions_for_spliceAnalysis.antisense.bb



exit


#determine splice-junctions using regtools

bedtools bamtobed -bed12 -i ${locTMP}spliced_reads_unpaired.bam |
  mawk -v OFS="\t" '{$5=1; print }' > ${locTMP}junctions.raw_unpaired.bed

mawk -v OFS="\t" -v locTMP=${locTMP} '
{
  n=split($11,splitEXON,/,/)
  split($12,splitINTRON,/,/)
  if(n==2){
    print $1, $2+splitEXON[1]-1, $2+splitINTRON[2]+1, $4, $5, $6 
  }else{
    print $0 > locTMP"junctions.complex.bed"
  }
}' ${locTMP}junctions.raw_unpaired.bed | 
  bedtools sort -i - | \
    bedtools groupby -i - \
        -g 1,2,3,6 \
        -c 4,5 -o first,sum > ${locTMP}junctions.collapsed_unpaired.bed

printf 'table interact
"interaction between two regions"
    ( 
    string chrom;        "Chromosome (or contig, scaffold, etc.). For interchromosomal, use 2 records" 
    uint chromStart;     "Start position of lower region. For interchromosomal, set to chromStart of this region" 
    uint chromEnd;       "End position of upper region. For interchromosomal, set to chromEnd of this region"
    string name;         "Name of item, for display.  Usually 'sourceName/targetName/exp' or empty"
    uint score;          "Score (0-1000)"
    double value;        "Strength of interaction or other data value. Typically basis for score"
    string exp;          "Experiment name (metadata for filtering). Use . if not applicable"
    string color;        "Item color.  Specified as r,g,b or hexadecimal #RRGGBB or html color name, as in //www.w3.org/TR/css3-color/#html4. Use 0 and spectrum setting to shade by score"
    string sourceChrom;  "Chromosome of source region (directional) or lower region. For non-directional interchromosomal, chrom of this region."
    uint sourceStart;    "Start position in chromosome of source/lower/this region"
    uint sourceEnd;      "End position in chromosome of source/lower/this region"
    string sourceName;   "Identifier of source/lower/this region"
    string sourceStrand; "Orientation of source/lower/this region: + or -.  Use . if not applicable"
    string targetChrom;  "Chromosome of target region (directional) or upper region. For non-directional interchromosomal, chrom of other region"
    uint targetStart;    "Start position in chromosome of target/upper/this region"
    uint targetEnd;      "End position in chromosome of target/upper/this region"
    string targetName;   "Identifier of target/upper/this region"
    string targetStrand; "Orientation of target/upper/this region: + or -.  Use . if not applicable"

    )
' > ${locTMP}interact.as

for STRAND in sense antisense; do
  if [[ $STRAND == sense ]]; then
    STRANDsign="+"
  else
    STRANDsign="-"
  fi
  awk -v STRANDsign=$STRANDsign -v CUTOFF=$CUTOFF '
    BEGIN {OFS="\t"} 
    {
        chrom      = $1
        chromStart = $2
        chromEnd   = $3
        junctionID = $5
        readCount  = $6
        strandStr  = $4

        if(readCount <= CUTOFF) next
        if(strandStr != STRANDsign )next

        # splice sites (donor = intron start, acceptor = intron end)
        sourceStart = chromStart
        sourceEnd   = chromStart + 1
        targetStart = chromEnd - 1
        targetEnd   = chromEnd

        # scale score
        score = (readCount > 1000 ? 1000 : readCount)

        # name for interact item
        name = chrom ":" sourceStart "-" sourceEnd "_" chrom ":" targetStart "-" targetEnd

        # color by strand
        if (strandStr == "+")  color = "0,0,255"
        else if (strandStr=="-") color = "255,0,0"
        else color="128,128,128"
        color=0

        print chrom, sourceStart, targetEnd, name, score, readCount, "regtools_SJ", color, \
              chrom, sourceStart, sourceEnd, "exon_donor", strandStr, \
              chrom, targetStart, targetEnd, "exon_acceptor", strandStr
    }' ${locTMP}junctions.collapsed_unpaired.bed \
    | sort -k1,1 -k2,2n > ${locTMP}splice_junctions_unpaired.${STRAND}.interact


  #add proper score 
    #convert to bed covering the 3 nt upstream of the splice junction
    mawk -v OFS="\t" -v US=$US -v STRAND=$STRAND '
    {
      if(STRAND == "sense") {
        print $1, $2-US, $2+1, $4, 0, $13
      }else{
        print $1, $3-1, $3+US, $4, 0, $13
      }
    }' ${locTMP}splice_junctions_unpaired.${STRAND}.interact | 
  bedtools map -a stdin -b ${locTMP}unique_mappers.${STRAND}_unpaired.bedGraph -c 4 -o mean -g ${CHRsizes} > ${locTMP}usCount_unpaired.${STRAND}.txt

  awk -v OFS="\t" -v INFILE=${locTMP}usCount_unpaired.${STRAND}.txt -v US=$US -v STRAND=$STRAND -v CUTOFFus=$CUTOFFus '
    BEGIN{
      while((getline < INFILE) > 0) {
        key = $4
        if(STRAND == "sense") {
          usCount[key] = $7
        } else {
          usCount[key] = $7*-1
        }
      }
    }
    {
      key = $4
      if(key in usCount) {
        us = usCount[key]
      } else {
        us = 0
      }
      if(us > 0) {
        score = int($6 / us * 1000)
        if(score > 1000) score = 1000
      } else {
        score = 500
      }
      $5=score
      $4=$4"_usCOUNT="us
      if(us > CUTOFFus) {
        print
      } 
    }' ${locTMP}splice_junctions_unpaired.${STRAND}.interact > ${locTMP}splice_junctions_unpaired.scored.${STRAND}.interact


    bedToBigBed -as=${locTMP}interact.as -type=bed5+13 -tab ${locTMP}splice_junctions_unpaired.scored.${STRAND}.interact ${CHRsizes} ${OPENdir}cluster-analysis/splice_junctions_unpaired.${STRAND}.interact
done
exit

####################################################################################################
#paired analysis

#map reads to assembly using STAR 
STAR --runThreadN $THREADS \
    --genomeDir ${locTMP}STARindex \
    --readFilesIn ${locTMP}R1.fastq ${locTMP}R2.fastq \
    --outFileNamePrefix ${locTMP}STAR/  \
    --outFilterMultimapNmax 1 --winAnchorMultimapNmax 2000 \
    --alignEndsType Local --twopassMode Basic \
    --outSAMtype BAM Unsorted --outStd SAM 


#sort all mappers
currTHREADS=$(( $THREADS / 2 ))
MEMperThread=$(( ( $MEM - 20 ) / $currTHREADS))
samtools sort --output-fmt BAM -m ${MEMperThread}G -@ $currTHREADS -o ${locTMP}all_mappers.bam ${locTMP}STAR/Aligned.out.bam &


samtools view -h ${locTMP}all_mappers.bam | 
    awk '$0 ~ /^@/ || $0 ~ /NH:i:1(\s|$)/' |
    samtools view -b -o ${locTMP}unique_mappers.bam

samtools index -@ $THREADS ${locTMP}unique_mappers.bam


samtools view -h ${locTMP}unique_mappers.bam | awk '
BEGIN {OFS="\t"}
/^@/ {print; next}  # Print headers
$6 ~ /[0-9]+N/ {print}  # Print reads with N in CIGAR
' | samtools view -b -o ${locTMP}spliced_reads.bam

samtools index -@ $THREADS ${locTMP}spliced_reads.bam

cp ${locTMP}spliced_reads.bam* ${OPENdir}cluster-analysis/

samtools view -h -f 64 -q 255 ${locTMP}unique_mappers.bam | 
  bedtools bamtobed -i - -split | 
  mawk -v OFS="\t" '{if($6 == "+") { $6="-";}else $6="+"; print}' > ${locTMP}alignments_reoriented.R1.bed

# second mates only, uniquely aligned
samtools view -h -f 128 -q 255 ${locTMP}unique_mappers.bam  |
  bedtools bamtobed -i - -split > ${locTMP}alignments_reoriented.R2.bed

sort -k1,1 -k2,2n --parallel=$THREADS -S${MEM}g ${locTMP}alignments_reoriented.R1.bed ${locTMP}alignments_reoriented.R2.bed > ${locTMP}alignments_reoriented.sorted.bed
  

#create bw file
bedtools genomecov -i ${locTMP}alignments_reoriented.sorted.bed  -bg -g ${CHRsizes} -strand + > ${locTMP}unique_mappers.sense.bedGraph

#create bw file
bedtools genomecov -i ${locTMP}alignments_reoriented.sorted.bed  -bg -g ${CHRsizes} -strand - |
  mawk -v OFS="\t" '{$4=-$4; print}'  > ${locTMP}unique_mappers.antisense.bedGraph

bedGraphToBigWig ${locTMP}unique_mappers.sense.bedGraph ${CHRsizes} ${OPENdir}cluster-analysis/unique_mappers.sense.bw
bedGraphToBigWig ${locTMP}unique_mappers.antisense.bedGraph ${CHRsizes} ${OPENdir}cluster-analysis/unique_mappers.antisense.bw



#determine splice-junctions using regtools
regtools junctions extract -a 8 -m 50 -M 500000  -s FS ${locTMP}unique_mappers.bam > ${locTMP}junctions.raw.bed

mawk -v OFS="\t" -v locTMP=${locTMP} '
{
  n=split($11,splitEXON,/,/)
  split($12,splitINTRON,/,/)
  if(n==2){
    print $1, $2+splitEXON[1]-1, $2+splitINTRON[2]+1, $4, $5, $6 
  }else{
    print $0 > locTMP"junctions.complex.bed"
  }
}' ${locTMP}junctions.raw.bed | 
  bedtools sort -i - | \
    bedtools groupby -i - \
        -g 1,2,3,4,6 \
        -c 5 -o sum > ${locTMP}junctions.collapsed.bed

printf 'table interact
"interaction between two regions"
    ( 
    string chrom;        "Chromosome (or contig, scaffold, etc.). For interchromosomal, use 2 records" 
    uint chromStart;     "Start position of lower region. For interchromosomal, set to chromStart of this region" 
    uint chromEnd;       "End position of upper region. For interchromosomal, set to chromEnd of this region"
    string name;         "Name of item, for display.  Usually 'sourceName/targetName/exp' or empty"
    uint score;          "Score (0-1000)"
    double value;        "Strength of interaction or other data value. Typically basis for score"
    string exp;          "Experiment name (metadata for filtering). Use . if not applicable"
    string color;        "Item color.  Specified as r,g,b or hexadecimal #RRGGBB or html color name, as in //www.w3.org/TR/css3-color/#html4. Use 0 and spectrum setting to shade by score"
    string sourceChrom;  "Chromosome of source region (directional) or lower region. For non-directional interchromosomal, chrom of this region."
    uint sourceStart;    "Start position in chromosome of source/lower/this region"
    uint sourceEnd;      "End position in chromosome of source/lower/this region"
    string sourceName;   "Identifier of source/lower/this region"
    string sourceStrand; "Orientation of source/lower/this region: + or -.  Use . if not applicable"
    string targetChrom;  "Chromosome of target region (directional) or upper region. For non-directional interchromosomal, chrom of other region"
    uint targetStart;    "Start position in chromosome of target/upper/this region"
    uint targetEnd;      "End position in chromosome of target/upper/this region"
    string targetName;   "Identifier of target/upper/this region"
    string targetStrand; "Orientation of target/upper/this region: + or -.  Use . if not applicable"

    )
' > ${locTMP}interact.as

awk '
  BEGIN {OFS="\t"} 
  {
      chrom      = $1
      chromStart = $2
      chromEnd   = $3
      junctionID = $4
      readCount  = $6
      strandStr  = $5
      
      if(readCount < 3) next
      if(strandStr != "+" )next

      # splice sites (donor = intron start, acceptor = intron end)
      sourceStart = chromStart
      sourceEnd   = chromStart + 1
      targetStart = chromEnd - 1
      targetEnd   = chromEnd

      # scale score
      score = (readCount > 1000 ? 1000 : readCount)

      # name for interact item
      name = chrom ":" sourceStart "-" sourceEnd "_" chrom ":" targetStart "-" targetEnd

      # color by strand
      if (strandStr == "+")  color = "0,0,255"
      else if (strandStr=="-") color = "255,0,0"
      else color="128,128,128"
      color=0

      print chrom, sourceStart, targetEnd, name, score, readCount, "regtools_SJ", color, \
            chrom, sourceStart, sourceEnd, "exon_donor", strandStr, \
            chrom, targetStart, targetEnd, "exon_acceptor", strandStr
  }' ${locTMP}junctions.collapsed.bed \
  | sort -k1,1 -k2,2n > ${locTMP}splice_junctions.interact


#add proper score 
US=1
  #convert to bed covering the 3 nt upstream of the splice junction
  mawk -v OFS="\t" -v US=$US '
  {
    print $1, $2-US, $2+1, $4, 0, $13
  }' ${locTMP}splice_junctions.interact | 
bedtools map -a stdin -b ${locTMP}unique_mappers.sense.bedGraph -c 4 -o mean -g ${CHRsizes} > ${locTMP}usCount.txt

mawk -v OFS="\t" -v INFILE=${locTMP}usCount.txt -v US=$US '
  BEGIN{
    while((getline < INFILE) > 0) {
      key = $4
      usCount[key] = $7
    }
  }
  {
    key = $4
    if(key in usCount) {
      us = usCount[key]
    } else {
      us = 0
    }
    if(us > 0) {
      score = int($6 / us * 1000)
      if(score > 1000) score = 1000
    } else {
      score = 500
    }
    $5=score
    $4=$4"_usCOUNT="us
    if(us > 10){
      print
    } 
  }' ${locTMP}splice_junctions.interact > ${locTMP}splice_junctions.scored.interact


  bedToBigBed -as=${locTMP}interact.as -type=bed5+13 -tab ${locTMP}splice_junctions.scored.interact ${CHRsizes} ${OPENdir}cluster-analysis/splice_junctions.interact

exit

#####
exit

###############################################################################
# bam_to_interact.sh
# Usage:  bam_to_interact.sh  input.bam  output.interact  [min_support]
#
# • Extracts splice-junctions (N operations in the CIGAR) directly from BAM
# • Builds a read-coverage bedGraph (±split) with bedtools genomecov
# • For every donor site, finds the median coverage in a ±50 bp window
# • SCORE = junction_count / median_cov   (if median_cov==0 → use max_junction_count)
# • Scaled to the UCSC 0–1000 range and written in full interact format
#
# Dependencies: samtools, bedtools, gawk  (GNU awk ≥4 for asort)
###############################################################################

set -euo pipefail

bam=$1
out=$2
min_support=${3:-2}          # default ≥2 supporting reads
window=50                    # ±50 bp window for median
tmp_junc=$(mktemp)
tmp_bg=$(mktemp)

###############################################################################
# 1.  Coverage bedGraph  (split-aware)
###############################################################################
echo "• Generating bedGraph coverage..."
bedtools genomecov -split -bg -ibam "$bam" > "$tmp_bg"

###############################################################################
# 2.  Extract junctions with counts  (samtools + gawk, no “length” var)
###############################################################################
echo "• Extracting junctions (min_support=${min_support})..."
samtools view -F 4 "$bam" | \
gawk -v min_support="$min_support" -v OFS="\t" '
function add_junc(ch,s,e){
    key = ch ":" s ":" e
    ++cnt[key]
    chr[key]=ch; start[key]=s; stop[key]=e
}
{
    chrom  = $3
    refPos = $4            # 1-based leftmost position
    cigar  = $6

    # iterate through CIGAR M/I/D/N/S/H/=/X tokens
    while (match(cigar, /[0-9]+[MIDNSHP=X]/, m)) {
        tok  = m[0]
        n    = substr(tok, 1, length(tok)-1) + 0
        op   = substr(tok, length(tok), 1)
        cigar = substr(cigar, RSTART + RLENGTH)   # trim processed part

        if      (op ~ /^[MD=X]$/) { refPos += n }        # consume reference
        else if (op == "N")  {                           # junction!
            s = refPos
            e = refPos + n
            add_junc(chrom, s, e)
            refPos += n
        }                                                # I/S/H/P do nothing
    }
}
END{
    for (k in cnt)
        if (cnt[k] >= min_support)
            print chr[k], start[k], stop[k], cnt[k]
}' > "$tmp_junc"

###############################################################################
# 3.  Global max junction count  (for zero-coverage fallback)
###############################################################################
max_cnt=$(awk '($4>m){m=$4} END{print (m?m:1)}' "$tmp_junc")
echo "• Max junction count = $max_cnt"

###############################################################################
# 4.  Build interact file  (load bedGraph once, then stream junctions)
###############################################################################
echo "• Writing interact track → $out"
gawk  -v bg="$tmp_bg" \
      -v max_cnt="$max_cnt" \
      -v win="$window" \
      -v OFS="\t" '
BEGIN{
    # ------------------------------------------------------------------
    # Load coverage into covArr[chrom,position]  (0-based, half-open)
    # ------------------------------------------------------------------
    while ((getline < bg) > 0) {
        c = $1; s = $2; e = $3; cov = $4
        for (pos=s; pos<e; pos++)  covArr[c, pos] = cov
    }
}
{
    chrom = $1;  start = $2;  end = $3;  cnt = $4
    # -------------------------------------------------------------- MEDIAN
    n=0
    for (p = start - win; p <= start + win; p++) {
        key = chrom SUBSEP p
        if ((key) in covArr)
            vals[n++] = covArr[key]
    }
    if (n==0) {
        med=0
    } else {
        asort(vals)                               # gawk built-in
        med = (n % 2) ? vals[int(n/2)+1] \
                      : (vals[n/2] + vals[n/2+1]) / 2
    }
    denom = (med>0 ? med : max_cnt)
    score = int((cnt/denom) * 1000)
    if (score > 1000) score = 1000
    # ---------------------------------------------------------- INTERACT
    name = "junction_" chrom "_" start "_" end
    print chrom, start, end, name, score, cnt, ".", "255,0,0", \
          chrom, start, start+1, chrom ":" start, ".", \
          chrom, end-1, end, chrom ":" end, "."
}
' "$tmp_junc" > "$out"

echo "• Done!  Upload '$out' as a custom track in the UCSC Genome Browser."

#------------------------------------------------------------------------------
# 5.  Cleanup
#------------------------------------------------------------------------------
rm -f "$tmp_junc" "$tmp_bg"