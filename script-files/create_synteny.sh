#!/bin/bash

#SBATCH --cpus-per-task=20
#SBATCH --mem=40g
#SBATCH --partition=c
#SBATCH -e "%x.e.%j.txt"
#SBATCH -o "%x.o.%j.txt"
#SBATCH --qos=medium
#SBATCH --time=20:00:00


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
locTMP=${TMPdir}/synteny/
mkdir -p ${locTMP}

source ${SCRIPTdir}tools

OPENdir=${OPENdir}/synteny/
mkdir -p ${OPENdir}
#***********************************no commenting above************************************************************

###################################################################################################
#filter main chromosomes

if [[ ! -s ${locTMP}ref_main.fa ]]; then
  seqkit grep -r -p "^2L$|^2R$|^3L$|^3R$|^X$|mito|^4$" $refFASTA | seqkit grep -r -v -p "mapped" | seqkit fx2tab | mawk -v OFS="\t" '{ $1="chr"$1; print }' | seqkit tab2fx --line-width 0 > ${locTMP}ref_main.fa
  seqkit fx2tab $refFASTA | mawk -v OFS="\t" '{ $1="chr"$1; print }' | seqkit tab2fx --line-width 0 > ${locTMP}ref_main.fa
fi

nASSEMBLY=$(wc -l ${assemblyFILE} | cut -f 1 -d ' ')

#get dot into the output directory
if [[ ! -d ${OPENdir}dot ]]; then
  cd $OPENdir
  git clone https://github.com/dnanexus/dot.git
fi

ml  build-env/f2022
ml  r/4.5.1-gfbf-2023b
for currASSEMBLY in $(seq 1 $nASSEMBLY); do
  currASSEMBLYline=$(sed -n ${currASSEMBLY}p ${assemblyFILE})
  currASSEMBLYseq=$(echo $currASSEMBLYline |tr ' ' '\t' | cut -f 1)  #remove file extension
  currASSEMBLYname=$(echo $currASSEMBLYline |tr ' ' '\t' | cut -f 2)  #remove file extension
  currASSEMBLYflam=$(echo $currASSEMBLYline |tr ' ' '\t' | cut -f 3)  #remove file extension

  if [[  ! -s ${locTMP}${currASSEMBLYname}.fa ]]; then
    seqkit seq --line-width 0 ${currASSEMBLYseq}  > ${locTMP}${currASSEMBLYname}.fa
   fi

  #align genomes using minimap
minimap2 -x asm10 -t ${CORES} \
  --secondary=yes  -N 100 \
  ${locTMP}ref_main.fa ${locTMP}${currASSEMBLYname}.fa > ${locTMP}aligned-genomes.dm6_${currASSEMBLYname}.paf

  paf2dotplot.r  -o ${OPENdir}${currASSEMBLYname}_vs_dm6-all.dot -fs -q 10000 -m 1000 --min-ref-len 10000 ${locTMP}aligned-genomes.dm6_${currASSEMBLYname}.paf
  paf2dotplot.r  -o ${OPENdir}${currASSEMBLYname}_vs_dm6-all.dot.filtered -fs  -m 20000 -q 500000 --min-ref-len 10000 ${locTMP}aligned-genomes.dm6_${currASSEMBLYname}.paf

 minimap2 -x asm5 -t ${CORES}  ${locTMP}ref_main.fa ${locTMP}${currASSEMBLYname}.fa > ${locTMP}aligned-genomes.dm6_${currASSEMBLYname}.siomi.paf
 
  paf2dotplot.r  -o ${OPENdir}${currASSEMBLYname}_vs_dm6.siomi.dot -fs -q 10000 -m 1000  ${locTMP}aligned-genomes.dm6_${currASSEMBLYname}.siomi.paf
  paf2dotplot.r  -o ${OPENdir}${currASSEMBLYname}_vs_dm6.dot.siomi.filtered -fs  -m 20000 -q 500000 ${locTMP}aligned-genomes.dm6_${currASSEMBLYname}.siomi.paf

  minimap2 -ax asm5 -t ${CORES} --eqx ${locTMP}ref_main.fa ${locTMP}${currASSEMBLYname}.fa |
    samtools sort -O BAM -@ $CORES - > ${locTMP}aligned-genomes.dm6_${currASSEMBLYname}.bam

  cp  ${locTMP}aligned-genomes.dm6_${currASSEMBLYname}.bam  ${OPENdir}aligned-genomes.dm6_${currASSEMBLYname}.bam
  samtools index ${OPENdir}aligned-genomes.dm6_${currASSEMBLYname}.bam
  nucmer -p ${locTMP}${currASSEMBLYname}_vs_dm6 --maxmatch -c 100 -b 500 -l 50 ${locTMP}${currASSEMBLYname}.fa ${locTMP}ref_main.fa  --threads $CORES      # Whole genome alignment. Any other alignment can also be used.
  deltaFilter -m -i 90 -l 100 ${locTMP}${currASSEMBLYname}_vs_dm6.delta > ${locTMP}${currASSEMBLYname}_vs_dm6.filtered.delta     # Remove small and lower quality alignments
  showCoords -THrd ${locTMP}${currASSEMBLYname}_vs_dm6.filtered.delta > ${locTMP}${currASSEMBLYname}_vs_dm6.filtered.coords      # Convert alignment information to a .TSV 


  #create dot input files
  python ${SCRIPTdir}DotPrep.py --delta ${locTMP}${currASSEMBLYname}_vs_dm6.delta --out ${OPENdir}${currASSEMBLYname}_vs_dm6 --overview 10000

  # Running syri for finding structural rearrangements between A and B
  syri -c ${locTMP}aligned-genomes.dm6_${currASSEMBLYname}.bam -r ${locTMP}${currASSEMBLYname}.fa -q  ${locTMP}ref_main.fa -F B --prefix ${currASSEMBLYname}_vs_dm6. --dir ${locTMP} --nc $CORES --tdmaxolp 0.9 --tdgaplen 100000 --unic 100 --unip 0.9  --inc 100 --all

  printf '#file	name	tags\n'  | tr ' ' '\t' > ${locTMP}genomes.txt
  echo  ${locTMP}${currASSEMBLYname}.fa ${currASSEMBLYname} lw:1.5 | tr ' ' '\t' >> ${locTMP}genomes.txt
  echo ${locTMP}ref_main.fa dm6 lw:1.5 | tr ' ' '\t' >> ${locTMP}genomes.txt

  plotsr \
      --sr ${locTMP}${currASSEMBLYname}_vs_dm6.syri.out \
      --genomes ${locTMP}genomes.txt \
      -o ${OPENdir}${currASSEMBLYname}_vs_dm6.synteny.png

  plotsr \
      --sr ${locTMP}${currASSEMBLYname}_vs_dm6.syri.out \
      --genomes ${locTMP}genomes.txt --reg OSC:chrX_RagTag:21666842-24266842 \
      -o ${OPENdir}${currASSEMBLYname}_vs_dm6.synteny.zoom.png -s 100 --rtr


done

#chrX dot-plot



paftools delta2paf ${locTMP}OSC_r1.01_vs_dm6.filtered.delta | cut -f 1-12> ${locTMP}aligned-genomes.dm6_OSC_r1.01.fromDelta.paf

#determine uniq or multiple mappings
cut -f1,3,4 ${locTMP}aligned-genomes.dm6_OSC_r1.01.fromDelta.paf | sort | uniq -c > ${locTMP}counts.1.txt
cut -f6,8,9 ${locTMP}aligned-genomes.dm6_OSC_r1.01.fromDelta.paf | sort | uniq -c > ${locTMP}counts.2.txt
mawk -v OFS="\t" -v COUNTin_1=${locTMP}counts.1.txt -v COUNTin_2=${locTMP}counts.2.txt -v TMP=${locTMP} '
BEGIN{
  while( (getline < COUNTin_1) > 0 ) {
    COUNT1[$2"~"$3"~"$4]=$1
  }
  while( (getline < COUNTin_2) > 0 ) {
    COUNT2[$2"~"$3"~"$4]=$1
  }
  print "queryID\tqueryLen\tqueryStart\tqueryEnd\tstrand\trefID\trefLen\trefStart\trefEnd\tnumResidueMatches\tlenAln\tmapQ\tREP" > TMP "aligned-genomes.dm6_OSC_r1.01.repeatANN.paf"
  print "queryID\tqueryLen\tqueryStart\tqueryEnd\tstrand\trefID\trefLen\trefStart\trefEnd\tnumResidueMatches\tlenAln\tmapQ\tREP" > TMP "aligned-genomes.dm6_OSC_r1.01.chrX.repeatANN.paf"
}
{
  TAG="uniq"
  if($1"~"$3"~"$4 in COUNT1){
    if(COUNT1[$1"~"$3"~"$4]>1){
      TAG="multi"
    }
    if($6"~"$8"~"$9 in COUNT2){
      if(COUNT2[$6"~"$8"~"$9]>1){
        TAG="multi"
      }
    }
    print $6,$7,$8,$9,$5,$1,$2,$3,$4,$10,$11,$12,TAG 
    print $0, TAG
  }else{
    print $0, "error"
  }
  
}' ${locTMP}aligned-genomes.dm6_OSC_r1.01.fromDelta.paf > ${locTMP}aligned-genomes.dm6_OSC_r1.01.repeatANN.paf

paf2dotplot/paf2dotplot.repeatColoring.r -o ${OPENdir}OSC_r1.01_vs_dm6.full-genome.delta-filtered.dot  -fs -q 100000 -m 1000 --min-ref-len 100000 ${locTMP}aligned-genomes.dm6_OSC_r1.01.repeatANN.paf --plot-size 8

exit
#extract chrX for Figure 4
grep chrX ${locTMP}aligned-genomes.dm6_OSC_r1.01.repeatANN.paf | grep X_RagTag >> ${locTMP}aligned-genomes.dm6_OSC_r1.01.chrX.repeatANN.paf
echo X_RagTag 22666842 24081114 flam_OSC 0 + | tr ' ' '\t' > ${locTMP}flam_OSC.bed
echo chrX 21624796 22447808 flam_dm6 0 + | tr ' ' '\t' > ${locTMP}dm6_OSC.bed
paf2dotplot/paf2dotplot.repeatColoring.r -o ${OPENdir}OSC_r1.01_vs_dm6.chrX.dot -e ${locTMP}dm6_OSC.bed -E ${locTMP}flam_OSC.bed -fs -q 1000 -m 1000 --min-ref-len 1000 ${locTMP}aligned-genomes.dm6_OSC_r1.01.chrX.repeatANN.paf --plot-size 8

exit


###################################################################################################
#flamenco alignments
  BRE_ASSEMBLY=$(grep -n OSC_r1.01 ${assemblyFILE} | cut -f 1 -d ':')
  BRE_ASSEMBLYline=$(sed -n ${BRE_ASSEMBLY}p ${assemblyFILE})
  BRE_ASSEMBLYseq=$(echo $BRE_ASSEMBLYline |tr ' ' '\t' | cut -f 1)  #remove file extension
  BRE_ASSEMBLYname=$(echo $BRE_ASSEMBLYline |tr ' ' '\t' | cut -f 2)  #remove file extension
  BRE_ASSEMBLYflam=$(echo $BRE_ASSEMBLYline |tr ' ' '\t' | cut -f 3)  #remove file extension



echo chrX 21624796 22252390 flam_dm6 0 + | tr ' ' '\t' > ${locTMP}flam_dm6.bed
bedtools getfasta -name -fi ${locTMP}ref_main.fa -bed ${locTMP}flam_dm6.bed -fo ${locTMP}flam_dm6.fasta

echo X_RagTag 22666842 23402245 flam_OSC 0 + | tr ' ' '\t' > ${locTMP}flam_OSC.bed
bedtools getfasta -name -fi $BRE_ASSEMBLYseq -bed ${locTMP}flam_OSC.bed -fo ${locTMP}flam_OSC.fasta

cd $locTMP
nucmer -p ${locTMP}out --maxmatch -c 100 -b 500 -l 50 ${locTMP}flam_OSC.fasta ${locTMP}flam_dm6.fasta       # Whole genome alignment. Any other alignment can also be used.
deltaFilter -m -i 90 -l 100 ${locTMP}out.delta > ${locTMP}out.filtered.delta     # Remove small and lower quality alignments
showCoords -THrd ${locTMP}out.filtered.delta > ${locTMP}out.filtered.coords      # Convert alignment information to a .TSV format as required by SyRI

syri -c ${locTMP}out.filtered.coords -d ${locTMP}out.filtered.delta -r ${locTMP}flam_OSC.fasta -q ${locTMP}flam_dm6.fasta --nosnp --dir ${locTMP} --tdgaplen 10000

printf '#file	name	tags\n'  | tr ' ' '\t' > ${locTMP}genomes.txt
echo ${locTMP}flam_OSC.fasta OSC lw:1.5 | tr ' ' '\t' >> ${locTMP}genomes.txt
echo ${locTMP}flam_dm6.fasta dm6 lw:1.5 | tr ' ' '\t' >> ${locTMP}genomes.txt

plotsr \j
    --sr ${locTMP}syri.out \
    --genomes ${locTMP}genomes.txt  -s 5000 \
    -o ${OPENdir}output_plot.flamenco.nucmer2.pdf 



for HAP in hap2 ; do
  currASSEMBLY=$(grep -n $HAP ${assemblyFILE} | cut -f 1 -d ':')
  currASSEMBLYline=$(sed -n ${currASSEMBLY}p ${assemblyFILE})
  currASSEMBLYseq=$(echo $currASSEMBLYline |tr ' ' '\t' | cut -f 1)  #remove file extension
  currASSEMBLYname=$(echo $currASSEMBLYline |tr ' ' '\t' | cut -f 2)  #remove file extension
  currASSEMBLYflam=$(echo $currASSEMBLYline |tr ' ' '\t' | cut -f 3)  #remove file extension

  
  echo $currASSEMBLYflam | 
    mawk -v OFS="\t" -v HAP=$HAP '{split($1,splitCOORD,/:|-/); print splitCOORD[1],splitCOORD[2],splitCOORD[3],HAP}' > ${locTMP}${HAP}.bed
  bedtools getfasta -name -fi $currASSEMBLYseq -bed ${locTMP}${HAP}.bed -fo ${locTMP}${HAP}.fasta

  cd $locTMP
  nucmer -p ${locTMP}${HAP}_out --maxmatch -c 100 -b 500 -l 50 ${locTMP}flam_OSC.fasta ${locTMP}${HAP}.fasta      # Whole genome alignment. Any other alignment can also be used.
  deltaFilter -m -i 90 -l 100 ${locTMP}${HAP}_out.delta > ${locTMP}${HAP}_out.filtered.delta     # Remove small and lower quality alignments
  showCoords -THrd ${locTMP}${HAP}_out.filtered.delta > ${locTMP}${HAP}_out.filtered.coords      # Convert alignment information to a .TSV format as required by SyRI

  syri -c ${locTMP}${HAP}_out.filtered.coords -d ${locTMP}${HAP}_out.filtered.delta -r ${locTMP}flam_OSC.fasta -q ${locTMP}${HAP}.fasta --nosnp --dir ${locTMP} --tdgaplen 10000 --prefix ${HAP}_

  printf '#file	name	tags\n'  | tr ' ' '\t' > ${locTMP}genomes.txt
  echo ${locTMP}flam_OSC.fasta OSC lw:1.5 | tr ' ' '\t' >> ${locTMP}genomes.txt
  echo ${locTMP}${HAP}.fasta $HAP lw:1.5 | tr ' ' '\t' >> ${locTMP}genomes.txt

  plotsr \
      --sr ${locTMP}${HAP}_syri.out \
      --genomes ${locTMP}genomes.txt  -s 5000 \
      -o ${OPENdir}output_plot.flamenco.nucmer.Siomi_${HAP}.pdf 

done


exit


###################################################################################################

