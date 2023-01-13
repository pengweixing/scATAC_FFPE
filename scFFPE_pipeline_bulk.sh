#################################################
#  File Name:cut-tag.sh
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Fri 30 Oct 2020 10:55:04 AM UTC
#################################################

#### For HG38 and mm10 #####


source /home/pengweixing/bashrc
export LD_LIBRARY_PATH=/disk1/pengweixing/software/preseq/htslib/lib:$LD_LIBRARY_PATH
echo `ls *_L001_R1_001.fastq.gz` > sample.list
# 1.QC####
if [ ! -d "QC" ];then
mkdir QC
fi
for each_R1 in `cat sample.list`
do
	/disk1/pengweixing/software/FastQC/fastqc -o QC -t 10 $each_R1 &
	sleep 0
done
wait
/disk1/pengweixing/software/anaconda3/bin/multiqc ./QC -n all_QC_report -o ./ 

refrence_genome=$1

if [ ARG"$2" == ARG ];
then processor=3
else
processor=$2
fi

#2. debarcode
if [ ! -d "sc_Debarcode" ];then
	mkdir sc_Debarcode
fi
cp /disk1/pengweixing/pipeline/scFFPE/barcode.txt ./
for each_R1 in `cat sample.list`
do
each_R2=`echo $each_R1|sed s#R1#R2#g`
out_name=`echo $each_R1 |awk -F '_L001' '{print $1}'`
/disk1/pengweixing/software/anaconda3/bin/python /disk1/pengweixing/pipeline/scFFPE/fastq_to_10x.PE75.py \
  -r1 $each_R1 \
  -r2 $each_R2 -b barcode.txt -o ./sc_Debarcode/$out_name \
  -d /disk1/pengweixing/software/cellranger-atac-2.0.0_modified/lib/python/atac/barcodes/1M.barcodes.txt &
done
wait

#2.triming###

echo "trimming adaptor"
for each_R1 in `cat sample.list`
do
out=`echo $each_R1 |sed s#.fastq.gz##g`
python /disk1/pengweixing/software/cut-tag/cut75.single.py $each_R1 $out.cut75.fastq.gz &
done
wait

for each_R1 in `cat sample.list`
do
	out=`echo $each_R1 |sed s#.fastq.gz##g`
/disk1/pengweixing/software/anaconda3/bin/python /disk1/pengweixing/software/cut-tag/pyadapter_trim.FFPE.multiprocess.py -i $out.cut75.fastq.gz -p $processor  &
done
wait
echo "trimming adaptor finished "

#3.mapping###
echo "mapping is running"
cores=20
echo $refrence_genome

if [[ "$refrence_genome" == "mm" ]]
then
ref=/disk1/pengweixing/database/mm10/bowtie2/mm10
elif [[ "$refrence_genome" == "hs" ]]
then
ref=/disk1/pengweixing/database/hg38/index/hg38.fa
fi
if [ ! -d "Mapping" ];then
mkdir Mapping
fi
if [ ! -d "Library_complexity" ];then
mkdir Library_complexity
fi
currdir=`pwd`
for each_R1 in `cat sample.list`
do
#	each_R2=`echo $each_R1 |sed 's/R1/R2/'`
#	each_R2_trim=`echo $each_R2 |sed 's#.fastq.gz#.trim.fastq.gz#g'`
	each_R1_trim=`echo $each_R1 |sed 's#.fastq.gz#.cut75.fastq.trim.gz#g'`
	out_name=`echo $each_R1 |awk -F '_L001' '{print $1}'`
	bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${cores} -x ${ref} -U $each_R1_trim -S $currdir/Mapping/${out_name}_bowtie2.sam &> $currdir/Mapping/${out_name}_bowtie2.txt
	echo "mapping finished"
	samtools view --threads 10 -S -b $currdir/Mapping/${out_name}_bowtie2.sam >  $currdir/Mapping/${out_name}_bowtie2.bam
	rm $currdir/Mapping/${out_name}_bowtie2.sam
	samtools sort -m 10G --threads 10 $currdir/Mapping/${out_name}_bowtie2.bam -o  $currdir/Mapping/${out_name}_bowtie2.sort.bam
	##### Library complexity ######
	reads_number=`samtools flagstat $currdir/Mapping/${out_name}_bowtie2.sort.bam |sed -n '1p' |awk '{print $1}'`
	echo $reads_number
	ss=`echo "scale=0;$reads_number/10" |bc`
	ee=`echo "scale=0;$reads_number*5000"|bc`
	if [ $ss -eq 0 ];
	then
	ss=1000
	fi
	/disk1/pengweixing/software/preseq/preseq-3.1.2/build/bin/preseq c_curve -s $ss -o $currdir/Library_complexity/${out_name}.curve.txt -B $currdir/Mapping/${out_name}_bowtie2.sort.bam
	/disk1/pengweixing/software/preseq/preseq-3.1.2/build/bin/preseq lc_extrap -s $ss -e $ee -o $currdir/Library_complexity/${out_name}.expect.txt -B $currdir/Mapping/${out_name}_bowtie2.sort.bam
	/disk1/pengweixing/software/anaconda3/bin/python /disk1/pengweixing/software/cut-tag/preseq_plot.py $currdir/Library_complexity/${out_name}.curve.txt $currdir/Library_complexity/${out_name}.expect.txt $currdir/Library_complexity/${out_name}_Lib_complexity.pdf ${out_name}

	###### duplicates removing####
	java -Xmx4G -jar /home/xingqichen/SOFTWARE/picard-tools-1.119/MarkDuplicates.jar INPUT=$currdir/Mapping/${out_name}_bowtie2.sort.bam OUTPUT=$currdir/Mapping/${out_name}_bowtie2.sort.rmdup.bam METRICS_FILE=$currdir/Mapping/${out_name}_Picard_Metrics_unfiltered_bam.txt VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=true &> $currdir/Mapping/${out_name}.Picard.log
	/disk1/pengweixing/software/sambamba/sambamba index -t 10 $currdir/Mapping/${out_name}_bowtie2.sort.bam
## Filte#r and keep the mapped read pairs
	samtools view -h -q 2 -F 0x4 $currdir/Mapping/${out_name}_bowtie2.sort.rmdup.bam |awk '$3!="chrM"' |samtools view -b - > $currdir/Mapping/${out_name}.q2.sort.rmdup.bam
	samtools flagstat $currdir/Mapping/${out_name}.q2.sort.rmdup.bam > $currdir/Mapping/${out_name}.q2.sort.rmdup.bam.stat
done

######6. peak calling #####
if [ ! -d "Peak_calling" ];then
mkdir Peak_calling
fi
for each_R1 in `cat sample.list`
do
out_name=`echo $each_R1 |awk -F '_L001' '{print $1}'`
bedtools genomecov -bg -ibam $currdir/Mapping/${out_name}.q2.sort.rmdup.bam > $currdir/Mapping/${out_name}_bowtie2.normalized.bedgraph
#bash /disk1/pengweixing/software/SEACR/SEACR_1.3.sh $currdir/Mapping/${out_name}_bowtie2.normalized.bedgraph 0.01 non stringent $currdir/Peak_calling/${out_name}_seacr_top0.01.peaks
bedtools bamtobed -i  $currdir/Mapping/${out_name}.q2.sort.rmdup.bam  >  $currdir/Mapping/${out_name}.q2.sort.rmdup.bed
/disk1/pengweixing/software/anaconda3/bin/sicer -t $currdir/Mapping/${out_name}.q2.sort.rmdup.bed -s hg38 --window_size 200 --gap_size 400 --output_directory $currdir/Peak_calling/
done

#Rscript /disk1/pengweixing/software/cut-tag/peak_stat.r sample.name 
/disk1/pengweixing/software/anaconda3/bin/python /disk1/pengweixing/pipeline/scFFPE/stat3.py sample.list $currdir > QCtable.txt
######7. Tss enrichment #####
if [ ! -d "bigwig" ];then
mkdir bigwig
fi
if [ ! -d "TSS_enrich" ];then
mkdir TSS_enrich
fi

echo $refrence_genome

if [[ "$refrence_genome" == "mm" ]]
then
refanno=/disk1/pengweixing/database/mm10/mm10.refGene.gtf
elif [[ "$refrence_genome" == "hs" ]]
then
refanno=/disk1/pengweixing/database/hg38/gencode.v29.primary_assembly.annotation_UCSC_names.gtf
fi

for each_R1 in `cat sample.list`
do
out_name=`echo $each_R1 |awk -F '_L001' '{print $1}'`
echo $out_name
/disk1/pengweixing/software/sambamba/sambamba index -t 10 $currdir/Mapping/${out_name}.q2.sort.rmdup.bam
cores=20
/disk1/pengweixing/software/miniconda3/bin//bamCoverage --numberOfProcessors $cores -b $currdir/Mapping/${out_name}.q2.sort.rmdup.bam -o $currdir/bigwig/${out_name}.q2.bw   --normalizeUsing CPM
/disk1/pengweixing/software/miniconda3/bin/computeMatrix scale-regions -S  ./bigwig/${out_name}.q2.bw \
                              -R $refanno \
                              --beforeRegionStartLength 3000 \
                              --regionBodyLength 2000 \
                              --afterRegionStartLength 3000 \
			      --binSize 10 \
			      --missingDataAsZero \
			      --sortRegions descend \
                              --skipZeros -o ./bigwig/${out_name}.matrix_gene.mat.gz -p $cores
#plotProfile --matrixFile ./bigwig/${out_name}.matrix_gene.mat.gz  --outFileName ./TSS_enrich/${out_name}.tss.enrich.pdf 
/disk1/pengweixing/software/miniconda3/bin/plotHeatmap --matrixFile ./bigwig/${out_name}.matrix_gene.mat.gz --sortRegions no --outFileName ./TSS_enrich/${out_name}.scale-regions.enrich.pdf  --colorMap Blues  --heatmapHeight 12 --legendLocation upper-left 
done

for each_R1 in `cat sample.list`
do
out_name=`echo $each_R1 |awk -F '_L001' '{print $1}'`
echo $out_name
cores=20
/disk1/pengweixing/software/miniconda3/bin/computeMatrix reference-point -S  ./bigwig/${out_name}.q2.bw \
			      -R $refanno \
			      --beforeRegionStartLength 3000 \
			      --afterRegionStartLength 3000 \
			      --binSize 10 \
			      --missingDataAsZero \
			      --sortRegions descend \
			      --skipZeros -o ./bigwig/${out_name}.matrix_gene.reference-point.mat.gz -p $cores
/disk1/pengweixing/software/miniconda3/bin/plotHeatmap --matrixFile ./bigwig/${out_name}.matrix_gene.reference-point.mat.gz --outFileName ./TSS_enrich/${out_name}.tss.reference-point.heat.enrich.pdf --sortRegions no  --colorMap Blues  --heatmapHeight 12 --legendLocation upper-left 
done
rm *.cut75.fastq.gz
