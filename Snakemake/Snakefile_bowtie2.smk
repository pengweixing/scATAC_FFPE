import os
with open('sample.group','r') as f:
    sampledict = {}
    all_line = [line.strip() for line in f.readlines()]
    ALL_Groups = set()
    for a in all_line:
        group = a.split()[0]
        if '_' in group:
            group = group.replace('_','.')
        ALL_Groups.add(group)
        tn5 = a.split()[1]
        if not group in sampledict:
            sampledict[group] = [tn5]
        else:
            sampledict[group].append(tn5)

CWD=os.getcwd()

rule all:
    input:
        expand("02.mapping/{groupname}.sort.markdup.merge.bam",groupname = ALL_Groups),
        expand("02.mapping/{groupname}.sort.markdup.merge.stat",groupname = ALL_Groups),
        expand("05.Peakcalling/{groupname}_peaks.filterBL.bed",groupname = ALL_Groups),
        expand("04.TSS/{groupname}.tss.reference-point.heat.enrich.pdf",groupname = ALL_Groups),
        expand("06.QC/{groupname}.singlecell.qc.txt",groupname = ALL_Groups),
        expand("06.QC/{groupname}_fragments_barplot.pdf",groupname = ALL_Groups),
        expand("06.QC/{groupname}_violinplot.pdf",groupname = ALL_Groups),
        expand("06.QC/{groupname}_Unique_FRiP.scatterplot.pdf",groupname = ALL_Groups),
        expand("07.singlecell/{groupname}.FRiP.selected.singlecell.fragments.gz",groupname = ALL_Groups),
        expand("07.singlecell/{groupname}.FRiT.selected.singlecell.fragments.gz",groupname = ALL_Groups),
        expand("07.singlecell/{groupname}.cell_peak.matrix",groupname = ALL_Groups),
        expand("07.singlecell/{groupname}.FRiP.selected.singlecell.bam",groupname = ALL_Groups),
        expand("07.singlecell/{groupname}.FRiT.selected.singlecell.bam",groupname = ALL_Groups),
        
   #     expand("All.QC.txt",groupname = ALL_Groups)
    output:
        "All.QC.txt"
    params:
        getbulkqc=config['scripts']['getbulkqc'],
        python_main=config['scripts']['python'],
    shell:
        "{params.python_main} {params.getbulkqc} {CWD}"

rule bwa_map:
    message: "do mapping"
    input:
        ref=config['ref'],
        r1="01.debarcoding/{sample}_S1_L001_R1_001.fastq.gz",
        r2="01.debarcoding/{sample}_S1_L001_R2_001.fastq.gz"
    
    threads: 20
    log: "log/{sample}.mem.log"

    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}"

    output:
        bam="02.mapping/{sample}.sort.bam"
    shell:  
        "bowtie2 --very-sensitive --no-discordant --no-mixed -p {threads} -X 1000 -x /disk1/pengweixing/database/mm10/index/mm10 -1 {input.r1} -2 {input.r2} |samtools view -@ 5 - -b |samtools sort -m 2G -@ 5 - -o {output.bam} >& {log}"
rule index:
    message: "indexing bam file"
    input:
        "02.mapping/{sample}.sort.bam"
    log: "log/{sample}.index.log"    
    threads: 20
    output:
        "02.mapping/{sample}.sort.bam.bai"
    shell:
        "samtools index -@ 5 {input}"

rule markdup:
    message: "mark duplicates"
    input:
        sortbam="02.mapping/{sample}.sort.bam",
        bai="02.mapping/{sample}.sort.bam.bai"
    log:
        "log/{sample}.markdup.log"
    output:
        markbam=temp("02.mapping/{sample}.sort.markdup.bam")
    threads: 3
    params:
        MarDup=config['scripts']['markdup'],
        python_main=config['scripts']['python']
    shell:
        "{params.python_main} {params.MarDup} -i {input.sortbam} -o {output.markbam} -n {wildcards.sample} -t pair >& {log}"

rule Merge:
    input:
         lambda wildcards: ["02.mapping/{sample}.sort.markdup.bam".format(sample=bam) for bam in sampledict[wildcards.groupname]]
    output:
        bam="02.mapping/{groupname}.sort.markdup.merge.bam",
        bai="02.mapping/{groupname}.sort.markdup.merge.bam.bai"
    log:
        "log/{groupname}.markdup.merge.log"
    threads: 10
    shell:
        "samtools merge -@ {threads} {output.bam} {input} >& {log} && samtools index -@ {threads} {output.bam}"

rule flagstat:
    message: "flagstat Merged bam file"
    input:
        "02.mapping/{groupname}.sort.markdup.merge.bam"
    threads: 10
    output:
        "02.mapping/{groupname}.sort.markdup.merge.stat"
    shell:
        "sambamba flagstat -t {threads} {input} > {output}"

rule RmDup_and_Filter_LQ:
    message: "remove duplicates and filter low mapping quality and extract proper paired reads"
    input:
        bam="02.mapping/{groupname}.sort.markdup.merge.bam",
        bai="02.mapping/{groupname}.sort.markdup.merge.bam.bai"
    output:
        bam="02.mapping/{groupname}.sort.rmdup.merge.pair.bam",
        flagstat="02.mapping/{groupname}.sort.rmdup.merge.pair.bam.flagstat",
        bai="02.mapping/{groupname}.sort.rmdup.merge.pair.bam.bai"
    threads: 5
    params:
        Quality=config['quality']
    shell:
        "samtools view -f 0x2 -q {params.Quality} -F 3072 {input.bam} -b > {output.bam} && samtools index -@ 5 {output.bam}"
        " && sambamba flagstat -t 5 {output.bam} > {output.flagstat}"

rule Bigwig:
    message: "generate bigwig"
    input:
        bam="02.mapping/{groupname}.sort.rmdup.merge.pair.bam",
        bai="02.mapping/{groupname}.sort.rmdup.merge.pair.bam.bai"
    output:
        "03.bigwig/{groupname}.sort.rmdup.merge.pair.bw"
    threads: 40
    params:
        bamCoverage=config['scripts']['bamCoverage'],
    shell:
        "{params.bamCoverage} --bam {input.bam} --numberOfProcessors {threads} --outFileName {output} --normalizeUsing CPM"

rule TSS_enrichment:
    message: "generate TSS enrichment"
    input:
        "03.bigwig/{groupname}.sort.rmdup.merge.pair.bw"
    output:
        matgz="04.TSS/{groupname}.matrix_gene.mat.gz",
        tsspdf="04.TSS/{groupname}.tss.reference-point.heat.enrich.pdf"
    threads: 20
    params:
        anno=config['refanno'],
        computeMatrix=config['scripts']['computeMatrix'],
        plotHeatmap=config['scripts']['plotHeatmap']
    shell:
        "{params.computeMatrix} reference-point -S {input}"
                               " -R {params.anno}"
                               " --beforeRegionStartLength 3000"
                               " --afterRegionStartLength 3000"
                               " --binSize 10"
                               " --missingDataAsZero"
                               " --sortRegions descend" 
                               " --skipZeros -o {output.matgz} -p {threads} && "
        "{params.plotHeatmap} --matrixFile {output.matgz}"
        " --outFileName {output.tsspdf} --sortRegions no"
        " --colorMap Blues  --heatmapHeight 12 --legendLocation upper-left"
 
rule Sortbyname:
    input:
        bam="02.mapping/{groupname}.sort.rmdup.merge.pair.bam",
        bai="02.mapping/{groupname}.sort.rmdup.merge.pair.bam.bai"
    threads: 5
    output:
        bam=temp("02.mapping/{groupname}.sortbyname.rmdup.merge.pair.bam"),
    shell:
        "sambamba sort -n -t 5 {input.bam} -o {output.bam}"


rule BamToFragments:
    input:
        bam="02.mapping/{groupname}.sortbyname.rmdup.merge.pair.bam",
    log:
        "log/{groupname}.bamtofragments.log"
    output:
        frag=temp("02.mapping/{groupname}.rmdup.merge.pair.fragments"),
        fragments=temp("02.mapping/{groupname}.rmdup.merge.pair.sort.fragments"),
        fragments2="02.mapping/{groupname}.rmdup.merge.pair.sort.fragments.gz",
        bam=temp("02.mapping/{groupname}.rmdup.merge.pair.singlcell.bam"),
        bam2="02.mapping/{groupname}.rmdup.merge.pair.singlcell.sort.bam"
    params:
        BamtoFrag=config['scripts']['BamtoFrag'],
        python_main=config['scripts']['python']
    shell:
        "{params.python_main} {params.BamtoFrag} -i {input.bam} -of {output.frag} -ob {output.bam} && bedtools sort -i {output.frag} > {output.fragments}"
        " && sambamba sort {output.bam} -o {output.bam2} && bgzip -c {output.fragments} > {output.fragments2}"

rule Peakcalling:
    message: "run peak calling"
    input:
        "02.mapping/{groupname}.rmdup.merge.pair.sort.fragments" 
    output:
        peaks="05.Peakcalling/{groupname}_peaks.narrowPeak",
        peakBL="05.Peakcalling/{groupname}_peaks.filterBL.bed"
    threads: 1
    params:
        aa="--nomodel  --shift 0 -q 0.01",
        gsize=config['gisze'],
        blacklist=config['blacklist']
    shell:
        "macs2  callpeak -t {input} -f BED"
        " -g {params.gsize} --outdir 05.Peakcalling/ -n {wildcards.groupname} {params.aa} && bedtools intersect "
        "-a 05.Peakcalling/{wildcards.groupname}_peaks.narrowPeak -b {params.blacklist} -v > 05.Peakcalling/{wildcards.groupname}_peaks.filterBL.bed && "
        " bedtools intersect -a 05.Peakcalling/{wildcards.groupname}_summits.bed -b {params.blacklist} -v > 05.Peakcalling/{wildcards.groupname}_summits.filterBL.bed &&"
        " sort -k5,5nr 05.Peakcalling/{wildcards.groupname}_summits.filterBL.bed > 05.Peakcalling/{wildcards.groupname}_summits.filterBL.sort.bed"

rule Get_FRiP_bam:
    message: "get the bam file within peak region"
    input:
        bam="02.mapping/{groupname}.sort.rmdup.merge.pair.bam",
        bai="02.mapping/{groupname}.sort.rmdup.merge.pair.bam.bai",
        peakBL="05.Peakcalling/{groupname}_peaks.filterBL.bed"
    output:
        bam=temp("02.mapping/{groupname}.sort.rmdup.merge.pair.withinpeaks.bam"),
        bai=temp("02.mapping/{groupname}.sort.rmdup.merge.pair.withinpeaks.bam.bai")
    shell:
        "samtools view -L {input.peakBL} {input.bam} -b > {output.bam} && samtools index {output.bam}"

rule Get_FRiT_bam:
    message: "get the bam file within tss region"
    input:
        bam="02.mapping/{groupname}.sort.rmdup.merge.pair.bam",
        bai="02.mapping/{groupname}.sort.rmdup.merge.pair.bam.bai",
        peakBL=config['tssbed']
    output:
        bam="02.mapping/{groupname}.sort.rmdup.merge.pair.withintss.bam",
        bai="02.mapping/{groupname}.sort.rmdup.merge.pair.withintss.bam.bai"
    shell:
        "samtools view -L {input.peakBL} {input.bam} -b > {output.bam} && samtools index {output.bam}"

rule Get_SingleCellQC:
    message: "get single cell qc statistics"
    input:
        bam="02.mapping/{groupname}.sort.markdup.merge.bam",
        bai="02.mapping/{groupname}.sort.markdup.merge.bam.bai",
        peakbam="02.mapping/{groupname}.sort.rmdup.merge.pair.withinpeaks.bam",
        tssbam="02.mapping/{groupname}.sort.rmdup.merge.pair.withintss.bam"
    log: 
        "log/{groupname}.GetSingleCellQC.log"
    output:
        txt="06.QC/{groupname}.singlecell.qc.txt"
    params:
        getQC=config['scripts']['getscQC'],
        python_main=config['scripts']['python'],
        Quality=config['quality']
    shell:
         "{params.python_main} {params.getQC} -b1 {input.bam} -b2 {input.peakbam} -b3 {input.tssbam} -o {output.txt} -q {params.Quality} >& {log}"

rule Plot_QC:
    message: "plot the single cell qc data with barplot and scatterplot"
    input:
        "06.QC/{groupname}.singlecell.qc.txt"
    output:
        bar="06.QC/{groupname}_fragments_barplot.pdf",
        violin="06.QC/{groupname}_violinplot.pdf",
        scatterplot="06.QC/{groupname}_Unique_FRiP.scatterplot.pdf",
        scatterplot2="06.QC/{groupname}_Unique_FRiT.scatterplot.pdf",
        cellnames="06.QC/{groupname}_FRiP_selected.cell.names",
        cellnames2="06.QC/{groupname}_FRiT_selected.cell.names"
    params:
        Rscript=config['scripts']['Rscript'],
        plotqc=config['scripts']['plotqc'],
    shell:
        "{params.Rscript} {params.plotqc} -i {input} -p 5 -t 10 -n 500 -o 06.QC/{wildcards.groupname} "

rule Extract_Cell:
    message: "extract the possible cells based on qc table"
    input:
        cellnames="06.QC/{groupname}_FRiP_selected.cell.names",
        cellnames2="06.QC/{groupname}_FRiT_selected.cell.names",
        fragments="02.mapping/{groupname}.rmdup.merge.pair.sort.fragments.gz",
        bamcell="02.mapping/{groupname}.rmdup.merge.pair.singlcell.sort.bam"
    output:
        selectedbam="07.singlecell/{groupname}.FRiP.selected.singlecell.bam",
        selectedfrags=temp("07.singlecell/{groupname}.FRiP.selected.singlecell.fragments"),
        selectedfrags2="07.singlecell/{groupname}.FRiP.selected.singlecell.fragments.gz",
        
        selectedbam_t="07.singlecell/{groupname}.FRiT.selected.singlecell.bam",
        selectedfrags_t=temp("07.singlecell/{groupname}.FRiT.selected.singlecell.fragments"),
        selectedfrags2_t="07.singlecell/{groupname}.FRiT.selected.singlecell.fragments.gz",
    params:
        extractcell=config['scripts']['extractcell'],
        python_main=config['scripts']['python'],
        Quality=config['quality']
    shell:
        "{params.python_main} {params.extractcell} -n {input.cellnames}  -f {input.fragments} "
        "-b {input.bamcell} -ob {output.selectedbam} -of {output.selectedfrags} && bgzip -c {output.selectedfrags} > {output.selectedfrags2}"
        " && sambamba index -t 5 {output.selectedbam} && "
        " {params.python_main} {params.extractcell} -n {input.cellnames2}  -f {input.fragments} "
        "-b {input.bamcell} -ob {output.selectedbam_t} -of {output.selectedfrags_t} && bgzip -c {output.selectedfrags_t} > {output.selectedfrags2_t}"
        " && sambamba index -t 5 {output.selectedbam_t}"

rule Cell_Peak_Matrix:
    message: "generate the cell*peak matrix"
    input:
        peaks="05.Peakcalling/{groupname}_peaks.filterBL.bed",
        selectcell="06.QC/{groupname}_FRiP_selected.cell.names",
        bamcell="07.singlecell/{groupname}.FRiP.selected.singlecell.bam",
    output:
        matrix="07.singlecell/{groupname}.cell_peak.matrix",
    params:
        extractcell=config['scripts']['cellpeak'],
        python_main=config['scripts']['python'],

    shell:
         "{params.python_main} {params.extractcell} -p {input.peaks} "
         " -n {input.selectcell} -ib {input.bamcell} -o {output.matrix}"


rule Cell_TSS_Matrix:
    message: "generate the cell*peak matrix"
    input:
        peakBL=config['tssbed'],
        selectcell="06.QC/{groupname}_FRiT_selected.cell.names",
        bamcell="07.singlecell/{groupname}.FRiT.selected.singlecell.bam",
    output:
        matrix="07.singlecell/{groupname}.cell_tss.matrix",
    params:
        extractcell=config['scripts']['cellpeak'],
        python_main=config['scripts']['python'],

    shell:
         "{params.python_main} {params.extractcell} -p {input.peaks} "
         " -n {input.selectcell} -ib {input.bamcell} -o {output.matrix}"