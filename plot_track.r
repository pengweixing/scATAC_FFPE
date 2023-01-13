#################################################
#  File Name:plot_track.r
#  Author: Pengwei.Xing
#  Mail: xingwei421@qq.com,pengwei.xing@igp.uu.se,xpw1992@gmail.com
#  Created Time: Fri Sep 30 15:50:09 2022
#################################################
#  chr	start	end	cell_ID_in each_cell_type	cell_type
"""
chr1	3000817	3000938	1958	B_cells
chr1	3000828	3000938	1958	B_cells
chr1	3000890	3000938	1958	B_cells
chr1	3000892	3000938	1958	B_cells
chr1	3000893	3000938	1958	B_cells
chr1	3001258	3001350	1444	B_cells
chr1	3004174	3004322	331	T_cells
chr1	3004303	3004366	155	T_cells
chr1	3004487	3004541	746	T_cells
"""
library(Gviz)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
data = read.table("fragments_cell_type.gviz.txt")
data_long = read.table("../04.archR_long/fragments_cell_type.gviz.txt")

aaa <- function(data2){
    size=0.3
    dat = data2[data2$Type=="B_cells",]$ID
    end = data2[data2$Type=="B_cells",]$end
    start = data2[data2$Type=="B_cells",]$start

    dat2 = data2[data2$Type=="T_cells",]$ID
    end2 = data2[data2$Type=="T_cells",]$end
    start2 = data2[data2$Type=="T_cells",]$start

    dat3 = data2[data2$Type=="other",]$ID
    end3 = data2[data2$Type=="other",]$end
    start3 = data2[data2$Type=="other",]$start

    dtrack <- DataTrack(data = dat, start = start,col="black",
                    end = end, chromosome = mychr, genome = 'mm10',
                    name = "B_cells",cex=size)
    dtrack2 <- DataTrack(data = dat2, start = start2,col="black",
                    end = end2, chromosome = mychr, genome = 'mm10',
                    name = "T_cells",cex=size)
    dtrack3 <- DataTrack(data = dat3, start = start3,col="black",
                    end = end3, chromosome = mychr, genome = 'mm10',
                    name = "other",cex=size)
    return(list(gtrack, grtrack,dtrack, dtrack2, dtrack3))
}
gtrack <- GenomeAxisTrack(col="darkgray",cex=0.8)
colnames(data)=c("chr","start","end","ID","Type")
colnames(data_long)=c("chr","start","end","ID","Type")
### B cell
mychr = 'chr7'
left= 126402155
right = 126427365

data2 = data[data$chr==mychr&data$start>left&data$end<right,]
data_long2 = data_long[data_long$chr==mychr&data_long$start>left&data_long$end<right,]
grtrack <- GeneRegionTrack(txdb, genome = 'mm10',stacking="dense",fill="#006064",
                           col="#006064",
                           transcriptAnnotation="symbol",
                           chromosome = mychr, name = "Gene")
short = aaa(data2)
long = aaa(data_long2)
pdf(file="Bcell.gviz.track.pdf",height=5,width=10)
plotTracks(short,chromosome=mychr,from=left, to=right)
dev.off()

pdf(file="../04.archR_long/Bcell.gviz.track.pdf",height=5,width=10)
plotTracks(long,chromosome=mychr,from=left, to=right)
dev.off()
