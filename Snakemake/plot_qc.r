#################################################
#  File Name:get_qc.r
#  Author: xingpengwei
#  Mail: xingwei421@qq.com,pengwei.xing@igp.uu.se,xpw1992@gmail.com
#  Created Time: Sat Jan  8 18:13:16 2022
#################################################

library(ggplot2)
library(scales)
library(getopt)

spec <- matrix(
  c(
    "FRiP", "p", 2, "integer", "FRiP value",
    "FRiT", "t", 2, "integer", "FRiT value",
    "input", "i", 2, "character", "Input: single cell table",
    "number", "n", 2, "character", "the unique fragments number",
    "output", "o", 2, "character", "The output name!",
    "help",   "h", 0, "logical",  "Rscript plot_qc.r -p 5 -t 10 -i Spleen.singlecell.qc.txt -o spleen"),
  byrow = TRUE, ncol = 5 
)

# 使用getopt方法
opt <- getopt(spec = spec)

# opt实际上就是一个列表，直接使用$来索引到对应的参数的值
number =  as.numeric(opt$number)
FRiP_value <- as.numeric(opt$FRiP)
FRiT_value <- as.numeric(opt$FRiT)
data <- read.table(opt$input,header=T)

data$ifcell <- "No"
data$ifcell[which(data$Final >= number & data$Final <= 100000)] <- "Yes"
data_cell <- data[which(data$Final >= number & data$Final <= 100000),]
total_raw <- paste0("The total raw fragments:", prettyNum(sum(data$total),big.mark = ","))
total_mapped_value = sum(data$total)

selected_cells <- paste0(opt$output,"_Frags",number,"_selected.cell.names")
write.table(x = data_cell$barcode ,file = selected_cells, quote = F,row.names = F,col.names = F)

total_mapped <- paste0("The total mapped fragments:", prettyNum(total_mapped_value,big.mark = ","))
total_passed <- paste0("The total final fragments:", prettyNum(sum(data$Final),big.mark = ","))
total_cells <- length(which(data$Final>=number))
total_cells2 <- paste0("The total single cells:",total_cells)
median_value <- median(data_cell$Final)
median_frag <- paste0("The median fragments of cell:",median_value)
data2 <- data[data$Final>0,]
####Bar plot 
p <- ggplot(data2, aes(x=Final,fill=ifcell)) + geom_histogram(color='white',binwidth = 0.05)
p <- p +   scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)),limits = c(50,100000))
p <- p + theme_bw()+xlab("Number of fragments")+ylab("Frequency")
p <- p + geom_vline(xintercept = c(number,100000),linetype = 'longdash',colour="#03406A")
p <- p + scale_fill_manual(values=c("#65A5D1","#03406A"))
p <- p + annotate("text", x = number, y = max(na.omit(layer_data(p)$count))*0.7,
   label = paste0(total_raw,"\n",total_mapped,"\n",total_passed,"\n",total_cells2,"\n",median_frag), hjust=0)
filename1 <- paste0(opt$output,"_fragments_barplot.pdf")
pdf(file=filename1,width=8,height=6)
p
dev.off()

dup_rate <- data_cell$Dups
pass_filter_frag <- data_cell$Final
firp_value <- data_cell$FRiP
data2 <- as.data.frame(cbind(dup_rate,pass_filter_frag,firp_value))
data2$firp_value = data2$firp_value

p_dup_rate <- ggplot(data2, aes(x='Duplication rate',y=dup_rate)) +
  geom_violin(trim=FALSE,fill='#51cfba', color="black")+
  geom_boxplot(width=0.05,fill="#FFAB73")+
  theme_classic()+
labs(y = "Percentage (%)",x="")+
theme(axis.text.x = element_text(color="black",size=15),axis.text.y = element_text(color="black"))+
scale_fill_brewer(palette="Blues")

p_frags <- ggplot(data2, aes(x='Fragments per cell',y=pass_filter_frag)) + 
  geom_violin(trim=FALSE,fill='#5ED0BD', color="black")+
  geom_boxplot(width=0.05,fill="#FFAB73")+
#geom_jitter(shape=16, size=1,position=position_jitter(0.2))+
  theme_classic()+
labs(y = "The fragments number (%)",x="")+
scale_fill_brewer(palette="Blues")+
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)))+
theme(axis.text.x = element_text(color="black",size=15),axis.text.y = element_text(color="black",size=15))

p_frip <- ggplot(data2, aes(x='FRiP per cell',y=firp_value)) +
  geom_violin(trim=FALSE,fill='#5ED0BD', color="black")+
  geom_boxplot(width=0.05,fill="#FFAB73")+
#geom_jitter(shape=16, size=1,position=position_jitter(0.2))+
  theme_classic()+
labs(y = "Percentage (%)",x="")+
scale_fill_brewer(palette="Blues")+
#scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
 #             labels = trans_format("log10", math_format(10^.x)))+
theme(axis.text.x = element_text(color="black",size=15),axis.text.y = element_text(color="black",size=15))

library(ggpubr)
filename2 = paste0(opt$output,"_violinplot.pdf")
pdf(file=filename2,width=8,height=6)
ggarrange(p_dup_rate, p_frags, p_frip,
          labels = c("A", "B","C"),
          ncol = 3, nrow = 1)

dev.off()

####scatter plot : unique fragments vs FRiP
data2 <- data[data$Final>200&data$Final<100000,]
data2$Group <- "unselected"
data2$Group[which((data2$Final > number) & (data2$FRiP > FRiP_value))] = 'selected'

data2$Group <- factor(data2$Group,levels = c('selected','unselected'))
num_cell2 <- nrow(data[data$Final>number&data$Final<10000&data$FRiP>FRiP_value,])
p <- ggplot(data=data2)+geom_point(aes(x=Final,y=FRiP, color = Group),size=1)+scale_color_manual(values=c("#2f5688","#CC0000"))+theme_bw()
p <- p +geom_vline(xintercept = number,col="red")+geom_hline(yintercept = FRiP_value,col="red")
filename3 <- paste0(opt$output,"_Unique",number,"_FRiP",FRiP_value,".scatterplot.pdf")
pdf(file=filename3,width=6,height=6)
p + annotate("text", x = 5000, y = 30,label = paste0("Selected cells = ",num_cell2,"\n" ,"FRiP=5","\n","Unique Frags = 500"), hjust=0)
dev.off()

filename4 <- paste0(opt$output,"_Unique",number,"_FRiP",FRiP_value,"_selected.cell.names")
data_out <- data[data$Final>number&data$Final<100000&data$FRiP>FRiP_value,]
write.table(x = data_out$barcode ,file = filename4, quote = F,row.names = F,col.names = F)
filename6 <- paste0(opt$output,"_Unique",number,"_FRiP",FRiP_value,"_selected.cell.csv")
write.csv(x = data_out ,file = filename6, quote = F,row.names = F)

######scatter plot : unique fragments vs FRiT
data2 <- data[data$Final>200&data$Final<100000,]
data2$Group <- "unselected"
data2$Group[which((data2$Final > number) & (data2$FRiT > FRiT_value))] = 'selected'

data2$Group <- factor(data2$Group,levels = c('selected','unselected'))
num_cell2 <- nrow(data[data$Final>number&data$Final<100000&data$FRiT>FRiT_value,])
p <- ggplot(data=data2)+geom_point(aes(x=Final,y=FRiT, color = Group),size=1)+scale_color_manual(values=c("#2f5688","#CC0000"))+theme_bw()
p <- p +geom_vline(xintercept = number,col="red")+geom_hline(yintercept = FRiT_value,col="red")
filename3 <- paste0(opt$output,"_Unique",number,"_FRiT",FRiT_value,".scatterplot.pdf")
pdf(file=filename3,width=6,height=6)
p + annotate("text", x = 5000, y = 30,label = paste0("Selected cells = ",num_cell2,"\n" ,"FRiT=",FRiT_value,"\n","Unique Frags = 500"), hjust=0)
dev.off()

filename5 <- paste0(opt$output,"_Unique",number,"_FRiT",FRiT_value,"_selected.cell.names")
data_out <- data[data$Final>number&data$Final<100000&data$FRiT>FRiT_value,]
write.table(x = data_out$barcode ,file = filename5, quote = F,row.names = F,col.names = F)

filename6 <- paste0(opt$output,"_Unique",number,"_FRiT",FRiT_value,"_selected.cell.csv")
write.csv(x = data_out ,file = filename6, quote = F,row.names = F)