#################################################
#  File Name:get_qc.r
#  Author: xingpengwei
#  Mail: xingwei421@qq.com,pengwei.xing@igp.uu.se,xpw1992@gmail.com
#  Created Time: Sat Jan  8 18:13:16 2022
#################################################

library(ggplot2)
library(scales)
args <- commandArgs(T)
data <- read.csv(args[1])
number <- args[3]
number = as.numeric(number)
data$ifcell="No"
data$ifcell[which(data$passed_filters>=number&data$passed_filters<=100000)] = "Yes"
data_cell = data[which(data$passed_filters>=number&data$passed_filters<=100000),]
total_raw = paste0("The total raw fragments:", prettyNum(sum(data$total),big.mark = ","))
total_mapped_value = sum(data$total)-sum(data$unmapped)
total_mapped = paste0("The total mapped fragments:", prettyNum(total_mapped_value,big.mark = ","))
total_passed = paste0("The total final fragments:", prettyNum(sum(data$passed_filters),big.mark = ","))
total_cells = length(which(data$passed_filters>=number))
total_cells2 = paste0("The total single cells:",total_cells)
median_value = median(data_cell$passed_filters)
median_frag = paste0("The median fragments of cell:",median_value)
data2 = data[data$passed_filters>0,]

p <- ggplot(data2, aes(x=passed_filters,fill=ifcell)) + geom_histogram(color='white',binwidth = 0.05)
p <- p +   scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)),limits = c(50,100000))
p <- p + theme_bw()+xlab("Number of fragments")+ylab("Frequency")
p <- p + geom_vline(xintercept = c(number,100000),linetype = 'longdash',colour="#03406A")
p <- p + scale_fill_manual(values=c("#65A5D1","#03406A"))
p <- p + annotate("text", x = number, y = max(na.omit(layer_data(p)$count))*0.7,
   label = paste0(total_raw,"\n",total_mapped,"\n",total_passed,"\n",total_cells2,"\n",median_frag), hjust=0)
filename1 <- paste0(args[2],"_fragments_barplot.pdf")
pdf(file=filename1,width=8,height=6)
p
dev.off()
dup_rate = data_cell$duplicate/(data_cell$total-data_cell$unmapped)
pass_filter_frag = data_cell$passed_filters
firp_value = data_cell$peak_region_fragments/data_cell$passed_filters
data2 = as.data.frame(cbind(dup_rate,pass_filter_frag,firp_value))
data2$dup_rate = data2$dup_rate*100
data2$firp_value = data2$firp_value*100

p_dup_rate <- ggplot(data2, aes(x='Duplication rate',y=dup_rate)) +
  geom_violin(trim=FALSE,fill='#5ED0BD', color="black")+
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
filename2 = paste0(args[2],"_violinplot.pdf")
pdf(file=filename2,width=8,height=6)
ggarrange(p_dup_rate, p_frags, p_frip,
          labels = c("A", "B","C"),
          ncol = 3, nrow = 1)

dev.off()
