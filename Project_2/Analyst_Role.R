setwd("C:/Users/vishw/Desktop/Vishwa/BF528/PROJECT 2")       #to set the directory path

diff <- read.table('gene_exp.diff.csv', header = TRUE)           #to read the table from working directory
order <- diff[order(diff$q_value) , ]                            #to sort the data 
View(order)

top10 <- head(order[ , c(3,8,9,10,12,13)], 10)                   #to find the top 10 from sorted data
View(top10)

histogram_1 <- hist(diff$log2.fold_change., breaks = 75, xlab = 'Log2 fold change', col = 'orange', main = 'Histogram with Log2 fold change of all genes')
#hist function to plot histogram of log2 fold change 

significant_diff <- subset(x = diff, significant == 'yes')   #making subset with significant genes 
View(significant_diff)

histogram_2 <- hist(significant_diff$log2.fold_change., breaks = 65, xlab = 'Log2 fold change', col= 'orange', main = 'Histogram with log2 fold change of significant genes')

p_value <- subset(diff, diff$p_value<0.01)   #subset the diff to find genes with p_value<0.01 i.e. 2374
View(p_value)

#subset of significant_diff with p_value < 0.01
sig_pvalue <- subset(significant_diff, significant_diff$p_value<0.01)  #2166 

#subset of diff with q_value<0.01
q_value <- subset(diff, diff$q_value<0.01)   #1144
View(q_value)

#subset of p value with log2 fold change greater than 0
up_pvalue <- subset(p_value, log2.fold_change. > 0)   #1189
View(up_pvalue)

#subset of p value with log2 fold change less than 0
down_pvalue <- subset(p_value, log2.fold_change. <0)  #1185

#subset of significant genes based for up and down regulated based on log2 fold change
up_sig <- subset(significant_diff, significant_diff$log2.fold_change.>0)     #1091
View(up_sig)
down_sig <- subset(significant_diff, significant_diff$log2.fold_change. < 0) #1075
View(down_sig)

write (up_sig$gene, 'up_regulatedgenes.csv', sep = ',')              #writing csv files for significant up and down genes
write(down_sig$gene, 'down_regulatedgenes.csv', sep = ',')
write(top10, 'Top10_differentiallyEG.csv')
