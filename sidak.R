library(tidyverse)
library(data.table)
args <- commandArgs(T)

work_dir = args[1]
prefix = args[2]
pc_num = args[3]

# calculate the lambdaGC
calculate_GC <- function(work_dir, prefix, pc_num) {
  output = vector("list", pc_num)
  for(i in 1:pc_num) {
    filename = paste0(work_dir, prefix, i, ".assoc.linear")
    df = as.data.frame(fread(filename)) %>% mutate(chisq = STAT^2)
    output[[i]] = data.frame(PC_index = i, lambdaGC = median(df$chisq)/qchisq(0.5,1))
  }
  write_tsv(bind_rows(output), paste0(work_dir, prefix, "all.lambdaGC.txt"))
}

calculate_GC(args[1], args[2], args[3])

# combine all GWAS results
combine_results = matrix(nrow=0, ncol=12)
for (i in 1:pc_num) {
  filename = paste0(work_dir, prefix, i, ".assoc.linear") 
  df = as.data.frame(fread(filename))
  lambda = read.table(paste0(work_dir, prefix, "all.lambdaGC.txt"), header=T)
  df = df %>% mutate(chisq = STAT^2, 
  		     chisq_adjusted = chisq/lambda[[i, 2]], 
    	             P_adjusted = exp(1)^pchisq(chisq_adjusted, 1, lower.tail = FALSE, log.p = TRUE))
  write_tsv(df, filename)
  combine_results = rbind(combine_results, df)
}


# select most significant P values for each SNP
num_snps=nrow(df)
combine_results = combine_results %>% arrange(CHR, BP, P_adjusted)
new_combine_results = combine_results[seq(1, as.integer(pc_num)*as.integer(num_snps), as.integer(pc_num)), ]
new_df = new_combine_results %>% mutate(P_concat = 1-(1-P_adjusted)^as.integer(pc_num))
new_df = new_df %>% mutate(P_concat = ifelse(P_concat == 0, P_adjusted*as.integer(pc_num), P_concat))
write_tsv(new_df, paste0(work_dir, prefix, "all.assoc.linear"))
