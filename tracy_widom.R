library(LEA)
args <- commandArgs(T)
sample_num = args[1]
file_dir = args[2]
filename = args[3]

EUR = new("pcaProject")
EUR@n = as.integer(sample_num)

EUR@projDir = file_dir
EUR@eigenvalue.file = filename
EUR_tw = tracy.widom(EUR)
# number of significant PCs: 
cat(pc_num = nrow(subset(EUR_tw, pvalues < 0.001)), "\n")
