args<- commandArgs(trailingOnly=TRUE)
print(args)

top = as.numeric(args[1])
sdLim = as.numeric(args[2])
result_path = args[3]

pca = read.delim(paste0(result_path, "/laser.RefPC.coord"))

for(i in 1:top)
        {
        m = which(pca[,2+i] > sdLim*sd(pca[,2+i]) | pca[,2+i] < -sdLim*sd(pca[,2+i]))
        if(length(m) > 0)
                {
                if(!exists("samples")){
                        samples = cbind(as.character(pca[m,2]),i)
                        }
                else{
                        samples = rbind(samples,cbind(as.character(pca[m,2]),i))
                        }
                }
        }

if(exists("samples"))
        {
        write.table(samples, paste0(result_path, "/remove.txt"), quote=F, sep="\t", row.names=F, col.names=F)
        }


