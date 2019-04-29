args <- commandArgs(trailingOnly = TRUE)

if (!require(data.table)) {
  install.packages("data.table")
  library(data.table)
}

if (!require(Matrix)) {
  install.packages("Matrix")
  library(Matrix)
}

file<-args[1]
out<-args[2]

f<-as.data.frame(fread(file,h=T,sep='\t'))

sm = sparseMatrix(as.numeric(factor(f$gene)), as.numeric(factor(f$cell)), x=f$count)

rownames(sm)<-levels(factor(f$gene))
colnames(sm)<-levels(factor(f$cell))

write.table(as.matrix(sm),file=paste0(out,'.txt'),sep='\t',col.names=T,row.names=T,quote=F)

saveRDS(sm,file=paste0(out,'.RDS'))
