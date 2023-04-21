Args <- commandArgs(trailingOnly=TRUE)

pfile = as.character(Args[1])
phe = as.character(Args[2])
pheno_names_file = as.character(Args[3])

pheno_names = paste0(read.table(pheno_names_file)$V1,collapse = ',')

if(!dir.exists('./1.plink_temp')) dir.create('./1.plink_temp')
setwd('./1.plink_temp')

task_name='1.plink_pre'
qsub_file=paste0(task_name,'_qsub.txt')
jobname=task_name  
logname=paste(task_name,'.log',sep='') 
errfile=paste(task_name,'.err',sep='') 

sink(qsub_file)

cat('#!/bin/bash\n')
cat('\n')
cat(sprintf("#PBS -N %s\n",jobname))
cat('\n')
cat(sprintf("#PBS -o %s\n",logname))
cat('\n')
cat(sprintf("#PBS -e %s\n",errfile))
cat('\n')
cat(sprintf("#PBS -q %s\n",'CGB'))
cat('\n')
cat('#PBS -m bea')
cat('\n')
cat('#PBS -l nodes=1:ppn=1')
cat('\n')

cat(paste('cd',getwd(),'\n'))
out=paste0(task_name)
cat('source ~/.bashrc\n')
cat('\n')


cat(sprintf("plink2 --glm allow-no-covars skip-invalid-pheno \\
            --read-freq ../0.freq.afreq \\
            --out %s \\
            --pfile %s \\
            --pheno %s \\
            --pheno-name %s", out,pfile,phe,pheno_names))
sink()
# com=paste('qsub',qsub_file)
# system(com)
