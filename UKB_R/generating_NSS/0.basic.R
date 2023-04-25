Args <- commandArgs(trailingOnly=TRUE)

pfile = as.character(Args[1])
# phe = as.character(Args[2])
# pheno_names_file = as.character(Args[3])


if(!dir.exists('./0.basic')) dir.create('./0.basic')
setwd('./0.basic')

task_name='0.freq'
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
cat('#PBS -l nodes=1:ppn=2')
cat('\n')

cat(paste('cd',getwd(),'\n'))
out=paste0(task_name)
cat('source ~/.bashrc\n')
cat('\n')


cat(sprintf("plink2 \\
      --pfile  %s \\
      --freq --out %s \\", pfile,out))
sink()
# com=paste('qsub',qsub_file)
# system(com)


task_name='0.count'
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
cat('#PBS -l nodes=1:ppn=2')
cat('\n')

cat(paste('cd',getwd(),'\n'))
out=paste0(task_name)
cat('source ~/.bashrc\n')
cat('\n')


cat(sprintf("plink2 \\
      --pfile  %s \\
      --geno-count --out %s \\", pfile, out))
sink()