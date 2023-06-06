for (i in 1:4) {
  task_name=paste('clump-',i,sep='')
  qsub_file=paste0(task_name,'.txt')
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
  cat('#PBS -m bea')
  cat('\n')
  cat('#PBS -l nodes=1:ppn=2')
  cat('\n')
  
  dir<-getwd()
  cat(paste('cd',dir,'\n'))
  out=task_name
  cat('source ~/.bashrc\n')
  cat('\n')
  # cat('conda init \n')
  cat('conda activate \n')
  # cat('conda activate r4.1\n')
  # cat('\n')
  # cat('export OPENBLAS_NUM_THREADS=1 \n')
  
  cat(sprintf("R CMD BATCH --vanilla '--args %d' clump_pbs.R %s.Rout \n",i))
  sink()
  com=paste('qsub',qsub_file)
  system(com)
  # wait the task to initialize
}