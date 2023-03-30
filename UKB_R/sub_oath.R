
oath_dir = ''
data_file = '/public3/project_users/chengb/hjc/projects/UKBioCoin/matrix/14M'
phe_names = paste0('X',read.table("/public3/project_users/chengb/hjc/projects/UKBioCoin/pheno_set_0320_3_162")$V1,'.0.0')
PC_nums = c(5,10)
out_dir = '/public3/project_users/chengb/hjc/projects/UKBioCoin/oath_results/'

for(phe_name in phe_names)
{
  for(PC_num in PC_nums){
    PCs = paste0('PC',1:PC_num,collapse = ',')
    covariates = paste0(PCs)
    
    task_name=paste('oath-',phe_name,'-',PC_num,'PC',sep='')
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
    cat(sprintf("#PBS -q %s\n",'CGB'))
    cat('\n')
    cat('#PBS -m bea')
    cat('\n')
    cat('#PBS -l nodes=1:ppn=1')
    cat('\n')
    
    dir<-getwd()
    cat(paste('cd',dir,'\n'))
    out=paste0(out_dir,task_name)
    cat('source ~/.bashrc\n')
    cat('\n')
    
    cat(sprintf("%s --file %s --phe %s --covar %s --out %s \n",oath_dir, data_file,phe_name,covariates,out))
    
    sink()
    com=paste('qsub',qsub_file)
    system(com)
    # wait the task to initialize
  }
}
