task_prefix = 'oath'
data_file = '/public3/project_users/chengb/hjc/projects/UKBioCoin/matrix/14M'
phe_names = c('X2443.0.0')
PC_nums = c(0,5,10)
covars = c('X31.0.0')
out_dir = '/public3/project_users/chengb/hjc/projects/UKBioCoin/oath_results/'

for(phe_name in phe_names)
{
  for(PC_num in PC_nums){
    if (PC_num==0){
      covariates = paste0(covars,collapse = ',')
    }else{
      PCs = paste0('PC',1:PC_num,collapse = ',')
      covariates = paste0(PCs,',',covars)
    }
    
    
    
    task_name=paste(task_prefix,'-',phe_name,'-',PC_num,'PC',sep='')
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
    
    cat(sprintf("/public3/project_users/chengb/hjc/projects/UKBioCoin/Cpp/UKBioCoin --file %s --phe %s --covar %s --out %s \n",data_file,phe_name,covariates,out))
    
    sink()
    com=paste('qsub',qsub_file)
    system(com)
  }
}
  
