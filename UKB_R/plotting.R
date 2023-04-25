library(ComplexHeatmap)
p=3
q=3

library(circlize)
col_fun = colorRamp2(c(-1, 0, 1), c('#EE4431', "white", '#FC8002'))

nsnp = 30
nc = 15
nv = 5

pdf('./figures/sample_cov_xy.pdf',width = 4,height = 5)
sam_cov_xy = matrix(c(runif(nsnp*nc,0.3,0.6),runif(nsnp*nv,0.3,1)),nsnp,(nv+nc))
Heatmap(sam_cov_xy, cluster_rows = F, cluster_columns = F, col = c('white','#FC8002'),
        column_split = factor(c(rep('Covariates(Non-PC)',nc),rep('PC',nv)), levels = c('Covariates(Non-PC)','PC')),
        column_title_gp = gpar(fontsize = 9),
        row_title_gp = gpar(fontsize = 9),
        row_title = 'SNP',
        cluster_row_slices = FALSE, 
        cluster_column_slices = FALSE,
        rect_gp = gpar(col = "white", lwd = 2),
        show_row_names = F,
        show_column_names = F,
        width = unit(6, "cm"), height = unit(9, "cm"),
        show_heatmap_legend = FALSE,border = F)
dev.off()

pdf('./figures/sample_cov_yy.pdf',width = 4,height = 4)
sam_cov_yy = cov_yy[c(1:nc,140:(139+nv)),c(1:nc,140:(139+nv))]
sam_cov_yy = sam_cov_yy/max(sample_cov_yy)
sam_cov_yy[(nc+1):(nc+nv),(nc+1):(nc+nv)] = diag(nv)
Heatmap(sam_cov_yy, cluster_rows = F, cluster_columns = F, col = c('white','#369F2D'),
        column_split = factor(c(rep('Covariates(Non-PC)',nc),rep('PC',nv)), levels = c('Covariates(Non-PC)','PC')),
        row_split = factor(c(rep('Covariates(Non-PC)',nc),rep('PC',nv)), levels = c('Covariates(Non-PC)','PC')),
        column_title_gp = gpar(fontsize = 9),
        row_title_gp = gpar(fontsize = 9),
        cluster_row_slices = FALSE, 
        cluster_column_slices = FALSE,
        rect_gp = gpar(col = "white", lwd = 2),
        show_row_names = F,
        show_column_names = F,
        width = unit(6, "cm"), height = unit(6, "cm"),
        show_heatmap_legend = FALSE,border = F)
dev.off()

pdf('./figures/sample_var_x.pdf',width = 2,height = 4)
sam_var_x = runif(nsnp)
sam_var_x = sam_var_x/max(sam_var_x)
Heatmap(sam_var_x, cluster_rows = F, cluster_columns = F, col = c('#FAC7B3','#EE4431'),
        column_title_gp = gpar(fontsize = 9),
        row_title_gp = gpar(fontsize = 9),
        row_title = 'SNP',
        cluster_row_slices = FALSE, 
        cluster_column_slices = FALSE,
        rect_gp = gpar(col = "white", lwd = 2),
        show_row_names = F,
        show_column_names = F,
        width = unit(9/nsnp, "cm"), height = unit(9, "cm"),
        show_heatmap_legend = FALSE,border = F)
dev.off()

pdf('./figures/sample_a.pdf',width = 2,height = 4)
Heatmap(sam_cov_xy[,(p+1)], cluster_rows = F, cluster_columns = F, col = c('white','#FC8002'),
        cluster_row_slices = FALSE, 
        cluster_column_slices = FALSE,
        rect_gp = gpar(col = "white", lwd = 2),
        show_row_names = F,
        show_column_names = F,
        width = unit(6, "cm"), height = unit(9, "cm"),
        show_heatmap_legend = FALSE,border = F)
dev.off()


p=3
q=3
Omega = matrix(0,1+p+q,1+p+q)
Omega[1,1] = sam_var_x[1]
Omega[1,2:(1+p+q)] = sam_cov_xy[1,c(1:p,(nc+nv-2):(nc+nv))]
Omega[2:(1+p+q),1] = sam_cov_xy[1,c(1:q,(nc+nv-2):(nc+nv))]
Omega[2:(1+p+q),2:(1+p+q)] = cov_yy[c(1:q,(nc+nv-2):(nc+nv)),c(1:q,(nc+nv-2):(nc+nv))]


library(circlize)
col_fun_var_x = colorRamp2(c(min(sam_var_x),max(sam_var_x)), c('#FAC7B3','#EE4431'))
col_fun_cov_xy = colorRamp2(c(min(sam_cov_xy),max(sam_cov_xy)), c('white','#FC8002'))
col_fun_cov_yy = colorRamp2(c(min(sam_cov_yy),max(sam_cov_yy)), c('white','#369F2D'))

cnames = c('SNP i','Sex','BMI','Age','PC1','PC2','PC3')
pdf('./figures/sample_Omega.pdf',width = 4,height = 4)
Heatmap(Omega, cluster_rows = F, cluster_columns = F,
        column_title_gp = gpar(fontsize = 9),
        row_title_gp = gpar(fontsize = 9),
        rect_gp = gpar(type="none"),
        row_split = factor(c("SNP",rep('Covariate',3),rep('PC',3)), levels = c('SNP','Covariate','PC')),
        column_split = factor(c("SNP",rep('Covariate',3),rep('PC',3)), levels = c('SNP','Covariate','PC')),
        cluster_row_slices = FALSE, 
        cluster_column_slices = FALSE,
        column_names_side = 'top',
        row_names_side = 'right',
        column_names_rot = 45,
        row_title = NULL,
        column_title = NULL,
        show_heatmap_legend = FALSE,border = F,row_labels = cnames, column_labels = cnames,
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(i==1&j==1){
            grid.rect(x = x, y = y, width = width, height = height, 
                      gp = gpar(col = "white", fill = col_fun_var_x(Omega[i,j])))
          }
          else if((i==1&j!=1)|(i!=1&j==1)){
            grid.rect(x = x, y = y, width = width, height = height, 
                      gp = gpar(col = "white", fill = col_fun_cov_xy(Omega[i,j])))
          }else{
            grid.rect(x = x, y = y, width = width, height = height, 
                      gp = gpar(col = "white", fill = col_fun_cov_yy(Omega[i,j])))
          }
        },
        width = unit(18/nsnp*(p+q+1), "cm"), height = unit(18/nsnp*(p+q+1), "cm"),
)
dev.off()

