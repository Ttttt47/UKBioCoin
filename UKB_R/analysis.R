library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(circlize)

cov_yy = as.matrix(read.table('./14M_cov_yy.table'))


# 绘制热图
ggplot(cov_yy_df, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank())
  
heatmap_object = Heatmap(cov_yy, col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
        , cluster_rows = TRUE
        , cluster_columns = TRUE,
        show_column_names =  F,
        show_row_names = F)
heatmap_object = draw(heatmap_object)
row_order(heatmap_object)
col_order(heatmap_object)

# 获取聚类结果
row_order <- heatmap_object@row_order
col_order <- heatmap_object@column_order

# 根据聚类结果提取关系较强的区块
strong_relation_rows <- cov_yy[row_order, ]
strong_relation_cols <- strong_relation_rows[, col_order]

# 显示提取的区块
print(strong_relation_cols)


##### clustering
# 对协方差矩阵进行层次聚类
row_dist <- dist(cov_yy)
col_dist <- dist(t(cov_yy))
row_hclust <- hclust(row_dist)
col_hclust <- hclust(col_dist)


# 选择一个高度值
height_threshold <- 2

# 在树状图中，选择具有较强关系的行和列
strong_relation_row_clusters <- cutree(row_hclust, h = height_threshold)

# 绘制行和列的树状图
par(mfrow=c(1,2))
plot(row_hclust, main = "Row Dendrogram", sub = "", xlab = "", labels = F)
summary(factor(strong_relation_row_clusters))
hist(strong_relation_row_clusters,breaks = 200)

#  related
phe = 'X102.0.0'
id = strong_relation_row_clusters[phe]
id = 3
summary(factor(strong_relation_row_clusters))[id]

# 找到具有较强关系的行和列的索引
strong_relation_rows <- which(strong_relation_row_clusters == id)
strong_relation_rows
# 提取具有较强关系的行和列
cov_yy[strong_relation_rows, phe]
strong_relation_block <- cov_yy[strong_relation_rows, strong_relation_rows]
# 显示提取的区块
print(strong_relation_block)

