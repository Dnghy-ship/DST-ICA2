# 读取metadata
metadata <-read.delim("ALS_GSE234297_pb_series_matrix.txt")
metadata <- as.data.frame(t(metadata))
write.csv(metadata,"ALS_GSE234297_pb_series_matrix.csv")
