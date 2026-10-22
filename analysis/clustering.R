library(cluster)
for(name in c('all','angust')) {
 x=as.matrix(read.table(paste0('analysis/data/',name,'_final.lfmm_imputed.lfmm')))
 p=colMeans(x)/2;keep=p>0 & p<1 & apply(x,2,sd)>0;x=x[,keep];p=p[keep]
 y=sweep(sweep(x,2,2*p),2,sqrt(2*p*(1-p)),'/');d=dist(y);h=hclust(d,method='ward.D2');meta=read.csv(paste0('analysis/results/',name,'_PCA.csv'))
 write.csv(data.frame(left=h$merge[,1],right=h$merge[,2],height=h$height),paste0('analysis/results/',name,'_ward_merge.csv'),row.names=FALSE)
 write.csv(data.frame(index=h$order,sample=meta$sample[h$order]),paste0('analysis/results/',name,'_ward_order.csv'),row.names=FALSE)
 ss=do.call(rbind,lapply(2:6,function(k)data.frame(K=k,silhouette=mean(silhouette(cutree(h,k),d)[,3]))))
 write.csv(ss,paste0('analysis/results/',name,'_ward_silhouette.csv'),row.names=FALSE)
 k=ss$K[which.max(ss$silhouette)];meta$ward_cluster=cutree(h,k);write.csv(meta,paste0('analysis/results/',name,'_ward_clusters.csv'),row.names=FALSE)
 print(name);print(ss)
}
