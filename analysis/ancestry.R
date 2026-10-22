.libPaths(c('analysis/Rlib',.libPaths()))
library(LEA)
root=normalizePath('analysis');setwd(file.path(root,'data'))
for (name in c('all','angust')) {
 input=paste0(name,'.lfmm');maxK=ifelse(name=='all',6,5)
 p=snmf(input,K=1:maxK,repetitions=5,entropy=TRUE,percentage=.05,CPU=2,alpha=10,iterations=1000,tolerance=1e-5,seed=20261009,project='new')
 cv=do.call(rbind,lapply(1:maxK,function(k)data.frame(K=k,run=1:5,cross_entropy=as.numeric(cross.entropy(p,K=k)))))
 write.csv(cv,file.path(root,'results',paste0(name,'_cross_entropy.csv')),row.names=FALSE)
 agg=aggregate(cross_entropy~K,cv,mean);bestK=agg$K[which.min(agg$cross_entropy)]
 for(k in 1:maxK){best=which.min(cross.entropy(p,K=k));write.csv(Q(p,K=k,run=best),file.path(root,'results',paste0(name,'_Q_K',k,'.csv')),row.names=FALSE)}
 # Hold out 5% of observed calls to assess imputation accuracy.
 x=as.matrix(read.table(input));set.seed(41009);observed=which(x!=9);hide=sample(observed,round(.05*length(observed)));xval=x;xval[hide]=9
 valfile=paste0(name,'_validation.lfmm');write.lfmm(xval,valfile)
 pv=snmf(valfile,K=bestK,repetitions=1,entropy=FALSE,CPU=2,alpha=10,iterations=1000,tolerance=1e-5,seed=20261010,project='new')
 impute(pv,valfile,method='mode',K=bestK,run=1)
 pred=read.lfmm(paste0(valfile,'_imputed.lfmm'))
 modegeno=apply(xval,2,function(v){tt=table(v[v!=9]);as.numeric(names(tt)[which.max(tt)])})
 colidx=col(x)[hide];true=x[hide];est=pred[hide];baseline=modegeno[colidx]
 metrics=data.frame(K=bestK,hidden=length(hide),accuracy=mean(true==est),mode_accuracy=mean(true==baseline),MAE=mean(abs(true-est)),mode_MAE=mean(abs(true-baseline)))
 write.csv(metrics,file.path(root,'results',paste0(name,'_imputation_validation.csv')),row.names=FALSE)
 write.csv(table(truth=true,prediction=est),file.path(root,'results',paste0(name,'_imputation_confusion.csv')),row.names=FALSE)
 # Fit the selected K using all observed calls, then impute missing calls.
 finalfile=paste0(name,'_final.lfmm');write.lfmm(x,finalfile)
 pf=snmf(finalfile,K=bestK,repetitions=1,entropy=FALSE,CPU=2,alpha=10,iterations=1000,tolerance=1e-5,seed=20261011,project='new')
 impute(pf,finalfile,method='mode',K=bestK,run=1)
 write.csv(Q(pf,K=bestK,run=1),file.path(root,'results',paste0(name,'_Q_final.csv')),row.names=FALSE)
 cat('\nFINISHED ',name,' K=',bestK,' validation accuracy=',metrics$accuracy,' baseline=',metrics$mode_accuracy,'\n')
}
writeLines(capture.output(sessionInfo()),file.path(root,'results','R_sessionInfo.txt'))
