args<-commandArgs(FALSE); script<-sub('^--file=','',args[grepl('^--file=',args)]); out<-normalizePath(file.path(dirname(script),'..')); root<-normalizePath(file.path(out,'../..'))
# Figure contract: quantify calibration adaptation and threshold-dependent net benefit.
# Quantitative grid; R-only figures; 180 x 125 mm; PDF and 300-dpi PNG.
library(ggplot2); library(patchwork)
auc<-function(y,p){a<-sum(y);b<-length(y)-a;(sum(rank(p)[y==1])-a*(a+1)/2)/(a*b)}
met<-function(y,lp){p<-plogis(lp); c(AUC=auc(y,p),Brier=mean((y-p)^2),O_E=sum(y)/sum(p),Intercept=unname(coef(glm(y~1+offset(lp),family=binomial()))),Slope=unname(coef(glm(y~lp,family=binomial()))[2]))}
ths<-seq(.05,.50,.01)
nb<-function(y,p) sapply(ths,function(t) mean((p>=t)*(y-(1-y)*t/(1-t))))
allraw<-allsummary<-dcalist<-list()
for(sc in c('survivor24')){
 d<-read.csv(file.path(out,'data',paste0(sc,'_external_predictions.csv')))
 stopifnot(!anyDuplicated(d$pid),all(d$death30%in%0:1),all(is.finite(d$lp_with)),all(is.finite(d$lp_without)))
 y<-d$death30;n<-nrow(d); set.seed(if(sc=='primary')2026102901 else 2026102902)
 draws<-vector('list',1000)
 for(b in 1:1000){
  tr<-sort(c(sample(which(y==0),floor(sum(y==0)/2)),sample(which(y==1),floor(sum(y==1)/2))))
  te<-setdiff(seq_len(n),tr)
  rr<-list()
  for(md in c('without','with')){
   lp<-d[[paste0('lp_',md)]]; f<-glm(y[tr]~lp[tr],family=binomial()); stopifnot(f$converged,coef(f)[2]>0)
   newlp<-unname(coef(f)[1]+coef(f)[2]*lp[te]); old<-met(y[te],lp[te]); new<-met(y[te],newlp)
   for(st in c('Frozen','Recalibrated')) rr[[paste(md,st)]]<-data.frame(scenario=sc,split=b,model=md,state=st,n_train=length(tr),events_train=sum(y[tr]),n_test=length(te),events_test=sum(y[te]),alpha=coef(f)[1],beta=coef(f)[2],t(if(st=='Frozen')old else new),row.names=NULL)
  };draws[[b]]<-do.call(rbind,rr)
 }
 raw<-do.call(rbind,draws);allraw[[sc]]<-raw
 summ<-list()
 for(md in c('without','with')) for(st in c('Frozen','Recalibrated')) for(mt in c('AUC','Brier','O_E','Intercept','Slope')){
  z<-raw[raw$model==md & raw$state==st,mt]
  summ[[paste(md,st,mt)]]<-data.frame(scenario=sc,model=md,state=st,metric=mt,n=n,events=sum(y),n_test=raw$n_test[1],events_test=raw$events_test[1],median=median(z),lower=quantile(z,.025),upper=quantile(z,.975),interval='2.5th-97.5th split percentiles; not a confidence interval',row.names=NULL)
 };allsummary[[sc]]<-do.call(rbind,summ)
 p0<-d$p_without;p1<-d$p_with
 point0<-nb(y,p0);point1<-nb(y,p1)
 set.seed(if(sc=='primary')2026102903 else 2026102904)
 boot<-replicate(1000,{ix<-sample.int(n,n,TRUE);nb(y[ix],p1[ix])-nb(y[ix],p0[ix])})
 dc<-data.frame(scenario=sc,n=n,events=sum(y),threshold=ths,without=point0,with=point1,all=mean(y)-(1-mean(y))*ths/(1-ths),none=0,difference=point1-point0,lower=apply(boot,1,quantile,.025),upper=apply(boot,1,quantile,.975))
 dcalist[[sc]]<-dc
 cat(sc,'completed\n');flush.console()
}
r<-do.call(rbind,allraw);s<-do.call(rbind,allsummary);dc<-do.call(rbind,dcalist)
write.csv(r,file.path(out,'results/recalibration_splits.csv'),row.names=FALSE)
write.csv(s,file.path(out,'results/recalibration_summary.csv'),row.names=FALSE)
write.csv(dc,file.path(out,'results/decision_curves.csv'),row.names=FALSE)
# Paired improvement from the same held-out patients, across splits.
pairs<-merge(r[r$state=='Frozen',],r[r$state=='Recalibrated',],by=c('scenario','split','model'),suffixes=c('_frozen','_updated'))
pairs$Brier_improvement<-pairs$Brier_frozen-pairs$Brier_updated
write.csv(pairs,file.path(out,'results/paired_recalibration.csv'),row.names=FALSE)
