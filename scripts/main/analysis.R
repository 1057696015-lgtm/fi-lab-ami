# Frozen MIMIC-IV -> eICU validation; no eICU outcome used in prediction fitting.
args <- commandArgs(trailingOnly=FALSE)
script <- sub('^--file=', '', args[grepl('^--file=',args)])
out <- normalizePath(file.path(dirname(script),'..'))
dest <- file.path(out,'results')
options(digits=15)
M <- read.csv(file.path(out,'data/M_design.csv'), check.names=FALSE)
E <- read.csv(file.path(out,'data/E_design.csv'), check.names=FALSE)
cov <- setdiff(names(M), c('fi01','death30','followup_days','id','pid'))
stopifnot(identical(names(M),names(E)), !anyNA(M[,c(cov,'fi01','death30')]), !anyNA(E[,c(cov,'fi01','death30')]))
auc <- function(y,p) { n1<-sum(y); n0<-length(y)-n1; (sum(rank(round(p,12),ties.method='average')[y==1])-n1*(n1+1)/2)/(n0*n1) }
fit <- function(x,y) {
  f<-glm.fit(x=x,y=y,family=binomial(),control=glm.control(epsilon=1e-10,maxit=100))
  stopifnot(f$converged,all(is.finite(f$coefficients)))
  f
}
metrics <- function(y,lp) {
  p<-plogis(lp)
  cali<-glm.fit(x=matrix(1,length(y),1),y=y,offset=lp,family=binomial())$coefficients[1]
  cals<-fit(cbind(1,lp),y)$coefficients[2]
  c(AUC=auc(y,p), Brier=mean((p-y)^2), O_E=sum(y)/sum(p), Calibration_intercept=unname(cali),Calibration_slope=unname(cals))
}
paired <- function(y,lp0,lp1) {
  a<-metrics(y,lp0);b<-metrics(y,lp1)
  c(setNames(a,paste0('Without.',names(a))),setNames(b,paste0('With.',names(b))),Delta_AUC=unname(b['AUC']-a['AUC']),Brier_improvement=unname(a['Brier']-b['Brier']))
}
res<-coefs<-internal<-tests<-list()
for (scenario in c('survivor24')) {
  m<-if(scenario=='primary')M else M[M$followup_days>1,]
  e<-if(scenario=='primary')E else E[E$followup_days>1,]
  x0<-as.matrix(m[,cov]); x1<-cbind(x0,fi01=m$fi01)
  z0<-as.matrix(e[,cov]); z1<-cbind(z0,fi01=e$fi01)
  f0<-fit(x0,m$death30);f1<-fit(x1,m$death30)
  b0<-f0$coefficients;b1<-f1$coefficients
  lp0<-drop(z0%*%b0);lp1<-drop(z1%*%b1)
  point<-paired(e$death30,lp0,lp1)
  set.seed(if(scenario=='primary')2026102801 else 2026102802)
  B<-1000L
  boot<-t(replicate(B,{
    ix<-sample.int(nrow(e),nrow(e),replace=TRUE)
    paired(e$death30[ix],lp0[ix],lp1[ix])
  }))
  stopifnot(all(is.finite(boot)))
  write.csv(data.frame(replicate=seq_len(B),boot,check.names=FALSE),file.path(dest,paste0(scenario,'_bootstrap_draws.csv')),row.names=FALSE)
  res[[scenario]]<-data.frame(scenario=scenario,n_development=nrow(m),events_development=sum(m$death30),n_validation=nrow(e),events_validation=sum(e$death30),metric=names(point),estimate=point,lower=apply(boot,2,quantile,.025),upper=apply(boot,2,quantile,.975),bootstrap_valid=B,row.names=NULL)
  r0<-pROC::roc(e$death30,round(plogis(lp0),12),levels=c(0,1),direction='<',quiet=TRUE)
  r1<-pROC::roc(e$death30,round(plogis(lp1),12),levels=c(0,1),direction='<',quiet=TRUE)
  dt<-pROC::roc.test(r1,r0,paired=TRUE,method='delong',alternative='two.sided')
  tests[[scenario]]<-data.frame(scenario=scenario,test='Paired two-sided DeLong test of external AUCs',p_value=dt$p.value,statistic=unname(dt$statistic))
  for (lab in c('without','with')) {
    b<-if(lab=='without')b0 else b1
    coefs[[paste(scenario,lab)]]<-data.frame(scenario=scenario,model=lab,term=names(b),coefficient=unname(b))
  }
  write.csv(data.frame(id=e$id,pid=e$pid,death30=e$death30,sofa=e$sofa,fi01=e$fi01,lp_without=lp0,lp_with=lp1,p_without=plogis(lp0),p_with=plogis(lp1)),file.path(out,'data',paste0(scenario,'_external_predictions.csv')),row.names=FALSE)
  saveRDS(list(without=f0,with=f1,predictor_columns=cov,scenario=scenario),file.path(out,'data',paste0(scenario,'_frozen_models.rds')))
  # Internal bootstrap optimism: MIMIC-IV refitting only; not an external CI.
  set.seed(if(scenario=='primary')2026102803 else 2026102804)
  app<-paired(m$death30,drop(x0%*%b0),drop(x1%*%b1))
  op<-t(replicate(200,{
    ix<-sample.int(nrow(m),nrow(m),replace=TRUE)
    bb0<-fit(x0[ix,,drop=FALSE],m$death30[ix])$coefficients
    bb1<-fit(x1[ix,,drop=FALSE],m$death30[ix])$coefficients
    paired(m$death30[ix],drop(x0[ix,,drop=FALSE]%*%bb0),drop(x1[ix,,drop=FALSE]%*%bb1))-paired(m$death30,drop(x0%*%bb0),drop(x1%*%bb1))
  }))
  internal[[scenario]]<-data.frame(scenario=scenario,n=nrow(m),events=sum(m$death30),metric=names(app),apparent=app,mean_optimism=colMeans(op),optimism_corrected=app-colMeans(op),bootstrap_valid=200,row.names=NULL)
  cat('\n',scenario,'\n');print(res[[scenario]][res[[scenario]]$metric%in%c('Without.AUC','With.AUC','Delta_AUC','Brier_improvement','Without.Calibration_slope','With.Calibration_slope'),]);flush.console()
}
write.csv(do.call(rbind,res),file.path(dest,'external_validation.csv'),row.names=FALSE)
write.csv(do.call(rbind,coefs),file.path(dest,'frozen_coefficients.csv'),row.names=FALSE)
write.csv(do.call(rbind,internal),file.path(dest,'internal_validation.csv'),row.names=FALSE)
write.csv(do.call(rbind,tests),file.path(dest,'external_auc_tests.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(dest,'R_sessionInfo.txt'))
