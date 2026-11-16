args<-commandArgs(FALSE);s<-sub('^--file=','',args[grepl('^--file=',args)]);out<-normalizePath(file.path(dirname(s),'..'))
flows<-tabs<-list()
for(k in c('M','E')){
 d<-read.csv(file.path(out,'data',paste0(k,'_design.csv')),check.names=FALSE)
 flows[[k]]<-data.frame(cohort=k,before=nrow(d),excluded_deaths=sum(d$followup_days<=1 & d$death30==1),excluded_other_discharges=sum(d$followup_days<=1 & d$death30==0),included=sum(d$followup_days>1),events=sum(d$death30[d$followup_days>1]))
 d<-d[d$followup_days>1,];stopifnot(all(d$followup_days>1));write.csv(d,file.path(out,'data',paste0(k,'_survivor24.csv')),row.names=FALSE)
 vars<-setdiff(names(d),c('const','id','pid','followup_days'))
 tabs[[k]]<-do.call(rbind,lapply(vars,function(v){x<-d[[v]];data.frame(cohort=k,variable=v,n=length(x),mean=mean(x),sd=sd(x),sum=sum(x),percent=100*mean(x))}))
}
write.csv(do.call(rbind,flows),file.path(out,'results/exclusion_audit.csv'),row.names=FALSE)
write.csv(do.call(rbind,tabs),file.path(out,'results/baseline.csv'),row.names=FALSE)
