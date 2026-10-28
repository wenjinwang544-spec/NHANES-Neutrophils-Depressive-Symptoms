# BEGIN MANAGED ENTRY POINT -- fails before study code without explicit workspace.
.neutro_source_root <- Sys.getenv("NEUTRO_SOURCE_ROOT")
if (!nzchar(.neutro_source_root)) {
 .entry <- grep("^--file=",commandArgs(FALSE),value=TRUE)
 if(length(.entry)!=1L)stop("Use Rscript with --workspace or the managed Python launcher")
 .cursor <- dirname(normalizePath(gsub("~+~"," ",sub("^--file=","",.entry),fixed=TRUE),mustWork=TRUE))
 repeat {if(file.exists(file.path(.cursor,"scripts/lib/safety.R"))){.neutro_source_root<-.cursor;break};.parent<-dirname(.cursor);if(.parent==.cursor)stop("Cannot locate managed source root");.cursor<-.parent}
}
base::source(file.path(.neutro_source_root,"scripts/lib/safety.R"),local=.GlobalEnv)
neutro_setup(.neutro_source_root)
if("--check-only" %in% commandArgs(TRUE)){cat("SAFE_ENTRY_CHECK_ONLY\n");quit(status=0L)}
# END MANAGED ENTRY POINT
neutro_domain_entry('scripts/nhanes/30_run_G3b.R')
# Derived candidate; see provenance/script_selection.csv. Invoke via run_job.py.
options(warn=1,digits=17)
suppressPackageStartupMessages(library(survey));suppressPackageStartupMessages(library(metafor))
out<-'work/heterogeneity';SDfixed<-0.381754443
csv<-function(x,n){p<-file.path(out,n);stopifnot(!neutro_file_exists(p));neutro_write_csv(x,p,row.names=FALSE)}
neutro_cat('START',as.character(Sys.time()),'\n');print(sessionInfo());neutro_capture_output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
inputs<-c('work/analytic/analytic_v3.rds','work/cohort/C2_full_adult_design.rds','work/shapley/G3_2_specifications.csv','work/shapley/G3_3_cycle_estimates.csv');hash<-neutro_md5(inputs)
d<-as.data.frame(neutro_readRDS(inputs[1]));des<-neutro_readRDS(inputs[2]);b<-des$variables;ix<-match(b$SEQN,d$SEQN);b$included<-!is.na(ix)
for(v in setdiff(names(d),c('SEQN','SDMVPSU','SDMVSTRA','WTMEC_adj')))b[[v]]<-d[[v]][ix]
for(v in c('Neutrophil','Lymphocyte','Platelet'))b[[paste0('log_',v)]]<-log(b[[paste0(v,'_count')]])
extra<-c('log10_cotinine','alcohol4','education','pir','diabetes','hypertension','any_nsaid','any_corticosteroid','any_statin');b$cc<-b$included&complete.cases(b[,extra]);des$variables<-b
full<-subset(des,included);cc<-subset(des,cc);stopifnot(nrow(full$variables)==33174,nrow(cc$variables)==29360)
cov1<-c('age','sex','bmi','smoking','race_eth');cells<-c('log_Neutrophil','log_Lymphocyte','log_Platelet')
cycles<-sort(unique(d$SDDSRVYR));fits<-list();reg<-data.frame();result<-data.frame()
neutro_cat('Cycle factor omitted within every cycle: constant. All cycle coefficient CIs and P use t, df=degf(cycle subset design). SE_SD=SE_log_N*0.381754443.\n')
for(mod in c('Model1','Model3'))for(cy in cycles){dd<-if(mod=='Model1')subset(full,SDDSRVYR==cy)else subset(cc,SDDSRVYR==cy);f<-reformulate(c(cells,cov1,if(mod=='Model3')extra),'depression');id<-paste(mod,cy,sep='_');neutro_cat('MODEL_START',id,'N',nrow(dd$variables),'\n');print(f);ww<-character()
 m<-withCallingHandlers(svyglm(f,design=dd,family=quasibinomial()),warning=function(w){ww<<-c(ww,conditionMessage(w));neutro_cat('WARNING',conditionMessage(w),'\n');invokeRestart('muffleWarning')});fits[[id]]<-m
 beta<-unname(coef(m)['log_Neutrophil']);se<-unname(sqrt(vcov(m)['log_Neutrophil','log_Neutrophil']));df<-degf(dd);q<-qt(.975,df);z<-beta*SDfixed;ss<-se*SDfixed
 result<-rbind(result,data.frame(model=mod,cycle=cy,period=paste0(2005+2*(cy-4),'-',2006+2*(cy-4)),N=nrow(dd$variables),cases=sum(dd$variables$depression==1),beta_ln_unit=beta,SE_ln_unit=se,ln_OR_SD=z,SE_ln_OR_SD=ss,OR_SD=exp(z),lower=exp(z-q*ss),upper=exp(z+q*ss),t=z/ss,df=df,P=2*pt(-abs(z/ss),df),residual_df=m$df.residual))
 reg<-rbind(reg,data.frame(id=id,formula=paste(deparse(f),collapse=' '),N=nobs(m),converged=m$converged,warnings=paste(ww,collapse=' | ')));neutro_cat('MODEL_END',id,'df',df,'\n')}
old<-neutro_read_csv(inputs[4]);j<-match(paste(result$model,result$cycle),paste(old$model,old$cycle));check<-data.frame(model=result$model,cycle=result$cycle,beta_abs_difference=abs(result$beta_ln_unit-old$beta[j]),SE_abs_difference=abs(result$SE_ln_unit-old$SE[j]));stopifnot(max(check$beta_abs_difference)<1e-10,max(check$SE_abs_difference)<1e-10)
csv(result,'G3b_cycle_estimates.csv');csv(reg,'G3b_model_registry.csv');csv(check,'G3b_G3_coefficient_check.csv');neutro_saveRDS(fits,file.path(out,'G3b_cycle_models.rds'))
metafits<-list();meta<-loo<-data.frame()
combine<-function(z,mod,omit='none'){
 neutro_cat('META_START',mod,'omitted',omit,'k',nrow(z),'method DL; normal z and unmodified HKSJ t (knha)\n')
 a<-rma.uni(yi=z$ln_OR_SD,vi=z$SE_ln_OR_SD^2,method='DL',test='z');h<-rma.uni(yi=z$ln_OR_SD,vi=z$SE_ln_OR_SD^2,method='DL',test='knha');metafits[[paste(mod,omit,'DL',sep='_')]]<<-a;metafits[[paste(mod,omit,'HKSJ',sep='_')]]<<-h
 row<-function(x,method)data.frame(model=mod,omitted_period=omit,method=method,k=x$k,ln_OR_SD=as.numeric(x$b),SE=as.numeric(x$se),OR=exp(as.numeric(x$b)),lower=exp(x$ci.lb),upper=exp(x$ci.ub),tau2=x$tau2,I2=x$I2,df=if(method=='HKSJ')x$k-1 else Inf,P=x$pval)
 neutro_cat('META_END',mod,omit,'\n');rbind(row(a,'DL_normal'),row(h,'HKSJ'))}
for(mod in c('Model1','Model3')){z<-result[result$model==mod,];meta<-rbind(meta,combine(z,mod));for(cy in cycles)loo<-rbind(loo,combine(z[z$cycle!=cy,],mod,z$period[z$cycle==cy]))}
csv(meta,'G3b_meta_DL_HKSJ.csv');csv(loo,'G3b_leave_one_DL_HKSJ.csv');neutro_saveRDS(metafits,file.path(out,'G3b_meta_objects.rds'))
forest<-function(){par(mfrow=c(1,2),mar=c(4,10,3,2));for(mod in c('Model1','Model3')){z<-result[result$model==mod,];m<-meta[meta$model==mod,];or<-c(z$OR_SD,m$OR);lo<-c(z$lower,m$lower);hi<-c(z$upper,m$upper);yy<-9:1;plot(or,yy,xlim=range(c(lo,hi,1)),ylim=c(.5,9.5),yaxt='n',xlab='OR per 1 SD',ylab='',main=mod,pch=c(rep(16,7),18,17));segments(lo,yy,hi,yy);axis(2,at=yy,labels=c(paste0(z$period,' (N=',z$N,')'),'DL normal','HKSJ'),las=2,cex.axis=.85);abline(v=1,lty=2)}}
neutro_png(file.path(out,'G3b_cycle_forest.png'),width=14,height=7,units='in',res=300);forest();neutro_dev_off();neutro_pdf(file.path(out,'G3b_cycle_forest.pdf'),width=14,height=7);forest();neutro_dev_off()
# Existing specification model estimates only; no specification refitting.
s<-neutro_read_csv(inputs[3]);a<-s[s$restriction=='age20',];n<-s[s$restriction=='none',];j<-match(paste(a$mask,a$exposure),paste(n$mask,n$exposure));stopifnot(all(a$N==n$N[j]),all(abs(a$ln_OR_SD-n$ln_OR_SD[j])<1e-12),all(abs(a$lower-n$lower[j])<1e-12),all(abs(a$upper-n$upper[j])<1e-12))
u<-s[s$restriction!='age20',];stopifnot(nrow(u)==192);u<-u[order(u$ln_OR_SD,u$id),];u$rank<-seq_len(nrow(u));csv(u,'G3b_specifications_192.csv')
stat<-function(z,label)data.frame(set=label,n=nrow(z),median_ln_OR_SD=median(z$ln_OR_SD),min_ln_OR_SD=min(z$ln_OR_SD),max_ln_OR_SD=max(z$ln_OR_SD),median_OR=median(z$OR),min_OR=min(z$OR),max_OR=max(z$OR),CI_excludes1=sum(z$lower>1|z$upper<1),CI_excludes1_proportion=mean(z$lower>1|z$upper<1),CI_contains1=sum(z$lower<=1&z$upper>=1))
stats<-rbind(stat(s,'before_256'),stat(u,'after_192'));csv(stats,'G3b_specification_comparison.csv');inside<-u[u$lower<=1&u$upper>=1,];inside<-inside[order(inside$restriction,inside$exposure,inside$mask),];csv(inside,'G3b_CI_contains1.csv')
csv(data.frame(removed_restriction='age20',retained_restriction='none',removed_n=64,identical_N=TRUE,identical_beta=TRUE,identical_CI=TRUE),'G3b_dedup_check.csv')
plot_spec<-function(){layout(matrix(c(1,2),2,1),heights=c(2,1.4));par(mar=c(1,9,2,1));plot(u$rank,u$ln_OR_SD,xlim=c(.5,192.5),xaxs='i',xaxt='n',pch=16,cex=.45,xlab='',ylab='ln OR per 1 SD',ylim=range(log(c(u$lower,u$upper))));segments(u$rank,log(u$lower),u$rank,log(u$upper),col=ifelse(u$lower>1|u$upper<1,'#20689b','#999999'));points(u$rank,u$ln_OR_SD,pch=16,cex=.4);abline(h=0,lty=2);par(mar=c(4,9,0,1));mat<-rbind(t(sapply(1:5,function(g)as.integer(bitwAnd(u$mask,bitwShiftL(1,g-1))>0))),joint=as.integer(u$exposure=='joint'),stroke0=as.integer(u$restriction=='stroke0'),extreme=as.integer(u$restriction=='exclude_extreme'));image(1:192,1:nrow(mat),t(mat[nrow(mat):1,,drop=FALSE]),xlim=c(.5,192.5),xaxs='i',col=c('white','#20689b'),axes=FALSE,xlab='Specification rank',ylab='');axis(1);axis(2,at=1:nrow(mat),labels=rev(c('A Cotinine','B Alcohol','C Education/PIR','D Diabetes/HTN','E Medication','Joint counts','Stroke=0','Exclude extreme')),las=2,cex.axis=.8)}
neutro_png(file.path(out,'G3b_specification_curve_192.png'),width=12,height=8,units='in',res=300);plot_spec();neutro_dev_off();neutro_pdf(file.path(out,'G3b_specification_curve_192.pdf'),width=12,height=8);plot_spec();neutro_dev_off()
stopifnot(all(hash==neutro_md5(inputs)));csv(data.frame(path=inputs,md5_before=unname(hash),md5_after=unname(neutro_md5(inputs)),unchanged=TRUE),'G3b_input_integrity.csv')
neutro_writeLines(c('Cycle coefficient intervals: t distribution, df=degf(cycle subset design); cycle factor omitted as constant.','Fixed log neutrophil SD: 0.381754443; reported SE_ln_OR_SD is standard error on per-SD log OR scale.','Both meta intervals use DL tau-squared estimation. DL_normal: inverse-variance random-effects estimate with normal CI. HKSJ: metafor test=knha, unmodified Hartung-Knapp variance adjustment and t(k-1) interval; no ad hoc lower bound on adjustment.','Specification estimates read from G3 CSV without refitting. 64 age>=20 specifications removed after verifying identical N, estimate and CI to unrestricted counterparts; unrestricted versions retained.','All originals unchanged; no package installation.'),file.path(out,'G3b_methods.txt'))
print(result);print(meta);print(loo);print(stats);neutro_cat('COMPLETE',as.character(Sys.time()),'\n')
