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
neutro_domain_entry('scripts/nhanes/17_run_E1.R')
# Derived candidate; see provenance/script_selection.csv. Invoke via run_job.py.
options(warn=1)
suppressPackageStartupMessages(library(survey))
out <- neutro_normalize('work/descriptive')
csv <- function(x,n) { p<-file.path(out,n);stopifnot(!neutro_file_exists(p));neutro_write_csv(x,p,row.names=FALSE) }
inputs<-c('work/analytic/analytic_v3.rds','work/cohort/C2_full_adult_design.rds')
hash<-neutro_md5(inputs); print(sessionInfo())
d<-as.data.frame(neutro_readRDS(inputs[1])); ds<-neutro_readRDS(inputs[2]); stopifnot(nrow(d)==33174,!anyDuplicated(d$SEQN))
for(v in c('Neutrophil','Lymphocyte','Platelet'))d[[paste0('log_',v)]]<-log(d[[paste0(v,'_count')]])
d$SII<-(d$Neutrophil_count/1000)*(d$Platelet_count/1000)/(d$Lymphocyte_count/1000)
d$log_SII<-log(d$SII);d$PHQ10<-as.numeric(d$PHQ9_total>=10)
stopifnot(all(d$PHQ10==d$depression),all(is.finite(d$log_SII)))
d$sex<-factor(d$sex,levels=c(1,2),labels=c('Male','Female'))
d$race_eth<-factor(d$race_eth);d$education<-factor(d$education,levels=c('Less than high school','High school/GED','Some college','College or above'))
for(v in c('diabetes','hypertension','any_nsaid','any_corticosteroid','any_statin','antidep'))d[[v]]<-factor(d[[v]],levels=c(0,1),labels=c('No','Yes'))
base<-ds$variables;ix<-match(base$SEQN,d$SEQN);base$included<-!is.na(ix)
for(v in names(d))if(!v%in%c('SEQN','SDMVPSU','SDMVSTRA','WTMEC_adj'))base[[v]]<-d[[v]][ix]
base$cycle_factor<-factor(base$SDDSRVYR);m2vars<-c('log10_cotinine','alcohol4','education','pir');base$m2common<-base$included&complete.cases(base[,m2vars]);ds$variables<-base
full<-subset(ds,included);common<-subset(ds,m2common);SDN<-.381754443;SDS<-sqrt(as.numeric(coef(svyvar(~log_SII,full))))
csv(data.frame(N=33174,PHQ10_cases=sum(d$PHQ10),Model2_common_N=nrow(common$variables),Neutrophil_log_SD=SDN,SII_log_SD=SDS,SII_weighted_mean=as.numeric(coef(svymean(~SII,full))),log_SII_weighted_mean=as.numeric(coef(svymean(~log_SII,full))),SII_min=min(d$SII),SII_max=max(d$SII)),'parameters.csv')
cov1<-'age + sex + bmi + smoking + cycle_factor + race_eth';cov2<-paste(cov1,'+ log10_cotinine + alcohol4 + education + pir');cells<-c('log_Neutrophil','log_Lymphocyte','log_Platelet')
registry<-data.frame();fits<-list()
fit<-function(id,f,des,poisson=FALSE,kind='analysis'){
  neutro_cat('\nMODEL_START',id,'\n');print(f);ww<-character()
  m<-withCallingHandlers(svyglm(f,design=des,family=if(kind=='Table1_linear')gaussian()else if(poisson)quasipoisson(link='log')else quasibinomial()),warning=function(w){ww<<-c(ww,conditionMessage(w));neutro_cat('WARNING',conditionMessage(w),'\n');invokeRestart('muffleWarning')})
  print(summary(m));registry<<-rbind(registry,data.frame(model=id,kind=kind,N=nobs(m),formula=paste(deparse(f),collapse=' '),family=family(m)$family,design_df=degf(des),residual_df=m$df.residual,converged=m$converged,warnings=paste(ww,collapse=' | ')))
  fits[[id]]<<-m; m
}
est<-function(m,v,des,delta=1,design_t=FALSE){b<-as.numeric(sum(coef(m)*v));se<-sqrt(as.numeric(t(v)%*%vcov(m)%*%v));df<-if(design_t)degf(des)else m$df.residual;q<-qt(.975,df);data.frame(beta=b,SE=se,delta=delta,effect=exp(b*delta),lower=exp((b-q*se)*delta),upper=exp((b+q*se)*delta),P=2*pt(-abs(b/se),df),df=df)}
vec<-function(m,term){v<-setNames(rep(0,length(coef(m))),names(coef(m)));v[term]<-1;v}
sii<-joint<-constraints<-data.frame()
for(y in c('PHQ10','PHQ9_total')){
 p<-y=='PHQ9_total';m<-fit(paste0('E1_1_SII_',y),as.formula(paste(y,'~ log_SII +',cov1)),full,p)
 for(scale in c('per_ln','per_SD'))sii<-rbind(sii,cbind(data.frame(outcome=y,N=nobs(m),scale=scale),est(m,vec(m,'log_SII'),full,if(scale=='per_ln')1 else SDS)))
 m<-fit(paste0('E1_1_joint_',y),as.formula(paste(y,'~',paste(cells,collapse=' + '),'+',cov1)),full,p)
 for(t in cells)joint<-rbind(joint,cbind(data.frame(outcome=y,term=t),est(m,vec(m,t),full)))
 R<-rbind(vec(m,cells[1])-vec(m,cells[3]),vec(m,cells[1])+vec(m,cells[2]));rb<-R%*%coef(m);rv<-R%*%vcov(m)%*%t(R);F<-as.numeric(t(rb)%*%solve(rv,rb))/2;df<-degf(full)
 constraints<-rbind(constraints,data.frame(outcome=y,N=nobs(m),F=F,df1=2,df2=df,P=pf(F,2,df,lower.tail=FALSE)))
 neutro_cat('SII_CONSTRAINT: R1=betaN-betaP; R2=betaN+betaL; Wald F=(Rb) V^-1 (Rb)/2; denominator df=degf(design)\n')
}
csv(sii,'E1_1_SII.csv');csv(joint,'E1_1_joint_coefficients.csv');csv(constraints,'E1_1_constraints.csv')

# Table1: all nonmissing domains are obtained using subset of the C2-based design.
# Continuous P: survey gaussian linear regression x~PHQ10, default model residual t df.
# Categorical P: svychisq(~category+PHQ10,statistic='F'), Rao-Scott second-order adjustment.
# Continuous SMD=(mean1-mean0)/sqrt((weighted variance1+weighted variance0)/2); report absolute.
# Categorical SMD=sqrt(delta' inverse((multinomial covariance0+covariance1)/2) delta), drop first category.
# Cotinine: geometric mean=10^weighted mean(log10_cotinine), quantiles on concentration scale;
# P and SMD computed on log10_cotinine scale. Percentages use variable-specific nonmissing denominator.
neutro_cat('\nTABLE1_METHODS\n');neutro_cat('Continuous: svyglm Gaussian x~PHQ10, residual-df t P; abs pooled weighted-variance SMD. Categorical: Rao-Scott second-order F P; pooled multinomial covariance generalized SMD. Cotinine geometric mean and Q25/Q75 in ng/mL, P/SMD on log10 scale.\n')
cont<-c('age','bmi','log10_cotinine','pir','PHQ9_total','Neutrophil_count','Lymphocyte_count','Platelet_count','WBC_count','SII')
catvars<-c('sex','race_eth','smoking','education','alcohol4','diabetes','hypertension','any_nsaid','any_corticosteroid','any_statin','antidep')
tab<-meta<-data.frame()
for(v in c(cont,catvars)){
 des<-full;des$variables$.x<-des$variables[[v]];ok<-subset(des,!is.na(.x));g0<-subset(ok,PHQ10==0);g1<-subset(ok,PHQ10==1)
 if(v%in%cont){m<-fit(paste0('Table1_',v),.x~PHQ10,ok,kind='Table1_linear');cs<-coef(summary(m));pv<-cs['PHQ10',ncol(cs)]
  mu0<-as.numeric(coef(svymean(~.x,g0)));mu1<-as.numeric(coef(svymean(~.x,g1)));var0<-as.numeric(coef(svyvar(~.x,g0)));var1<-as.numeric(coef(svyvar(~.x,g1)));smd<-abs(mu1-mu0)/sqrt((var0+var1)/2)
 }else{ok$variables$.x<-factor(ok$variables$.x);g0<-subset(ok,PHQ10==0);g1<-subset(ok,PHQ10==1);test<-svychisq(~.x+PHQ10,ok,statistic='F');pv<-test$p.value;neutro_cat('RAO_SCOTT',v,'\n');print(test)
  pr<-function(z)as.numeric(coef(svymean(~.x,z)));p0<-pr(g0);p1<-pr(g1);dd<-(p1-p0)[-1];C<-((diag(p0)-tcrossprod(p0))+(diag(p1)-tcrossprod(p1)))/2;C<-C[-1,-1,drop=FALSE];smd<-sqrt(as.numeric(t(dd)%*%solve(C,dd)))
 }
 meta<-rbind(meta,data.frame(variable=v,P=pv,SMD=smd,P_method=if(v%in%cont)'survey Gaussian x~PHQ10; residual-df t'else 'Rao-Scott second-order F',SMD_method=if(v%in%cont)'absolute standardized mean difference; pooled survey variance'else 'multinomial pooled covariance'))
 for(group in c('Overall','PHQ10_0','PHQ10_1')){
  gd<-if(group=='Overall')des else if(group=='PHQ10_0')subset(des,PHQ10==0)else subset(des,PHQ10==1);nd<-subset(gd,!is.na(.x));N<-nrow(nd$variables);missing<-nrow(gd$variables)-N
  if(v%in%cont){mu<-as.numeric(coef(svymean(~.x,nd)));sd<-sqrt(as.numeric(coef(svyvar(~.x,nd))));qq<-if(v=='log10_cotinine')as.numeric(coef(svyquantile(~.x,nd,quantiles=c(.25,.75),ci=FALSE)))else c(NA,NA)
   tab<-rbind(tab,data.frame(variable=v,category='',group=group,N_nonmissing=N,N_missing=missing,n=NA,weighted_percent=NA,mean=mu,SD=sd,geometric_mean=if(v=='log10_cotinine')10^mu else NA,Q25=if(v=='log10_cotinine')10^qq[1]else NA,Q75=if(v=='log10_cotinine')10^qq[2]else NA,P=pv,SMD=smd))
  }else{pp<-as.numeric(coef(svymean(~.x,nd)));nn<-table(factor(nd$variables$.x,levels=levels(ok$variables$.x)));tab<-rbind(tab,data.frame(variable=v,category=names(nn),group=group,N_nonmissing=N,N_missing=missing,n=as.numeric(nn),weighted_percent=100*pp,mean=NA,SD=NA,geometric_mean=NA,Q25=NA,Q75=NA,P=pv,SMD=smd))}
 }
}
csv(tab,'E1_3_Table1.csv');csv(meta,'E1_3_methods.csv')
csv(registry,'model_registry_part1.csv')
neutro_saveRDS(list(design=ds,registry=registry,fits=fits),file.path(out,'E1_part1_state.rds'))
stopifnot(all(hash==neutro_md5(inputs)));csv(data.frame(path=inputs,md5_before=unname(hash),md5_after=unname(neutro_md5(inputs))),'source_hashes_part1.csv')
neutro_cat('PART1_COMPLETE\n')
