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
neutro_domain_entry('scripts/nhanes/dependency_run_B6.R')
options(warn=1)
suppressPackageStartupMessages(library(survey))
out<-neutro_normalize('work/inflammation');csv<-function(x,n){p<-file.path(out,n);stopifnot(!neutro_file_exists(p));neutro_write_csv(x,p,row.names=FALSE)}
paths<-c('work/cohort/supplement_input.rds','work/alcohol/nhanes_final_enriched_with_alcohol4_B4.rds');hash<-neutro_md5(paths);original<-neutro_readRDS(paths[1]);d<-neutro_readRDS(paths[2]);stopifnot(nrow(d)==33174,all(vapply(names(original),function(v)identical(original[[v]],d[[v]]),logical(1))))
neutro_cat('START',as.character(Sys.time()),'\n');print(sessionInfo());neutro_cat('Contrast test = t, df = degf(design). B6-4 M3 and E1 use the same age>=20 E1 complete-case domain as B6-3.\n')
d$cycle_factor<-factor(d$SDDSRVYR);ev<-c('alcohol4','education','pir','diabetes','hypertension','any_nsaid','any_corticosteroid','any_statin');d$E1_complete<-complete.cases(d[,ev])
ds<-svydesign(id=~SDMVPSU,strata=~SDMVSTRA,weights=~WTMEC_adj,nest=TRUE,data=d);adult<-subset(ds,age>=20);common<-subset(adult,E1_complete)
csv(data.frame(sample=c('original','age_ge20','age_ge20_E1_common'),N=c(nrow(d),nrow(adult$variables),nrow(common$variables)),depression_cases=c(sum(d$depression==1),sum(adult$variables$depression==1),sum(common$variables$depression==1)),fewer_than_33174=33174-c(nrow(d),nrow(adult$variables),nrow(common$variables))),'B6_1_sample_counts.csv')
excluded<-do.call(rbind,lapply(4:10,function(cy){x<-d[d$SDDSRVYR==cy&d$age<20,];data.frame(cycle=paste(1997+2*cy,1998+2*cy,sep='-'),N=nrow(x),depression_cases=sum(x$depression==1))}));csv(excluded,'B6_1_excluded_by_cycle.csv');stopifnot(all(d$age[d$age<20]%in%18:19))
base<-'log_Neutrophil + log_Lymphocyte + log_Platelet + age + sex + bmi + smoking + cycle_factor + race_eth';cells<-c('log_Neutrophil','log_Lymphocyte','log_Platelet');extras<-c(M3='',E2='+ alcohol4 + education + pir',E1='+ alcohol4 + education + pir + diabetes + hypertension + any_nsaid + any_corticosteroid + any_statin')
registry<-coeffs<-cts<-data.frame();fits<-list()
fit<-function(id,des,extra='',outcome='depression',pois=FALSE,contrast=FALSE){f<-as.formula(paste(outcome,'~',base,extra));ww<-character();neutro_cat('\nMODEL_START',id,'N=',nrow(des$variables),'\n');print(f)
 m<-withCallingHandlers(svyglm(f,design=des,family=if(pois)quasipoisson(link='log')else quasibinomial()),warning=function(w){ww<<-c(ww,conditionMessage(w));neutro_cat('WARNING',conditionMessage(w),'\n');invokeRestart('muffleWarning')});print(summary(m));ci<-confint(m);print(ci);cs<-coef(summary(m));b<-coef(m)[cells]
 coeffs<<-rbind(coeffs,data.frame(model=id,N=nobs(m),term=cells,beta=unname(b),SE=unname(sqrt(diag(vcov(m)))[cells]),effect=unname(exp(b)),lower=unname(exp(ci[cells,1])),upper=unname(exp(ci[cells,2])),P=unname(cs[cells,ncol(cs)])))
 if(contrast)for(j in 2:3){v<-setNames(rep(0,length(coef(m))),names(coef(m)));v[cells[1]]<-1;v[cells[j]]<--1;c<-svycontrast(m,v);bb<-as.numeric(coef(c));se<-as.numeric(SE(c));ct<-data.frame(model=id,contrast=paste(cells[1],cells[j],sep=' - '),beta=bb,SE=se,t=bb/se,df=degf(des),P=2*pt(-abs(bb/se),degf(des)));cts<<-rbind(cts,ct);print(ct)}
 registry<<-rbind(registry,data.frame(model=id,N=nobs(m),formula=paste(deparse(f),collapse=' '),family=family(m)$family,design_df=degf(des),residual_df=m$df.residual,converged=m$converged,warning_count=length(ww),warnings=paste(ww,collapse=' | ')));fits[[id]]<<-m;neutro_cat('MODEL_END',id,'\n')}
fit('B6_2_age20_M3',adult,contrast=TRUE)
for(id in names(extras))fit(paste0('B6_3_common_',id),common,extras[id])
lsd<-sqrt(as.numeric(coef(svyvar(~log_Neutrophil,common))));scale<-do.call(rbind,lapply(names(extras),function(id){m<-fits[[paste0('B6_3_common_',id)]];b<-coef(m)['log_Neutrophil'];ci<-confint(m)['log_Neutrophil',];data.frame(model=id,N=nobs(m),weighted_log_SD=lsd,OR_per_SD=unname(exp(b*lsd)),lower=unname(exp(ci[1]*lsd)),upper=unname(exp(ci[2]*lsd)))}));csv(scale,'B6_3_effect_per_SD.csv')
fit('B6_4_common_M3_PHQtotal',common,outcome='PHQ9_total',pois=TRUE);fit('B6_4_common_E1_PHQtotal',common,extras['E1'],outcome='PHQ9_total',pois=TRUE)
csv(coeffs,'B6_coefficients.csv');csv(cts,'B6_t_contrasts.csv');csv(registry,'B6_model_registry.csv');csv(data.frame(path=paths,md5_before=unname(hash),md5_after=unname(neutro_md5(paths))),'B6_source_hashes.csv')
stopifnot(nrow(registry)==6,all(hash==neutro_md5(paths)),all(registry$converged),all(registry$warning_count==0),all(registry$N[-1]==nrow(common$variables)),nrow(cts)==2)
# Independently verify the t-contrast formula from coefficients and covariance.
m<-fits$B6_2_age20_M3;for(j in 2:3){se<-sqrt(vcov(m)[cells[1],cells[1]]+vcov(m)[cells[j],cells[j]]-2*vcov(m)[cells[1],cells[j]]);pv<-2*pt(-abs((coef(m)[cells[1]]-coef(m)[cells[j]])/se),degf(adult));stopifnot(abs(pv-cts$P[j-1])<1e-12)}
csv(data.frame(check=c('original_columns_identical_to_B4','source_hashes_unchanged','six_specified_models','common_sample_consistent','t_contrasts_verified'),passed=TRUE),'B6_validation.csv');csv(data.frame(status='COMPLETE_WAITING_USER',models_run=6),'B6_status.csv');neutro_cat('END COMPLETE_WAITING_USER\n')
