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
neutro_domain_entry('scripts/nhanes/16_run_D2_3.R')
# Derived candidate; see provenance/script_selection.csv. Invoke via run_job.py.
options(warn=1);suppressPackageStartupMessages(library(survey))
out<-neutro_normalize('work/analytic');csv<-function(x,n){p<-file.path(out,n);stopifnot(!neutro_file_exists(p));neutro_write_csv(x,p,row.names=FALSE)}
neutro_cat('NEW INSTRUCTION; no previous approvals retained.\n');print(sessionInfo());inputs<-c('work/clinical/analytic_v2_antidepressant_pending.rds','work/cohort/C2_full_adult_design.rds');sources<-data.frame(path=inputs,md5=unname(neutro_md5(inputs)));d<-neutro_readRDS(inputs[1]);ds<-neutro_readRDS(inputs[2])
base<-ds$variables;ix<-match(base$SEQN,d$SEQN);base$D2_included<-!is.na(ix);for(v in c('alcohol4','log10_cotinine','LBXCOT','CRP_mgL','log_CRP','education','pir','diabetes','hypertension','any_nsaid','any_corticosteroid','any_statin','race_eth','sex','smoking'))base[[v]]<-d[[v]][ix];base$cycle_factor<-factor(base$SDDSRVYR);ev<-c('log10_cotinine','alcohol4','education','pir','diabetes','hypertension','any_nsaid','any_corticosteroid','any_statin');base$common<-base$D2_included&complete.cases(base[,ev]);base$PHQ15<-as.numeric(base$PHQ9_total>=15);base$cotinine3<-factor(ifelse(is.na(base$LBXCOT),NA,ifelse(base$LBXCOT<.015,'below_LOD',ifelse(base$LBXCOT<3,'LOD_to_lt3','ge3'))),levels=c('below_LOD','LOD_to_lt3','ge3'));ds$variables<-base;full<-subset(ds,D2_included);common<-subset(ds,common);stopifnot(nrow(common$variables)==29360)
covars<-'log_Neutrophil + log_Lymphocyte + log_Platelet + age + sex + bmi + smoking + cycle_factor + race_eth';cells<-c('log_Neutrophil','log_Lymphocyte','log_Platelet');fits<-list();registry<-cf<-ct<-data.frame();extra<-c(Model1='',Model2='+ log10_cotinine + alcohol4 + education + pir',Model3='+ log10_cotinine + alcohol4 + education + pir + diabetes + hypertension + any_nsaid + any_corticosteroid + any_statin')
fit<-function(id,domain,extra='',outcome='depression',pois=FALSE){f<-as.formula(paste(outcome,'~',covars,extra));ww<-character();neutro_cat('\nMODEL_START',id,'\n');print(f);m<-withCallingHandlers(svyglm(f,design=domain,family=if(pois)quasipoisson(link='log')else quasibinomial()),warning=function(w){ww<<-c(ww,conditionMessage(w));neutro_cat('WARNING',conditionMessage(w),'\n');invokeRestart('muffleWarning')});print(summary(m));ci<-confint(m);cs<-coef(summary(m));b<-coef(m)[cells];cf<<-rbind(cf,data.frame(model=id,N=nobs(m),term=cells,beta=unname(b),SE=unname(sqrt(diag(vcov(m)))[cells]),effect=unname(exp(b)),lower=unname(exp(ci[cells,1])),upper=unname(exp(ci[cells,2])),P=unname(cs[cells,ncol(cs)])))
 for(j in 2:3){v<-setNames(rep(0,length(coef(m))),names(coef(m)));v[cells[1]]<-1;v[cells[j]]<--1;c<-svycontrast(m,v);bb<-as.numeric(coef(c));se<-as.numeric(SE(c));r<-data.frame(model=id,contrast=paste(cells[1],cells[j],sep=' - '),beta=bb,SE=se,df=degf(domain),P=2*pt(-abs(bb/se),degf(domain)));ct<<-rbind(ct,r);print(r)}
 registry<<-rbind(registry,data.frame(model=id,N=nobs(m),formula=paste(deparse(f),collapse=' '),family=family(m)$family,design_df=degf(domain),residual_df=m$df.residual,converged=m$converged,warning_count=length(ww),warnings=paste(ww,collapse=' | ')));fits[[id]]<<-m;neutro_cat('MODEL_END',id,'\n')}


for(id in names(extra))fit(paste0('D2_3a_',id),common,extra[id])
for(id in names(extra)){fit(paste0('D2_3b_PHQtotal_',id),common,extra[id],outcome='PHQ9_total',pois=TRUE);fit(paste0('D2_3b_PHQ15_',id),common,extra[id],outcome='PHQ15')}
crp<-subset(ds,common&SDDSRVYR%in%c(4,5,6,9,10)&!is.na(log_CRP));fit('D2_3c_Model3',crp,extra['Model3']);fit('D2_3c_Model3_logCRP',crp,paste(extra['Model3'],'+ log_CRP'))
fit('D2_3d_Model2_cotinine3',common,'+ cotinine3 + alcohol4 + education + pir')
lsd<-sqrt(as.numeric(coef(svyvar(~log_Neutrophil,common))));qq<-as.numeric(coef(svyquantile(~Neutrophil_count,common,quantiles=c(.25,.75),ci=FALSE)));sc<-data.frame();for(id in names(extra)){m<-fits[[paste0('D2_3a_',id)]];b<-coef(m)['log_Neutrophil'];ci<-confint(m)['log_Neutrophil',];delta<-c(lsd,log(qq[2]/qq[1]));sc<-rbind(sc,data.frame(model=id,scale=c('per_weighted_SD','Q25_to_Q75'),OR=unname(exp(b*delta)),lower=unname(exp(ci[1]*delta)),upper=unname(exp(ci[2]*delta))))};csv(sc,'D2_3a_effect_scales.csv');csv(data.frame(N=29360,log_SD=lsd,Q25=qq[1],Q75=qq[2]),'D2_3a_scale_parameters.csv');csv(as.data.frame(table(common$variables$cotinine3)),'D2_3d_cotinine3_counts.csv')
csv(cf,'D2_coefficients.csv');csv(ct,'D2_t_contrasts.csv');csv(registry,'D2_model_registry.csv');sources<-unique(sources);stopifnot(all(sources$md5==neutro_md5(sources$path)),nrow(registry)==12);csv(sources,'D2_source_hashes.csv');csv(data.frame(status='D2_3_COMPLETE_D2_1_D2_2_CONFIRMATION_PENDING',models_run=12),'D2_status.csv');neutro_cat('END MI_CONFIRMATION_PENDING\n')
