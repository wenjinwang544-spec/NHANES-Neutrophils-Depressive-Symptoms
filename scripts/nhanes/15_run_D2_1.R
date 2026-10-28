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
neutro_domain_entry('scripts/nhanes/15_run_D2_1.R')
# Derived candidate; see provenance/script_selection.csv. Invoke via run_job.py.
options(warn=1);suppressPackageStartupMessages({library(survey);library(haven)})
out<-neutro_normalize('work/analytic');csv<-function(x,n){p<-file.path(out,n);stopifnot(!neutro_file_exists(p));neutro_write_csv(x,p,row.names=FALSE)}
neutro_cat('START',as.character(Sys.time()),'\n');print(sessionInfo());neutro_cat('All five candidate drugs included per user. All designs from C2 full base then subset; contrasts t df=degf(design). MI awaits confirmation. D2-3 CRP paired models use common complete cases for Model3 and logCRP.\n')
inputs<-c('work/clinical/analytic_v2_antidepressant_pending.rds','work/cohort/C2_full_adult_design.rds');sources<-data.frame(path=inputs,md5=unname(neutro_md5(inputs)));d<-neutro_readRDS(inputs[1]);ds<-neutro_readRDS(inputs[2]);manifest<-neutro_read_csv('work/cohort/C3_input_manifest.csv');d$antidep<-0L;counts<-data.frame()
for(cy in 4:10){yr<-paste(1997+2*cy,1998+2*cy,sep='-');readmod<-function(mod){p<-manifest$path[manifest$cycle==yr&manifest$module==mod];stopifnot(length(p)==1);sources<<-rbind(sources,data.frame(path=p,md5=unname(neutro_md5(p))));neutro_read_xpt(p)};dict<-readmod('rxq_drug');rx<-readmod('rxq_rx');cn<-grep('^RXD.*CN',names(dict),value=TRUE);hit<-Reduce(`|`,lapply(cn,function(v)grepl('ANTIDEPRESSANT',toupper(as.character(dict[[v]])),fixed=TRUE)));ids<-unique(dict$RXDDRGID[hit]);users<-unique(rx$SEQN[rx$RXDDRGID%in%ids]);ix<-d$SDDSRVYR==cy;d$antidep[ix]<-as.integer(d$SEQN[ix]%in%users);counts<-rbind(counts,data.frame(cycle=yr,N=sum(ix),antidep_users=sum(d$antidep[ix])))}
csv(counts,'D2_1_antidep_counts.csv');p<-file.path(out,'analytic_v3.rds');stopifnot(!neutro_file_exists(p));neutro_saveRDS(d,p)
base<-ds$variables;ix<-match(base$SEQN,d$SEQN);base$D2_included<-!is.na(ix);for(v in c('antidep','alcohol4','log10_cotinine','LBXCOT','CRP_mgL','log_CRP','education','pir','diabetes','hypertension','any_nsaid','any_corticosteroid','any_statin','race_eth','sex','smoking'))base[[v]]<-d[[v]][ix];base$cycle_factor<-factor(base$SDDSRVYR);ev<-c('log10_cotinine','alcohol4','education','pir','diabetes','hypertension','any_nsaid','any_corticosteroid','any_statin');base$common<-base$D2_included&complete.cases(base[,ev]);base$PHQ15<-as.numeric(base$PHQ9_total>=15);base$cotinine3<-factor(ifelse(is.na(base$LBXCOT),NA,ifelse(base$LBXCOT<.015,'below_LOD',ifelse(base$LBXCOT<3,'LOD_to_lt3','ge3'))),levels=c('below_LOD','LOD_to_lt3','ge3'));ds$variables<-base;full<-subset(ds,D2_included);common<-subset(ds,common);stopifnot(nrow(common$variables)==29360)
covars<-'log_Neutrophil + log_Lymphocyte + log_Platelet + age + sex + bmi + smoking + cycle_factor + race_eth';cells<-c('log_Neutrophil','log_Lymphocyte','log_Platelet');fits<-list();registry<-cf<-ct<-data.frame();extra<-c(Model1='',Model2='+ log10_cotinine + alcohol4 + education + pir',Model3='+ log10_cotinine + alcohol4 + education + pir + diabetes + hypertension + any_nsaid + any_corticosteroid + any_statin')
fit<-function(id,domain,extra='',outcome='depression',pois=FALSE){f<-as.formula(paste(outcome,'~',covars,extra));ww<-character();neutro_cat('\nMODEL_START',id,'\n');print(f);m<-withCallingHandlers(svyglm(f,design=domain,family=if(pois)quasipoisson(link='log')else quasibinomial()),warning=function(w){ww<<-c(ww,conditionMessage(w));neutro_cat('WARNING',conditionMessage(w),'\n');invokeRestart('muffleWarning')});print(summary(m));ci<-confint(m);cs<-coef(summary(m));b<-coef(m)[cells];cf<<-rbind(cf,data.frame(model=id,N=nobs(m),term=cells,beta=unname(b),SE=unname(sqrt(diag(vcov(m)))[cells]),effect=unname(exp(b)),lower=unname(exp(ci[cells,1])),upper=unname(exp(ci[cells,2])),P=unname(cs[cells,ncol(cs)])))
 for(j in 2:3){v<-setNames(rep(0,length(coef(m))),names(coef(m)));v[cells[1]]<-1;v[cells[j]]<--1;c<-svycontrast(m,v);bb<-as.numeric(coef(c));se<-as.numeric(SE(c));r<-data.frame(model=id,contrast=paste(cells[1],cells[j],sep=' - '),beta=bb,SE=se,df=degf(domain),P=2*pt(-abs(bb/se),degf(domain)));ct<<-rbind(ct,r);print(r)}
 registry<<-rbind(registry,data.frame(model=id,N=nobs(m),formula=paste(deparse(f),collapse=' '),family=family(m)$family,design_df=degf(domain),residual_df=m$df.residual,converged=m$converged,warning_count=length(ww),warnings=paste(ww,collapse=' | ')));fits[[id]]<<-m;neutro_cat('MODEL_END',id,'\n')}

fit('D2_1_M3_antidep',full,'+ antidep');fit('D2_1_M3_antidep_excluded',subset(full,antidep==0));fit('D2_1_Model3_antidep',common,paste(extra['Model3'],'+ antidep'))

csv(cf,'D2_1_coefficients.csv');csv(ct,'D2_1_t_contrasts.csv');csv(registry,'D2_1_model_registry.csv');csv(sources,'D2_1_source_hashes.csv');stopifnot(nrow(registry)==3,all(sources$md5==neutro_md5(sources$path)));neutro_cat('D2_1_COMPLETE\n')
