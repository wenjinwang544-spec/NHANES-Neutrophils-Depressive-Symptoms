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
neutro_domain_entry('scripts/nhanes/32_run_H10_v3.R')
# Derived candidate; see provenance/script_selection.csv. Invoke via run_job.py.
options(warn=1,digits=17)
out <- 'work/attenuation'
stopifnot(dir.exists(out), !neutro_file_exists(file.path(out,'H10_attenuation_summary.csv')))
inputs <- c('work/analytic/analytic_v3.rds','work/cohort/C2_full_adult_design.rds','work/shapley/G3_1_subset_models.csv','work/shapley/G3_1_bootstrap_betas.csv','work/shapley/G3_bootstrap_model_registry.csv')
sha256 <- function(p) system2('shasum', c('-a','256', shQuote(p)), stdout=TRUE)
meta_before <- neutro_file_info(inputs); hash_before <- vapply(inputs,function(p)sub(' .*','',sha256(p)),character(1))
logfile <- file.path(out,'H10_run.log'); logcon <- neutro_file(logfile,open='a'); neutro_sink(logcon,split=TRUE); neutro_sink(logcon,type='message')
neutro_cat('Started:',format(Sys.time(),tz='UTC'),'UTC\nR:',R.version.string,'\nSeed: 20261014\nBootstrap repetitions: 1000\n')
neutro_cat('Method: paired PSU resampling within SDMVSTRA, same PSU multiplier applied to all models.\n')
neutro_cat('Complete-case sample: G3 cc definition; expected N=29360.\n')
neutro_cat('Model 1 = cells + age + sex + bmi + smoking + cycle factor + race_eth.\nModel 2 = Model 1 + log10_cotinine + alcohol4 + education + pir.\nModel 3 = Model 2 + diabetes + hypertension + any_nsaid + any_corticosteroid + any_statin.\n')
neutro_cat('Each replicate fits only these three models using glm.fit and the G3 survey design matrices; no Shapley subset fits.\n')
suppressPackageStartupMessages(library(survey))
d <- as.data.frame(neutro_readRDS(inputs[1])); des <- neutro_readRDS(inputs[2]); b <- des$variables
ix <- match(b$SEQN,d$SEQN); b$included <- !is.na(ix)
for(v in setdiff(names(d),c('SEQN','SDMVPSU','SDMVSTRA','WTMEC_adj'))) b[[v]] <- d[[v]][ix]
b$cycle_factor <- factor(b$SDDSRVYR)
groups <- list(A='log10_cotinine',B='alcohol4',C=c('education','pir'),D=c('diabetes','hypertension'),E=c('any_nsaid','any_corticosteroid','any_statin'))
b$cc <- b$included & complete.cases(b[,unlist(groups)])
des$variables <- b; full <- subset(des,included); cc <- subset(des,cc)
stopifnot(nrow(full$variables)==33174,nrow(cc$variables)==29360)
cov1 <- c('age','sex','bmi','smoking','cycle_factor','race_eth'); cells <- c('log_Neutrophil','log_Lymphocyte','log_Platelet')
forms <- list('Model 1'=reformulate(c(cells,cov1),'depression'),'Model 2'=reformulate(c(cells,cov1,unlist(groups[c('A','B','C')])),'depression'),'Model 3'=reformulate(c(cells,cov1,unlist(groups)),'depression'))
model_m <- c('Model 1'=0L,'Model 2'=7L,'Model 3'=31L)
X <- lapply(forms,function(f)model.matrix(f,data=cc$variables)); y <- cc$variables$depression; w <- weights(cc)
pointtab <- neutro_read_csv(inputs[3],check.names=FALSE); point <- setNames(pointtab$beta[match(model_m,pointtab$mask)],names(model_m))
stopifnot(all(is.finite(point))); neutro_cat('Original sample coefficients from saved G3 subset output:\n'); print(point)
psu <- interaction(b$SDMVSTRA,b$SDMVPSU,drop=TRUE); pid <- as.integer(psu)
pairs <- unique(data.frame(pid=pid,stratum=b$SDMVSTRA)); strata <- split(pairs$pid,pairs$stratum)
set.seed(20261014); mult <- matrix(0,nlevels(psu),1000)
for(r in 1:1000) for(ids in strata) mult[ids,r] <- tabulate(match(sample(ids,length(ids),replace=TRUE),ids),length(ids))
take <- match(cc$variables$SEQN,b$SEQN); stopifnot(!anyNA(take),length(take)==29360)
boot <- matrix(NA_real_,1000,3,dimnames=list(NULL,names(forms))); conv <- matrix(FALSE,1000,3,dimnames=list(NULL,names(forms))); warns <- matrix('',1000,3,dimnames=list(NULL,names(forms))); errs <- matrix('',1000,3,dimnames=list(NULL,names(forms)))
for(r in 1:1000){wr<-w*mult[pid[take],r];keep<-wr>0
 for(j in seq_along(forms)){ww<-character();er<-'';fit<-tryCatch(withCallingHandlers(glm.fit(x=X[[j]][keep,,drop=FALSE],y=y[keep],weights=wr[keep]/mean(wr[keep]),family=quasibinomial()),warning=function(e){ww<<-c(ww,conditionMessage(e));invokeRestart('muffleWarning')}),error=function(e){er<<-conditionMessage(e);NULL});warns[r,j]<-paste(unique(ww),collapse=' | ');errs[r,j]<-er
  if(!is.null(fit)){val<-fit$coefficients['log_Neutrophil'];boot[r,j]<-unname(val);conv[r,j]<-isTRUE(fit$converged)&&is.finite(val)}}
 if(r%%100==0)neutro_cat('Completed replicate',r,'/1000 at',format(Sys.time(),tz='UTC'),'\n')
}
valid<-apply(conv,1,all)&apply(is.finite(boot),1,all);neutro_cat('Valid paired replicates:',sum(valid),'/1000\n');neutro_cat('Model-level converged finite fits:\n');print(colSums(conv));neutro_cat('Warning counts:\n');print(colSums(matrix(nzchar(warns),nrow=1000,ncol=3)));neutro_cat('Error counts:\n');print(colSums(matrix(nzchar(errs),nrow=1000,ncol=3)));if(!any(valid))stop('No successful paired replicates')
q7<-function(z)as.numeric(quantile(z,c(.025,.975),na.rm=TRUE,type=7,names=FALSE))
summ<-list();distributions<-list()
for(comp in c('Model 2','Model 3')){dp<-point['Model 1']-point[comp];ap<-100*dp/point['Model 1'];db<-boot[valid,'Model 1']-boot[valid,comp];ab<-100*db/boot[valid,'Model 1'];ci_d<-q7(db);ci_a<-q7(ab)
 summ[[comp]]<-data.frame(comparison=paste('Model 1 minus',comp),N=29360,beta_Model1=unname(point['Model 1']),beta_comparator=unname(point[comp]),coefficient_difference_point=unname(dp),coefficient_difference_q2_5=ci_d[1],coefficient_difference_q97_5=ci_d[2],attenuation_percent_point=unname(ap),attenuation_percent_q2_5=ci_a[1],attenuation_percent_q97_5=ci_a[2],valid_replicates=sum(valid))
 distributions[[comp]]<-data.frame(comparison=paste('Model 1 minus',comp),minimum=min(ab),q2_5=quantile(ab,.025,type=7,names=FALSE),q25=quantile(ab,.25,type=7,names=FALSE),median=quantile(ab,.5,type=7,names=FALSE),q75=quantile(ab,.75,type=7,names=FALSE),q97_5=quantile(ab,.975,type=7,names=FALSE),maximum=max(ab),Model1_beta_le_zero=sum(boot[valid,'Model 1']<=0),attenuation_gt_100=sum(ab>100),attenuation_lt_0=sum(ab<0),valid_replicates=sum(valid))}
summary<-do.call(rbind,summ);dist<-do.call(rbind,distributions)
neutro_write_csv(summary,file.path(out,'H10_attenuation_summary.csv'),row.names=FALSE);neutro_write_csv(dist,file.path(out,'H10_attenuation_distribution.csv'),row.names=FALSE)
bootout<-data.frame(replicate=1:1000,Model1_beta=boot[,1],Model2_beta=boot[,2],Model3_beta=boot[,3],all_three_success=valid,Model1_converged=conv[,1],Model2_converged=conv[,2],Model3_converged=conv[,3],Model1_warnings=warns[,1],Model2_warnings=warns[,2],Model3_warnings=warns[,3],Model1_error=errs[,1],Model2_error=errs[,2],Model3_error=errs[,3])
neutro_write_csv(bootout,file.path(out,'H10_bootstrap_coefficients.csv'),row.names=FALSE)
meta_after<-neutro_file_info(inputs);hash_after<-vapply(inputs,function(p)sub(' .*','',sha256(p)),character(1));integrity<-data.frame(path=inputs,size_before=meta_before$size,size_after=meta_after$size,mtime_before=as.character(meta_before$mtime),mtime_after=as.character(meta_after$mtime),sha256_before=hash_before,sha256_after=hash_after,unchanged=meta_before$size==meta_after$size&meta_before$mtime==meta_after$mtime&hash_before==hash_after);neutro_write_csv(integrity,file.path(out,'H10_input_integrity.csv'),row.names=FALSE);stopifnot(all(integrity$unchanged))
neutro_cat('Input files unchanged:',nrow(integrity),'/',nrow(integrity),'\nSummary:\n');print(summary,row.names=FALSE);neutro_cat('Distribution diagnostics:\n');print(dist,row.names=FALSE);neutro_cat('Completed:',format(Sys.time(),tz='UTC'),'UTC\n');neutro_sink(type='message');neutro_sink();close(logcon)
