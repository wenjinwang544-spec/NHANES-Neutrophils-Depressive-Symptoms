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
neutro_domain_entry('scripts/nhanes/29_run_G3_remaining.R')
# Derived candidate; see provenance/script_selection.csv. Invoke via run_job.py.
options(warn=1,digits=17)
suppressPackageStartupMessages(library(survey))
out<-'work/shapley';SDfixed<-0.381754443
csv<-function(x,n){p<-file.path(out,n);stopifnot(!neutro_file_exists(p));neutro_write_csv(x,p,row.names=FALSE,na='NA')}
neutro_cat('START',as.character(Sys.time()),'\n');print(sessionInfo());
inputs<-c('work/analytic/analytic_v3.rds','work/cohort/C2_full_adult_design.rds');hash<-neutro_md5(inputs)
d<-as.data.frame(neutro_readRDS(inputs[1]));des<-neutro_readRDS(inputs[2]);b<-des$variables;ix<-match(b$SEQN,d$SEQN);b$included<-!is.na(ix)
for(v in setdiff(names(d),c('SEQN','SDMVPSU','SDMVSTRA','WTMEC_adj')))b[[v]]<-d[[v]][ix]
b$cycle_factor<-factor(b$SDDSRVYR)
groups<-list(A='log10_cotinine',B='alcohol4',C=c('education','pir'),D=c('diabetes','hypertension'),E=c('any_nsaid','any_corticosteroid','any_statin'))
b$cc<-b$included&complete.cases(b[,unlist(groups)])
b$ln_NLR<-log(b$Neutrophil_count/b$Lymphocyte_count);b$ln_PLR<-log(b$Platelet_count/b$Lymphocyte_count);b$ln_SII<-log((b$Neutrophil_count/1000)*(b$Platelet_count/1000)/(b$Lymphocyte_count/1000))
des$variables<-b;full<-subset(des,included);cc<-subset(des,cc);stopifnot(nrow(full$variables)==33174,nrow(cc$variables)==29360)
cov1<-c('age','sex','bmi','smoking','cycle_factor','race_eth');cells<-c('log_Neutrophil','log_Lymphocyte','log_Platelet')
sets<-lapply(0:31,function(mask)which(as.logical(intToBits(mask)[1:5])))
setname<-function(i)if(length(sets[[i]]))paste(names(groups)[sets[[i]]],collapse='')else 'none'
formula_for<-function(i,joint=TRUE,cycle=TRUE){v<-c(if(joint)cells else cells[1],if(cycle)cov1 else setdiff(cov1,'cycle_factor'),unlist(groups[sets[[i]]]));reformulate(v,'depression')}
registry<-list();fits<-list();counter<-0
fit<-function(id,f,dd){neutro_cat('MODEL_START',id,'N',nrow(dd$variables),'\n');print(f);ww<-character();err<-''
 m<-tryCatch(withCallingHandlers(svyglm(f,design=dd,family=quasibinomial()),warning=function(w){ww<<-c(ww,conditionMessage(w));neutro_cat('WARNING',id,conditionMessage(w),'\n');invokeRestart('muffleWarning')}),error=function(e){err<<-conditionMessage(e);neutro_cat('ERROR',err,'\n');NULL})
 counter<<-counter+1;registry[[counter]]<<-data.frame(id=id,replicate=0,N=nrow(dd$variables),formula=paste(deparse(f),collapse=' '),converged=if(is.null(m))FALSE else m$converged,warnings=paste(ww,collapse=' | '),error=err)
 if(!is.null(m))fits[[id]]<<-m;neutro_cat('MODEL_END',id,'\n');m}
effect<-function(m,delta=SDfixed){if(is.null(m))return(data.frame(beta=NA,SE=NA,ln_OR_SD=NA,OR=NA,lower=NA,upper=NA,P=NA,df=NA));beta<-coef(m)['log_Neutrophil'];se<-sqrt(vcov(m)['log_Neutrophil','log_Neutrophil']);ci<-confint(m)['log_Neutrophil',];data.frame(beta=unname(beta),SE=unname(se),ln_OR_SD=unname(beta*delta),OR=unname(exp(beta*delta)),lower=unname(exp(ci[1]*delta)),upper=unname(exp(ci[2]*delta)),P=coef(summary(m))['log_Neutrophil',4],df=m$df.residual)}
neutro_cat('Coefficient CI/P: t with model residual df, consistent with preceding D2/G1. Joint Wald: denominator degf(domain design). Bootstrap ordinary paired within-stratum PSU resampling from full base, then subset; percentile CI; 200 repetitions; seed 20261014.\n')
# Model1 full sample constraints, F(2,degf(design)).
m<-fit('G3_4_Model1_joint_full',formula_for(1),full);cn<-names(coef(m));constraint<-function(name,rows){L<-matrix(0,2,length(cn),dimnames=list(NULL,cn));for(j in 1:2)L[j,names(rows[[j]])]<-rows[[j]];theta<-as.vector(L%*%coef(m));V<-L%*%vcov(m)%*%t(L);F<-as.numeric(t(theta)%*%solve(V,theta)/2);df<-degf(full);data.frame(index=name,F=F,df1=2,df2=df,P=pf(F,2,df,lower.tail=FALSE))}
wt<-rbind(constraint('NLR',list(c(log_Neutrophil=1,log_Lymphocyte=1),c(log_Platelet=1))),constraint('PLR',list(c(log_Platelet=1,log_Lymphocyte=1),c(log_Neutrophil=1))))
indices<-data.frame();for(v in c('ln_NLR','ln_PLR','ln_SII')){sd<-sqrt(as.numeric(coef(svyvar(reformulate(v),full))));mi<-fit(paste0('G3_4_',v),reformulate(c(v,cov1),'depression'),full);beta<-coef(mi)[v];se<-sqrt(vcov(mi)[v,v]);ci<-confint(mi)[v,];indices<-rbind(indices,data.frame(index=v,N=nrow(full$variables),weighted_SD=sd,beta=unname(beta),SE=unname(se),OR_SD=unname(exp(beta*sd)),lower=unname(exp(ci[1]*sd)),upper=unname(exp(ci[2]*sd)),P=coef(summary(mi))[v,4],df=mi$df.residual))}
csv(wt,'G3_4_constraints.csv');csv(indices,'G3_4_indices.csv');csv(do.call(rbind,registry),'G3_remaining_model_registry.csv');neutro_saveRDS(fits,file.path(out,'G3_constraint_models.rds'));print(wt);print(indices);neutro_cat('G3_4_COMPLETE\n')
