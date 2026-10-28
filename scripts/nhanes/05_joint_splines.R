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
neutro_domain_entry('scripts/nhanes/05_joint_splines.R')
out<-file.path(getwd(),'work/spline');
suppressPackageStartupMessages({library(survey);library(Hmisc)})
stopifnot(!neutro_file_exists(file.path(out,'candidate_curves.csv')));neutro_capture_output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
d<-neutro_readRDS('work/analytic/analytic_v3.rds');ds<-neutro_readRDS('work/cohort/C2_full_adult_design.rds');b<-ds$variables;ix<-match(b$SEQN,d$SEQN)
for(v in setdiff(names(d),c('SEQN','SDMVPSU','SDMVSTRA','WTMEC_adj','WTMEC2YR')))b[[v]]<-d[[v]][ix]
b$included<-!is.na(ix);b$sex<-factor(b$sex,levels=c(1,2),labels=c('Male','Female'));b$race_eth<-factor(b$race_eth);b$education<-factor(b$education,levels=c('Less than high school','High school/GED','Some college','College or above'));b$cycle_factor<-factor(b$SDDSRVYR);b$PHQ10<-as.numeric(b$PHQ9_total>=10)
for(v in c('Neutrophil','Lymphocyte','Platelet'))b[[paste0('log_',v)]]<-log(b[[paste0(v,'_count')]])
b$cc2<-b$included&complete.cases(b[,c('log10_cotinine','alcohol4','education','pir')]);b$cc3<-b$included&complete.cases(b[,c('log10_cotinine','alcohol4','education','pir','diabetes','hypertension','any_nsaid','any_corticosteroid','any_statin')]);ds$variables<-b
knots<-neutro_read_csv(file.path(repo,'supplementary/validated/shared_knots.csv'))$knot_log
check<-as.numeric(coef(svyquantile(~log_Neutrophil,subset(ds,included),c(.05,.35,.65,.95),ci=FALSE,qrule='math')));stopifnot(max(abs(knots-check))<1e-12)
basis<-matrix(NA_real_,nrow(b),3);basis[b$included,]<-rcspline.eval(b$log_Neutrophil[b$included],knots=knots,inclx=TRUE);b[,c('rcs_linear','rcs_nl1','rcs_nl2')]<-basis;ds$variables<-b
grid<-sort(unique(neutro_read_csv(file.path(repo,'supplementary/validated/candidate_curves.csv'))$Neutrophil_cells_uL));terms<-c('rcs_linear','rcs_nl1','rcs_nl2');B<-rcspline.eval(log(grid),knots=knots,inclx=TRUE);B0<-as.numeric(rcspline.eval(log(4000),knots=knots,inclx=TRUE));delta<-sweep(B,2,B0,'-')
ids<-c('Model1_full','Model1_common','Model2_common','Model3_complete_case');formbase<-'rcs_linear + rcs_nl1 + rcs_nl2 + log_Lymphocyte + log_Platelet + age + sex + bmi + smoking + cycle_factor + race_eth';extra<-c('','','+ log10_cotinine + alcohol4 + education + pir','+ log10_cotinine + alcohol4 + education + pir + diabetes + hypertension + any_nsaid + any_corticosteroid + any_statin')
curves<-tests<-reg<-data.frame();fits<-list()
for(j in 1:4)for(y in c('PHQ10','PHQ9_total')){
 des<-if(j==1)subset(ds,included)else if(j%in%2:3)subset(ds,cc2)else subset(ds,cc3);f<-as.formula(paste(y,'~',formbase,extra[j]));warn<-character();neutro_cat('FIT',ids[j],y,nrow(des$variables),'\n')
 m<-withCallingHandlers(svyglm(f,design=des,family=if(y=='PHQ10')quasibinomial()else quasipoisson(link='log')),warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart('muffleWarning')});stopifnot(m$converged);wt<-regTermTest(m,~rcs_nl1+rcs_nl2,method='Wald',df=degf(des));z<-as.numeric(delta%*%coef(m)[terms]);V<-vcov(m)[terms,terms];se<-sqrt(pmax(0,rowSums((delta%*%V)*delta)));q<-qt(.975,degf(des))
 curves<-rbind(curves,data.frame(model=ids[j],outcome=y,N=nobs(m),Neutrophil_cells_uL=grid,Neutrophil_10e9_L=grid/1000,reference_cells_uL=4000,link_difference=z,SE=se,effect=exp(z),lower=exp(z-q*se),upper=exp(z+q*se),df=degf(des)))
 tests<-rbind(tests,data.frame(model=ids[j],outcome=y,N=nobs(m),Wald_F=as.numeric(wt$Ftest),df1=wt$df,df2=degf(des),P_nonlinear=as.numeric(wt$p)))
 reg<-rbind(reg,data.frame(model=ids[j],outcome=y,N=nobs(m),formula=paste(deparse(f),collapse=' '),design_df=degf(des),residual_df=m$df.residual,converged=m$converged,warnings=paste(warn,collapse=' | ')));fits[[paste(ids[j],y,sep='_')]]<-m
}
neutro_write_csv(curves,file.path(out,'candidate_curves.csv'),row.names=FALSE);neutro_write_csv(tests,file.path(out,'candidate_nonlinearity.csv'),row.names=FALSE);neutro_write_csv(reg,file.path(out,'candidate_model_registry.csv'),row.names=FALSE);neutro_saveRDS(fits,file.path(out,'candidate_models.rds'))
