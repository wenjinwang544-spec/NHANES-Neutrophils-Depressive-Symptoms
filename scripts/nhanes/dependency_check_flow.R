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
neutro_domain_entry('scripts/nhanes/dependency_check_flow.R')
options(warn=1);out<-neutro_normalize('work/export');csv<-function(x,n){p<-file.path(out,n);stopifnot(!neutro_file_exists(p));neutro_write_csv(x,p,row.names=FALSE)}
input<-'work/cohort/C2_full_adult_design.rds';h<-neutro_md5(input);b<-neutro_readRDS(input)$variables
q<-as.matrix(b[paste0('DPQ0',1:9,'0')]);q[!q%in%0:3]<-NA;phq<-complete.cases(q)
cbc<-complete.cases(b[,c('LBDNENO','LBDLYMNO','LBXPLTSI')])&b$LBDNENO>0&b$LBDLYMNO>0&b$LBXPLTSI>0;cbc[is.na(cbc)]<-FALSE
smoking_ok<-b$SMQ020%in%2|(b$SMQ020%in%1&b$SMQ040%in%1:3);bmi_ok<-!is.na(b$BMXBMI);stage<-phq&cbc
inc<-stage&smoking_ok&bmi_ok;early<-stage&!smoking_ok&b$RIDAGEYR%in%c(18,19)&b$SDDSRVYR%in%4:7
check<-data.frame(step=c('age_ge18_base','PHQ9_incomplete','PHQ9_complete','CBC_missing_or_nonpositive','CBC_complete_positive','BMI_only_missing','smoking_only_missing','BMI_and_smoking_missing','early_18_19_smoking_missing','final_included'),computed=c(nrow(b),sum(!phq),sum(phq),sum(phq&!cbc),sum(stage),sum(stage&!bmi_ok&smoking_ok),sum(stage&bmi_ok&!smoking_ok),sum(stage&!bmi_ok&!smoking_ok),sum(early),sum(inc)),requested=c(42143,5884,36259,1511,34748,316,1234,24,1240,33174));check$difference<-check$computed-check$requested;csv(check,'Figure1_flow_validation.csv');print(check)
stopifnot(all(check$difference==0),all(b$SEQN[inc]%in%neutro_readRDS('work/analytic/analytic_v3.rds')$SEQN),h==neutro_md5(input))
csv(data.frame(SEQN=b$SEQN,PHQ9_complete=phq,CBC_complete_positive=cbc,BMI_complete=bmi_ok,smoking_defined=smoking_ok,included=inc,early_18_19_smoking_missing=early),'Figure1_participant_flags.csv')
csv(data.frame(path=input,md5=unname(h)),'flow_input_hash.csv');neutro_cat('FLOW_GATE_PASS\n')
