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
neutro_domain_entry('scripts/nhanes/22_export_remaining.R')
# Derived candidate; see provenance/script_selection.csv. Invoke via run_job.py.
out<-neutro_normalize('work/export');csv<-function(x,n){p<-file.path(out,n);stopifnot(!neutro_file_exists(p));neutro_write_csv(x,p,row.names=FALSE);neutro_cat('EXPORT',n,nrow(x),'\n')};prov<-function(x,s,f,d){data.frame(source_script_and_output=paste(s,f,paste0('export=',d),sep=' | '),x,check.names=FALSE)}
for(kind in c('coefficients','contrasts')){
 src<-paste0('work/clinical/D1_',if(kind=='coefficients')'coefficients'else 't_contrasts','.csv');x<-neutro_read_csv(src)
 for(group in c('CRP','age20')){sel<-if(group=='CRP')grepl('^CRP_',x$model)else grepl('^adult_',x$model);z<-x[sel,];stopifnot(nrow(z)>0);dest<-paste0('Supplement_',group,'_',kind,'.csv');csv(prov(z,'work/clinical/run_D1.R',src,dest),dest)}
 src<-paste0('work/analytic/D2_',if(kind=='coefficients')'coefficients'else 't_contrasts','.csv');x<-neutro_read_csv(src);z<-subset(x,grepl('D2_3d',model));stopifnot(nrow(z)>0);dest<-paste0('Supplement_cotinine_categories_',kind,'.csv');csv(prov(z,'work/analytic/run_D2_3.R',src,dest),dest)
}
src<-'work/export/Figure1_flow_validation.csv';csv(prov(neutro_read_csv(src),'work/export/check_flow.R',src,'Figure1_plot_data.csv'),'Figure1_plot_data.csv')
