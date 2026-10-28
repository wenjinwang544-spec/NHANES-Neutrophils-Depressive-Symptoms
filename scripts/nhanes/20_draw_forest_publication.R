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
neutro_domain_entry('scripts/nhanes/20_draw_forest_publication.R')
# Derived candidate; see provenance/script_selection.csv. Invoke via run_job.py.
out<-neutro_normalize('work/export');d<-neutro_read_csv(file.path(out,'Figure2_plot_data.csv'));draw<-function(){layout(matrix(1:6,nrow=3,byrow=TRUE),widths=c(4,1.65));for(y in c('depression','PHQ9_total','PHQ15')){z<-d[d$outcome==y,];z<-z[match(c('Model1','Model2','Model3'),z$model),];par(mar=c(3.4,5,2.6,.5));plot(NA,xlim=range(c(.98,z$SD_lower_MI,z$SD_upper_MI))*c(.99,1.01),ylim=c(.5,3.5),xlab='',ylab='',yaxt='n',main=switch(y,depression='PHQ-9 >=10: OR per 1 SD',PHQ9_total='PHQ-9 total score: RR per 1 SD',PHQ15='PHQ-9 >=15: OR per 1 SD'),bty='l');axis(2,at=3:1,labels=c('Model 1','Model 2','Model 3'),las=1);abline(v=1,lty=2,col='grey55');segments(z$SD_lower_MI,3:1,z$SD_upper_MI,3:1,lwd=2,col='#176B91');points(z$SD_OR_MI,3:1,pch=15,cex=1.1,col='#176B91');par(mar=c(3.4,0,2.6,1));plot.new();plot.window(xlim=c(0,1),ylim=c(.5,3.5));text(.01,3:1,sprintf('%.3f (%.3f-%.3f)',z$SD_OR_MI,z$SD_lower_MI,z$SD_upper_MI),adj=0,cex=.84);text(.01,3.5,'Estimate (95% CI)',adj=0,cex=.8)}}
for(type in c('png','pdf')){p<-file.path(out,paste0('Figure2_forest_publication.',type));stopifnot(!neutro_file_exists(p));if(type=='png')neutro_png(p,width=10,height=8,units='in',res=300,type='quartz',bg='white')else neutro_pdf(p,width=10,height=8,family='Helvetica',useDingbats=FALSE);draw();neutro_dev_off()};neutro_cat('FOREST_FINAL_COMPLETE\n')
