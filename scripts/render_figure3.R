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
neutro_domain_entry('scripts/render_figure3.R')
# Render only: uses validated Stage 3 coordinates; no model fitting.
a<-commandArgs(TRUE);input<-file.path(.neutro$src,'supplementary/validated');j<-match('--curve-input',a);if(!is.na(j)){stopifnot(j<length(a));input<-neutro_normalize(a[j+1],mustWork=TRUE)};out<-file.path(.neutro$out,'figures');neutro_dir_create(out,recursive=TRUE,showWarnings=FALSE)
curves<-neutro_read_csv(file.path(input,'candidate_curves.csv'));old<-curves
ids<-c('Model1_full','Model1_common','Model2_common','Model3_complete_case')
stopifnot(nrow(curves)==2456,all(curves$reference_cells_uL==4000))
cols<-c('#0072B2','#777777','#E69F00','#009E73');ltys<-c(1,2,1,1);labels<-c('Model 1, full (N=33,174)','Model 1, complete cases (N=29,415)','Model 2, complete cases (N=29,415)','Model 3, complete cases (N=29,360)');rug<-neutro_read_csv(file.path(input,'Figure3_distribution.csv'))$count_10e9_L
plotfig<-function(compare=FALSE){par(mfrow=c(1,2),mar=c(4.7,4.7,3.3,1),family='Helvetica');for(y in c('PHQ10','PHQ9_total')){z<-subset(curves,outcome==y);zz<-subset(old,outcome==y);yr<-range(z$lower,z$upper,if(compare)c(zz$lower,zz$upper));plot(NA,xlim=range(z$Neutrophil_10e9_L),ylim=yr,xlab=expression('Neutrophils ('*10^9*'/L)'),ylab=if(y=='PHQ10')'OR versus 4.0'else'RR versus 4.0',main=if(y=='PHQ10')expression('PHQ-9 ' >= 10)else'PHQ-9 total score',bty='l');abline(h=1,lty=2,col='grey55');abline(v=4,lty=3,col='grey65');for(j in c(1,4)){q<-z[z$model==ids[j],];polygon(c(q$Neutrophil_10e9_L,rev(q$Neutrophil_10e9_L)),c(q$lower,rev(q$upper)),col=adjustcolor(cols[j],alpha.f=.12),border=NA)};for(j in 1:4){q<-z[z$model==ids[j],];lines(q$Neutrophil_10e9_L,q$effect,col=cols[j],lty=ltys[j],lwd=2);if(compare){q<-zz[zz$model==ids[j],];lines(q$Neutrophil_10e9_L,q$effect,col=adjustcolor(cols[j],alpha.f=.7),lty=3,lwd=1.5)}};rug(rug,side=1,ticksize=.03);legend('topleft',legend=labels,col=cols,lty=ltys,lwd=2,bty='n',cex=.62);if(compare)mtext('Solid/dashed: candidate; dotted: archived',side=3,line=.3,cex=.7)}}

for(ext in c('png','pdf')){f<-file.path(out,paste0('Figure3_candidate_jointcounts.',ext));if(neutro_file_exists(f))stop('Refusing to overwrite output: ',f);if(ext=='png')neutro_png(f,width=3000,height=1800,res=300)else neutro_pdf(f,width=10,height=6,family='Helvetica');plotfig(FALSE);neutro_dev_off()}
