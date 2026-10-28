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
neutro_domain_entry('scripts/nhanes/19_draw_figures.R')
# Derived candidate; see provenance/script_selection.csv. Invoke via run_job.py.
out<-neutro_normalize('work/export');neutro_cat('START_FIGURES\n');cols<-c('#176B91','#C45D35');new_device<-function(name,type,w,h){p<-file.path(out,paste0(name,'.',type));stopifnot(!neutro_file_exists(p));if(type=='png')neutro_png(p,width=w,height=h,units='in',res=300,type='quartz',bg='white')else neutro_pdf(p,width=w,height=h,family='Helvetica',useDingbats=FALSE)}
flow<-neutro_read_csv(file.path(out,'Figure1_plot_data.csv'));n<-setNames(flow$computed,flow$step)
fig1<-function(){par(mar=c(.2,.2,.2,.2));plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));boxx<-function(x0,y0,x1,y1,label,fill='#EDF5F8',size=.9){rect(x0,y0,x1,y1,col=fill,border='#507A8D',lwd=1.1);text((x0+x1)/2,(y0+y1)/2,label,cex=size)}
 ys<-c(.89,.68,.47,.20);h<-.09
 boxx(.12,ys[1],.67,ys[1]+h,paste0('NHANES 2005-2018\nAge >=18: N = ',format(n['age_ge18_base'],big.mark=',')))
 boxx(.12,ys[2],.67,ys[2]+h,paste0('Complete PHQ-9\nN = ',format(n['PHQ9_complete'],big.mark=',')))
 boxx(.12,ys[3],.67,ys[3]+h,paste0('Three cell counts complete and >0\nN = ',format(n['CBC_complete_positive'],big.mark=',')))
 boxx(.12,ys[4],.67,ys[4]+h,paste0('BMI and smoking available\nAnalysis N = ',format(n['final_included'],big.mark=',')),fill='#DCEEF5')
 for(i in 1:3)arrows(.395,ys[i],.395,ys[i+1]+h,length=.08,lwd=1.1)
 boxx(.73,.79,.99,.88,paste0('PHQ-9 incomplete\nn = ',format(n['PHQ9_incomplete'],big.mark=',')),fill='#FCF3EA',size=.78);arrows(.395,.82,.72,.82,length=.06)
 boxx(.73,.58,.99,.67,paste0('Cell counts missing\nor nonpositive\nn = ',format(n['CBC_missing_or_nonpositive'],big.mark=',')),fill='#FCF3EA',size=.73);arrows(.395,.61,.72,.61,length=.06)
 boxx(.70,.31,.99,.45,paste0('BMI only missing: ',n['BMI_only_missing'],'\nSmoking only missing: ',format(n['smoking_only_missing'],big.mark=','),'\nBoth missing: ',n['BMI_and_smoking_missing']),fill='#FCF3EA',size=.71);arrows(.395,.38,.69,.38,length=.06)
 text(.50,.105,paste0('2005-2012, age 18-19: ',format(n['early_18_19_smoking_missing'],big.mark=','),' excluded with smoking unavailable\nSMQ020 not applicable'),cex=.8)
}

for(t in c('png','pdf')){new_device('Figure1_flow',t,8,10);fig1();neutro_dev_off()}
