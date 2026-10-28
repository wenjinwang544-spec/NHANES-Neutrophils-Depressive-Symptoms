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
neutro_domain_entry('scripts/nhanes/03_add_public_modules.R')
options(digits=17);suppressPackageStartupMessages({library(haven);library(dplyr);library(tidyr)})
out<-file.path(getwd(),'work/cohort');raw<-'data/nhanes';cycles<-c('2005-2006','2007-2008','2009-2010','2011-2012','2013-2014','2015-2016','2017-2018');suffix<-letters[4:10]
log_event<-function(...)neutro_cat(format(Sys.time(),'%Y-%m-%dT%H:%M:%S%z'),paste(...,collapse=' '),'\n',file=file.path(out,'C3_public_modules_run.log'),append=TRUE)
wn<-function(x,n){p<-file.path(out,n);stopifnot(!neutro_file_exists(p));neutro_write_csv(x,p,row.names=FALSE,na='NA');log_event('WRITE',n)}
sn<-function(x,n){p<-file.path(out,n);stopifnot(!neutro_file_exists(p));neutro_saveRDS(x,p);log_event('SAVE_RDS',n)}
meta<-list();manifest<-list();read_mod<-function(i,folder,mod,suff=TRUE){name<-if(suff)paste0(mod,'_',suffix[i],'.xpt')else paste0(mod,'.xpt');fs<-neutro_list_files(file.path(raw,cycles[i],folder),full.names=TRUE);f<-fs[tolower(basename(fs))==name];stopifnot(length(f)==1,!grepl('^P_',basename(f),ignore.case=TRUE));log_event('READ_XPT',f);d<-neutro_read_xpt(f);manifest[[length(manifest)+1]]<<-data.frame(cycle=cycles[i],module=mod,path=f,rows=nrow(d));for(v in names(d))meta[[length(meta)+1]]<<-data.frame(cycle=cycles[i],module=mod,variable=v,label=if(is.null(attr(d[[v]],'label')))NA_character_ else attr(d[[v]],'label'));d}
x<-neutro_readRDS(file.path(out,'C1_rebuilt_aligned_v3.rds'));base<-neutro_readRDS(file.path(out,'C2_adult_base.rds'));joined<-list();cotstats<-list();crpstats<-list();drugnames<-list();drugcounts<-list();drugperson<-list();alcreasons<-list();alcrows<-list();slvars<-list();rates<-list()
statrow<-function(d,i,pop){v<-d$CRP_mgL;data.frame(cycle=cycles[i],population=pop,N=nrow(d),nonmissing=sum(!is.na(v)),missing=sum(is.na(v)),missing_pct=mean(is.na(v))*100,median_mgL=if(any(!is.na(v)))median(v,na.rm=TRUE)else NA_real_,maximum_mgL=if(any(!is.na(v)))max(v,na.rm=TRUE)else NA_real_)}
for(i in seq_along(cycles)){
 z<-x[x$cycle==cycles[i],];stopifnot(!anyDuplicated(z$SEQN))
 cotmod<-if(i%in%c(2,3,4))'cotnal'else'cot';cot<-read_mod(i,'Laboratory',cotmod);stopifnot('LBXCOT'%in%names(cot));lc<-intersect(names(cot),'LBDCOTLC')
 if(length(lc)){below<-cot[[lc]]==1;method<-lc}else{below<-cot$LBXCOT==0.011;method<-'no LC variable; LBXCOT=0.011 per CDC codebook'}
 cot$COT_below_LOD<-ifelse(is.na(cot$LBXCOT),NA,below);cot$COT_LOD_ngmL<-0.015
 z<-left_join(z,cot[,c('SEQN','LBXCOT','COT_below_LOD','COT_LOD_ngmL')],by='SEQN')
 for(pop in c('raw_XPT','C1')){d<-if(pop=='raw_XPT')cot else z;den<-sum(!is.na(d$COT_below_LOD));cotstats[[length(cotstats)+1]]<-data.frame(cycle=cycles[i],population=pop,N=nrow(d),LBXCOT_nonmissing=sum(!is.na(d$LBXCOT)),missing_pct=mean(is.na(d$LBXCOT))*100,LOD_flag=method,LOD_ngmL=0.015,below_LOD=sum(d$COT_below_LOD,na.rm=TRUE),flag_denominator=den,below_LOD_pct=100*sum(d$COT_below_LOD,na.rm=TRUE)/den)}
 if(i%in%c(1,2,3,6,7)){
  mod<-if(i<=3)'crp'else'hscrp';crp<-read_mod(i,'Laboratory',mod);v<-if(i<=3)'LBXCRP'else'LBXHSCRP';label<-attr(crp[[v]],'label');stopifnot(grepl(if(i<=3)'mg/dL'else'mg/L',label,fixed=TRUE));crp$CRP_mgL<-as.numeric(crp[[v]])*if(i<=3)10 else 1
  z<-left_join(z,crp[,c('SEQN','CRP_mgL')],by='SEQN');crpstats[[length(crpstats)+1]]<-statrow(crp,i,'raw_XPT')
 }else z$CRP_mgL<-NA_real_
 crpstats[[length(crpstats)+1]]<-statrow(z,i,'C1')
 slq<-read_mod(i,'Questionnaire','slq');duration<-intersect(c('SLD010H','SLD012','SLD013'),names(slq));clock<-intersect(c('SLQ300','SLQ310','SLQ320','SLQ330'),names(slq));sleep<-c(duration,clock)
 for(v in sleep)slvars[[length(slvars)+1]]<-data.frame(cycle=cycles[i],variable=v,label=attr(slq[[v]],'label'),unit=if(v%in%duration)'hours'else'HH:MM',raw_special_77_99=if(is.numeric(slq[[v]]))sum(slq[[v]]%in%c(77,99))else sum(slq[[v]]%in%c('77777','99999','77:77','99:99')))
 z<-left_join(z,slq[,c('SEQN',sleep)],by='SEQN')
 alq<-read_mod(i,'Questionnaire','alq');alcoholvars<-setdiff(names(alq),'SEQN');for(v in alcoholvars)if(v%in%names(z))z[[v]]<-NULL;z<-left_join(z,alq,by='SEQN')
 missingz<-z[is.na(z$alcohol_drinks),];get<-function(v)if(v%in%names(missingz))as.numeric(missingz[[v]])else rep(NA_real_,nrow(missingz))
 # mutually exclusive reporting categories; no recoding applied.
 never<-get('ALQ111')==2;pastnone<-if(i==7)get('ALQ121')==0 else get('ALQ120Q')==0;low_lifetime<-get('ALQ110')==2
 relevant<-intersect(c('ALQ101','ALQ110','ALQ111','ALQ120Q','ALQ120U','ALQ121','ALQ130'),names(missingz));refused<-rep(FALSE,nrow(missingz));unknown<-refused
 for(v in relevant){vals<-as.numeric(missingz[[v]]);refcodes<-if(v%in%c('ALQ101','ALQ110','ALQ111','ALQ120U'))7 else if(v=='ALQ121')77 else 777;unkcodes<-if(v%in%c('ALQ101','ALQ110','ALQ111','ALQ120U'))9 else if(v=='ALQ121')99 else 999;refused<-refused|vals%in%refcodes;unknown<-unknown|vals%in%unkcodes}
 never[is.na(never)]<-FALSE;pastnone[is.na(pastnone)]<-FALSE;low_lifetime[is.na(low_lifetime)]<-FALSE
 reason<-case_when(never~'never_drink_explicit_ALQ111_2',pastnone~'no_drink_past_year_frequency_0',low_lifetime~'lifetime_less_than12_skip_not_proven_never',refused&unknown~'refused_and_dontknow',refused~'refused',unknown~'dontknow',i<=3&missingz$age<20~'age18_19_not_in_exact_ALQ_file',TRUE~'other_missing_or_nonresponse')
 alcrows[[i]]<-data.frame(cycle=cycles[i],SEQN=missingz$SEQN,age=missingz$age,reason=reason,ALQ101=get('ALQ101'),ALQ110=get('ALQ110'),ALQ111=get('ALQ111'),ALQ120Q=get('ALQ120Q'),ALQ121=get('ALQ121'),ALQ130=get('ALQ130'))
 counts<-as.data.frame(table(reason));counts$cycle<-cycles[i];names(counts)[1:2]<-c('reason','N');alcreasons[[i]]<-counts
 rx<-read_mod(i,'Questionnaire','rxq_rx');dict<-read_mod(i,'Questionnaire','rxq_drug',FALSE);cn<-grep('^RXD.*CN',names(dict),value=TRUE);hit<-Reduce(`|`,lapply(cn,function(v)grepl('ANTIDEPRESSANT',toupper(as.character(dict[[v]])),fixed=TRUE)));ids<-unique(dict$RXDDRGID[hit]);rx$candidate_antidepressant<-rx$RXDDRGID%in%ids
 matched<-rx[rx$candidate_antidepressant,];namesused<-sort(unique(as.character(matched$RXDDRUG)));drugnames[[i]]<-data.frame(cycle=cycles[i],drug_name=namesused)
 users<-unique(matched$SEQN);drugcounts[[i]]<-data.frame(cycle=cycles[i],raw_RXQ_users=length(users),adult_base_users=sum(base$SEQN[base$cycle==cycles[i]]%in%users),C1_users=sum(z$SEQN%in%users),C1_N=nrow(z),unmatched_nonmissing_drug_ID_records=sum(!is.na(rx$RXDDRGID)&!rx$RXDDRGID%in%dict$RXDDRGID),pending_confirmation=TRUE)
 drugperson[[i]]<-data.frame(cycle=cycles[i],SEQN=users)
 joined[[i]]<-z
}
merged<-bind_rows(joined);stopifnot(nrow(merged)==nrow(x),!anyDuplicated(merged$SEQN),setequal(merged$SEQN,x$SEQN))
newvars<-c(sprintf('DPQ%03d',seq(10,90,10)),'LBXCOT','COT_below_LOD','CRP_mgL','SLD010H','SLD012','SLD013','SLQ300','SLQ310','SLQ320','SLQ330',unique(unlist(lapply(joined,function(z)grep('^ALQ',names(z),value=TRUE)))))
for(i in seq_along(cycles)){z<-merged[merged$cycle==cycles[i],];for(v in newvars){u<-z[[v]];rates[[length(rates)+1]]<-data.frame(cycle=cycles[i],variable=v,N=nrow(z),missing_N=sum(is.na(u)),missing_pct=mean(is.na(u))*100,raw_special_77_99_N=if(is.numeric(u))sum(u%in%c(77,99))else sum(u%in%c('77777','99999','77:77','99:99')),definition='raw_NA except DPQ7/9 recoded; no sleep/alcohol recoding applied')}}
wn(bind_rows(manifest),'C3_input_manifest.csv');wn(bind_rows(meta),'C3_variable_labels_all.csv');wn(bind_rows(slvars),'C3_sleep_duration_clock_metadata.csv');wn(bind_rows(cotstats),'C3_cotinine_LOD_summary.csv');wn(bind_rows(crpstats),'C3_CRP_summary.csv');wn(bind_rows(drugnames),'C3_antidepressant_names_by_cycle_pending.csv');wn(data.frame(drug_name=sort(unique(unlist(lapply(drugnames,function(x)x$drug_name))))),'C3_antidepressant_unique_names_pending.csv');wn(bind_rows(drugcounts),'C3_antidepressant_users_pending.csv');wn(bind_rows(drugperson),'C3_antidepressant_candidate_users_pending.csv');wn(bind_rows(alcreasons),'C3_alcohol_missing_reason_counts.csv');wn(bind_rows(alcrows),'C3_alcohol_missing_participant_details.csv');wn(bind_rows(rates),'C4_new_variable_missingness_by_cycle.csv');wn(data.frame(dataset=c('C1','C1_plus_confirmed_new_variables'),N=c(nrow(x),nrow(merged)),unique_SEQN=c(length(unique(x$SEQN)),length(unique(merged$SEQN)))),'C4_merge_N.csv');sn(merged,'C3_C1_new_variables_LOD_confirmed_pending_alcohol_antidepressants.rds')
log_event('C3_C4_FINISHED','no new models','alcohol not recoded','antidepressants not merged','COTNAL_E LOD pending')
print(bind_rows(cotstats));print(bind_rows(crpstats));print(bind_rows(drugcounts));print(bind_rows(alcreasons));print(sort(unique(unlist(lapply(drugnames,function(x)x$drug_name)))));print(nrow(merged));
