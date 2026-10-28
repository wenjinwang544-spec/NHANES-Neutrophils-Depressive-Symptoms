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
neutro_domain_entry('scripts/nhanes/01_build_cohort.R')
options(stringsAsFactors=FALSE, digits=17)
suppressPackageStartupMessages({library(haven);library(dplyr);library(tidyr)})
project <- neutro_normalize(getwd()); out <- file.path(project,'work/cohort'); raw <- 'data/nhanes'
logfile<-file.path(out,'run_log.txt')
log_event<-function(...)neutro_cat(format(Sys.time(),'%Y-%m-%dT%H:%M:%S%z'),paste(...,collapse=' '),'\n',file=logfile,append=TRUE)
write_new<-function(x,name){p<-file.path(out,name);if(neutro_file_exists(p))stop('禁止覆盖: ',p);neutro_write_csv(x,p,row.names=FALSE,na='NA');log_event('WRITE',p,nrow(x))}
save_new<-function(x,name){p<-file.path(out,name);if(neutro_file_exists(p))stop('禁止覆盖: ',p);neutro_saveRDS(x,p);log_event('WRITE_RDS',p,nrow(x))}
cycles<-c('2005-2006','2007-2008','2009-2010','2011-2012','2013-2014','2015-2016','2017-2018'); suffix<-letters[4:10]
manifest<-list();labels<-list()
read_mod<-function(i,folder,mod,cols=NULL){
 expected<-paste0(mod,'_',suffix[i],'.xpt');dir<-file.path(raw,cycles[i],folder)
 f<-neutro_list_files(dir,full.names=TRUE);f<-f[tolower(basename(f))==expected]
 if(length(f)!=1)stop('精确模块文件匹配失败: ',expected)
 stopifnot(!grepl('^P_',basename(f),ignore.case=TRUE));log_event('READ_XPT',f)
 d<-neutro_read_xpt(f); manifest[[length(manifest)+1]]<<-data.frame(cycle=cycles[i],module=mod,path=f,rows=nrow(d),columns=ncol(d),size_bytes=neutro_file_info(f)$size)
 for(v in names(d))labels[[length(labels)+1]]<<-data.frame(cycle=cycles[i],module=mod,variable=v,label=if(is.null(attr(d[[v]],'label')))NA_character_ else attr(d[[v]],'label'))
 if(!is.null(cols)){d<-d[,intersect(names(d),c('SEQN',cols)),drop=FALSE];for(v in setdiff(cols,names(d)))d[[v]]<-NA_real_}
 if(mod!='rxq_rx' && anyDuplicated(d$SEQN))stop('SEQN 非唯一: ',f)
 d
}
nsaid<-c('ASPIRIN','IBUPROFEN','NAPROXEN','CELECOXIB','DICLOFENAC','MELOXICAM','INDOMETHACIN','KETOPROFEN','PIROXICAM','SULINDAC','ETODOLAC','KETOROLAC','FLURBIPROFEN','OXAPROZIN','TOLMETIN')
cort<-c('PREDNISONE','PREDNISOLONE','METHYLPREDNISOLONE','DEXAMETHASONE','HYDROCORTISONE','TRIAMCINOLONE','BUDESONIDE','BETAMETHASONE','CORTISONE','FLUDROCORTISONE','FLUTICASONE','BECLOMETHASONE','MOMETASONE')
statin<-c('ATORVASTATIN','SIMVASTATIN','ROSUVASTATIN','PRAVASTATIN','LOVASTATIN','FLUVASTATIN','PITAVASTATIN','CERIVASTATIN')
phq<-sprintf('DPQ%03d',seq(10,90,10));parts<-list();drug_matches<-list()
for(i in seq_along(cycles)){
 d<-read_mod(i,'Demographics','demo',c('SDDSRVYR','RIAGENDR','RIDAGEYR','WTMEC2YR','SDMVPSU','SDMVSTRA','RIDRETH1','RIDRETH3','DMDEDUC2','INDFMPIR'))
 mods<-list(c('Examination','bmx','BMXBMI'),c('Laboratory','cbc','LBXWBCSI','LBDNENO','LBDLYMNO','LBDNEUT','LBDLYMP','LBXPLTSI'),c('Questionnaire','smq','SMQ020','SMQ040'),c('Questionnaire','dpq',phq),c('Questionnaire','mcq','MCQ160F'),c('Questionnaire','diq','DIQ010'),c('Questionnaire','bpq','BPQ020'),c('Questionnaire','alq','ALQ130'))
 for(m in mods)d<-left_join(d,read_mod(i,m[1],m[2],m[-c(1,2)]),by='SEQN')
 rx<-read_mod(i,'Questionnaire','rxq_rx','RXDDRUG');upper<-toupper(as.character(rx$RXDDRUG))
 hit<-function(k)as.integer(Reduce(`|`,lapply(k,function(kw)grepl(kw,upper,fixed=TRUE))))
 rx$any_nsaid<-hit(nsaid);rx$any_corticosteroid<-hit(cort);rx$any_statin<-hit(statin)
 drug_matches[[i]]<-rx %>% filter(any_nsaid==1|any_corticosteroid==1|any_statin==1) %>% distinct(RXDDRUG,any_nsaid,any_corticosteroid,any_statin) %>% mutate(cycle=cycles[i])
 rx<-rx %>% group_by(SEQN) %>% summarise(across(starts_with('any_'),~as.integer(any(.==1,na.rm=TRUE))),.groups='drop')
 d<-left_join(d,rx,by='SEQN') %>% mutate(across(starts_with('any_'),~replace_na(.,0L)),cycle=cycles[i])
 parts[[i]]<-d
}
a<-bind_rows(parts);stopifnot(!anyDuplicated(a$SEQN))
a<-a %>% mutate(across(all_of(phq),~ifelse(. %in% c(7,9),NA_real_,as.numeric(.))))
if(any(!is.na(as.matrix(a[,phq])) & !(as.matrix(a[,phq])%in%0:3)))stop('PHQ 非0–3/缺失编码')
a$PHQ9_total<-rowSums(a[,phq]);a$PHQ9_SO<-rowSums(a[,c('DPQ030','DPQ040','DPQ050','DPQ080')]);a$PHQ9_CA<-rowSums(a[,c('DPQ010','DPQ020','DPQ060','DPQ070','DPQ090')])
a<-a %>% mutate(age=as.numeric(RIDAGEYR),sex=as.numeric(RIAGENDR),bmi=as.numeric(BMXBMI),smoking=case_when(SMQ020==2~3,SMQ020==1&SMQ040==3~2,SMQ020==1&SMQ040%in%c(1,2)~1,TRUE~NA_real_),depression=ifelse(!is.na(PHQ9_total)&PHQ9_total>=10,1,0),Neutrophil_count=LBDNENO*1000,Lymphocyte_count=LBDLYMNO*1000,Platelet_count=LBXPLTSI*1000,WBC_count=LBXWBCSI*1000,WTMEC_adj=WTMEC2YR/length(cycles),stroke=case_when(MCQ160F==1~1,MCQ160F==2~0,TRUE~NA_real_),diabetes=case_when(DIQ010==1~1,DIQ010%in%c(2,3)~0,TRUE~NA_real_),hypertension=case_when(BPQ020==1~1,BPQ020==2~0,TRUE~NA_real_),race_var=coalesce(RIDRETH3,RIDRETH1),race_eth=case_when(race_var==3~'Non-Hispanic White',race_var==4~'Non-Hispanic Black',race_var%in%c(1,2)~'Hispanic',race_var%in%c(5,6,7)~'Other/Multiracial',TRUE~NA_character_),education=case_when(DMDEDUC2%in%c(1,2)~'Less than high school',DMDEDUC2==3~'High school/GED',DMDEDUC2==4~'Some college',DMDEDUC2==5~'College or above',TRUE~NA_character_),pir=ifelse(INDFMPIR>=0&INDFMPIR<=5,as.numeric(INDFMPIR),NA_real_),alcohol_drinks=case_when(ALQ130%in%c(777,999)~NA_real_,ALQ130==0~0,TRUE~as.numeric(ALQ130)))
eligible<-a %>% filter(age>=18,!is.na(PHQ9_total))
cbc<-eligible %>% filter(Neutrophil_count>0,Lymphocyte_count>0,Platelet_count>0)
complete<-cbc %>% filter(!is.na(bmi),!is.na(smoking))
# Old script also checks age/sex, finite SII, and design metadata; record these separately without changing Figure 1 counts.
complete$legacy_extra_valid<-with(complete,!is.na(age)&!is.na(sex)&is.finite(log(Neutrophil_count*Platelet_count/Lymphocyte_count))&!is.na(WTMEC2YR)&!is.na(SDMVPSU)&!is.na(SDMVSTRA))
flow<-data.frame(step=c('all_DEMO','adult_complete_PHQ9','positive_complete_Neutrophil_Lymphocyte_Platelet','complete_BMI_and_smoking','legacy_extra_metadata_invalid','depression_ge10','stroke_missing'),computed=c(nrow(a),nrow(eligible),nrow(cbc),nrow(complete),sum(!complete$legacy_extra_valid),sum(complete$depression==1),sum(is.na(complete$stroke))),reference=c(NA,36259,34748,33174,NA,2902,756))
flow$difference<-flow$computed-flow$reference
write_new(flow,'C1_flow_counts_v3.csv');write_new(bind_rows(manifest),'C1_input_manifest_v3.csv');write_new(bind_rows(labels),'C1_input_variable_labels_v3.csv');write_new(bind_rows(drug_matches),'C1_legacy_medication_keyword_matches_v3.csv')

complete$smoking<-factor(c("Current","Former","Never")[complete$smoking],levels=c("Never","Current","Former"))
save_new(complete,"C1_rebuilt_aligned_v3.rds");save_new(a,"C1_all_DEMO_reconstructed_v3.rds")
stopifnot(nrow(complete)==33174,sum(complete$depression)==2902,all(is.na(flow$difference)|flow$difference==0))
