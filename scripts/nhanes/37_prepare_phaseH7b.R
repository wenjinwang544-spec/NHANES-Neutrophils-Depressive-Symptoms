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
neutro_domain_entry('scripts/nhanes/37_prepare_phaseH7b.R')
# Derived candidate; see provenance/script_selection.csv. Invoke via run_job.py.
options(warn=1)
suppressPackageStartupMessages(library(haven))
out <- 'work/sleep'
read_paths <- character() # Actual XPT inputs, not a historical download folder
source_analytic <- 'work/analytic/analytic_v3.rds'
stopifnot(neutro_file_exists(source_analytic))
dat <- as.data.frame(neutro_readRDS(source_analytic))
stopifnot(nrow(dat)==33174L, !anyDuplicated(dat$SEQN))
cycle_map <- data.frame(SDDSRVYR=4:10, year=c(2005,2007,2009,2011,2013,2015,2017),
  suffix=LETTERS[4:10], cycle=c('2005-2006','2007-2008','2009-2010','2011-2012','2013-2014','2015-2016','2017-2018'))
`%||%` <- function(x,y) if(is.null(x)) y else x
special <- c(7,9,77,99,777,999,7777,9999,77777,99999)
recode_special <- function(x) { if(is.numeric(x)) x[x %in% special] <- NA_real_; x }
recode_sleep <- function(x) { if(is.numeric(x)) x[x %in% c(77,99)] <- NA_real_; x }
yn <- function(x) { x <- recode_special(as.numeric(x)); ifelse(x==1,1,ifelse(x==2,0,NA_real_)) }
needed <- c('asthma','emphysema','bronchitis','copd','arthritis','cancer','thyroid','liver','kidney','sleep_hours')
merged <- data.frame(SEQN=dat$SEQN)
cycle_avail <- list(); label_rows <- list(); value_rows <- list()
for(i in seq_len(nrow(cycle_map))) {
  cm <- cycle_map[i,]
  ids <- dat$SEQN[dat$SDDSRVYR==cm$SDDSRVYR]
  getx <- function(stem) { folder<-file.path('data/nhanes',cm$cycle,'Questionnaire');fs<-neutro_list_files(folder,full.names=TRUE);f<-fs[tolower(basename(fs))==tolower(paste0(stem,'_',cm$suffix,'.xpt'))];stopifnot(length(f)==1);read_paths <<- unique(c(read_paths,neutro_normalize(f,mustWork=TRUE)));neutro_read_xpt(f) }
  mcq <- getx('MCQ'); kiq <- getx('KIQ_U'); paq <- getx('PAQ'); slq <- getx('SLQ')
  item_map <- list(asthma=c('MCQ010'), emphysema=c('MCQ160G'), bronchitis=c('MCQ160K'),
    copd=c('MCQ160O'), arthritis=c('MCQ160A'), cancer=c('MCQ220'), thyroid=c('MCQ160M'), liver=c('MCQ160L'), kidney=c('KIQ022'))
  rec <- data.frame(SEQN=ids)
  for(v in names(item_map)) {
    src <- item_map[[v]][1]
    tab <- if(v=='kidney') kiq else mcq
    has <- src %in% names(tab)
    if(has) { vv <- yn(tab[[src]][match(ids,tab$SEQN)]); rec[[v]] <- vv
      label_rows[[length(label_rows)+1]] <- data.frame(cycle=cm$cycle,variable=src,label=attr(tab[[src]],'label')%||%NA_character_,available=TRUE)
      value_rows[[length(value_rows)+1]] <- data.frame(cycle=cm$cycle,variable=src,observed=sum(!is.na(vv)),missing=sum(is.na(vv)))
    } else { rec[[v]] <- NA_real_; label_rows[[length(label_rows)+1]] <- data.frame(cycle=cm$cycle,variable=src,label=NA_character_,available=FALSE) }
  }
  sl_var <- if(cm$SDDSRVYR<=8) 'SLD010H' else 'SLD012'
  sl_has <- sl_var %in% names(slq)
  if(sl_has) { sl <- recode_sleep(as.numeric(slq[[sl_var]][match(ids,slq$SEQN)])); rec$sleep_hours <- sl
    label_rows[[length(label_rows)+1]] <- data.frame(cycle=cm$cycle,variable=sl_var,label=attr(slq[[sl_var]],'label')%||%NA_character_,available=TRUE)
    value_rows[[length(value_rows)+1]] <- data.frame(cycle=cm$cycle,variable=sl_var,observed=sum(!is.na(sl)),missing=sum(is.na(sl)))
  } else rec$sleep_hours <- NA_real_
  merge_idx <- match(merged$SEQN,rec$SEQN)
  for(v in setdiff(names(rec),'SEQN')) merged[[v]][dat$SDDSRVYR==cm$SDDSRVYR] <- rec[[v]][match(merged$SEQN[dat$SDDSRVYR==cm$SDDSRVYR],rec$SEQN)]
  pa_vars <- grep('^PAQ(605|620|635|650|665)$|^PAQ(100|180)$',names(paq),value=TRUE)
  cycle_avail[[i]] <- data.frame(cycle=cm$cycle,MCQ_file=paste0('MCQ_',cm$suffix,'.xpt'),KIQ_file=paste0('KIQ_U_',cm$suffix,'.xpt'),
    PAQ_file=paste0('PAQ_',cm$suffix,'.xpt'),SLQ_file=paste0('SLQ_',cm$suffix,'.xpt'),
    MCQ160O_available='MCQ160O'%in%names(mcq),PAQ_core_items=paste(pa_vars,collapse=';'),SLQ_sleep_field=if(sl_has)sl_var else 'not found')
}
merged$chronic_respiratory <- ifelse(merged$asthma==1|merged$emphysema==1|merged$bronchitis==1,1,
  ifelse(!is.na(merged$asthma)&!is.na(merged$emphysema)&!is.na(merged$bronchitis),0,NA_real_))
# Candidate that excludes asthma, for a possible separate asthma + COPD/airway construction.
merged$other_respiratory <- ifelse(merged$emphysema==1|merged$bronchitis==1|merged$copd==1,1,
  ifelse(!is.na(merged$emphysema)&!is.na(merged$bronchitis)&!is.na(merged$copd),0,NA_real_))
merged$activity_any_uniform <- NA_real_ # no past-week harmonized construction for 2005-06 PAQ_D
merged$cycle <- dat$SDDSRVYR
merged$depression <- dat$depression
merged$log_Neutrophil <- dat$log_Neutrophil
merged$log_Lymphocyte <- dat$log_Lymphocyte
merged$log_Platelet <- dat$log_Platelet
merged$WTMEC_adj <- dat$WTMEC_adj
neutro_saveRDS(merged,file.path(out,'nhanes_h7b_merged.rds'))
neutro_write_csv(do.call(rbind,cycle_avail),file.path(out,'cycle_file_variable_availability.csv'),row.names=FALSE)
neutro_write_csv(do.call(rbind,label_rows),file.path(out,'variable_labels_by_cycle.csv'),row.names=FALSE)
# Official CDC codebook URL crosswalk.
urls <- do.call(rbind,lapply(seq_len(nrow(cycle_map)),function(i){cm<-cycle_map[i,]; base<-sprintf('https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/%d/DataFiles/',cm$year);data.frame(cycle=cm$cycle,
 MCQ=paste0(base,'MCQ_',cm$suffix,'.htm'),KIQ_U=paste0(base,'KIQ_U_',cm$suffix,'.htm'),PAQ=paste0(base,'PAQ_',cm$suffix,'.htm'),SLQ=paste0(base,'SLQ_',cm$suffix,'.htm'))}))
neutro_write_csv(urls,file.path(out,'CDC_codebook_urls_by_cycle.csv'),row.names=FALSE)
neutro_write_csv(do.call(rbind,value_rows),file.path(out,'raw_variable_counts_by_cycle.csv'),row.names=FALSE)
analysis_vars <- c('sleep_hours','asthma','emphysema','bronchitis','chronic_respiratory','arthritis','cancer','thyroid','liver','kidney')
miss <- do.call(rbind,lapply(analysis_vars,function(v){x<-merged[[v]];data.frame(variable=v,N=length(x),missing_N=sum(is.na(x)),missing_pct=100*mean(is.na(x)),over_10pct=mean(is.na(x))>.10)}))
neutro_write_csv(miss,file.path(out,'variable_missingness_N33174.csv'),row.names=FALSE)
# Download manifest (SHA-256 of every new XPT file).
files<-sort(unique(read_paths));manifest<-neutro_input_manifest(files,expected=nrow(cycle_map)*4L)
neutro_write_csv(manifest,file.path(out,'raw_download_manifest.csv'),row.names=FALSE)
neutro_cat('N=',nrow(merged),'\n');print(do.call(rbind,cycle_avail));print(miss);print(data.frame(file=basename(files),size=neutro_file_info(files)$size))
