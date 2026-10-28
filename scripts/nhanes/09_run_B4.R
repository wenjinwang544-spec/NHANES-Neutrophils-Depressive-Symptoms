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
neutro_domain_entry('scripts/nhanes/09_run_B4.R')
# Derived candidate; see provenance/script_selection.csv. Invoke via run_job.py.
options(warn=1)
suppressPackageStartupMessages({library(survey);library(haven)})
out<-neutro_normalize("work/alcohol");input<-neutro_normalize("work/cohort/supplement_input.rds")
csv<-function(x,n){p<-file.path(out,n);stopifnot(!neutro_file_exists(p));neutro_write_csv(x,p,row.names=FALSE,na="")}
neutro_cat("START ",format(Sys.time(),tz="Asia/Shanghai"),"\nINPUT ",input,"\nMD5 ",unname(neutro_md5(input)),"\n")
print(sessionInfo());neutro_cat("M3 unchanged; contrast=t df=degf(design); only 5 prespecified models; all domains via subset(design).\n")
d<-neutro_readRDS(input);stopifnot(nrow(d)==33174,sum(d$depression==1)==2902)
# B4-0 is completed before any construction or model.
dpqvars<-sprintf("DPQ%03d",seq(10,90,10));early<-sourcefiles<-headers<-data.frame()
for(cy in 4:6){yr<-paste(1997+2*cy,1998+2*cy,sep="-");root<-file.path("data/nhanes",yr)
 paths<-c(DEMO=file.path(root,"Demographics",paste0("demo_",letters[cy],".xpt")),DPQ=file.path(root,"Questionnaire",paste0("dpq_",letters[cy],".xpt")))
 demo<-neutro_read_xpt(paths[["DEMO"]]);dpq<-neutro_read_xpt(paths[["DPQ"]]);stopifnot(!anyDuplicated(demo$SEQN),!anyDuplicated(dpq$SEQN),all(dpqvars%in%names(dpq)))
 sourcefiles<-rbind(sourcefiles,data.frame(cycle=yr,type=names(paths),path=unname(paths),md5=unname(neutro_md5(paths))))
 for(v in names(dpq)){lab<-attr(dpq[[v]],"label");if(is.null(lab))lab<-"";headers<-rbind(headers,data.frame(cycle=yr,variable=v,label=lab))}
 a<-demo[demo$RIDAGEYR%in%c(18,19),];ix<-match(a$SEQN,dpq$SEQN);q<-dpq[ix,dpqvars];valid<-sapply(q,function(z)!is.na(z)&z%in%0:3)
 early<-rbind(early,data.frame(cycle=yr,DEMO_age18_19_N=nrow(a),DPQ_records_matched=sum(!is.na(ix)),DPQ_no_record=sum(is.na(ix)),DPQ_nine_valid_0_to_3=sum(rowSums(valid)==9),DPQ_matched_nine_incomplete=sum(!is.na(ix)&rowSums(valid)<9),DPQ_matched_all_nine_NA=sum(!is.na(ix)&rowSums(!is.na(q))==0),DPQ_matched_partial_or_invalid=sum(!is.na(ix)&rowSums(valid)<9 & rowSums(!is.na(q))>0)))
 person<-data.frame(SEQN=as.numeric(a$SEQN),age=as.numeric(a$RIDAGEYR),DPQ_record_present=!is.na(ix),nine_valid=rowSums(valid)==9,valid_items=rowSums(valid));csv(person,paste0("B4_0_age18_19_",yr,".csv"))
}
csv(early,"B4_0_DEMO_DPQ_age18_19.csv");csv(headers,"B4_0_DPQ_header_labels.csv")
age<-do.call(rbind,lapply(4:10,function(cy){x<-d[d$SDDSRVYR==cy,];data.frame(SDDSRVYR=cy,cycle=paste(1997+2*cy,1998+2*cy,sep="-"),N=nrow(x),minimum_age=min(x$age))}));csv(age,"B4_0_RDS_minimum_age.csv")
csv(data.frame(cycle=age$cycle,eligible_age_minimum=12,public_release_age_minimum=18,url=paste0("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/",seq(2005,2017,2),"/DataFiles/DPQ_",LETTERS[4:10],".htm")),"B4_0_DPQ_official_age_sources.csv")
neutro_cat("B4_0_COMPLETE\n");print(early);print(age)
# Construct only user-authorized alcohol4 and lifetime skip flag, in a new copy.
vs<-c("ALQ101","ALQ110","ALQ111","ALQ120Q","ALQ121","ALQ130");alqs<-list()
for(cy in 4:10){yr<-paste(1997+2*cy,1998+2*cy,sep="-");p<-file.path("data/nhanes",yr,"Questionnaire",paste0("alq_",letters[cy],".xpt"));x<-neutro_read_xpt(p);stopifnot(!anyDuplicated(x$SEQN));sourcefiles<-rbind(sourcefiles,data.frame(cycle=yr,type="ALQ",path=p,md5=unname(neutro_md5(p))))
 z<-data.frame(SEQN=as.numeric(x$SEQN),SDDSRVYR=cy);for(v in vs)z[[v]]<-if(v%in%names(x))as.numeric(x[[v]])else NA_real_;alqs[[yr]]<-z
}
alq<-do.call(rbind,alqs);ix<-match(paste(d$SDDSRVYR,d$SEQN),paste(alq$SDDSRVYR,alq$SEQN));stopifnot(!anyNA(ix));a<-alq[ix,]
old<-d$SDDSRVYR<=9;lt12<-old & a$ALQ101%in%2 & a$ALQ110%in%2
ref<-(!old & a$ALQ111%in%2)|(old & a$ALQ120Q%in%0)|(!old & a$ALQ121%in%0)|lt12
valid<-!is.na(a$ALQ130)&a$ALQ130>=1&!(a$ALQ130%in%c(777,999))
stopifnot(sum(lt12)==4271,sum(ref)==10911,sum(ref&valid)==0)
cat4<-rep(NA_character_,nrow(d));cat4[ref]<-"不饮酒/极少饮酒";cat4[valid & a$ALQ130<=1]<-"ALQ130<=1";cat4[valid & a$ALQ130==2]<-"ALQ130=2";cat4[valid & a$ALQ130>=3]<-"ALQ130>=3"
aug<-d;aug$alcohol4<-factor(cat4,levels=c("不饮酒/极少饮酒","ALQ130<=1","ALQ130=2","ALQ130>=3"));aug$lifetime_less12_skip<-as.logical(lt12)
stopifnot(sum(is.na(aug$alcohol4))==57,all(vapply(names(d),function(v)identical(d[[v]],aug[[v]]),logical(1))))
p<-file.path(out,"nhanes_final_enriched_with_alcohol4_B4.rds");stopifnot(!neutro_file_exists(p));neutro_saveRDS(aug,p)
aug$cycle_factor<-factor(aug$SDDSRVYR,levels=4:10)
ds<-svydesign(id=~SDMVPSU,strata=~SDMVSTRA,weights=~WTMEC_adj,nest=TRUE,data=aug)
desc<-do.call(rbind,lapply(c(levels(aug$alcohol4),"Missing"),function(cat){ix<-if(cat=="Missing")is.na(aug$alcohol4)else !is.na(aug$alcohol4)&aug$alcohol4==cat
 flagds<-update(ds,category_indicator=as.numeric(ix));prop<-svymean(~category_indicator,flagds)
 domain<-subset(ds,ix);pr<-svymean(~depression,domain)
 data.frame(alcohol4=cat,N=sum(ix),weighted_proportion_all_N33174=as.numeric(coef(prop)),depression_N=sum(aug$depression[ix]==1),depression_proportion_unweighted=mean(aug$depression[ix]),depression_proportion_weighted=as.numeric(coef(pr)))}))
csv(desc,"B4_1_alcohol4_distribution_depression.csv")
covars<-"age + sex + bmi + smoking + cycle_factor + race_eth";cells<-c("log_Neutrophil","log_Lymphocyte","log_Platelet")
registry<-data.frame();warnrows<-data.frame(model=character(),warning=character());results<-list()
fit<-function(id,des,extra,contrasts=TRUE){f<-as.formula(paste("depression ~",paste(cells,collapse=" + "),"+",covars,extra));warn<-character()
 neutro_cat("MODEL_START ",id," N=",nrow(des$variables)," formula=",paste(deparse(f),collapse=" "),"\n",sep="");print(des$call)
 m<-withCallingHandlers(svyglm(f,design=des,family=quasibinomial()),warning=function(w){warn<<-c(warn,conditionMessage(w));neutro_cat("MODEL_WARNING ",id," ",conditionMessage(w),"\n");invokeRestart("muffleWarning")})
 print(summary(m));print(confint(m));registry<<-rbind(registry,data.frame(model=id,N=nobs(m),formula=paste(deparse(f),collapse=" "),family="quasibinomial",df_residual=m$df.residual,df_design=degf(des),converged=isTRUE(m$converged),iterations=m$iter,warning_count=length(warn)))
 if(length(warn))warnrows<<-rbind(warnrows,data.frame(model=id,warning=warn))
 cs<-coef(summary(m))[cells,,drop=FALSE];ci<-confint(m)[cells,,drop=FALSE]
 cf<-data.frame(model=id,N=nobs(m),term=cells,beta=unname(coef(m)[cells]),SE=unname(sqrt(diag(vcov(m)))[cells]),OR=unname(exp(coef(m)[cells])),CI_lower=unname(exp(ci[,1])),CI_upper=unname(exp(ci[,2])),P_raw=unname(cs[,ncol(cs)]))
 ct<-data.frame()
 if(contrasts){one<-function(expr,label){o<-svycontrast(m,expr);b<-as.numeric(coef(o));se<-as.numeric(SE(o));t<-b/se;df<-degf(des);data.frame(model=id,N=nobs(m),comparison=label,log_difference=b,SE=se,t=t,df=df,P_raw=2*pt(-abs(t),df))};ct<-rbind(one(quote(log_Neutrophil-log_Lymphocyte),"Neutrophil-Lymphocyte"),one(quote(log_Neutrophil-log_Platelet),"Neutrophil-Platelet"));print(ct)}
 results[[id]]<<-list(coefficients=cf,contrasts=ct);neutro_cat("MODEL_END ",id," CONVERGED=",isTRUE(m$converged)," WARNING_COUNT=",length(warn),"\n",sep="")
}
sub4<-subset(ds,!is.na(alcohol4));stopifnot(nrow(sub4$variables)==33117)
fit("B4_2_same_sample_M3",sub4,"");fit("B4_2_M3_alcohol4",sub4,"+ alcohol4")
evars<-c("alcohol4","education","pir","diabetes","hypertension","any_nsaid","any_corticosteroid","any_statin")
ds<-update(ds,E1_complete=complete.cases(aug[,evars]));subE1<-subset(ds,E1_complete)
fit("B4_3_same_sample_M3",subE1,"");fit("B4_3_E1",subE1,"+ alcohol4 + education + pir + diabetes + hypertension + any_nsaid + any_corticosteroid + any_statin")
# Five-category sensitivity reference remains the non/rare group excluding lifetime<12.
sens<-as.character(aug$alcohol4);sens[aug$lifetime_less12_skip]<-"终生少于12杯跳问"
alcohol5<-factor(sens,levels=c(levels(aug$alcohol4),"终生少于12杯跳问"))
ds5<-update(ds,alcohol5=alcohol5);sub5<-subset(ds5,!is.na(alcohol5));stopifnot(nrow(sub5$variables)==33117)
fit("B4_4_lifetime_skip_separate",sub5,"+ alcohol5",FALSE)
csv(do.call(rbind,lapply(results,`[[`,"coefficients")),"B4_joint_coefficients.csv");csv(do.call(rbind,lapply(results,`[[`,"contrasts")),"B4_t_contrasts.csv")
csv(registry,"B4_model_registry.csv");csv(warnrows,"B4_model_warnings.csv");csv(sourcefiles,"B4_source_files.csv")
for(task in c("B4_2","B4_3")){ids<-if(task=="B4_2")c("B4_2_same_sample_M3","B4_2_M3_alcohol4")else c("B4_3_same_sample_M3","B4_3_E1")
 metrics<-c("N",unlist(lapply(cells,function(x)paste(x,c("OR","CI_lower","CI_upper","P_raw"),sep="_"))),"Neutrophil-Lymphocyte_P","Neutrophil-Platelet_P")
 tb<-data.frame(metric=metrics)
 for(id in ids){r<-results[[id]];tb[[id]]<-c(r$coefficients$N[1],unlist(lapply(cells,function(x)as.numeric(r$coefficients[r$coefficients$term==x,c("OR","CI_lower","CI_upper","P_raw")]))),r$contrasts$P_raw)}
 csv(tb,paste0(task,"_paired_summary.csv"))
}
csv(data.frame(status="COMPLETE_WAITING_USER",models_run=nrow(registry),alcohol4_nonmissing=33117,alcohol4_missing=57,new_RDS=p),"B4_status.csv")
neutro_cat("END COMPLETE_WAITING_USER MODELS_RUN=",nrow(registry),"\n")
