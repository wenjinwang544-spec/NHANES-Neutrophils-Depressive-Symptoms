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
neutro_domain_entry('scripts/nhanes/dependency_run_B5.R')
options(warn=1)
suppressPackageStartupMessages({library(survey);library(haven)})
out<-neutro_normalize('work/education'); csv<-function(x,n){p<-file.path(out,n);stopifnot(!neutro_file_exists(p));neutro_write_csv(x,p,row.names=FALSE)}
p<-c('work/cohort/supplement_input.rds','work/alcohol/nhanes_final_enriched_with_alcohol4_B4.rds');hash<-neutro_md5(p);d0<-neutro_readRDS(p[1]);d<-neutro_readRDS(p[2]);stopifnot(all(vapply(names(d0),function(v)identical(d0[[v]],d[[v]]),logical(1))))
neutro_cat('START',as.character(Sys.time()),'\n');print(sessionInfo());neutro_cat('Contrast test = t; df = degf(design). E2 is layer 3 on same E1 common sample, reused without duplicate fit.\n')
audit<-persons<-data.frame();sources<-data.frame(path=p,md5=unname(hash))
for(cy in 4:7){yr<-paste(1997+2*cy,1998+2*cy,sep='-');root<-file.path('data/nhanes',yr);types<-c('demo','dpq','cbc','bmx','smq');dirs<-c('Demographics','Questionnaire','Laboratory','Examination','Questionnaire');paths<-file.path(root,dirs,paste0(types,'_',letters[cy],'.xpt'));x<-lapply(paths,read_xpt);names(x)<-types;sources<-rbind(sources,data.frame(path=paths,md5=unname(neutro_md5(paths))))
 a<-x$demo[x$demo$RIDAGEYR%in%18:19,];pull<-function(t,v)as.numeric(x[[t]][[v]][match(a$SEQN,x[[t]]$SEQN)])
 q<-sapply(sprintf('DPQ%03d',seq(10,90,10)),function(v)pull('dpq',v));s1<-rowSums(!is.na(q)&q%in%0:3)==9
 n<-if('LBDNENO'%in%names(x$cbc))pull('cbc','LBDNENO')else pull('cbc','LBDNEUT')*pull('cbc','LBXWBCSI')/100
 l<-if('LBDLYMNO'%in%names(x$cbc))pull('cbc','LBDLYMNO')else pull('cbc','LBDLYMP')*pull('cbc','LBXWBCSI')/100
 pl<-pull('cbc','LBXPLTSI');bmi<-pull('bmx','BMXBMI');s020<-pull('smq','SMQ020');s040<-pull('smq','SMQ040');s2<-s1&complete.cases(n,l,pl)&n>0&l>0&pl>0;s3<-s2&!is.na(bmi);s4<-s3&!is.na(s020)&!is.na(s040)
 smoking_valid<-s020%in%2|(s020%in%1&s040%in%1:3);inrds<-a$SEQN%in%d0$SEQN[d0$SDDSRVYR==cy]
 counts<-c(nrow(a),sum(s1),sum(s2),sum(s3),sum(s4),sum(inrds));drop<-head(counts,-1)-tail(counts,-1);steps<-c('PHQ9九项完整','三计数完整且>0','BMI非缺失','SMQ020和SMQ040均非缺失','进入原RDS')
 audit<-rbind(audit,data.frame(cycle=yr,DEMO_N=counts[1],PHQ_complete=counts[2],counts_positive=counts[3],BMI_nonmissing=counts[4],both_SMQ_nonmissing=counts[5],RDS_N=counts[6],largest_drop_step=steps[which.max(drop)],largest_drop=max(drop),SMQ020_nonmissing_all=sum(!is.na(s020)),SMQ040_nonmissing_all=sum(!is.na(s040)),smoking_constructible_after_BMI=sum(s3&smoking_valid)))
 persons<-rbind(persons,data.frame(cycle=yr,SEQN=a$SEQN,age=a$RIDAGEYR,PHQ_complete=s1,counts_positive=s2,BMI_nonmissing=s3,SMQ020=s020,SMQ040=s040,both_SMQ_nonmissing=s4,smoking_constructible=smoking_valid,in_original_RDS=inrds))
}
csv(audit,'B5_0_flow.csv');csv(persons,'B5_0_person_audit.csv');print(audit)
csv(data.frame(cycle=paste(seq(2005,2017,2),seq(2006,2018,2),sep='-'),SMQ_file_min_age=12,SMQ020_min_age=c(rep(20,4),rep(18,3)),url=paste0('https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/',seq(2005,2017,2),'/DataFiles/SMQ_',LETTERS[4:10],'.htm')),'B5_0_official_SMQ_ages.csv')
d$cycle_factor<-factor(d$SDDSRVYR);ev<-c('alcohol4','education','pir','diabetes','hypertension','any_nsaid','any_corticosteroid','any_statin');d$common<-complete.cases(d[,ev]);d$PHQ15<-as.numeric(d$PHQ9_total>=15)
ds<-svydesign(id=~SDMVPSU,strata=~SDMVSTRA,weights=~WTMEC_adj,nest=TRUE,data=d);des<-subset(ds,common);stopifnot(nrow(des$variables)==29669)
base<-'log_Neutrophil + log_Lymphocyte + log_Platelet + age + sex + bmi + smoking + cycle_factor + race_eth';cells<-c('log_Neutrophil','log_Lymphocyte','log_Platelet');fits<-list();coeffs<-cts<-registry<-data.frame()
fit<-function(id,extra='',outcome='depression',pois=FALSE){f<-as.formula(paste(outcome,'~',base,extra));ww<-character();neutro_cat('\nMODEL_START',id,'\n');print(f)
 m<-withCallingHandlers(svyglm(f,design=des,family=if(pois)quasipoisson(link='log')else quasibinomial()),warning=function(w){ww<<-c(ww,conditionMessage(w));neutro_cat('WARNING',conditionMessage(w),'\n');invokeRestart('muffleWarning')});print(summary(m));ci<-confint(m);print(ci);cs<-coef(summary(m));b<-coef(m)[cells]
 coeffs<<-rbind(coeffs,data.frame(model=id,N=nobs(m),term=cells,beta=unname(b),SE=sqrt(diag(vcov(m)))[cells],effect=exp(b),lower=exp(ci[cells,1]),upper=exp(ci[cells,2]),P=cs[cells,ncol(cs)]))
 for(j in 2:3){v<-setNames(rep(0,length(coef(m))),names(coef(m)));v[cells[1]]<-1;v[cells[j]]<--1;c<-svycontrast(m,v);bb<-as.numeric(coef(c));se<-as.numeric(SE(c));ct<-data.frame(model=id,contrast=paste(cells[1],cells[j],sep=' - '),beta=bb,SE=se,t=bb/se,df=degf(des),P=2*pt(-abs(bb/se),degf(des)));cts<<-rbind(cts,ct);print(ct)}
 registry<<-rbind(registry,data.frame(model=id,N=nobs(m),formula=paste(deparse(f),collapse=' '),family=family(m)$family,design_df=degf(des),residual_df=m$df.residual,converged=m$converged,warning_count=length(ww),warnings=paste(ww,collapse=' | ')));fits[[id]]<<-m;neutro_cat('MODEL_END',id,'\n')}
extra<-c('','+ alcohol4','+ alcohol4 + education + pir','+ alcohol4 + education + pir + diabetes + hypertension','+ alcohol4 + education + pir + diabetes + hypertension + any_nsaid + any_corticosteroid + any_statin');ids<-c('M3','M3_alcohol4','E2','E2_diabetes_hypertension','E1')
for(i in seq_along(ids))fit(ids[i],extra[i])
b0<-coef(fits$M3)['log_Neutrophil'];changes<-data.frame(model=ids,beta=sapply(fits[ids],function(m)coef(m)['log_Neutrophil']));changes$percent_change_vs_M3<-100*(changes$beta-b0)/b0;csv(changes,'B5_1_logOR_percent_change.csv')
# Weighted SD and quartiles on the same E1 complete-case domain.
lsd<-sqrt(as.numeric(coef(svyvar(~log_Neutrophil,des))));qq<-as.numeric(coef(svyquantile(~Neutrophil_count,des,quantiles=c(.25,.75),ci=FALSE)));deltas<-c(1,lsd,log(5000/4000),log(qq[2]/qq[1]));scales<-c('per_1_ln','per_1_weighted_SD','4000_to_5000','weighted_Q25_to_Q75');sc<-data.frame()
for(id in c('M3','E2','E1')){m<-fits[[id]];b<-coef(m)['log_Neutrophil'];ci<-confint(m)['log_Neutrophil',];sc<-rbind(sc,data.frame(model=id,N=nobs(m),scale=scales,delta_log=deltas,OR=exp(b*deltas),lower=exp(ci[1]*deltas),upper=exp(ci[2]*deltas)))}
csv(sc,'B5_2_effect_scales.csv');csv(data.frame(N=29669,weighted_log_SD=lsd,weighted_count_Q25=qq[1],weighted_count_Q75=qq[2]),'B5_2_scale_parameters.csv')
fit('PHQ15_M3',outcome='PHQ15');fit('PHQ15_E1',extra[5],outcome='PHQ15');fit('PHQtotal_M3',outcome='PHQ9_total',pois=TRUE);fit('PHQtotal_E1',extra[5],outcome='PHQ9_total',pois=TRUE)
csv(coeffs,'B5_coefficients.csv');csv(cts,'B5_t_contrasts.csv');csv(registry,'B5_model_registry.csv');csv(sources,'B5_sources.csv')
stopifnot(nrow(registry)==9,all(registry$N==29669),all(cts$df==109),all(neutro_md5(sources$path)==sources$md5))
# Validate E1 and common-sample M3 against B4 without fitting additional models.
b4<-neutro_read_csv('work/alcohol/B4_joint_coefficients.csv');for(pair in list(c('M3','B4_3_same_sample_M3'),c('E1','B4_3_E1'))){a<-coeffs[coeffs$model==pair[1],];b<-b4[b4$model==pair[2],];stopifnot(max(abs(a$effect-b$OR))<1e-10)}
csv(data.frame(check=c('9_prespecified_models','all_N_29669','contrast_df_109','source_hashes_unchanged','B4_M3_E1_reproduced'),passed=TRUE),'B5_validation.csv')
csv(data.frame(status='COMPLETE_WAITING_USER',models_run=9,common_sample_N=29669),'B5_status.csv')
neutro_cat('END COMPLETE_WAITING_USER\n')
