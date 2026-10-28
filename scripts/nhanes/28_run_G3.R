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
neutro_domain_entry('scripts/nhanes/28_run_G3.R')
# Derived candidate; see provenance/script_selection.csv. Invoke via run_job.py.
options(warn=1,digits=17)
suppressPackageStartupMessages(library(survey))
out<-'work/shapley';SDfixed<-0.381754443
csv<-function(x,n){p<-file.path(out,n);stopifnot(!neutro_file_exists(p));neutro_write_csv(x,p,row.names=FALSE,na='NA')}
neutro_cat('START',as.character(Sys.time()),'\n');print(sessionInfo());neutro_capture_output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
inputs<-c('work/analytic/analytic_v3.rds','work/cohort/C2_full_adult_design.rds');hash<-neutro_md5(inputs)
d<-as.data.frame(neutro_readRDS(inputs[1]));des<-neutro_readRDS(inputs[2]);b<-des$variables;ix<-match(b$SEQN,d$SEQN);b$included<-!is.na(ix)
for(v in setdiff(names(d),c('SEQN','SDMVPSU','SDMVSTRA','WTMEC_adj')))b[[v]]<-d[[v]][ix]
b$cycle_factor<-factor(b$SDDSRVYR)
groups<-list(A='log10_cotinine',B='alcohol4',C=c('education','pir'),D=c('diabetes','hypertension'),E=c('any_nsaid','any_corticosteroid','any_statin'))
b$cc<-b$included&complete.cases(b[,unlist(groups)])
b$ln_NLR<-log(b$Neutrophil_count/b$Lymphocyte_count);b$ln_PLR<-log(b$Platelet_count/b$Lymphocyte_count);b$ln_SII<-log((b$Neutrophil_count/1000)*(b$Platelet_count/1000)/(b$Lymphocyte_count/1000))
des$variables<-b;full<-subset(des,included);cc<-subset(des,cc);stopifnot(nrow(full$variables)==33174,nrow(cc$variables)==29360)
cov1<-c('age','sex','bmi','smoking','cycle_factor','race_eth');cells<-c('log_Neutrophil','log_Lymphocyte','log_Platelet')
sets<-lapply(0:31,function(mask)which(as.logical(intToBits(mask)[1:5])))
setname<-function(i)if(length(sets[[i]]))paste(names(groups)[sets[[i]]],collapse='')else 'none'
formula_for<-function(i,joint=TRUE,cycle=TRUE){v<-c(if(joint)cells else cells[1],if(cycle)cov1 else setdiff(cov1,'cycle_factor'),unlist(groups[sets[[i]]]));reformulate(v,'depression')}
registry<-list();fits<-list();counter<-0
fit<-function(id,f,dd){neutro_cat('MODEL_START',id,'N',nrow(dd$variables),'\n');print(f);ww<-character();err<-''
 m<-tryCatch(withCallingHandlers(svyglm(f,design=dd,family=quasibinomial()),warning=function(w){ww<<-c(ww,conditionMessage(w));neutro_cat('WARNING',id,conditionMessage(w),'\n');invokeRestart('muffleWarning')}),error=function(e){err<<-conditionMessage(e);neutro_cat('ERROR',err,'\n');NULL})
 counter<<-counter+1;registry[[counter]]<<-data.frame(id=id,replicate=0,N=nrow(dd$variables),formula=paste(deparse(f),collapse=' '),converged=if(is.null(m))FALSE else m$converged,warnings=paste(ww,collapse=' | '),error=err)
 if(!is.null(m))fits[[id]]<<-m;neutro_cat('MODEL_END',id,'\n');m}
effect<-function(m,delta=SDfixed){if(is.null(m))return(data.frame(beta=NA,SE=NA,ln_OR_SD=NA,OR=NA,lower=NA,upper=NA,P=NA,df=NA));beta<-coef(m)['log_Neutrophil'];se<-sqrt(vcov(m)['log_Neutrophil','log_Neutrophil']);ci<-confint(m)['log_Neutrophil',];data.frame(beta=unname(beta),SE=unname(se),ln_OR_SD=unname(beta*delta),OR=unname(exp(beta*delta)),lower=unname(exp(ci[1]*delta)),upper=unname(exp(ci[2]*delta)),P=coef(summary(m))['log_Neutrophil',4],df=m$df.residual)}
neutro_cat('Coefficient CI/P: t with model residual df, consistent with preceding D2/G1. Joint Wald: denominator degf(domain design). Bootstrap ordinary paired within-stratum PSU resampling from full base, then subset; percentile CI; 200 repetitions; seed 20261014.\n')
subsets<-data.frame();xlist<-list();betas<-numeric(32)
for(i in 1:32){id<-paste0('G3_subset_',i-1);m<-fit(id,formula_for(i),cc);stopifnot(!is.null(m));betas[i]<-coef(m)['log_Neutrophil'];xlist[[i]]<-model.matrix(m);subsets<-rbind(subsets,cbind(data.frame(mask=i-1,groups=setname(i),N=nrow(cc$variables)),effect(m,1)))}
# Shapley weights equal marginal averages over all 120 permutations.
shapley<-function(beta){s<-numeric(5);for(g in 1:5)for(mask in 0:31)if(bitwAnd(mask,bitwShiftL(1,g-1))==0){k<-length(sets[[mask+1]]);s[g]<-s[g]+factorial(k)*factorial(4-k)/factorial(5)*(beta[bitwOr(mask,bitwShiftL(1,g-1))+1]-beta[mask+1])};s}
sh<-shapley(betas);total<-betas[32]-betas[1];stopifnot(abs(sum(sh)-total)<1e-12)
psu<-interaction(b$SDMVSTRA,b$SDMVPSU,drop=TRUE);pid<-as.integer(psu);pairs<-unique(data.frame(pid=pid,stratum=b$SDMVSTRA));strata<-split(pairs$pid,pairs$stratum)
set.seed(20261014);mult<-matrix(0,nlevels(psu),200)
for(r in 1:200)for(ids in strata)mult[ids,r]<-tabulate(match(sample(ids,length(ids),replace=TRUE),ids),length(ids))
take<-match(cc$variables$SEQN,b$SEQN);w<-weights(cc);y<-cc$variables$depression
boot_one<-function(r){wr<-w*mult[pid[take],r];keep<-wr>0;bb<-numeric(32);rr<-list()
 for(i in 1:32){ww<-character();er<-'';m<-tryCatch(withCallingHandlers(glm.fit(xlist[[i]][keep,,drop=FALSE],y[keep],weights=wr[keep]/mean(wr[keep]),family=quasibinomial()),warning=function(e){ww<<-c(ww,conditionMessage(e));invokeRestart('muffleWarning')}),error=function(e){er<<-conditionMessage(e);NULL});bb[i]<-if(is.null(m))NA else m$coefficients['log_Neutrophil'];rr[[i]]<-data.frame(id=paste0('G3_subset_',i-1),replicate=r,N=sum(keep),converged=if(is.null(m))FALSE else m$converged,warnings=paste(ww,collapse=' | '),error=er)}
 ss<-shapley(bb);tt<-bb[32]-bb[1];list(values=data.frame(replicate=r,group=names(groups),Shapley=ss,total_change=tt,percent=100*ss/tt),registry=do.call(rbind,rr),betas=data.frame(replicate=r,mask=0:31,beta=bb))}
boot<-list();for(start in seq(1,200,by=10)){ids<-start:min(start+9,200);res<-parallel::mclapply(ids,boot_one,mc.cores=4,mc.set.seed=FALSE);boot<-c(boot,res);neutro_cat('BOOTSTRAP',max(ids),'/200',as.character(Sys.time()),'\n')}
bv<-do.call(rbind,lapply(boot,`[[`,'values'));br<-do.call(rbind,lapply(boot,`[[`,'registry'));bb<-do.call(rbind,lapply(boot,`[[`,'betas'))
ci<-t(sapply(names(groups),function(g)quantile(bv$percent[bv$group==g],c(.025,.975),na.rm=TRUE)))
shout<-data.frame(group=names(groups),variables=sapply(groups,paste,collapse='+'),baseline_beta=betas[1],full_beta=betas[32],total_change=total,Shapley=sh,percent=100*sh/total,percent_lower=ci[,1],percent_upper=ci[,2]);stopifnot(abs(sum(shout$percent)-100)<1e-8)
csv(subsets,'G3_1_subset_models.csv');csv(shout,'G3_1_Shapley.csv');csv(bv,'G3_1_bootstrap_Shapley.csv');csv(bb,'G3_1_bootstrap_betas.csv');csv(br,'G3_bootstrap_model_registry.csv')
# 256 specifications, reuse the 32 identical unrestricted joint fits.
restrictions<-c('none','stroke0','age20','exclude_extreme');spec<-data.frame()
for(restr in restrictions){dd<-switch(restr,none=cc,stroke0=subset(cc,stroke==0),age20=subset(cc,age>=20),exclude_extreme=subset(cc,Neutrophil_count<=9600&Lymphocyte_count<=5000))
 for(joint in c(TRUE,FALSE))for(i in 1:32){id<-paste('G3_spec',restr,if(joint)'joint'else'neut_only',i-1,sep='_');if(restr=='none'&&joint){m<-fits[[paste0('G3_subset_',i-1)]];neutro_cat('REUSED',id,paste0('G3_subset_',i-1),'\n')}else m<-fit(id,formula_for(i,joint),dd)
 spec<-rbind(spec,cbind(data.frame(id=id,mask=i-1,groups=setname(i),exposure=if(joint)'joint'else'neut_only',restriction=restr,N=nrow(dd$variables)),effect(m)))} }
spec$rank<-rank(spec$ln_OR_SD,ties.method='first');spec<-spec[order(spec$rank),];csv(spec,'G3_2_specifications.csv')
summary<-data.frame(specifications=nrow(spec),median_ln_OR_SD=median(spec$ln_OR_SD),min_ln_OR_SD=min(spec$ln_OR_SD),max_ln_OR_SD=max(spec$ln_OR_SD),median_OR=median(spec$OR),min_OR=min(spec$OR),max_OR=max(spec$OR),CI_excludes1=sum(spec$lower>1|spec$upper<1),CI_excludes1_proportion=mean(spec$lower>1|spec$upper<1));csv(summary,'G3_2_summary.csv')
plot_spec<-function(){layout(matrix(c(1,2),2,1),heights=c(2,1.4));par(mar=c(1,4,2,1));plot(spec$rank,spec$ln_OR_SD,pch=16,cex=.45,xlab='',ylab='ln OR per 1 SD',ylim=range(log(c(spec$lower,spec$upper))));segments(spec$rank,log(spec$lower),spec$rank,log(spec$upper),col=ifelse(spec$lower>1,'#20689b','#999999'));points(spec$rank,spec$ln_OR_SD,pch=16,cex=.4);abline(h=0,lty=2);par(mar=c(4,8,0,1));mat<-rbind(sapply(1:5,function(g)as.integer(bitwAnd(spec$mask,bitwShiftL(1,g-1))>0)) |> t(),joint=as.integer(spec$exposure=='joint'),stroke0=as.integer(spec$restriction=='stroke0'),age20=as.integer(spec$restriction=='age20'),extreme=as.integer(spec$restriction=='exclude_extreme'));image(1:256,1:nrow(mat),t(mat[nrow(mat):1,,drop=FALSE]),col=c('white','#20689b'),axes=FALSE,xlab='Specification rank',ylab='');axis(1);axis(2,at=1:nrow(mat),labels=rev(c('A Cotinine','B Alcohol','C Education/PIR','D Diabetes/HTN','E Medication','Joint counts','Stroke=0','Age>=20','Exclude extreme')),las=2,cex.axis=.8)}
neutro_png(file.path(out,'G3_specification_curve.png'),width=12,height=8,units='in',res=300);plot_spec();neutro_dev_off();neutro_pdf(file.path(out,'G3_specification_curve.pdf'),width=12,height=8);plot_spec();neutro_dev_off()
# Cycle is constant in a within-cycle domain; remove its unidentifiable factor term only.
neutro_cat('Within-cycle Model1/3: cycle factor is constant; omit unidentifiable cycle factor. Other covariates unchanged.\n')
cycles<-sort(unique(d$SDDSRVYR));cyres<-data.frame()
for(mod in c('Model1','Model3'))for(cy in cycles){dd<-if(mod=='Model1')subset(full,SDDSRVYR==cy)else subset(cc,SDDSRVYR==cy);m<-fit(paste('G3_cycle',mod,cy,sep='_'),formula_for(if(mod=='Model1')1 else 32,cycle=FALSE),dd);cyres<-rbind(cyres,cbind(data.frame(model=mod,cycle=cy,N=nrow(dd$variables),cases=sum(dd$variables$depression==1),design_df=degf(dd)),effect(m)))}
# DerSimonian-Laird random effects: inverse-variance weights; normal pooled CI.
dl<-function(z){yi<-z$beta*SDfixed;vi<-(z$SE*SDfixed)^2;w<-1/vi;k<-length(yi);mu0<-sum(w*yi)/sum(w);Q<-sum(w*(yi-mu0)^2);C<-sum(w)-sum(w^2)/sum(w);tau2<-max(0,(Q-(k-1))/C);wr<-1/(vi+tau2);mu<-sum(wr*yi)/sum(wr);se<-sqrt(1/sum(wr));data.frame(k=k,ln_OR_SD=mu,SE=se,OR=exp(mu),lower=exp(mu-qnorm(.975)*se),upper=exp(mu+qnorm(.975)*se),tau2=tau2,I2=if(Q>0)max(0,(Q-(k-1))/Q)*100 else 0,Q=Q)}
meta<-loo<-data.frame();for(mod in c('Model1','Model3')){z<-cyres[cyres$model==mod,];meta<-rbind(meta,cbind(data.frame(model=mod),dl(z)));for(cy in cycles)loo<-rbind(loo,cbind(data.frame(model=mod,omitted_cycle=cy),dl(z[z$cycle!=cy,])))}
csv(cyres,'G3_3_cycle_estimates.csv');csv(meta,'G3_3_meta_DL.csv');csv(loo,'G3_3_leave_one_cycle.csv')

csv(do.call(rbind,registry),'G3_model_registry.csv');stopifnot(all(hash==neutro_md5(inputs)));neutro_cat('SHAPLEY_SPECIFICATION_CYCLE_COMPLETE\n')
