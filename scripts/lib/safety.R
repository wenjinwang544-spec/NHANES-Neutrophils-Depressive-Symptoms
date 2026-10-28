# Managed I/O only. Scientific formulas and statistical functions are not modified.
.neutro <- new.env(parent=emptyenv())
neutro_fail <- function(...)stop(paste0('NEUTRO_SAFETY: ',paste(...,collapse=' ')),call.=FALSE)
neutro_overlap <- function(a,b)a==b||startsWith(a,paste0(b,'/'))||startsWith(b,paste0(a,'/'))
neutro_real <- function(p){p<-path.expand(p);if(!grepl('^/',p))p<-file.path(getwd(),p);if(file.exists(p)||dir.exists(p))return(normalizePath(p,mustWork=TRUE));parent<-dirname(p);if(parent==p)neutro_fail('Cannot resolve path');file.path(neutro_real(parent),basename(p))}
neutro_setup <- function(source_root){
 a<-commandArgs(TRUE);getarg<-function(k){i<-match(k,a);if(is.na(i)||i==length(a))return('');a[i+1]}
 w<-Sys.getenv('NEUTRO_WORKSPACE');if(!nzchar(w))w<-getarg('--workspace');if(!nzchar(w))neutro_fail('Explicit initialized --workspace required; direct source without it is forbidden')
 # jsonlite is an explicit runtime dependency, used only for the path manifest.
 if(!requireNamespace('jsonlite',quietly=TRUE))neutro_fail('jsonlite is required for the safety configuration')
 c<-jsonlite::fromJSON(file.path(w,'run_config.json'),simplifyVector=FALSE)
 src<-normalizePath(source_root,mustWork=TRUE);w<-normalizePath(w,mustWork=TRUE);out<-normalizePath(file.path(w,'output'),mustWork=TRUE);inp<-normalizePath(file.path(w,'inputs'),mustWork=TRUE)
 if(c$schema!=1||neutro_real(c$workspace)!=w||neutro_real(c$source_root)!=src||neutro_real(c$output_root)!=out||neutro_real(c$input_root)!=inp||length(c$protected_roots)<1)neutro_fail('Invalid initialized layout')
 if(any(vapply(c(file.path(w,'output'),file.path(w,'inputs')),function(p){v<-Sys.readlink(p);!is.na(v)&&nzchar(v)},logical(1))))neutro_fail('Root symlink forbidden')
 protected<-vapply(c$protected_roots,neutro_real,character(1));external<-if(length(c$external_inputs))vapply(c$external_inputs,neutro_real,character(1))else character()
 policy_file<-file.path(src,'config/storage_policy.json')
 if(file.exists(policy_file)){policy<-jsonlite::fromJSON(policy_file);if(!identical(policy$schema,1L)||!isTRUE(policy$protect_source_parent))neutro_fail('Invalid archive storage policy');protected<-unique(c(protected,dirname(src)));for(r in protected)if(neutro_overlap(inp,r))neutro_fail('Inputs overlap archive storage parent')}

 for(r in c(inp,src,protected,external))if(neutro_overlap(out,r))neutro_fail('Output overlaps input/source/protected root')
 for(item in c$imports){f<-normalizePath(file.path(inp,'gwas_cache',item$name),mustWork=TRUE);if(!startsWith(f,paste0(inp,'/')))neutro_fail('Imported cache escaped input root');if(unname(tools::sha256sum(f))!=item$sha256)neutro_fail('Imported cache hash changed')}

 .neutro$out<-out;.neutro$inp<-inp;.neutro$src<-src;.neutro$protected<-c(inp,src,protected,external);.neutro$config<-c;.neutro$owned<-character();.neutro$devices<-integer();.neutro$connections<-character()
 setwd(out);options(device=function(...)neutro_fail('Implicit/default graphics device forbidden; no Rplots.pdf fallback'))
 invisible(c)
}
neutro_input <- function(p){if(inherits(p,'connection'))return(p);if(!is.character(p))return(p);vapply(p,function(x){if(startsWith(x,paste0(.neutro$out,'/data/')))x<-substring(x,nchar(.neutro$out)+2L);if(startsWith(x,'data/nhanes')){r<-.neutro$config$external_inputs$nhanes;if(is.null(r))neutro_fail('NHANES input root not configured');return(file.path(r,sub('^data/nhanes/?','',x)))};if(startsWith(x,'data/gwas/'))return(file.path(.neutro$inp,'gwas_cache',sub('^data/gwas/','',x)));if(startsWith(x,'data/'))return(file.path(.neutro$inp,sub('^data/','',x)));if(x%in%c('work/imputation/imputations_standardized.rds','work/imputation/imputation_scaling.csv')&&!base::file.exists(x)){return(file.path(.neutro$inp,'gwas_cache',basename(x)))};x},character(1),USE.NAMES=FALSE)}
neutro_output <- function(p,append=FALSE){
 if(inherits(p,'connection'))neutro_fail('Unmanaged output connection forbidden')
 if(length(p)!=1||!is.character(p)||!nzchar(p)||p%in%c('','|','stdout','stderr'))neutro_fail('Explicit output filename required')
 if(grepl('%[0-9]*d',p))neutro_fail('Unbounded graphics filename template forbidden')
 # Reject lexical traversal before canonical resolution, even when its endpoint would fall inside.
 if(any(strsplit(p,'/',fixed=TRUE)[[1]]=='..'))neutro_fail('Parent traversal forbidden')
 q<-neutro_real(p);out<-.neutro$out
 if(!startsWith(q,paste0(out,'/')))neutro_fail('Output escapes explicit output root: ',p)
 if(startsWith(q,paste0(out,'/data/')))neutro_fail('Logical data namespace is read-only')
 for(r in .neutro$protected)if(neutro_overlap(q,r))neutro_fail('Output overlaps input/source/protected root')
 raw<-if(grepl('^/',p))p else file.path(getwd(),p);if(!startsWith(raw,paste0(out,'/')))neutro_fail('Lexical output root mismatch');parts<-strsplit(substring(raw,nchar(out)+2L),'/',fixed=TRUE)[[1]];cur<-out
 for(part in parts){cur<-file.path(cur,part);link<-Sys.readlink(cur);if(!is.na(link)&&nzchar(link))neutro_fail('Symlink in output path forbidden')}
 if(append&&!grepl('[.](log|txt)$',tolower(q)))neutro_fail('Append is restricted to managed text logs')
 if(file.exists(q)||dir.exists(q)){if(!(append&&q%in%.neutro$owned))neutro_fail('Refusing existing output: ',q)}
 dir.create(dirname(q),recursive=TRUE,showWarnings=FALSE);.neutro$owned<-unique(c(.neutro$owned,q));q
}
neutro_write_csv <- function(x,file,...){utils::write.csv(x,file=neutro_output(file),...)}
neutro_write_table <- function(x,file,...){utils::write.table(x,file=neutro_output(file),...)}
neutro_saveRDS <- function(object,file,...){base::saveRDS(object,file=neutro_output(file),...)}
neutro_save <- function(...,file){base::save(...,file=neutro_output(file))}
neutro_writeLines <- function(text,con,...){base::writeLines(text,con=neutro_output(con),...)}
neutro_cat <- function(...,file='',append=FALSE){if(identical(file,''))return(base::cat(...));base::cat(...,file=neutro_output(file,append),append=append)}
neutro_capture_output <- function(...,file=NULL,append=FALSE,type=c('output','message'),split=FALSE){if(is.null(file))return(utils::capture.output(...));utils::capture.output(...,file=neutro_output(file,append),append=append,type=type,split=split)}
neutro_sink <- function(file=NULL,append=FALSE,type=c('output','message'),split=FALSE){
 if(is.null(file))return(base::sink(type=type,split=split))
 if(inherits(file,'connection')){id<-as.character(as.integer(file));known<-.neutro$connections[id];if(is.na(known)||!identical(summary(file)$description,unname(known)))neutro_fail('Unmanaged sink connection forbidden');return(base::sink(file=file,append=append,type=type,split=split))}
 base::sink(file=neutro_output(file,append),append=append,type=type,split=split)
}
neutro_file <- function(description='',open='',...){if(grepl('[wa+]',open)){p<-neutro_output(description,grepl('a',open));con<-base::file(p,open=open,...);.neutro$connections[as.character(as.integer(con))]<-p;return(con)};base::file(neutro_input(description),open=open,...)}
neutro_readRDS <- function(file,...)base::readRDS(neutro_input(file),...)
neutro_read_csv <- function(file,...)utils::read.csv(neutro_input(file),...)
neutro_readLines <- function(con,...)base::readLines(neutro_input(con),...)
neutro_read_xpt <- function(file,...)haven::read_xpt(neutro_input(file),...)
neutro_fread <- function(input,...)data.table::fread(neutro_input(input),...)
neutro_load <- function(file,...)base::load(neutro_input(file),...)
neutro_list_files <- function(path='.',...)base::list.files(neutro_input(path),...)
neutro_normalize <- function(path,...)base::normalizePath(neutro_input(path),...)
neutro_dir_create <- function(path,...){q<-neutro_real(path);if(q!=.neutro$out&&!startsWith(q,paste0(.neutro$out,'/')))neutro_fail('Directory escapes output root');if(dir.exists(q))return(invisible(FALSE));neutro_output(file.path(path,'.directory_guard'));base::dir.create(q,...)}
neutro_device <- function(factory,filename,...){
 p<-neutro_output(filename);old<-grDevices::dev.cur();ok<-FALSE
 on.exit({if(!ok&&grDevices::dev.cur()!=old)try(grDevices::dev.off(),silent=TRUE)},add=TRUE)
 withCallingHandlers(factory(p,...),warning=function(w)neutro_fail('Graphics initialization warning: ',conditionMessage(w)))
 now<-grDevices::dev.cur();if(now==old||now==1L)neutro_fail('Graphics device did not initialize')
 .neutro$devices<-c(.neutro$devices,unname(now));ok<-TRUE;invisible(now)
}
neutro_png <- function(filename,...)neutro_device(grDevices::png,filename,...)
neutro_pdf <- function(file,...)neutro_device(grDevices::pdf,file,...)
neutro_tiff <- function(filename,...)neutro_device(grDevices::tiff,filename,...)
neutro_jpeg <- function(filename,...)neutro_device(grDevices::jpeg,filename,...)
neutro_svg <- function(filename,...)neutro_device(grDevices::svg,filename,...)
neutro_dev_off <- function(which=grDevices::dev.cur()){if(!unname(which)%in%.neutro$devices)neutro_fail('Unmanaged graphics device close');grDevices::dev.off(which)}
neutro_ggsave <- function(filename,plot=ggplot2::last_plot(),...){dots<-list(...);if('path'%in%names(dots))neutro_fail('ggsave path override forbidden; use explicit filename');p<-neutro_output(filename);withCallingHandlers(ggplot2::ggsave(filename=p,plot=plot,...),warning=function(w)neutro_fail('ggsave warning: ',conditionMessage(w)));if(!base::file.exists(p)||base::file.info(p)$size<=0)neutro_fail('ggsave did not produce explicit output')}
neutro_source <- function(file,...){p<-neutro_real(file);if(!startsWith(p,paste0(.neutro$src,'/')))neutro_fail('Source outside repository');base::source(p,...)}
neutro_file_exists <- function(...)base::file.exists(neutro_input(c(...)))
neutro_file_info <- function(...)base::file.info(neutro_input(c(...)))
neutro_md5 <- function(files)tools::md5sum(neutro_input(files))
neutro_sha256 <- function(files)tools::sha256sum(neutro_input(files))
neutro_setwd <- function(dir){p<-neutro_real(dir);if(p!=.neutro$out&&!startsWith(p,paste0(.neutro$out,'/')))neutro_fail('CWD outside managed output');base::setwd(p)}
neutro_input_manifest <- function(files,expected=NULL){
 files<-unique(files);if(length(files)<1L)neutro_fail('Empty actual input manifest is forbidden')
 if(!is.null(expected)&&length(files)!=expected)neutro_fail('Actual input count differs from expected module count')
 p<-vapply(files,function(x)normalizePath(neutro_input(x),mustWork=TRUE),character(1));roots<-c(.neutro$inp,unlist(.neutro$config$external_inputs,use.names=FALSE))
 if(!all(vapply(p,function(x)any(vapply(roots,function(r)x==r||startsWith(x,paste0(r,'/')),logical(1))),logical(1))))neutro_fail('Manifest input is outside declared input roots')
 h<-unname(tools::sha256sum(p));if(anyNA(h)||any(!nzchar(h)))neutro_fail('Missing input SHA256')
 data.frame(file=basename(p),path=p,size_bytes=base::file.info(p)$size,sha256=h)
}
neutro_domain_entry <- function(entry,env=parent.frame()){
 # Preserve the approved Stage5 launcher behavior for direct R entrypoints too.
 b<-basename(entry);needed<-any(startsWith(b,c('06_run_B','07_run_B','08_run_B','09_run_B','dependency_run_B')))||identical(Sys.getenv('NEUTRO_FORCE_DOMAIN'),'1')
 if(!needed)return(invisible(FALSE))
 suppressPackageStartupMessages(library(survey))
 assign('full',neutro_readRDS('work/cohort/C2_full_adult_design.rds'),envir=env)
 base::sys.source(file.path(.neutro$src,'scripts/lib/domain_design.R'),envir=env)
 assign('svydesign',get('domain_svydesign',envir=env,inherits=FALSE),envir=env)
 invisible(TRUE)
}
neutro_tabix_table <- function(tabixFile,...)seqminer::tabix.read.table(tabixFile=neutro_input(tabixFile),...)
