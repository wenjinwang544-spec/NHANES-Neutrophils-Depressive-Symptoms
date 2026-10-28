"""No-model safety tests. Creates synthetic fixtures only under this archive's docs/validation."""
from pathlib import Path
import sys,os,tempfile,json,subprocess,hashlib,shutil
sys.dont_write_bytecode=True
Q=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(Q/'scripts/lib'))
from safety import create_workspace,load,output,SafetyError
V=Q/'docs/validation';V.mkdir(exist_ok=True)
T=Path(tempfile.mkdtemp(prefix='synthetic_',dir=V))
rows=[]
def check(name,fn,reject=False):
 try:fn();ok=not reject;detail='accepted'
 except Exception as e:ok=reject;detail=type(e).__name__+': '+str(e)
 detail=detail.replace(str(Q),'<ARCHIVE>')
 rows.append({'check':name,'passed':ok,'detail':detail})
 if not ok:raise AssertionError(name+' '+detail)
src=T/'source';src.mkdir();old=T/'protected';old.mkdir();(old/'sentinel').write_text('UNCHANGED');raw=T/'raw';raw.mkdir();toy=raw/'toy.txt';toy.write_text('synthetic')
w=T/'valid';check('synthetic fresh init',lambda:create_workspace(w,src,[old],raw,[('toy.txt',toy)]));c=load(w)
for name,path in [('protected child',old/'child'),('source child',src/'child'),('raw child',raw/'child'),('ancestor',T),('existing workspace',w)]:check(name,lambda p=path:create_workspace(p,src,[old],raw),True)
check('synthetic explicit output',lambda:output(c,'work/toy.csv'))
p=output(c,'work/existing.csv');p.write_text('OLD')
for name,path in [('overwrite','work/existing.csv'),('traversal','../bad.csv'),('protected output',str(old/'bad.csv')),('data namespace','data/gwas_cache/bad.csv')]:check(name,lambda p=path:output(c,p),True)
(w/'output/link').symlink_to(old,target_is_directory=True)
check('symlink escape',lambda:output(c,'link/bad.csv'),True)
# Actual archive policy rejects before creating the requested workspace.
check('actual archive parent auto-protected',lambda:create_workspace(T/'forbidden_archive_workspace',Q,[old]),True)
assert not(T/'forbidden_archive_workspace').exists()
# Test same relative policy with only synthetic roots, including positive initialization outside its parent.
ps=T/'storage_parent/source';(ps/'config').mkdir(parents=True);(ps/'config/storage_policy.json').write_bytes((Q/'config/storage_policy.json').read_bytes())
check('synthetic policy denies parent descendant',lambda:create_workspace(ps.parent/'bad',ps,[old]),True)
pw=T/'outside_synthetic_storage_parent';check('synthetic policy accepts disjoint workspace',lambda:create_workspace(pw,ps,[old]))
# Verify R also accepts the policy with a disjoint synthetic workspace, without study code.
r=subprocess.run([os.environ.get('RSCRIPT','Rscript'),'--vanilla','-e',"a<-commandArgs(TRUE);base::source(file.path(a[1],'scripts/lib/safety.R'));Sys.setenv(NEUTRO_WORKSPACE=a[2]);neutro_setup(a[3]);cat('POLICY_VALID_SYNTHETIC')",str(Q),str(pw),str(ps)],cwd=T,capture_output=True,text=True)
assert r.returncode==0 and 'POLICY_VALID_SYNTHETIC'in r.stdout,r.stderr
rows.append({'check':'R policy permits disjoint synthetic workspace','passed':True})
# Forged initialized config under actual research parent must be denied by both Python and direct R.
fw=T/'forged';fw.mkdir();(fw/'inputs').mkdir();(fw/'output').mkdir()
config={'schema':1,'workspace':str(fw),'source_root':str(Q),'output_root':str(fw/'output'),'input_root':str(fw/'inputs'),'protected_roots':[str(old)],'external_inputs':{},'imports':[]}
(fw/'run_config.json').write_text(json.dumps(config));check('Python load rejects omitted research protection',lambda:load(fw),True)
# Execute exact R safety source against synthetic config; no scientific entry/body is sourced.
rfile=T/'synthetic_guard_tests.R'
rfile.write_text('''a<-commandArgs(TRUE);base::source(file.path(a[1],'scripts/lib/safety.R'));Sys.setenv(NEUTRO_WORKSPACE=a[2]);neutro_setup(a[3]);z<-list();test<-function(n,e,reject=FALSE){err<-tryCatch({force(e);NULL},error=function(e)conditionMessage(e));z[[length(z)+1]]<<-list(check=n,passed=if(reject)!is.null(err)else is.null(err))}
test('R explicit synthetic CSV',neutro_write_csv(data.frame(value=1),'work/good.csv',row.names=FALSE))
test('R no overwrite',neutro_write_csv(data.frame(value=2),'work/good.csv',row.names=FALSE),TRUE)
test('R default plot blocked',plot(1:3),TRUE)
test('R initialization exception blocked',neutro_device(function(...)stop('injected'),'figures/error.png'),TRUE)
test('R no device opened blocked',neutro_device(function(...)invisible(NULL),'figures/no_device.png'),TRUE)
test('R explicit PDF device',{neutro_pdf('figures/toy.pdf');plot(1:3);neutro_dev_off()})
test('R unsafe input write blocked',neutro_writeLines('bad',file.path(.neutro$inp,'bad.txt')),TRUE)
cat(jsonlite::toJSON(z,auto_unbox=TRUE));if(any(!vapply(z,function(x)x$passed,logical(1))))quit(status=1)
''')
env=os.environ.copy();env.pop('NEUTRO_SOURCE_ROOT',None);env.pop('NEUTRO_WORKSPACE',None)
r=subprocess.run([os.environ.get('RSCRIPT','Rscript'),'--vanilla',str(rfile),str(Q),str(w),str(src)],cwd=T,env=env,capture_output=True,text=True)
(T/'R_synthetic.log').write_text((r.stdout+r.stderr).replace(str(Q),'<ARCHIVE>'));assert r.returncode==0,r.stderr
rows.extend(json.loads(r.stdout))
r=subprocess.run([os.environ.get('RSCRIPT','Rscript'),'--vanilla',str(Q/'scripts/nhanes/05_joint_splines.R'),'--workspace',str(fw),'--check-only'],cwd=T,env=env,capture_output=True,text=True)
(T/'direct_R_policy.log').write_text((r.stdout+r.stderr).replace(str(Q),'<ARCHIVE>'))
assert r.returncode!=0 and 'NEUTRO_SAFETY'in r.stderr,(r.returncode,r.stdout,r.stderr)
rows.append({'check':'direct R rejects forged workspace under research parent','passed':True,'detail':r.stderr.replace(str(Q),'<ARCHIVE>')})
# No study body is entered, even when direct entry is used without initialization.
for p in [*list((Q/'scripts/nhanes').glob('*.R')),*list((Q/'scripts/mr').glob('*.R')),Q/'scripts/render_figure3.R']:
 r=subprocess.run([os.environ.get('RSCRIPT','Rscript'),'--vanilla',str(p)],cwd=T,env=env,capture_output=True,text=True)
 assert r.returncode!=0 and 'NEUTRO_SAFETY'in r.stderr,p
 rows.append({'check':'uninitialized direct R rejected '+str(p.relative_to(Q)),'passed':True,'detail':'guard stops before study code'})
assert not list(T.rglob('Rplots.pdf')) and (old/'sentinel').read_text()=='UNCHANGED' and toy.read_text()=='synthetic'
rows.append({'check':'no Rplots; protected synthetic files unchanged','passed':True})
# Configuration values contain local paths, so preserve only a sanitized test record publicly.
# Fixture files are audit-only and excluded by .gitignore; they never contain study data.
with (V/(T.name+'_results.json')).open('x') as f:json.dump(rows,f,indent=2)
shutil.rmtree(T) # Only the synthetic directory created by this invocation; never study inputs or source.
print(json.dumps({'tests':len(rows),'passed':sum(r['passed']for r in rows),'models_run':0,'synthetic_fixtures_removed':True},indent=2))
