from pathlib import Path
import argparse,subprocess,os,json,sys
sys.dont_write_bytecode = True
sys.path.insert(0,str(Path(__file__).resolve().parent/'lib'))
from safety import load,output,real,SafetyError
p=argparse.ArgumentParser();p.add_argument('script');p.add_argument('--workspace',required=True);p.add_argument('--rscript',default='Rscript');p.add_argument('--full-domain',action='store_true');p.add_argument('--check-only',action='store_true');p.add_argument('--curve-input');a=p.parse_args()
repo=Path(__file__).resolve().parents[1];script=(repo/a.script).resolve()
if not script.is_relative_to(repo/'scripts')or script.suffix!='.R':raise SafetyError('Study script must be within source scripts/')
c=load(a.workspace)
if real(c['source_root'])!=repo:raise SafetyError('Source root differs from initialized source')
if a.check_only:
 print(json.dumps({'script':a.script,'path_contract':'validated','executes_study_code':False}));raise SystemExit(0)
key=script.stem;log=output(c,'logs/'+key+'.log');driver=output(c,'logs/'+key+'_driver.R')
# A new driver and log are exclusive creations; study code is never copied into an input directory.
pre='repo<-'+json.dumps(str(repo))+'\n'
pre+='base::source('+json.dumps(str(script))+',local=TRUE)\n'
with driver.open('x')as f:f.write(pre)
env=os.environ.copy();env.update(NEUTRO_WORKSPACE=str(real(a.workspace)),NEUTRO_SOURCE_ROOT=str(repo));env.pop('NEUTRO_FORCE_DOMAIN',None)
if a.full_domain:env['NEUTRO_FORCE_DOMAIN']='1'
with log.open('x')as f:r=subprocess.run([a.rscript,'--vanilla',str(driver),*(['--curve-input',a.curve_input]if a.curve_input else[])],cwd=c['output_root'],env=env,stdout=f,stderr=subprocess.STDOUT)
print(json.dumps({'job':key,'exit_code':r.returncode,'log':str(log)}));raise SystemExit(r.returncode)
