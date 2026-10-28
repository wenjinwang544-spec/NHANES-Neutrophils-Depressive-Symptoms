"""Path contract for managed analysis outputs. No marker grants write permission."""
from pathlib import Path
import os,json,shutil,hashlib
class SafetyError(RuntimeError):pass

def overlaps(a,b):return a==b or a in b.parents or b in a.parents

def real(p):return Path(p).expanduser().resolve()

def disjoint(output,roots):
 for r in roots:
  if overlaps(real(output),real(r)):raise SafetyError(f'Output/input/source/protected paths overlap: {output} and {r}')

def no_links(p,root):
 p=Path(p);root=real(root)
 if '..'in p.parts:raise SafetyError('Parent traversal forbidden')
 if real(p)!=root and root not in real(p).parents:raise SafetyError('Resolved path escapes output root')
 rel=p.absolute().relative_to(root)
 c=root
 for part in rel.parts:
  c=c/part
  if c.is_symlink():raise SafetyError('Symlinks are forbidden in the output tree')
 return p

def storage_roots(source):
 p=real(source)/"config/storage_policy.json"
 if not p.exists():return []
 policy=json.loads(p.read_text())
 if policy.get("schema")!=1 or policy.get("protect_source_parent") is not True:raise SafetyError("Invalid archive storage policy")
 return [real(source).parent]

def load(workspace):
 w=real(workspace);c=json.loads((w/'run_config.json').read_text())
 if c.get('schema')!=1 or real(c['workspace'])!=w:raise SafetyError('Invalid initialization config')
 out=real(c['output_root']);inp=real(c['input_root']);src=real(c['source_root']);prot=[real(x)for x in c['protected_roots']]+storage_roots(src)
 if out!=w/'output' or inp!=w/'inputs' or not prot:raise SafetyError('Invalid fixed layout or missing protected roots')
 disjoint(out,[inp,src,*prot]);disjoint(inp,[src,*prot])
 if (w/'output').is_symlink()or(w/'inputs').is_symlink():raise SafetyError('Managed roots may not be symlinks')
 if not out.is_dir()or not inp.is_dir():raise SafetyError('Missing initialized roots')
 for r in c.get('external_inputs',{}).values():disjoint(out,[r])
 for item in c.get('imports',[]):
  f=inp/'gwas_cache'/item['name']
  if not f.is_file()or f.is_symlink()or not f.resolve().is_relative_to(inp):raise SafetyError('Imported input escaped its input root')
  if hashlib.sha256(f.read_bytes()).hexdigest()!=item['sha256']:raise SafetyError('Imported input hash changed')
 return c

def output(c,path,append=False,owned=None):
 out=real(c['output_root']);p=Path(path)
 if not p.is_absolute():p=out/p
 # Normalize lexical '..' and resolve every existing parent before creating anything.
 rp=real(p);no_links(p,out)
 if (out/'data')==rp or (out/'data')in rp.parents:raise SafetyError('Logical data namespace is read-only')
 disjoint(rp,[c['input_root'],c['source_root'],*c['protected_roots'],*c.get('external_inputs',{}).values()])
 if rp.exists() and not(append and owned is not None and str(rp)in owned):raise SafetyError(f'Refusing existing output: {rp}')
 rp.parent.mkdir(parents=True,exist_ok=True);no_links(rp,out)
 return rp

def create_workspace(workspace,source,protected,raw=None,caches=()):
 w=real(workspace);src=real(source);prot=[real(x)for x in protected]+storage_roots(src)
 if not prot:raise SafetyError('At least one --protected-root is required')
 roots=[src,*prot]+([real(raw)]if raw else[])+[real(v)for _,v in caches]
 for r in roots:
  if not r.exists():raise SafetyError(f'Input/protected/source root missing: {r}')
 disjoint(w,roots)
 if w.exists():raise SafetyError('Initialization requires a nonexistent workspace directory')
 # Resolve symlinked ancestors before creation, then validate again.
 w.mkdir(parents=True);out=w/'output';inp=w/'inputs';out.mkdir();inp.mkdir()
 for sub in WORK_DIRS:(out/'work'/sub).mkdir(parents=True,exist_ok=True)
 (out/'logs').mkdir();(out/'figures').mkdir()
 external={}
 if raw:external['nhanes']=str(real(raw))
 imports=[]
 for name,value in caches:
  if Path(name).name!=name or name in ['', '.', '..']:raise SafetyError('Cache name must be a plain filename')
  srcfile=real(value)
  if not srcfile.is_file():raise SafetyError('Imported cache must be a regular file')
  dst=inp/'gwas_cache'/name;dst.parent.mkdir(exist_ok=True)
  with srcfile.open('rb')as r,dst.open('xb')as f:shutil.copyfileobj(r,f)
  imports.append({'name':name,'sha256':hashlib.sha256(dst.read_bytes()).hexdigest(),'bytes':dst.stat().st_size})
 c={'schema':1,'workspace':str(w),'source_root':str(src),'output_root':str(out),'input_root':str(inp),'protected_roots':[str(x)for x in prot],'external_inputs':external,'imports':imports}
 with(w/'run_config.json').open('x')as f:json.dump(c,f,indent=2)
 load(w)
 return c

WORK_DIRS=['cohort','clinical','analytic','imputation','mi_validation','symptoms_base','symptoms_extended','thresholds','alcohol','education','inflammation','descriptive','export','covariates','discrimination','shapley','heterogeneity','attenuation','race','sleep','spline','sex','sex_total','mr_export','mr_online','mr_refit','coloc']
