from pathlib import Path
import argparse,json,sys
sys.dont_write_bytecode = True
sys.path.insert(0,str(Path(__file__).resolve().parent/'lib'))
from safety import create_workspace
p=argparse.ArgumentParser();p.add_argument('--workspace',required=True);p.add_argument('--protected-root',action='append',required=True);p.add_argument('--raw-nhanes');p.add_argument('--import-cache',action='append',default=[]);a=p.parse_args()
caches=[]
for x in a.import_cache:
 if '='not in x:p.error('--import-cache expects filename=/absolute/source/file')
 caches.append(x.split('=',1))
c=create_workspace(a.workspace,Path(__file__).resolve().parents[1],a.protected_root,a.raw_nhanes,caches)
print(json.dumps({'initialized':True,'workspace':c['workspace'],'output_root':c['output_root'],'inputs_are_separate':True}))
