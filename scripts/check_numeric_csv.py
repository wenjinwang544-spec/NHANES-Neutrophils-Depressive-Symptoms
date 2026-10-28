import csv,math,argparse
p=argparse.ArgumentParser();p.add_argument('reference');p.add_argument('candidate');p.add_argument('--tolerance',type=float,default=1e-9);a=p.parse_args()
def read(p):
 with open(p,encoding='utf-8-sig',newline='')as f:return list(csv.DictReader(f))
x=read(a.reference);y=read(a.candidate);assert len(x)==len(y),'Row count differs';count=0;maximum=0
for i,(r,s)in enumerate(zip(x,y),2):
 assert r.keys()==s.keys(),'Column names differ'
 for k,v in r.items():
  if any(z in k.lower()for z in ['source','path','mtime']):continue
  try:u=float(v);w=float(s[k])
  except (ValueError,TypeError):assert v==s[k],(i,k,'non-numeric label differs');continue
  if math.isnan(u)and math.isnan(w):continue
  assert math.isfinite(u)==math.isfinite(w),(i,k,'missingness differs');d=abs(u-w);maximum=max(maximum,d);count+=1;assert d<=a.tolerance,(i,k,u,w)
print({'numeric_cells':count,'max_abs_difference':maximum,'passed':True})
