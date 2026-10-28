from pathlib import Path
from decimal import Decimal,ROUND_HALF_UP
import csv,argparse,sys
sys.dont_write_bytecode = True
sys.path.insert(0,str(Path(__file__).resolve().parent/"lib"))
from safety import load,output
p=argparse.ArgumentParser();p.add_argument('--workspace',required=True);a=p.parse_args();config=load(a.workspace);w=Path(config['output_root']);out=output(config,'tables_revised');out.mkdir(exist_ok=False)
def read(p):
 with p.open(encoding='utf-8-sig',newline='')as f:return list(csv.DictReader(f))
def fixed(x,n=3):return format(Decimal(str(x)).quantize(Decimal(1).scaleb(-n),rounding=ROUND_HALF_UP),f'.{n}f')
def pval(x):
 x=Decimal(str(x));return '<0.001'if x<Decimal('.001')else format(x.quantize(Decimal(1).scaleb(x.adjusted()-2),rounding=ROUND_HALF_UP),'f')
def write(name,rows):
 with (out/name).open('x',newline='')as f:z=csv.DictWriter(f,fieldnames=list(rows[0]));z.writeheader();z.writerows(rows)
scalars=read(w/'work/sex/interaction_hypotheses_comparison.csv');total=read(w/'work/sex_total/total_score_interactions.csv');glob=read(w/'work/sex/global_D1_tests.csv');tg=read(w/'work/sex_total/total_score_global_D1.csv')
rows=[]
for x in scalars+total:
 if x['type']!='full':continue
 rows.append(dict(model=x['model'],outcome=x.get('outcome','PHQ10'),N='33,174',beta=fixed(x['beta']),SE=fixed(x['SE']),CI=fixed(x['lower'])+'–'+fixed(x['upper']),P=pval(x['P'])))
write('S5_neutrophil_interaction.csv',rows)
rows=[dict(model=x['model'],outcome='PHQ10',N='33,174',F=fixed(x['F']),df1=x['df1'],df2=fixed(x['df2'],1),P=pval(x['P']))for x in glob];rows.append(dict(model='Model1',outcome='PHQ9_total',N='33,174',F=fixed(tg[0]['F.value']),df1=tg[0]['df1'],df2=fixed(tg[0]['df2'],1),P=pval(tg[0]['P(>F)'])));write('S5_omnibus.csv',rows)
x=read(w/'work/spline/candidate_nonlinearity.csv');write('S8_nonlinearity.csv',[dict(model=r['model'],outcome=r['outcome'],N=f"{int(r['N']):,}",F=fixed(r['Wald_F']),df1=r['df1'],df2=r['df2'],P=pval(r['P_nonlinear']))for r in x])
x=read(w/'work/spline/candidate_curves.csv');write('S8_reference_points.csv',[dict(model=r['model'],outcome=r['outcome'],N=f"{int(r['N']):,}",reference=fixed(r['reference_cells_uL'],0),count_10e9_L=fixed(r['Neutrophil_10e9_L'],1),effect=fixed(r['effect']),lower=fixed(r['lower']),upper=fixed(r['upper']))for r in x if float(r['Neutrophil_cells_uL'])in [2000,3000,4000,5000,7000,9000]])
print('Revised S5/S8 display CSVs exported using direct rounding')
