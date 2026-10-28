from pathlib import Path
import argparse,subprocess,sys,json
sys.dont_write_bytecode = True
sys.path.insert(0,str(Path(__file__).resolve().parent/'lib'))
from safety import load,output
PRE=['01_build_cohort.R','02_build_full_adult_design.R','03_add_public_modules.R','04_clinical_and_analyses.R','15_run_D2_1.R']
POST=['06_prepare_supplement_inputs.R','06_run_B1_2_to_6.R','dependency_run_B2_1_and_audit.R','07_run_B2_remaining.R','08_run_B3_and_ALQ_audit.R','09_run_B4.R','dependency_run_B5.R','dependency_run_B6.R','16_run_D2_3.R','17_run_E1.R','dependency_check_flow.R','21_export_final.R','22_export_remaining.R','19_draw_figures.R','20_draw_forest_publication.R','25_run_G1.R','27_run_G2.R','26_finalize_G2_v2.R','28_run_G3.R','29_run_G3_remaining.R','30_run_G3b.R','32_run_H10_v3.R','36_run_phaseH7.R','37_prepare_phaseH7b.R','38_run_phaseH7b.R','05_joint_splines.R','sex_interactions.R','sex_total_score.R']
p=argparse.ArgumentParser();p.add_argument('--workspace',required=True);p.add_argument('--rscript',default='Rscript');p.add_argument('--frozen-imputations',action='store_true');p.add_argument('--plan',action='store_true');a=p.parse_args();repo=Path(__file__).resolve().parents[1];c=load(a.workspace);w=Path(c['output_root'])
mi=['verify_primary_from_frozen_imputations.R']if a.frozen_imputations else['14_run_sequence.R']
plan=[*PRE,*mi,*POST,'export_revised_tables.py','render_figure3.R']
if a.plan:
 print(json.dumps({'plan':plan,'models_executed':False,'frozen_imputation_tier':a.frozen_imputations},indent=2));raise SystemExit(0)
def run(script,extra=()):subprocess.run([sys.executable,str(repo/'scripts/run_job.py'),'scripts/nhanes/'+script,'--workspace',a.workspace,'--rscript',a.rscript,*extra],check=True)
for s in PRE:run(s)
for s in mi:run(s)
# Frozen branch reads imported inputs directly; copy only newly generated validation exports into named downstream outputs.
if a.frozen_imputations:
 for name in ['Rubin_MI_CC_comparison.csv','imputation_scaling.csv']:
  src=w/'work/mi_validation'/name;dst=output(c,'work/imputation/'+name)
  with src.open('rb')as r,dst.open('xb')as f:f.write(r.read())
for s in POST:run(s)
subprocess.run([sys.executable,str(repo/'scripts/export_revised_tables.py'),'--workspace',a.workspace],check=True)
rug=output(c,'work/spline/Figure3_distribution.csv')
with (repo/'supplementary/validated/Figure3_distribution.csv').open('rb')as r,rug.open('xb')as f:f.write(r.read())
subprocess.run([sys.executable,str(repo/'scripts/run_job.py'),'scripts/render_figure3.R','--workspace',a.workspace,'--rscript',a.rscript,'--curve-input',str(w/'work/spline')],check=True)
print('NHANES steps including final S5/S8 display exports completed. MR upstream remains separate and UNVALIDATED.')
