# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Final permanent file seal after48 restorations and independent replay."""
import engine as run
import cavity_search as search
from pathlib import Path
import json
ROOT=run.ROOT
def main():
    result=json.loads((ROOT/'RESULT.json').read_text());assert result['status']=='PASS'
    assert json.loads((ROOT/'completion.json').read_text())['completed']==48
    frozen=search.frozen_checks();assert all(v['unchanged']for v in frozen.values())
    run.write(ROOT/'freeze_checks.json',dict(original_eight_frozen_trees=frozen,run_log_sha256=run.sha(ROOT/'run.log'),audit_log_sha256=run.sha(ROOT/'audit.log'),source_file_sha256={str(p.relative_to(ROOT)):run.sha(p)for p in sorted((ROOT/'source').glob('*.py'))},audit_status=result['status'],scope='Post-exit immutable file seal; no new mathematical acceptance claim.'))
    run.write(ROOT/'status.json',dict(stage='finished and permanently frozen',budget=64,optimized=48,audit_status=result['status'],material_improvements=result['material_candidate_count'],best_final_radius=result['best_final_radius'],best_matches_current146_full100=result['best_matches_current146_full100_geometry_and_inventories'],all_triple_oracle_calls=result['all_triple_oracle_calls'],old_frozen_trees=frozen))
    mf={str(p.relative_to(ROOT)):run.sha(p)for p in sorted(ROOT.rglob('*'))if p.is_file()and p.name!='MANIFEST_SHA256.json'};run.write(ROOT/'MANIFEST_SHA256.json',mf);print(json.dumps(dict(files=len(mf),status='FROZEN',best=result['best_final_radius'],material=result['material_candidate_count'],old_eight_unchanged=True)))
if __name__=='__main__':main()
