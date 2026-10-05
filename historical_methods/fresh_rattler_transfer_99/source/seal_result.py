# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Final seal after optimizer/audit exit; no prior frozen tree is modified."""
import engine as run
import rattler_search as search
from pathlib import Path
import json
ROOT=run.ROOT
def main():
    result=json.loads((ROOT/'RESULT.json').read_text());assert result['status']=='PASS'
    assert json.loads((ROOT/'completion.json').read_text())['completed']==96
    frozen=search.frozen_checks();assert all(v['unchanged']for v in frozen.values())
    run.write(ROOT/'freeze_checks.json',dict(original_six_frozen_trees=frozen,source_input_sha256={f'source{i:02d}':run.sha(ROOT/'inputs'/f'source{i:02d}.npy')for i in range(3)},run_log_sha256=run.sha(ROOT/'run.log'),audit_log_sha256=run.sha(ROOT/'audit.log'),source_file_sha256={str(p.relative_to(ROOT)):run.sha(p)for p in sorted((ROOT/'source').glob('*.py'))},audit_status=result['status'],scope='Post-exit file integrity seal, not a new mathematical proof.'))
    run.write(ROOT/'status.json',dict(stage='finished and permanently frozen',optimized=96,budget=96,audit_status=result['status'],material_improvements=result['material_candidate_count'],best_final_radius=result['best_final_radius'],best_known_full100_matches=result['best_known_full100_matches'],all_triple_oracle_calls=result['all_triple_oracle_calls'],frozen_trees=frozen))
    mf={str(p.relative_to(ROOT)):run.sha(p)for p in sorted(ROOT.rglob('*'))if p.is_file()and p.name!='MANIFEST_SHA256.json'};run.write(ROOT/'MANIFEST_SHA256.json',mf);print(json.dumps(dict(files=len(mf),status='FROZEN',best=result['best_final_radius'],material=result['material_candidate_count'],old_six_unchanged=True)))
if __name__=='__main__':main()
