# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Permanent file seal after87 metric restorations and full independent replay."""
import engine as run
import metric_search as search
import json
ROOT=run.ROOT

def main():
    result=json.loads((ROOT/'RESULT.json').read_text());assert result['status']=='PASS';assert result['optimized']==87
    completed=json.loads((ROOT/'completion.json').read_text());assert completed['completed']==87;assert all(x==0 for x in completed['worker_exitcodes'])
    frozen=search.frozen_checks();assert all(v['unchanged']for v in frozen.values());assert frozen['catalog_handoff']['owner_unchanged']
    m=json.loads((ROOT/'inputs/history_snapshot/MANIFEST_SHA256.json').read_text());assert all(run.sha(ROOT/'inputs/history_snapshot'/p)==h for p,h in m.items())
    run.write(ROOT/'freeze_checks.json',dict(old_ten_frozen_trees_and_catalog_handoff=frozen,run_log_sha256=run.sha(ROOT/'run.log'),audit_log_sha256=run.sha(ROOT/'audit.log'),source_file_sha256={str(p.relative_to(ROOT)):run.sha(p)for p in sorted((ROOT/'source').glob('*.py'))},audit_status=result['status'],scope='Post-exit immutable file seal; no new mathematical claim.'))
    run.write(ROOT/'status.json',dict(stage='finished and permanently frozen',budget=96,optimized=87,skipped=9,audit_status=result['status'],material_improvements=result['material_candidate_count'],best_radius=result['best_radius'],best_independent_index=result['best_independent_index'],independent_oracle_calls=result['independent_oracle_calls'],old_frozen_trees=frozen))
    manifest={str(p.relative_to(ROOT)):run.sha(p)for p in sorted(ROOT.rglob('*'))if p.is_file()and p!=ROOT/'MANIFEST_SHA256.json'};run.write(ROOT/'MANIFEST_SHA256.json',manifest)
    assert all(run.sha(ROOT/p)==h for p,h in manifest.items());print(json.dumps(dict(files=len(manifest),status='FROZEN',manifest_sha256=run.sha(ROOT/'MANIFEST_SHA256.json'),best=result['best_radius'],material=result['material_candidate_count'])))

if __name__=='__main__':main()
