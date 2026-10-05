# Curated historical source: publication does not claim complete runnable replay.
# See PUBLICATION_SCOPE.json and historical_methods/README.md; omitted dependencies are intentional.
#!/usr/bin/env python3
"""Seal after the audit process exits, including its completed output log."""
import engine as run
import json
from pathlib import Path
ROOT=run.ROOT
def main():
    report=json.loads((ROOT/'RESULT.json').read_text());assert report['status']=='PASS'
    assert json.loads((ROOT/'completion.json').read_text())['completed']==128
    frozen={}
    for name in ('rattler_contact_transfer','best147_boundary_transfer','next146_boundary_transfer'):
        manifest=json.loads((ROOT.parent/name/'MANIFEST_SHA256.json').read_text());frozen[name]=dict(count=len(manifest),all_unchanged=all(run.sha(ROOT.parent/name/p)==h for p,h in manifest.items()))
    assert all(v['all_unchanged'] for v in frozen.values())
    run.write(ROOT/'freeze_checks.json',dict(original_three_trees=frozen,baseline_file_sha256=run.sha(ROOT/'inputs/baseline.npy'),audit_completed_status=report['status'],audit_log_file_sha256=run.sha(ROOT/'audit.log'),run_log_file_sha256=run.sha(ROOT/'run.log'),source_file_sha256={str(p.relative_to(ROOT)):run.sha(p) for p in sorted((ROOT/'source').glob('*.py'))},scope='Final file-integrity seal after optimizer and independent audit exit; no further mathematical claims.'))
    run.write(ROOT/'status.json',dict(stage='finished and frozen',optimized=128,audit_status=report['status'],material_improvements=report['material_candidate_count'],best_final_radius=report['best_final_radius'],best_full100_baseline_match=report['best_matches_baseline_full100_geometry_and_inventories'],original_frozen_trees=frozen,budget=128))
    manifest={str(p.relative_to(ROOT)):run.sha(p) for p in sorted(ROOT.rglob('*')) if p.is_file() and p.name!='MANIFEST_SHA256.json'};run.write(ROOT/'MANIFEST_SHA256.json',manifest)
    print(json.dumps(dict(files=len(manifest),status='FROZEN',best=report['best_final_radius'],material=report['material_candidate_count'],original_three_unchanged=True)))
if __name__=='__main__':main()
