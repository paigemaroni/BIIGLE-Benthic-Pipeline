#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-.}"
python3 - "$ROOT" <<'PY2'
import csv, re, sys
from pathlib import Path
from collections import Counter, defaultdict
root=Path(sys.argv[1]).resolve()
sheets=root/'Sheets'; out=root/'Analyses'/'00_Quality_Control'; out.mkdir(parents=True,exist_ok=True)
req={
'combined_biigle_annotations.csv':['dive_id','label_name','label_hierarchy','filename','shape_name','annotation_id'],
'gl_latlong.csv':['FullID','Lat','Long'],
'depth_temp_frameid.csv':['FrameID','video_time','depth_m','temperature'],
'taxon_list.csv':['Layer_01','cpc_codes','kingdom','phylum','class','taxonomic_resolution','common_id_short','common_id_mid','common_id_full','releif','substrate','type','size'],
'vme_taxon_list.csv':['Layer_01','cpc_codes','kingdom','phylum','class','taxonomic_resolution','common_id_short','common_id_mid','common_id_full'],
}
def n(s): return str(s or '').replace('\ufeff','').strip().casefold()
def norm_id(s): return re.sub(r'[-_]+','-',str(s or '').strip().casefold())
report=[]
for fn,cols in req.items():
 p=sheets/fn
 if not p.exists(): report.append((fn,'ERROR','missing file')); continue
 with p.open(newline='',encoding='utf-8-sig') as f:
  r=csv.reader(f); h=next(r,[]); rows=list(r)
 nh={n(x) for x in h}; miss=[c for c in cols if n(c) not in nh]
 report.append((fn,'PASS' if not miss else 'ERROR',f'{len(rows):,} rows; missing required columns: {", ".join(miss) if miss else "none"}'))
# Join coverage
main=sheets/'combined_biigle_annotations.csv'
if main.exists():
 with main.open(newline='',encoding='utf-8-sig') as f: data=list(csv.DictReader(f))
 frames={r['filename'].strip() for r in data if r.get('filename','').strip()}; dives={norm_id(r.get('dive_id')) for r in data if r.get('dive_id')}
 lat=sheets/'gl_latlong.csv'
 if lat.exists():
  with lat.open(newline='',encoding='utf-8-sig') as f: ids={norm_id(r.get('FullID')) for r in csv.DictReader(f) if r.get('FullID')}
  unmatched=sorted(dives-ids); report.append(('JOIN dive_id ↔ FullID','PASS' if not unmatched else 'WARNING',f'{len(dives)-len(unmatched)}/{len(dives)} unique dives matched'))
  (out/'00_unmatched_dive_ids.txt').write_text('\n'.join(unmatched),encoding='utf-8')
 meta=sheets/'depth_temp_frameid.csv'
 if meta.exists():
  with meta.open(newline='',encoding='utf-8-sig') as f: mf={str(r.get('FrameID','')).strip().casefold() for r in csv.DictReader(f) if r.get('FrameID')}
  unmatched=sorted(x for x in frames if x.casefold() not in mf); report.append(('JOIN filename ↔ FrameID','PASS' if not unmatched else 'WARNING',f'{len(frames)-len(unmatched)}/{len(frames)} unique BIIGLE frames matched'))
  with (out/'00_unmatched_frame_ids.csv').open('w',newline='',encoding='utf-8') as h:
   w=csv.writer(h); w.writerow(['filename']); w.writerows([[x] for x in unmatched])
with (out/'00_input_validation_summary.csv').open('w',newline='',encoding='utf-8') as h:
 w=csv.writer(h); w.writerow(['check','status','details']); w.writerows(report)
print('BIIGLE input validation')
print('=======================')
for a,b,c in report: print(f'{b:7}  {a}: {c}')
PY2
