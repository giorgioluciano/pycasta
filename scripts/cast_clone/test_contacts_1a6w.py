"""Optional real-data regression test. Explicitly supply the two JSON paths."""
import argparse,json
from pathlib import Path
from cast_contact_audit import classify_contacts
p=argparse.ArgumentParser();p.add_argument('geometry',type=Path);p.add_argument('topology',type=Path)
a=p.parse_args()
r=classify_contacts(json.loads(a.geometry.read_text()),json.loads(a.topology.read_text()))
expected={'global_pre_alpha':{'removed_alpha_edge':22,'cross_region_contact':3},
          'region_raw_boundary':{'removed_alpha_edge':221}}
if r['counts']!=expected:raise AssertionError(('Counts changed',r['counts']))
if r['unresolved_contact_count'] or r['included_mouth_edge_failure_count']:raise AssertionError('New failures')
actual={tuple(x['edge']):x['regions'] for x in r['records'] if x['category']=='cross_region_contact'}
if actual!={(935,1764):[0,4],(1315,1327):[0,9],(1140,1163):[0,11]}:raise AssertionError(actual)
print('PASS: fixed 1a6w contact regression; not certification of original CAST.')
