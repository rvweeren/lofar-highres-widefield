import sys

import GRID_LRT

from cloudant.client import CouchDB
from GRID_LRT.application import submit
from GRID_LRT.auth.get_picas_credentials import picas_cred
from GRID_LRT.token import caToken
from GRID_LRT.token import TokenJsonBuilder
from GRID_LRT.token import TokenList
from GRID_LRT.token import TokenView
from GRID_LRT.Staging import srmlist

print('Using GRID_LRT from ' + GRID_LRT.__file__)

srms = sys.argv[1]
obsid = sys.argv[2]
cal_obsid = sys.argv[3]
tok_type = 'AAA_prefactor_' + obsid
jsonfile = '/project/sksp/Software/lofar-highres-widefield/testdir/test_with_GRID_LRT/lb_step4_1asec.json'

print('1) Initializing PiCaS credentials and connect to database.')
pc = picas_cred()
client = CouchDB(pc.user, pc.password, url='https://picas-lofar.grid.surfsara.nl:6984', connect=True)
db = client[pc.database]
print('Connected to database ' + pc.database)

print('2) Creating a list of paths to the files.')
s = srmlist.srmlist()
with open(srms, 'r') as f:
    for l in f.readlines():
        s.append(l.strip())

print('3) Slicing list into groups.')
# g = s.sbn_dict(pref='SB', suff='_')
# Temporary hack
g = s.sbn_dict(pref='t_', suff='MHz')
d = srmlist.slice_dicts(g, 999)
print(d)

print('4) Building token list.')
tl = TokenList(token_type=tok_type, database=db)
for k, v in d.items():
    tok = caToken(database=db, token_type=tok_type, token_id=tok_type + '_1asec' + str(cal_obsid) + '_SB' + str(k))
    with open('temp_srm.txt', 'w') as f:
        f.write('\n'.join(v))
    tok.build(TokenJsonBuilder(jsonfile))
    tok.save()
    tok.add_attachment(attachment_name='srm.txt', filename='temp_srm.txt')
    tok.add_attachment(attachment_name='step4_1asec_image.parset', filename='/project/sksp/Software/lofar-highres-widefield/testdir/test_with_GRID_LRT/step4_1asec_image.parset')
    tl.append(tok)

tl.save()

for tok in tl:
    tok['OBSID'] = obsid
    tok['PIPELINE_STEP'] = 'lb_1asec_cal1'
    tok['status'] = 'queued'
    tok.save()

print('5) Adding status views.')
view_cal = TokenView('step4_1asec_image', condition='doc.PIPELINE_STEP=="lb_1asec_cal1"', emit_values=('doc._id', 'doc.status'))
tl.add_token_views()
tl.add_view(view_cal)

tl.save()

print('6) Create and launch the jobs.')
j = submit.SpiderLauncher(numjobs=len(d.keys()), token_type=tok_type, wholenode=False, parameter_step=999, NCPU=24)

with j:
    url = j.launch()
print('Job ID: ' + str(url))
