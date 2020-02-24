import os
import re
import sys

from cloudant.client import CouchDB
from GRID_LRT.auth.get_picas_credentials import picas_cred
from GRID_LRT.token import caToken
from GRID_LRT.token import TokenJsonBuilder 
from GRID_LRT.token import TokenList
from GRID_LRT.token import TokenView

import GRID_LRT
print('Using GRID_LRT from ' + GRID_LRT.__file__)

srms = sys.argv[1]
obsid = sys.argv[2]
cal_obsid = sys.argv[3]
tok_type = 'AAA_prefactor_' + obsid
jsonfile = '/project/sksp/Software/lofar-highres-widefield/testdir/test_with_GRID_LRT/lb_step3_selfcal.json'

print('1) Initializing PiCaS credentials and connect to database.')
pc=picas_cred()       
client = CouchDB(pc.user,pc.password, url='https://picas-lofar.grid.surfsara.nl:6984',connect=True)
db=client[pc.database]  
print('Connected to database ' + pc.database)

print('2) Creating a list of paths to the files.')
from GRID_LRT.Staging import srmlist
s = srmlist.srmlist(check_OBSID=False)
with open(srms, 'r') as f:
    for l in f.readlines():
            s.append(l.strip())

print('3) Slicing list into groups.')
d = {}
sources = []
source = ''
for v in sorted(s):
    match = re.search('S\d{1,4}', v)
    if not match:
        raise ValueError('No sourcename extracted!')
    # First iteration, source is empty.
    if not source:
        source = match.group(0)
    # Check if we switch to a new source
    if match.group(0) != source:
        d[source] = sources
        sources = [v] 
    else:
        sources.append(v)
    source = match.group(0)
else:
    # No additional source to trigger the if statement, add the final source after the loop finishes.
    d[source] = sources
    del match, source, sources

print('4) Building token list.')
tl = TokenList(token_type=tok_type,database=db)
tokens = tl.list_view_tokens('step3_selfcal_cals')
token_ids = [token['_id'] for token in tokens]
#sys.exit()
for k,v in d.items():
    match = re.search('S\d{1,4}', v[0])
    if not match:
        raise ValueError('No sourcename extracted!')
    else:
        source = match.group(0)
    if (tok_type + '_scal_' + source + str(cal_obsid)) not in token_ids:
        tok = caToken(database = db, token_type = tok_type, token_id = tok_type + '_scal_' + source + str(cal_obsid))
        with open('temp_srm_{:s}.txt'.format(source), 'w') as f:
            f.write('\n'.join(v))
        tok.build(TokenJsonBuilder(jsonfile))
        tok.save()
        tok.add_attachment(attachment_name = 'srm.txt', filename = 'temp_srm_{:s}.txt'.format(source))
        tl.append(tok)
    else:
        continue
tl.add_attachment(attachment_name = 'step3_selfcal_calibrators.parset', filename = '/project/sksp/Software/lofar-highres-widefield/testdir/test_with_GRID_LRT/step3_selfcal_calibrators.parset')
tl.save()

for tok in tl:
    tok['OBSID'] = obsid
    tok['PIPELINE_STEP'] = 'lb_selfcal_cal1'
    tok['status'] = 'queued'
    tok.save()

print('5) Adding status views.')
view_cal = TokenView('step3_selfcal_cals', condition='doc.PIPELINE_STEP=="lb_selfcal_cal1"', emit_values=('doc._id', 'doc.status'))
tl.add_token_views()
tl.add_view(view_cal)

tl.save()

print('6) Create and launch the jobs.')
from GRID_LRT.application import submit
#j = submit.JdlLauncher(numjobs=len(d.keys()), token_type=tok_type, wholenodes=False, parameter_step=4, NCPU=2)
j = submit.SpiderLauncher(numjobs=len(d.keys()), token_type=tok_type, wholenode=False, parameter_step=64, NCPU=8)

with j:
    url=j.launch()
print('Job ID: ' + str(url))
