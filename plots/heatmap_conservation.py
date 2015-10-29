#!/usr/bin/env python

import manipulate_biosqldb
import pandas as pd
import numpy
# DataFrame.tojson(orient="split")


server, db = manipulate_biosqldb.load_db('chlamydia_03_15')

sql = 'SELECT * FROM orth_chlamydia_03_15.group_1111;'

data = numpy.array([list(i) for i in server.adaptor.execute_and_fetchall(sql,)])

locus_list = '"' + '","'.join(data[0:,1]) + '"'

sql2 = 'select locus_tag, organism from orthology_detail_chlamydia_03_15 where locus_tag in (%s)' % locus_list
print sql2
locus2organism = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))

print locus2organism

#print data

print data[0:,2:]

columns = data[0:,1]
rows = [i + " (%s)" % locus2organism[i] for i in columns]
frame = pd.DataFrame(data[0:,2:], index=rows, columns=columns)
frame = frame.astype(float)
frame = frame/100

for i in range(0, len(frame)):

    frame.ix[i, i] = None

print frame

with open('group_1523.json', 'w') as f:
    f.write(frame.to_json(orient="split"))
