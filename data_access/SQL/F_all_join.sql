SELECT
p.objid,
p.subfield,  
p.alpha, 
p.delta, 
p.R,
p.Bdered,
p.Vdered,
p.Rdered,
p.zdered,
s.e1, 
s.e2, 
s.de, 
s.a, 
s.b,
s.theta,
s.status,
z.*
FROM
RC1c_public2.PhotoObj AS p INNER JOIN RC1c_public2.Bpz AS z ON p.objid = z.objid
INNER JOIN RC1Stage.Shapes2 AS s ON p.objid = s.objid
WHERE
p.objid IS NOT NULL
# qa from dls courtesy perry
# AND NOT po.processflags & 0x20 !=0 
AND p.Rdered is NOT NULL
AND p.Bdered is NOT NULL
AND p.Vdered is NOT NULL
AND p.zdered is NOT NULL
AND p.FLAGSB < 4
AND p.FLAGSV < 4
AND p.FLAGSR < 4
AND p.FLAGSz < 4
AND p.R > 18 AND p.R < 24.5
# The R band probability that object is a point source `d.Dlsqc_prob`
AND p.dlsqc_prob<0.1
#AND z.z_b>0.3
#AND z.z_b<1.
AND p.excluded = 0
