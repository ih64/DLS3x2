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
p.FLAGSB,
p.FLAGSV,
p.FLAGSR,
p.FLAGSz,
s.e1, 
s.e2, 
s.de, 
s.a, 
s.b,
s.theta,
s.status,
z.z_b, 
z.odds
FROM
RC1c_public.PhotoObj AS p INNER JOIN RC1c_public.Bpz AS z ON p.objid = z.objid
INNER JOIN RC1Stage.Shapes2 AS s ON p.objid = s.objid
WHERE
p.objid IS NOT NULL
# qa from dls courtesy perry
# AND NOT po.processflags & 0x20 !=0 
AND p.Rdered is NOT NULL
AND p.Bdered is NOT NULL
AND p.Vdered is NOT NULL
AND p.zdered is NOT NULL
AND p.R>18 AND p.R<24.5
# The R band probability that object is a point source `d.Dlsqc_prob`
AND p.dlsqc_prob<0.1
AND z.z_b>0.1
AND z.z_b<1.
AND p.excluded = 0
