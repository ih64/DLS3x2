SELECT
po.subfield,  
po.alpha, 
po.delta, 
po.distancetoborder,
po.distancetocelledge,
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
z.z_b, 
z.odds
FROM
RC1Stage.PhotoObj as po,
RC1c_public2.PhotoObj AS p, 
RC1c_public2.Dlsqc AS d, 
RC1c_public2.Bpz AS z, 
RC1Stage.Shapes2 AS s 
WHERE
d.objid=s.objid
AND p.objid = s.objid
AND p.objid = z.objid
and p.objid = po.objid
AND p.objid IS NOT NULL
AND po.processflags<8
AND s.status = 1
AND p.Rdered is NOT NULL
AND p.Bdered is NOT NULL
AND p.Vdered is NOT NULL
AND p.zdered is NOT NULL
AND p.Rdered>18 AND p.Rdered<24.5
# The R band probability that object is a point source `d.Dlsqc_prob`
AND d.Dlsqc_prob<0.1
# Shape cut
#AND s.b>0.4
AND z.z_b>0.1
AND z.z_b<1.
# Ellipticity error cut
#AND s.de<0.25
AND p.subfield LIKE 'F2%'
