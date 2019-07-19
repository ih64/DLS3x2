SELECT
p.subfield,  
p.alpha, 
p.delta, 
p.Rdered,
p.Vdered,
p.Rdered, 
p.zdered,
s.e1, 
s.e2, 
s.de, 
s.a, 
s.b, 
s.status
FROM
RC1c_public.PhotoObj AS p, 
RC1c_public.Dlsqc as d,
RC1Stage.Shapes2 AS s 
WHERE
d.objid=s.objid
AND p.objid = s.objid
AND p.objid = d.objid
AND p.objid IS NOT NULL
AND p.FLAGSR < 4
#AND p.Rdered < 20
AND p.Rdered is NOT NULL
AND d.Dlsqc_prob > .3
# Shape cut
#AND s.b > 0.4
#AND s.de < .3
#AND s.status = 1
