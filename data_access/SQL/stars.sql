SELECT
p.subfield,  
p.alpha, 
p.delta, 
p.r, 
s.e1, 
s.e2, 
s.de, 
s.a, 
s.b, 
p.processflags, 
s.flux_radius,
s.status
FROM
RC1Stage.PhotoObjAll AS p, 
RC1c_public2.Dlsqc as d,
RC1Stage.Shapes2 AS s 
WHERE
d.objid=s.objid
AND p.objid = s.objid
AND p.objid = d.objid
AND p.objid IS NOT NULL
AND p.processflags<8
AND p.r is NOT NULL
AND d.Dlsqc_prob > .3

