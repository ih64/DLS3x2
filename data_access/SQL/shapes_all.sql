SELECT
s.id, 
s.x, 
s.y, 
s.e1, 
s.e2, 
s.de, 
s.flux_radius, 
s.a, 
s.b, 
s.theta, 
s.alpha, 
s.delta, 
s.subfield, 
p.r, 
z.z_b
FROM
RC1Stage.Shapes2 as s,
RC1c_public.Dlsqc AS d,
RC1Stage.PhotoObjAll AS p,
RC1c_public.Bpz AS z
WHERE
p.objid=s.objid
AND p.objid = s.objid
AND z.objid = d.objid
AND d.Dlsqc_prob<0.1
AND s.objid is NOT NULL
AND p.r is NOT NULL 
AND s.de < 99
