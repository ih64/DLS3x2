SELECT po.objid,po.Rdered, po.Bdered - po.Vdered, po.Vdered - po.Rdered, po.FLAGSB, po.FLAGSV, po.FLAGSR, po.FLAGSz, z.z_b, s.de, s.b, s.status, s.e1, s.e2
FROM
RC1c_public.PhotoObj AS po
LEFT JOIN RC1c_public.Bpz as z
ON po.objid = z.objid
LEFT JOIN RC1Stage.Shapes2 as s
ON po.objid = s.objid
WHERE
po.excluded = 0
AND po.dlsqc_prob < 0.1
AND po.R BETWEEN 18 AND 27
