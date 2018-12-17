SELECT b.*
FROM
RC1c_public2.PhotoObj AS p INNER JOIN RC1c_public2.Bpz AS z on p.objid = z.objid
INNER JOIN RC1Stage.Shapes2 AS s on p.objid = s.objid
INNER JOIN RC1c_public2.Probs as b on p.objid = b.objid
WHERE
p.objid IS NOT NULL
AND p.Rdered is NOT NULL
AND p.Bdered is NOT NULL
AND p.Vdered is NOT NULL
AND p.zdered is NOT NULL
AND p.R < 21
AND p.R > 18
AND p.FLAGSB < 4
AND p.FLAGSV < 4
AND p.FLAGSR < 4
AND p.FLAGSz < 4
AND p.excluded = 0
# The R band probability that object is a point source `d.Dlsqc_prob`
AND p.dlsqc_prob<0.1
