function pValue = significanceBinomial(p1,p2,N)
p = (p1+p2) / 2;
sigmaP1P2 = sqrt(p*(1-p)*2/N);
z = abs(p1-p2)/sigmaP1P2;
pValue = 2*(1 - normcdf(z, 0, 1));
