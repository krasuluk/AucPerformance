function pvalue = calpvalue(aucs, sigma)
l = [1, -1];
z = abs(diff(aucs)) / sqrt(l * sigma * l');
pvalue = 2 * (1 - normcdf(z, 0, 1));
