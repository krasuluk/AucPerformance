function [pHM,CI] = significanceHM(A,B,AUCs)
% By Lukas Krasula

n_met = size(A,1);

CorrA = corr(A','type','Kendall');
CorrB = corr(B','type','Kendall');

pHM = ones(n_met);
CI = ones(n_met,1);

for i=1:n_met-1
    
    [CI(i),SE1] = AUC_CI(size(A,2),size(B,2),AUCs(i));
    
    for j=i+1:n_met
        [CI(j),SE2] = AUC_CI(size(A,2),size(B,2),AUCs(j));
        
        load('Hanley_McNeil.mat');
        
        rA = (CorrA(i,j) + CorrB(i,j))/2;
        AA = (AUCs(i) + AUCs(j))/2;
        
        [~,rr] = min(abs(rA-rA_vec));
        [~,aa] = min(abs(AA-AA_vec));
        r = Table_HM(rr,aa);
        
        z = abs(AUCs(i) - AUCs(j)) / sqrt( SE1^2 + SE2^2 + 2*r*SE1*SE2 );
        pHM(i,j) = 1-normcdf(z);
        pHM(j,i) = pHM(i,j);
    end
end
