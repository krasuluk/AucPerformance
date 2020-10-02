function [results] = Metrics_performance(objScoDif, signif, doPlot)

% INPUT:    objScoDif   : differences of objective scores [M,N]
%                         M : number of metrics
%                         N : number of pairs
%           signif      : statistical outcome of paired comparison [1,N]
%                          0 : no difference
%                         -1 : first stimulus is worse
%                          1 : first stimulus is better
%           doPlot      : boolean indicating if graphs should be plotted
%
% OUTPUT:   results - structure with following fields
%
%           AUC_DS      : Area Under the Curve for Different/Similar ROC 
%                         analysis
%           pDS_DL      : p-values for AUC_DS from DeLong test
%           pDS_HM      : p-values for AUC_DS from Hanley and McNeil test
%           AUC_BW      : Area Under the Curve for Better/Worse ROC
%                         analysis
%           pBW_DL      : p-values for AUC_BW from DeLong test
%           pBW_HM      : p-values for AUC_BW from Hanley and McNeil test
%           CC_0        : Correct Classification @ DeltaOM = 0 for 
%                         Better/Worse ROC analysis
%           pCC0_b      : p-values for CC_0 from binomial test
%           pCC0_F      : p-values for CC_0 from Fisher's exact test
%           THR         : threshold for 95% probability that the stimuli
%                         are different

%
% REFERENCES
%
% The method is described in the paper:
% L. Krasula, K. Fliegel, P. Le Callet, M.Klima, "On the accuracy of 
% objective image and video quality models: New methodology for 
% performance evaluation", QoMEX 2016.
%
% When you use our method in your research, please, cite the above stated
% paper.
%
% Copyright (c) 2016
% Lukas Krasula <l.krasula@gmail.com>

% Permission to use, copy, modify, and/or distribute this software for any
% purpose with or without fee is hereby granted, provided that the above
% copyright notice and this permission notice appear in all copies.
%
% THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
% WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR
% ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
% WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
% ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
% OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
%
% This sotware also uses the code described in:
% X. Sun and W. Xu, "Fast Implementation of DeLong's Algorithm for
% Comparing the Areas Under Correlated Receiver Operating Characteristic
% Curves," IEEE Signal Processing Letters, vol. 21, no. 11, pp. 1389-1393, 
% 2014.
%
% and
%
% P. Hanhart, L. Krasula, P. Le Callet, T. Ebrahimi, "How to benchmark 
% objective quality metrics from pair comparison data", QoMEX 2016.

addpath('util');

M = size(objScoDif,1);

D = abs(objScoDif(:,signif ~= 0));
S = abs(objScoDif(:,signif == 0));
samples.spsizes = [size(D,2),size(S,2)];
samples.ratings = [D,S];

% calculate AUCs

[AUC_DS,C] = fastDeLong(samples);

% significance calculation

pDS_DL = ones(M);

for i=1:M-1
    for j=i+1:M
        pDS_DL(i,j) = calpvalue(AUC_DS([i,j]), C([i,j],[i,j]));
        pDS_DL(j,i) = pDS_DL(i,j);
    end
end


[pDS_HM,CI_DS] = significanceHM(S, D, AUC_DS);

THR = prctile(S',95);

%%%%%%%%%%%%%%%%%%%%%%% Better / Worse %%%%%%%%%%%%%%%%%%%%%%%%%%%

B = [objScoDif(:,signif == 1),-objScoDif(:,signif == -1)];
W = -B;
samples.ratings = [B,W];
samples.spsizes = [size(B,2),size(W,2)];

% calculate AUCs

[AUC_BW,C] = fastDeLong(samples);

% calculate correct classification for DeltaOM = 0

L = size(B,2) + size(W,2);
CC_0 = zeros(M,1);
for m=1:M
    CC_0(m) = (sum(B(m,:)>0) + sum(W(m,:)<0)) / L;
end

% significance calculation

pBW_DL = ones(M);
pCC0_b = ones(M);
pCC0_F = ones(M);

for i=1:M-1
    for j=i+1:M
        pBW_DL(i,j) = calpvalue(AUC_BW([i,j]), C([i,j],[i,j]));
        pBW_DL(j,i) = pBW_DL(i,j);
        
        pCC0_b(i,j) = significanceBinomial(CC_0(i), CC_0(j), L);
        pCC0_b(j,i) = pCC0_b(i,j);
        
        pCC0_F(i,j) = fexact(CC_0(i)*L, 2*L, CC_0(i)*L + CC_0(j)*L, L, 'tail', 'b')/2;
        pCC0_F(j,i) = pCC0_F(i,j);
    end
end

[pBW_HM,CI_BW] = significanceHM(B, W, AUC_BW);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Adding outputs to the structure

results.AUC_DS = AUC_DS;
results.pDS_DL = pDS_DL;
results.pDS_HM = pDS_HM;
results.AUC_BW = AUC_BW;
results.pBW_DL = pBW_DL;
results.pBW_HM = pBW_HM;
results.CC_0 = CC_0;
results.pCC0_b = pCC0_b;
results.pCC0_F = pCC0_F;
results.THR = THR;

%%%%%%%%%%%%%%%%%%%%%%%% Plot Results %%%%%%%%%%%%%%%%%%%%%%%%%%%

if(doPlot == 1)

% Using Benjamini-Hochberg procedure for multiple comparisons in plots
% (note: correlation between groups has to be positive)

plot_auc(results.pDS_HM,results.AUC_DS, CI_DS, 'AUC (-)','Different/Similar')
plot_cc(results.pCC0_F,results.CC_0,'C_0 (%)','Better/Worse')
plot_auc(results.pBW_HM,results.AUC_BW, CI_BW, 'AUC (-)','Better/Worse')

end
