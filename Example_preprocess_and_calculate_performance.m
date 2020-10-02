% Example of preprocessing for "Metrics_performance" function using
% Toyama database (Z. M. Parvez Sazzad, Y. Kawayoke, and Y. Horita.
% Image quality evaluation database. http://mict.eng.u-toyama.ac.jp/database_toyama/.)
%
% In case of any questions, contact Lukas Krasula (L.Krasula@gmail.com)

clear all; close all; clc;

load('util/Example_Toyama.mat')
% MOS       - column vector of MOS values
% SD        - column vector of respective standard deviations
% num_obs   - column vector of number of observers that evaluated each stimulus
% psnr      - column vector of PSNR values
% ssim      - column vector of SSIM values

N = size(MOS,1);
metrics = {'psnr','ssim'};

%% calculate significance of differences according to ITU-T Rec. J.149
% Note that for data from indirect scaling experiment, the procedure will
% be different. 
% The goal is to obtain a matrix "signif" as
%           signif =  0 : pair is not significantly different
%                    -1 : first stimulus is statistically significantly worse
%                     1 : first stimulus is statistically significantly better

d_mos = repmat(MOS,1,N) - repmat(MOS',N,1); % calculate MOS difference 
d_mos = d_mos(~tril(ones(N))); % take elements above diagonal (we are not interested in self differences) 
    
se_mos = sqrt( repmat(((SD.^2)./num_obs),1,N) + repmat(((SD.^2)./num_obs)',N,1));
se_mos = se_mos(~tril(ones(N)));
se_mos(se_mos == 0) = 1e-6; % to avoid division by zero;

z = d_mos./se_mos;
    
p_z = .5+erf(abs(z)/sqrt(2))/2;
signif = p_z;
signif(p_z<0.95) = 0;
signif(p_z>=0.95) = sign(d_mos(p_z>=0.95)); 

    
%% calculate the ObjScoDif (i.e. differences of objective scores)
objScoDif = [];
for m = 1:length(metrics)
               
    eval(['OM = ' metrics{m} ';']) % read the metric values
        
    d_om = repmat(OM,1,N) - repmat(OM',N,1);
    objScoDif(:,m) = d_om(~tril(ones(N)));
        
   SROCC(m) = corr(MOS,OM,'type','spearman');
   KROCC(m) = corr(MOS,OM,'type','kendall');
        
end

% transposition, since the function requires signif to be a row vector
% see help Metrics_performance
signif = signif';
objScoDif = objScoDif';

%% Performance evaluation
results = Metrics_performance(objScoDif,signif,1);

