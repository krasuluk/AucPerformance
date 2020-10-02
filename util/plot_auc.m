function plot_auc(P,A,CI,Value,Analysis)

n_metrics = length(A);
temp = 1;

for ii = 1:n_metrics-1
        for jj = ii+1:n_metrics
            pDI_vec(temp) = P(ii,jj);
            temp = temp+1;
        end
end

% Using Benjamini-Hochberg procedure for multiple comparisons 
% (note: correlation between groups has to be positive)

[P_sorted,rank_P] = sort(pDI_vec);
    for ii = length(P_sorted):-1:1
        Alpha = 0.05*(ii/length(P_sorted));
        if(P_sorted(ii) < Alpha)
            P_sorted(1:ii) = 1;
            break;
        else
            P_sorted(ii) = 0;
        end
    end
    pDI_vec(rank_P) = P_sorted;
    
    temp00 = 1;
    for ii = 1:n_metrics-1
        for jj = ii+1:n_metrics
            if(A(ii)>A(jj))
                Sig_DI(ii,jj) = pDI_vec(temp00);
                Sig_DI(jj,ii) = -pDI_vec(temp00);
            else
                Sig_DI(ii,jj) = -pDI_vec(temp00);
                Sig_DI(jj,ii) = pDI_vec(temp00);
            end
            temp00 = temp00 + 1;
        end
    end
    
    figure
    subplot(121)
    errorbar(A,CI,'.k');
    [ymin,posmin] = min(A);
       ymin = ymin-CI(posmin); 
       [ymax,posmax] = max(A);
       ymax = ymax+CI(posmax);
       ed = 0.1*(ymax-ymin);
       MAXy = ymax;
        MINy = ymin;
    if(MAXy == MINy)
        MAXy = ymax+0.1;
        MINy = ymin-0.1;
    end
    MAXx = n_metrics;
    MINx = 1;
    set(gca,'XLim', [MINx-0.1*(MAXx-MINx),MAXx+0.1*(MAXx-MINx)],'YLim', [MINy-0.1*(MAXy-MINy),MAXy+0.1*(MAXy-MINy)],'XTick',1:n_metrics)
    
       set(gcf,'Position',[100,100,350,170])
    ylabel(Value)
    xlabel('Metrics (-)')
    title(Analysis)
    subplot(122)
    imagesc(Sig_DI)
    colormap(gray)
    title('Significance')
    set(gca,'XTick',1:n_metrics,'YTick',1:n_metrics)
    