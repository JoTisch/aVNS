function radial_hist(beta, label, colorstr, pltflg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -Description-
%Plot radial histogram of beta
%
% -Input- 
%beta: stimulation locations during the heart beat (R-peak == 0°) 
%label: stimulation protocol label (1=non-syn, 2=sys-syn, 3=dia-sync) 
%colorstr: colors related to the labels 
%pltflg: 1: plot figure; 0: no 
%
% -Output- 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if pltflg == 1
    figure()
    for j = 3:-1:1 
        rad_angle = deg2rad(wrapTo360(beta(label==j))); 
        polarhistogram(rad_angle, 'BinWidth',pi/18,'FaceColor',colorstr(j),'FaceAlpha',.2)
        hold on
    end 

    % edit axes
    polarAxs = gca; 
    polarAxs.Color = 'none'; % no background
                
    polarAxs.ThetaAxis.FontSize = 12;
    polarAxs.LineWidth = 1;
    polarAxs.GridColor = 'k';
    polarAxs.GridAlpha = 0.5;
    
    % nifty way to dynamically create minor grid lines in 10 deg spacing, skipping
    % the major grid lines. No one will ever want to understand this.
    minorGrid = 10:10:350;
    minorGrid = minorGrid(logical(mod(minorGrid,polarAxs.ThetaTick(2))));
    polarAxs.ThetaAxis.MinorTickValues = minorGrid;
    polarAxs.ThetaMinorGrid = 'on';
    polarAxs.MinorGridColor = 'k';
    polarAxs.MinorGridAlpha = 0.5;
    % add legend 
    legend('systole-sync aVNS','diastole-sync aVNS','non-sync aVNS','FontSize',12)
end 

    
