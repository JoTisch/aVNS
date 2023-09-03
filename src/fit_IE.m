function [] = fit_IE(alpha, RR, flagR, flagS, theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Description: 
%
%Input: 
%
%Output: 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                    Shift data by theta degrees                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if theta > 0 
    %
    [alpha,s_idx] = sort(alpha); 
    alpha(1:sum(alpha < theta),1) = alpha(1:sum(alpha < theta),1) + 360; 
    % 
    for k = 1:5 
    RR{k,:} = RR{k,:}(s_idx); 
    end 
    flagS = flagS(s_idx); 
    flagR = flagR(s_idx); 
    B = ["Non-sync aVNS", "Systole-sync aVNS", "Diastole-sync aVNS"];
    [~, index] = ismember(flagS, B);
    [~, s] = sort(index);
    flagS = flagS(s)';
    flagR = flagR(s)'; 
    alpha = alpha(s); 
    for k = 1:5 
    RR{k,:} = RR{k,:}(s); 
    end 
end 
x_lim = [theta,360+theta]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            Linear approximation                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:4 
    sIfit{k,:}= fit(alpha,((RR{k+1,:} ./ RR{1,:}) -1 )* 100,'poly1','Exclude',or(or(flagS ~= "Systole-sync aVNS",flagR ~= "Inspiration"),alpha<200)); 
    sEfit{k,:}= fit(alpha,((RR{k+1,:} ./ RR{1,:}) -1 )* 100,'poly1','Exclude',or(or(flagS ~= "Systole-sync aVNS",flagR ~= "Expiration"),alpha<200)); 
    dIfit{k,:}= fit(alpha,((RR{k+1,:} ./ RR{1,:}) -1 )* 100,'poly1','Exclude',or(or(flagS ~= "Diastole-sync aVNS",flagR ~= "Inspiration"),alpha > 250)); 
    dEfit{k,:}= fit(alpha,((RR{k+1,:} ./ RR{1,:}) -1 )* 100,'poly1','Exclude',or(or(flagS ~= "Diastole-sync aVNS",flagR ~= "Expiration"),alpha > 250)); 
    nfit{k,:}= fit(alpha,((RR{k+1,:} ./ RR{1,:}) -1 )* 100,'poly1','Exclude',flagS ~= "Non-sync aVNS");
end 

xfit1 = linspace(250,360+theta,100); 
xfit2 = linspace(theta,250,100); 
xfit3 = linspace(theta,360,150); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               Plot                                      %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colorstr = ["r","k","b"];  
figure()
hp1 = uipanel('position',[0 .5 .5 .5]);
hp2 = uipanel('position',[0 0 .5 .5]);
hp3 = uipanel('position',[.5 .5 .5 .5]);
hp4 = uipanel('position',[0.5 0 .5 .5]);
hp = {hp1; hp2; hp3; hp4}; 
g = findgroups(flagS); 
    for k = 1:4
        %%%
        s = scatterhist(alpha,((RR{k+1,:} ./ RR{1,:}) -1 )* 100,'NBins',[30 20], ...
            'Color','kbr', 'LineStyle',{'-','-.',':'},'Marker','+od','Style','bar', ...
            'Parent',hp{k},'Group',flagS,'Direction','in');
        for j=1:length(unique(g))
         cla(s(2))
         cla(s(3))
         hold(s(2),'on')
         hold(s(3),'on')
         arrayfun(@(j)histogram(s(2),alpha(g==j),'BinWidth',11.7,'FaceColor',...
             colorstr(j),'Normalization','probability'),unique(g))
        end 
         boxplot(s(3),((RR{k+1,:} ./ RR{1,:}) -1 )* 100,flagS,'orientation','horizontal',...
             'label', repmat({''},length(unique(flagS)),1),'color','kbr');
        axis(s(2:3),'tight')
        axis(s(2:3),'on')
        set(s(2),'xtick',[],'Xcolor','w','box','off')    
        hold on 
        %
        plot(xfit1,sIfit{k,:}.p1*xfit1+sIfit{k,:}.p2, 'Color','r', ... 
           'LineStyle',':','LineWidth',2,...
           'DisplayName', sprintf("Systole/Inspiration: m=%4.2f", ...
           sIfit{k,:}.p1)); 
        plot(xfit1,sEfit{k,:}.p1*xfit1+sEfit{k,:}.p2, 'Color','r', ... 
           'LineStyle','--','LineWidth',2,...
           'DisplayName', sprintf("Systole/Expiration: m=%4.2f", ...
           sEfit{k,:}.p1)); 
        %
        plot(xfit2,dIfit{k,:}.p1*xfit2+dIfit{k,:}.p2, 'Color','b', ... 
           'LineStyle',':','LineWidth',2,...
           'DisplayName', sprintf("Diastole/Inspiration: m=%4.2f", ...
           dIfit{k,:}.p1)); 
        plot(xfit2,dEfit{k,:}.p1*xfit2+dEfit{k,:}.p2, 'Color','b', ... 
           'LineStyle','--','LineWidth',2,...
           'DisplayName', sprintf("Diastole/Expiration: m=%4.2f", ...
           dEfit{k,:}.p1));
        %
        h3 = plot(nfit{k,:});
        set(h3,'Color','m','LineStyle','--','LineWidth',2, ...
            'DisplayName',sprintf("Non-sync: m=%4.2f",nfit{k,:}.p1)) 
        %
        title(['$\Delta RR = \frac{RR_{\mathrm{i+',num2str(k),'}} - RR_{\mathrm{i}}}{RR_{\mathrm{i}}} \cdot 100$'],'Interpreter','latex')
        xlim(x_lim)
        ylim([-30,30])
        ylim(s(3),[0,4])
        xlabel('$deg (^\circ)$','interpreter','latex'); 
        ylabel('$\Delta RR \; \mathrm{(rel. units)}$','interpreter','latex');
    end 