function [] = RdR_fit(alpha, RR, groupflag, pltflag, theta, showfit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%
% Input: 
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
    groupflag = groupflag(s_idx); 
    B = ["Non-sync aVNS", "Systole-sync aVNS", "Diastole-sync aVNS"];
    [~, index] = ismember(groupflag, B);
    [~, s] = sort(index);
    groupflag = groupflag(s)';
    alpha = alpha(s); 
    for k = 1:5 
    RR{k,:} = RR{k,:}(s); 
    end 
end 
x_lim = [theta,360+theta]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               FIT                                       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:4 
    if any(groupflag == "Expiration")
    ifit{k,:}= fit(alpha,((RR{k+1,:} ./ RR{1,:}) -1 )* 100,'poly1','Exclude',groupflag ~= "Inspiration"); 
    efit{k,:}= fit(alpha,((RR{k+1,:} ./ RR{1,:}) -1 )* 100,'poly1', 'Exclude', groupflag ~= "Expiration"); 
    elseif any(groupflag == "Non-sync aVNS")
    sfit{k,:}= fit(alpha,((RR{k+1,:} ./ RR{1,:}) -1 )* 100,'poly1','Exclude',or(groupflag ~= "Systole-sync aVNS",alpha < 200)); 
    dfit{k,:}= fit(alpha,((RR{k+1,:} ./ RR{1,:}) -1 )* 100,'poly1','Exclude',or(groupflag ~="Diastole-sync aVNS", alpha > 250)); 
    nfit{k,:}= fit(alpha,((RR{k+1,:} ./ RR{1,:}) -1 )* 100,'poly1','Exclude',groupflag ~= "Non-sync aVNS");
    end
end 
xfit1 = linspace(250,360+theta,100); 
xfit2 = linspace(theta,250,100); 
xfit3 = linspace(theta,360,150); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               Plot                                      %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if pltflag == 1 
    if any(groupflag == "Expiration")
        colorstr = ["k","b"]; 
    else 
        colorstr = ["r","k","b"]; 
    end 
%       figure()
%     hp1 = uipanel('position',[0 .5 .5 .5]);
%     hp2 = uipanel('position',[0 0 .5 .5]);
%     hp3 = uipanel('position',[.5 .5 .5 .5]);
%     hp4 = uipanel('position',[.5 0 .5 .5]);
%     hp = {hp1; hp2; hp3; hp4}; 
    g = findgroups(groupflag); 
    for k = 1:4
        figure() % delete
        s = scatterhist(alpha,((RR{k+1,:} ./ RR{1,:}) -1 )* 100,'NBins',[30 20], ...
            'Color','kbr', 'LineStyle',{'-','-.',':'},'Marker','+od','Style','bar', ...
            'Group',groupflag,'Direction','in'); % 'parent', hp{k},
        %
        for j=1:length(unique(g))
         cla(s(2))
         cla(s(3))
         hold(s(2),'on')
         hold(s(3),'on')
         arrayfun(@(j)histogram(s(2),alpha(g==j),'BinWidth',11.7,'FaceColor',...
             colorstr(j),'Normalization','pdf'),unique(g))
        end 
        %
        boxplot(s(3),((RR{k+1,:} ./ RR{1,:}) -1 )* 100,groupflag,'orientation','horizontal',...
             'label', repmat({''},length(unique(groupflag)),1),'color','kbr');

        %
        axis(s(2:3),'tight')
        axis(s(2:3),'on')
        set(s(2),'xtick',[],'Xcolor','w','box','off')    
        hold on 
        %
    
        if any(groupflag == "Expiration") && showfit == 1 
        h1 = plot(xfit3, ifit{k,:}.p1*xfit3+ifit{k,:}.p2, ...
            'DisplayName',sprintf("Inspiration: m=%4.2f, int.=%4.2f",ifit{k,:}.p1,-ifit{k,:}.p2/ifit{k,:}.p1)); 
        h2 = plot(xfit3, efit{k,:}.p1*xfit3+efit{k,:}.p2, ...
            'DisplayName',sprintf("Expiration: m=%4.2f, int.=%4.2f",efit{k,:}.p1,-efit{k,:}.p2/-efit{k,:}.p1)); 
        set(h1, 'LineWidth',2)
        set(h2, 'LineWidth',2)
        elseif any(groupflag == "Non-sync aVNS") && showfit == 1 
        plot(xfit1,sfit{k,:}.p1*xfit1+sfit{k,:}.p2, ...
            'Color','r','LineStyle','--','LineWidth',2,'DisplayName',sprintf("Systole: m=%4.2f, int.=%4.2f",sfit{k,:}.p1,-sfit{k,:}.p2/sfit{k,:}.p1)); 
        plot(xfit2,dfit{k,:}.p1*xfit2+dfit{k,:}.p2, ...
            'Color','b','LineStyle','--','LineWidth',2,'DisplayName',sprintf("Diastole: m=%4.2f, int.=%4.2f",dfit{k,:}.p1,-dfit{k,:}.p2/dfit{k,:}.p1)); 
        h3 = plot(nfit{k,:});
        set(h3,'Color','m','LineStyle','--','LineWidth',2,'DisplayName',sprintf("Non-sync: m=%4.2f",nfit{k,:}.p1))
        end 
        %
        legend('off')
        xlim(x_lim)
        ylim([-30,30])
        ylim(s(3),[0,4])
        xlabel('$\alpha (^\circ)$','interpreter','latex'); 
        ylabel(['$\Delta RR_{\mathrm{i+',num2str(k),'}} \, (\%)$'],'interpreter','latex');
    end 
end 
end


 



