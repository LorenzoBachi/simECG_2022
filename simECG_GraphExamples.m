% Code used to generate the ECG examples in the paper "ECG Modeling for
% Simulation of Arrhythmias in Time-Varying Conditions" (2023).
% 
% Copyright (c), Cristina Pérez, University of Zaragoza
% last revision: 08 2022

clc
%1) Load examples (each row is a single-lead ECG signal of 15-s duration)
fs = 1000;
examples = zeros(15*fs,4);
examples = xxx %    COMPLEATE
t = (0:15*fs-1)./fs;

%2) Representation
figure
figu = tiledlayout(5,11);
figu.Padding = 'compact';
figu.TileSpacing = 'Tight';

%--> Single-lead ECG signals
for ii =1:4
    ax(ii) = nexttile((ii*11+1)-11,[1 10]);
    % gridlines ---------------------------
    g_y = [-2:0.1:2]; % user defined grid Y [start:spaces:end] %in mV
    g_x = [0:0.04:15]; % user defined grid X [start:spaces:end] %in seconds
    
    for i=1:length(g_x)
        if mod(g_x(i),0.2) == 0
            plot([g_x(i) g_x(i)],[g_y(1) g_y(end)],'-','Color',[1, 0, 0, 0.5],'LineWidth',0.5) %y grid lines
        else
            plot([g_x(i) g_x(i)],[g_y(1) g_y(end)],'-','Color',[1, 0, 0, 0.5],'LineWidth',0.2) %y grid lines
        end
        hold on
    end
    for i=1:length(g_y)
        if mod(g_y(i),0.5) == 0
            plot([g_x(1) g_x(end)],[g_y(i) g_y(i)],'-','Color',[1, 0, 0, 0.5],'LineWidth',0.5) %x grid lines
        else
            plot([g_x(1) g_x(end)],[g_y(i) g_y(i)],'-','Color',[1, 0, 0, 0.5],'LineWidth',0.2) %x grid lines
        end
        hold on
    end
    % ---------
    plot(t, examples(:,ii),'k')
    set(gca,'fontsize',9,'FontName','Times','box','off','xtick',(0:1:15),'ytick',(-2:1:2))
end
xlabel('Time (s)','FontName', 'Times','fontsize',9)
linkaxes([ax],'xy')

%--> Boxes
for ii = 1:4
    ab(ii) = nexttile(ii*11,[1 1]);
    rectangle('Position',[0 0 0.5 0.5])
    text(0.1,0.7,strcat('#',num2str(ii)),'FontName', 'Times') %Example numbers
    ylim([-0.5 1]), xlim([0 1])
    set(gca,'fontsize',9,'FontName', 'Times','visible','off')
end
nexttile(11,[1 1]);
text(0,1.2,{'Check box'; 'if unrealistic'},'FontName', 'Times')

ab(6) = nexttile(45,[1 11]);
rectangle('Position',[0 0 2 2])
x1 = repmat(0.01,1,3)';
y1 = [0.5 1 1.5]';
x2 = repmat(1.99,1,3)';
y2 = y1;
for ii = 1:3
    line([x1(ii) x2(ii)],[y1(ii) y2(ii)],'Color','k','LineStyle','-')
    text(0.01,y1(3-ii+1)+0.3,strcat('#',num2str(ii),':'),'FontName', 'Times')%Example numbers
end
text(0.01,0.3,strcat('#',num2str(4),':'),'FontName', 'Times')%Example numbers

xlim([0 2]), ylim([0 2])
set(gca,'fontsize',9,'FontName', 'Times','visible','off')


%--> Page number
sgtitle('Page 1','FontName', 'Times')


%3) Save pdf automatically
orient(gcf,'landscape')
print(gcf,'-dpdf', '-fillpage', 'simu1.pdf','-r1500')





