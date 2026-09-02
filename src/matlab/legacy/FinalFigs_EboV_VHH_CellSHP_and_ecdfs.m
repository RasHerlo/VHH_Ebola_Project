%% EboV_VHH_FinalFigures

% Written by Rasmus Herlo, 220318

% In absence of former figures, this the figures are regenerated here
%{
- Overview of wholecells with AP2, VSV and VHHs (G10 and G84)
    - ZOOM-INs
        - VSV-internalized (+VHH)
        - VSV-int. (wo/VHH)
        - VSV on surface(w/AP2-RFP)
        - VSV on surface (wo/AP2-RFP)
- ecdf-graphs of:
    - VHH effect on Cell-Association (previously "precipitation")
    - VHH effect on AP2-association
    - VHH effect on internalization
%}


%% CellShape WholeCell Overviews

%{

Some prerun cellshapes can be found in:
G84:
G:\LLSM_VSV_VHH\TEMPORARYDATA\CS_G84
G10:
- G:\LLSM_VSV_VHH\TEMPORARYDATA\CS3_G10_Hidde1

If new needed, use:
'C:\Users\Ras\Dropbox\WORK\WORK POST-DOC NEW YORK\Code Repository\Code Repositories\CellClusterSurfShapes2_RH.m'

%}


%% G84 VERSION: Loading Elements, File-Overviews and Data-Collections

% DataOverview:
rtdir = 'G:\LLSM_VSV_VHH\TEMPORARYDATA\CS_G84';
load(fullfile(rtdir,'data488master.mat'))
load(fullfile(rtdir,'M488.mat'))
NumsX = unique(M488.ExpNum);


for dt = 1:numel(data)
    rtdirX = data(dt).source;
    rtdirX = strsplit(rtdirX,filesep);
    rtdirX = strjoin(rtdirX(1:end-2),filesep);
    % CellSHP
    CellSHP = load(fullfile(rtdirX,'CellSHP.mat'));
    open(fullfile(rtdirX,'CellSHP_fig.fig'))
    
    % Change Cluster Properties
    Chl = get(gca, 'Children'); % Chl(1) = scatter data
    Sz = Chl(1).SizeData;
    Chl(1).SizeData = Sz * 1.5; % Sizes of Red Clusters
    % Make AlphaShape gray and Transparent
    Chl(2).FaceColor = [0.75 0.75 0.75];
    set(gca,'ZDir','reverse')
    
    Idx = find(M488.ExpNum == NumsX(dt));
    % plot the scatter, play with sizes and cutoff
    n = 2;
    while n == 2
        Pms = inputdlg({'Upper percentile:','Lower percentile:','Marker Size:','Marker Edge Size:'},'Input',[1 35],{'0.85','0.20','0.2','0.002'});
        hold on
        [~,Idx2] = mink(M488.A_tot(1,Idx),round(numel(M488.A_tot(1,Idx))*str2double(Pms{1})));
        [~,Idx3] = mink(M488.A_tot(1,Idx),round(numel(M488.A_tot(1,Idx))*str2double(Pms{2})));
        Idx4 = setdiff(Idx(Idx2),Idx(Idx3));
        % sc3 = scatter3(M488.x_tot(1,Idx3),M488.y_tot(1,Idx3),M488.z_tot(1,Idx3),M488.A_tot(1,Idx3)*str2double(Pms{2}),...
        %    'g','filled','MarkerEdgeColor','b');
        for h = 1:numel(M488.A_tot(3,Idx4))
            if M488.A_tot(1,Idx4(h))>40
                scatter3(M488.x_tot(1,Idx4(h)),M488.y_tot(1,Idx4(h)),M488.z_tot(1,Idx4(h)),M488.A_tot(1,Idx4(h))*str2double(Pms{3}),...
                    'g','filled','MarkerEdgeColor','b','LineWidth',M488.A_tot(1,Idx4(h))*str2double(Pms{4}));
            else
                scatter3(M488.x_tot(1,Idx4(h)),M488.y_tot(1,Idx4(h)),M488.z_tot(1,Idx4(h)),M488.A_tot(1,Idx4(h))*str2double(Pms{3}),...
                    'g','filled','MarkerEdgeColor','k','LineWidth',M488.A_tot(1,Idx4(h))*str2double(Pms{4}));
            end
            if mod(h,150)==0
                disp(['Plotted ' num2str(h) ' out of  ' num2str(numel(M488.A_tot(3,Idx4))) 'scatters...'])
            end
        end
        b = gca; legend(b,'off');
        title(['Ex#' num2str(NumsX(dt)) ', Perc: [' Pms{1} ' ' Pms{2} '], MSz = ' Pms{3} ', ESz = ' Pms{4}])
        n = menu('Want to keep the look of it?','Yes','No'); 
    end
    % Store Parameters
    Xstrc(dt).Perc = str2double(Pms{1});
    Xstrc(dt).MSz = str2double(Pms{2});
    Xstrc(dt).ESz = str2double(Pms{3});
    [Xstrc(dt).az,Xstrc(dt).el] = view;
end

%% Specific for EX7

xlim([50 450])
zlim([0 100])
view([-30.5513,50.7614])
title('')
set(gcf,'Color','w')
xticklabels([])
yticklabels([])
zticklabels([])
grid on

%% G10 version


% DataOverview:
rtdir = 'G:\LLSM_VSV_VHH\TEMPORARYDATA\CS3_G10_Hidde1';
load(fullfile(rtdir,'data488master.mat'))
load(fullfile(rtdir,'M488.mat'))
NumsX = unique(M488.ExpNum);


for dt = 1:numel(data)
    rtdirX = data(dt).source;
    rtdirX = strsplit(rtdirX,filesep);
    rtdirX = strjoin(rtdirX(1:end-2),filesep);
    rtdirX = ['G' rtdirX(2:end)]; % different drive
    % CellSHP
    CellSHP = load(fullfile(rtdirX,'CellSHP.mat'));
    
%     open(fullfile(rtdirX,'CellSHPFig.fig'))
%     % Change Cluster Properties
%     Chl = get(gca, 'Children'); % Chl(1) = scatter data
%     Sz = Chl(1).SizeData;
%     Chl(1).SizeData = Sz * 1.5; % Sizes of Red Clusters
%     % Make AlphaShape gray and Transparent
%     Chl(2).FaceColor = [0.75 0.75 0.75];
%     set(gca,'ZDir','reverse')
    
    Idx = find(M488.ExpNum == NumsX(dt));
    % plot the scatter, play with sizes and cutoff
    n = 2;
    while n == 2
        open(fullfile(rtdirX,'CellSHPFig.fig'))
        % Change Cluster Properties
        Chl = get(gca, 'Children'); % Chl(1) = scatter data
        Sz = Chl(1).SizeData;
        Chl(1).SizeData = Sz * 1.5; % Sizes of Red Clusters
        % Make AlphaShape gray and Transparent
        Chl(2).FaceColor = [0.75 0.75 0.75];
        set(gca,'ZDir','reverse')
        
        Pms = inputdlg({'Upper percentile:','Lower percentile:','Marker Size:','Marker Edge Size:'},'Input',[1 35],{'0.85','0.20','0.2','0.002'});
        hold on
        [~,Idx2] = mink(M488.A_tot(1,Idx),round(numel(M488.A_tot(1,Idx))*str2double(Pms{1})));
        [~,Idx3] = mink(M488.A_tot(1,Idx),round(numel(M488.A_tot(1,Idx))*str2double(Pms{2})));
        Idx4 = setdiff(Idx(Idx2),Idx(Idx3));
        % sc3 = scatter3(M488.x_tot(1,Idx3),M488.y_tot(1,Idx3),M488.z_tot(1,Idx3),M488.A_tot(1,Idx3)*str2double(Pms{2}),...
        %    'g','filled','MarkerEdgeColor','b');
        for h = 1:numel(M488.A_tot(3,Idx4))
            if M488.A_tot(1,Idx4(h))>40
                scatter3(M488.x_tot(1,Idx4(h)),M488.y_tot(1,Idx4(h)),M488.z_tot(1,Idx4(h)),M488.A_tot(1,Idx4(h))*str2double(Pms{3}),...
                    'g','filled','MarkerEdgeColor','b','LineWidth',M488.A_tot(1,Idx4(h))*str2double(Pms{4}));
            else
                scatter3(M488.x_tot(1,Idx4(h)),M488.y_tot(1,Idx4(h)),M488.z_tot(1,Idx4(h)),M488.A_tot(1,Idx4(h))*str2double(Pms{3}),...
                    'g','filled','MarkerEdgeColor','k','LineWidth',M488.A_tot(1,Idx4(h))*str2double(Pms{4}));
            end
            if mod(h,150)==0
                disp(['Plotted ' num2str(h) ' out of  ' num2str(numel(M488.A_tot(3,Idx4))) 'scatters...'])
            end
        end
        b = gca; legend(b,'off');
        title(['Ex#' num2str(NumsX(dt)) ', Perc: [' Pms{1} ' ' Pms{2} '], MSz = ' Pms{3} ', ESz = ' Pms{4}])
        n = menu('Want to keep the look of it?','Yes','No'); 
    end
    % Store Parameters
    Xstrc(dt).Perc = str2double(Pms{1});
    Xstrc(dt).MSz = str2double(Pms{2});
    Xstrc(dt).ESz = str2double(Pms{3});
    [Xstrc(dt).az,Xstrc(dt).el] = view;
end

%% Specific for G10-Ex11

xlim([0 500])
zlim([0 80])
view([ -174.0470,41.8301])
title('')
set(gcf,'Color','w')
xticklabels([])
yticklabels([])
zticklabels([])
grid on

% zoom:
xlim([200 300])
ylim([350 600])

% Zoom 2:
xlim([30 300])
ylim([50 400])

%% ECDFS STARTING HERE

%% Precipitation G10 and G84

rtdir = 'G:\LLSM_VSV_VHH\FIGURE FILES\CFDs';

% Precipitation G10
open(fullfile(rtdir,'G10_Precipitation.fig'))
Lines = get(gca,'Children');
xCB = Lines(1).XData;
yCB = Lines(1).YData; % CB = cellbound
xPr = Lines(2).XData; 
yPr = Lines(2).YData; % Pr = Precipitated

% make 20 brackets and take ratios within those
for t = 1:40
    CB(t) = mean(yCB(xCB>(t-1)*25 & xCB<=t*25));
    Pr(t) = mean(yPr(xPr>(t-1)*25 & xPr<=t*25));
end

Ptot = CB./(CB+10*Pr);
Pnorm = Ptot/Ptot(1);
save('Pnorm.mat','Pnorm')

yyaxis right
plot(1:1000/numel(Pnorm):1000,Pnorm,'r--')
ylim([0.75 1.05])
xlim([0 600])
xlabel('Relative VHH-binding (a.u.)')
leg = legend({'Not bound VSV','Bound VSV','Rel. Binding Prob.'});
set(leg,'Box','Off')
ylabel('Relative Binding Probability')
haxes = gca;
haxes(2).YColor = [0 0 0];
yyaxis left
ylabel('Cumulative Frequency')
title('VHH G10 effect on cell-binding')

% Precipitation G84
open(fullfile(rtdir,'G84_Precipitation.fig'))
Lines = get(gca,'Children');
xCB = Lines(1).XData;
yCB = Lines(1).YData; % CB = cellbound
xPr = Lines(2).XData; 
yPr = Lines(2).YData; % Pr = Precipitated

% make 20 brackets and take ratios within those
for t = 1:40
    CB(t) = mean(yCB(xCB>(t-1)*25 & xCB<=t*25));
    Pr(t) = mean(yPr(xPr>(t-1)*25 & xPr<=t*25));
end

Ptot = CB./(CB+10*Pr);
Pnorm = Ptot/Ptot(1);
save('Pnorm.mat','Pnorm')

yyaxis right
plot(1:1000/numel(Pnorm):1000,Pnorm,'r--')
ylim([0.75 1.05])
xlim([0 600])
xlabel('Relative VHH-binding (a.u.)')
leg = legend({'Not bound VSV','Bound VSV','Rel. Binding Prob.'});
set(leg,'Box','Off')
ylabel('Relative Binding Probability')
haxes = gca;
haxes.YColor = [0 0 0];
yyaxis left
ylabel('Cumulative Frequency')
title('VHH G84 effect on cell-binding')

%% AP2-association, G10 and G84

rtdir = 'G:\LLSM_VSV_VHH\FIGURE FILES\CFDs';

% Precipitation G10
open(fullfile(rtdir,'G10_AP2association.fig'))
Lines = get(gca,'Children');
xCB = Lines(1).XData;
yCB = Lines(1).YData; % CB = associated
xPr = Lines(2).XData; 
yPr = Lines(2).YData; % Pr = not associated

for t = 1:40
    CB(t) = mean(yCB(xCB>(t-1)*25 & xCB<=t*25));
    Pr(t) = mean(yPr(xPr>(t-1)*25 & xPr<=t*25));
end

CBtot = CB*numel(xCB);
Prtot = Pr*numel(xPr);

AP2ass = CBtot./(CBtot+Prtot);
AP2assP = AP2ass/AP2ass(1);

yyaxis right
plot(1:1000/numel(AP2assP):1000,AP2assP,'k--')
ylim([0.8 1.0])
yticks(0.8:0.05:1.0)
xlim([0 600])
xlabel('Relative VHH-binding (a.u.)')
leg = legend({'Not AP2-assoc. VSV','AP2-assoc. VSV','Rel. Assoc. Prob.'});
set(leg,'Box','Off')
ylabel('Relative AP2-association Probability')
haxes = gca;
haxes.YColor = [0 0 0];
yyaxis left
ylabel('Cumulative Frequency')
title('VHH G10 effect on AP2-association')

% Precipitation G84
open(fullfile(rtdir,'G84_AP2association.fig'))
Lines = get(gca,'Children');
xCB = Lines(1).XData;
yCB = Lines(1).YData; % CB = associated
xPr = Lines(2).XData; 
yPr = Lines(2).YData; % Pr = not associated

for t = 1:40
    CB(t) = mean(yCB(xCB>(t-1)*25 & xCB<=t*25));
    Pr(t) = mean(yPr(xPr>(t-1)*25 & xPr<=t*25));
end

CBtot = CB*numel(xCB);
Prtot = Pr*numel(xPr);

AP2ass = CBtot./(CBtot+Prtot);
AP2assP = AP2ass/AP2ass(1);
save('AP2assP.mat','AP2assP')

yyaxis right
plot(1:1000/numel(AP2assP):1000,AP2assP,'k--')
ylim([0.8 1.0])
yticks(0.8:0.05:1.0)
xlim([0 600])
xlabel('Relative VHH-binding (a.u.)')
leg = legend({'Not AP2-assoc. VSV','AP2-assoc. VSV','Rel. Assoc. Prob.'});
set(leg,'Box','Off')
ylabel('Relative AP2-association Probability')
haxes = gca;
haxes.YColor = [0 0 0];
yyaxis left
ylabel('Cumulative Frequency')
title('VHH G84 effect on AP2-association')

%% Internalization, G10 and G84

rtdir = 'G:\LLSM_VSV_VHH\FIGURE FILES\CFDs';

% Internalization G10
open(fullfile(rtdir,'G10_Internalization.fig'))
Lines = get(gca,'Children');
xCB = Lines(1).XData;
yCB = Lines(1).YData; % CB = surface
xPr = Lines(2).XData; 
yPr = Lines(2).YData; % Pr = internalized

for t = 1:40
    CB(t) = mean(yCB(xCB>(t-1)*25 & xCB<=t*25));
    Pr(t) = mean(yPr(xPr>(t-1)*25 & xPr<=t*25));
end

CBtot = CB*numel(xCB);
Prtot = Pr*numel(xPr);

Int = Prtot./(CBtot+Prtot);
IntP = (Int-0.9*min(Int))/(Int(1)-0.9*min(Int));
save('IntP.mat','IntP')

yyaxis right
plot(1:1000/numel(IntP):1000,IntP,'k--')
ylim([0.6 1.10])
yticks(0.6:0.05:1.1)
xlim([0 600])
xlabel('Relative VHH-binding (a.u.)')
leg = legend({'Surface VSV','Internalized VSV','Rel. Intern. Prob.'});
set(leg,'Box','Off')
ylabel('Relative Internalization Probability')
haxes = gca;
haxes.YColor = [0 0 0];
yyaxis left
ylabel('Cumulative Frequency')
title('VHH G10 effect on Internalization')

% G84 and internalization

rtdir = 'G:\LLSM_VSV_VHH\FIGURE FILES\CFDs';

% Internalization G84
open(fullfile(rtdir,'G84_Internalization.fig'))
Lines = get(gca,'Children');
xCB = Lines(1).XData;
yCB = Lines(1).YData; % CB = surface
xPr = Lines(2).XData; 
yPr = Lines(2).YData; % Pr = internalized

for t = 1:40
    CB(t) = mean(yCB(xCB>(t-1)*25 & xCB<=t*25));
    Pr(t) = mean(yPr(xPr>(t-1)*25 & xPr<=t*25));
end

CBtot = CB*numel(xCB);
Prtot = Pr*numel(xPr);

Int = Prtot./(CBtot+Prtot);
IntP = (Int-0.9*min(Int))/(Int(1)-0.9*min(Int));
save('IntP.mat','IntP')

yyaxis right
plot(1:1000/numel(IntP):1000,IntP,'k--')
ylim([0.6 1.10])
yticks(0.6:0.05:1.1)
xlim([0 600])
xlabel('Relative VHH-binding (a.u.)')
leg = legend({'Surface VSV','Internalized VSV','Rel. Intern. Prob.'});
set(leg,'Box','Off')
ylabel('Relative Internalization Probability')
haxes = gca;
haxes.YColor = [0 0 0];
yyaxis left
ylabel('Cumulative Frequency')
title('VHH G84 effect on Internalization')


