%% FIGURES for Cells

% Files used for making figures chronologically:

% Figures found:
% G:\LLSM_VSV_VHH\2021 Processing In Vivo\New DataProc VHH Ebola Cellular Data.ppt
% G:\LLSM_VSV_VHH\FIGURE FILES

%% Calibration of undisturbed AP2s: (images and histograms) – NoMOI cells: AP2-pit = 30a.u. 

rtdir = 'G:\LLSM_VSV_VHH\SUM52\CS1_SUM52_NoMOI';
cd(rtdir)

load('dataM560.mat')
load('M560.mat')

% Check Histogram
figure,
histogram(M560.A_tot,10:10:1000)
set(gcf, 'Color','w')
xlim([10 1000])
title('NoMOI cells')
xlabel('AP2-RFP fluorescence (a.u.)')
ylabel('Counts (#)')

figure,
histogram(M560.A_tot,15:5:150)
set(gcf, 'Color','w')
xlim([15 150])
title('NoMOI cells')
xlabel('AP2-RFP fluorescence (a.u.)')
ylabel('Counts (#)')

% Check spatial distributions and histograms for each cell
% (Use M560.ExpNum)
Exps = unique(M560.ExpNum);

figure,
for f = 1:numel(Exps)
    subplot(1,numel(Exps),f)
    Idx = find(M560.ExpNum==Exps(f));
    Idx2 = find(M560.A_tot(Idx)>10 & M560.A_tot(Idx)<100);
    Idx3 = Idx(Idx2);
    scatter3(M560.x_tot(Idx3)',M560.y_tot(Idx3)',M560.z_tot(Idx3)',M560.A_tot(Idx3)','r','filled')
    title(['Ex' num2str(Exps(f)) ' 10:100']) 
end
set(gcf,'Color','w')

% Check BleedThrough/Cross-Activation Values in Green-Channel
% pretty sure GFP is 2nd row, but check both 2nd and 3rd

IdxA = find(M560.A_tot(1,:)>10 & M560.A_tot(1,:)<200); % suitable RFP-pits 

figure, histogram(M560.A_tot(2,IdxA),1:5:200)
title('BleedThrough in NoMOI')
set(gcf,'Color','w')

figure, scatter(M560.A_tot(1,IdxA),M560.A_tot(2,IdxA))
xlabel('560-fluorescence')
ylabel('488-fluorescence')
title('BleedThrough in NoMOI')
set(gcf,'Color','w')

%% Association between EboV NoVHH and AP2 – at static state (fr1)

rtdir = 'G:\LLSM_VSV_VHH\SUM51 for hidde summary\SUM51 for hidde summary\CS4_NoVHH';
cd(rtdir)

% OBS: ExpNum = [10,11,12] instead of [15,16,18]
% FN = dir(rtdir); FN = FN([4 5 7]); numeration done with this

load('M488.mat')
Exps = unique(M488.ExpNum);

figure('Position',[34 593 1857 420]),
for f = 1:numel(Exps)
    subplot(1,numel(Exps),f)
    Idx = find(M488.ExpNum==Exps(f));
    Idx2 = find(M488.A_tot(2,Idx)>10 & M488.A_tot(2,Idx)<100); % suitable RFP-pits
    Idx3 = Idx(Idx2); % 
    scatter3(M488.x_tot(2,Idx3)',M488.y_tot(2,Idx3)',M488.z_tot(2,Idx3)',0.5*M488.A_tot(2,Idx3)','r','filled')
    IdxB = find(M488.A_tot(1,:)>200 & M488.A_tot(1,:)<3000); % MISTAKE!!
    % COMBINED ALL 
    hold on
    scatter3(M488.x_tot(1,IdxB)',M488.y_tot(1,IdxB)',M488.z_tot(1,IdxB)',0.01*M488.A_tot(1,IdxB)','g','filled')
    hold off
    title(['Ex' num2str(Exps(f)) ' 10:100']) 
end
set(gcf,'Color','w')
% OBSERVATION: Small Associations of Virus with AP2-RFP pits, check sizes

IdxA = find(M488.A_tot(2,:)>10 & M488.A_tot(2,:)<200); % suitable RFP-pits 

figure, histogram(M488.A_tot(1,IdxA),1:100:10000)
title('Associated GFP-EboV to AP2-RFP pits')
set(gcf,'Color','w')

% MASTER 560 for RFP, and M488 for GFP

load('M488.mat')
load('M560.mat')

Exps = unique(M488.ExpNum);

GFPmin = 200; GFPmax = 6000;
RFPmin = 10; RFPmax = 100;

figure('Position',[34 593 1857 420]),
for f = 1:numel(Exps)
    subplot(1,numel(Exps),f)
    % 488
    Idx = find(M488.ExpNum==Exps(f)); % 488
    Idx2 = find(M488.A_tot(1,Idx)>GFPmin & M488.A_tot(1,Idx)<GFPmax); % suitable RFP-pits
    Idx3 = Idx(Idx2); % 
    scatter3(M488.x_tot(1,Idx3)',M488.y_tot(1,Idx3)',M488.z_tot(1,Idx3)',0.01*M488.A_tot(1,Idx3)','g','filled')
    % 560
    IdxB = find(M560.ExpNum==Exps(f));
    IdxB2 = find(M560.A_tot(1,IdxB)>RFPmin & M560.A_tot(1,IdxB)<RFPmax);
    IdxB3 = IdxB(IdxB2); % 560
    hold on
    scatter3(M560.x_tot(1,IdxB3)',M560.y_tot(1,IdxB3)',M560.z_tot(1,IdxB3)',0.5*M560.A_tot(1,IdxB3)','r','filled')
    hold off
    title(['FalseEx' num2str(Exps(f)) 'GFP = ' num2str(GFPmin) ':' num2str(GFPmax) ', RFP = ' num2str(RFPmin) ':' num2str(RFPmax)]) 
end
set(gcf,'Color','w')

%% AP2-bound VSV

% M488
IdxA = find(M488.A_tot(2,:)>10 & M488.A_tot(2,:)<200); % suitable RFP-pits 

figure('Position',[681 567 623 532]),
scatter3(M488.x_tot(2,IdxA)',M488.y_tot(2,IdxA)',M488.z_tot(2,IdxA)',0.5*M488.A_tot(2,IdxA)','r','filled')
hold on
scatter3(M488.x_tot(1,IdxA)',M488.y_tot(1,IdxA)',M488.z_tot(1,IdxA)',0.01*M488.A_tot(1,IdxA)','g','filled')
set(gcf,'Color','w')
title('M488-based: RFP-selected: Association')

% 560
IdxB = find(M560.A_tot(2,:)>200 & M560.A_tot(2,:)<6000); % suitable GFP-units

figure('Position',[681 567 623 532]),
scatter3(M560.x_tot(1,IdxB)',M560.y_tot(1,IdxB)',M560.z_tot(1,IdxB)',0.5*M560.A_tot(1,IdxB)','r','filled')
hold on
scatter3(M560.x_tot(2,IdxB)',M560.y_tot(2,IdxB)',M560.z_tot(2,IdxB)',0.01*M560.A_tot(2,IdxB)','g','filled')
set(gcf,'Color','w')
title('M560-based: GFP-selected: Association')

%% M488-based, GFP-selected, 1FOV vis. Tot-quantification

% Plot first cell with four colors:
% 488 - not ass. (green)
% 560 - not ass. (red)
% 488 - ass. (cyan)
% 560 - ass. (mag)

rtdir = 'G:\LLSM_VSV_VHH\SUM51 for hidde summary\SUM51 for hidde summary\CS4_NoVHH';
cd(rtdir)

load('M488.mat')
load('M560.mat')

Exps = unique(M488.ExpNum);

Idx488gfp = find(M488.A_tot(1,:)>200 & M488.A_tot(1,:)<6000); % 488gfp
Idx488rfp = find(M488.A_tot(2,:)>10 & M488.A_tot(2,:)<100); % 488rfp

Idx560rfp = find(M560.A_tot(1,:)>10 & M560.A_tot(1,:)<100); % 560rfp
Idx560gfp = find(M560.A_tot(2,:)>200 & M560.A_tot(2,:)<6000); % 560gfp

I4GRa = intersect(Idx488gfp, Idx488rfp);
I4G = setdiff(Idx488gfp, Idx488rfp);

I4RGa = intersect(Idx560rfp, Idx560gfp);
I4R = setdiff(Idx560rfp, Idx560gfp);

% plot FOV1

IdxFOV1 = find(M488.ExpNum==Exps(1)); % FOV1

I4GRa1 = intersect(I4GRa, IdxFOV1); % 488ass in FOV 1
I4G1 = intersect(I4G, IdxFOV1); % 488gfp in FOV 1
I4RGa1 = intersect(I4RGa, IdxFOV1); % 560ass in FOV 1
I4R1 = intersect(I4R, IdxFOV1); % 560rfp in FOV 1

figure('Position',[681 567 623 532]),
scatter3(M560.x_tot(1,I4R1)',M560.y_tot(1,I4R1)',M560.z_tot(1,I4R1)',0.5*M560.A_tot(1,I4R1)','r','filled')
hold on
scatter3(M488.x_tot(1,I4G1)',M488.y_tot(1,I4G1)',M488.z_tot(1,I4G1)',0.01*M488.A_tot(1,I4G1)','g','filled')
scatter3(M560.x_tot(1,I4RGa1)',M560.y_tot(1,I4RGa1)',M560.z_tot(1,I4RGa1)',0.5*M560.A_tot(1,I4RGa1)','m','filled')
scatter3(M488.x_tot(1,I4GRa1)',M488.y_tot(1,I4GRa1)',M488.z_tot(1,I4GRa1)',0.01*M488.A_tot(1,I4GRa1)','c','filled')
hold off
set(gcf,'Color','w')
set(gca,'Zdir','reverse')
title('Single Cell - No VHH added')
view(88.6565, -46.2734)
grid off
set(gca,'xtick',[]), set(gca,'ytick',[]), set(gca,'ztick',[]) 
leg = legend('AP2 - alone','VSV EboV - alone','AP2 w/EboV','VSV EboV in AP2');
set(leg,'Box','off')

% Quantify all, base on M488 w/RFP-ass. or base on M560 w/GFP-ass.
pie([numel(setdiff(Idx488gfp,Idx488rfp)),numel(intersect(Idx488gfp,Idx488rfp))],{'VSV EboV alone','VSV in AP2-pits'})
ax = gca;
newColors = [0, 1, 0; 0, 1, 1];
ax.Colormap = newColors; 
set(gcf,'Color','w')
title('VSV G-EboV')

pie([numel(setdiff(Idx560rfp,Idx560gfp)),numel(intersect(Idx560rfp,Idx560gfp))],{'AP2 alone','AP2 w/VSV EboV'})
ax = gca;
newColors = [1, 0, 0; 1, 0, 1];
ax.Colormap = newColors; 
set(gcf,'Color','w')
title('AP2 pits')

%% G10 equivalent (frame 1, Pre-mixed, AP2-association)

rtdir = 'G:\LLSM_VSV_VHH\TEMPORARYDATA\CS3_G10_Hidde1';
% Original stacks damaged with single lacking files
% See DataProc_211102_VHH_EboV_data_in_cells.m
cd(rtdir)

load('M488.mat') % [488,560,642]
load('M560.mat') % [560,488,642]
Exps = unique(M488.ExpNum);

% Plot all 5 560-Cells with 488-GFP-VHH attached
Idx488gfp = find(M488.A_tot(1,:)>200 & M488.A_tot(1,:)<6000); % 488gfp

% FIRST CHECK THE CUTOFF VALUE FOR VHH-VALUES
figure, histogram(M488.A_tot(3,Idx488gfp),-500:10:2000)
set(gcf,'Color','w'),title('G10 VHH-values associated with VSV EboV')

% use VHH cut-off at 80

figure,
for g = 1:numel(Exps)
    subplot(1,numel(Exps),g)
    % Idx Exp
    ExpIdx = find(M488.ExpNum==Exps(g));
    % 488
    Idx488gfp = find(M488.A_tot(1,ExpIdx)>200 & M488.A_tot(1,ExpIdx)<6000); % 488gfp
    I488Exp = ExpIdx(Idx488gfp);
    % Modify VHH-values (everything below 80 == 0)
    VHHvals = M488.A_tot(3,I488Exp);
    VHHvals(VHHvals<80) = 0.1;
    % 560
    Idx560rfp = find(M560.A_tot(1,ExpIdx)>10 & M560.A_tot(1,ExpIdx)<100); % 560rfp
    I560Exp = ExpIdx(Idx560rfp);
    % plot them
    scatter3(M560.x_tot(1,I560Exp)',M560.y_tot(1,I560Exp)',M560.z_tot(1,I560Exp)',0.5*M560.A_tot(1,I560Exp)','r','filled')
    hold on
    scatter3(M488.x_tot(1,I488Exp)',M488.y_tot(1,I488Exp)',M488.z_tot(1,I488Exp)',0.01*M488.A_tot(1,I488Exp)','g','filled')
    scatter3(M488.x_tot(3,I488Exp)',M488.y_tot(3,I488Exp)',M488.z_tot(3,I488Exp)',0.1*VHHvals','b','filled')
    hold off
    % OBS: preventing negative VHH-values
    set(gca,'Zdir','reverse')
    view(98.0882,-64.9092)
    grid off
    title(['Ex' num2str(Exps(g))]) 
end
set(gcf,'Color','w')
    
%% G10: Histogram and pie-diagram (using only Ex10 as example, as there is very little precipitation outside cell)

ExpIdx488 = find(M488.ExpNum==Exps(1));
ExpIdx560 = find(M560.ExpNum==Exps(1));

% Pie Diagram
Idx488gfp = find(M488.A_tot(1,ExpIdx488)>200 & M488.A_tot(1,ExpIdx488)<6000); % 488gfp
Idx488rfp = find(M488.A_tot(2,ExpIdx488)>10 & M488.A_tot(2,ExpIdx488)<100); % 488rfp

Idx560rfp = find(M560.A_tot(1,ExpIdx560)>10 & M560.A_tot(1,ExpIdx560)<100); % 560rfp
Idx560gfp = find(M560.A_tot(2,ExpIdx560)>200 & M560.A_tot(2,ExpIdx560)<6000); % 560gfp

I4GRa = intersect(Idx488gfp, Idx488rfp);
I4G = setdiff(Idx488gfp, Idx488rfp);

I4RGa = intersect(Idx560rfp, Idx560gfp);
I4R = setdiff(Idx560rfp, Idx560gfp);

% Quantify all, base on M488 w/RFP-ass. or base on M560 w/GFP-ass.
figure,
pie([numel(setdiff(Idx488gfp,Idx488rfp)),numel(intersect(Idx488gfp,Idx488rfp))],{'VSV EboV alone','VSV in AP2-pits'})
ax = gca;
newColors = [0, 1, 0; 0, 1, 1];
ax.Colormap = newColors; 
set(gcf,'Color','w')
title('VSV G-EboV G10-VHH added')

% Histograms
IdxAss = ExpIdx488(I4GRa);
IdxNAss = ExpIdx488(I4G);

% Normalize Histogram -100:600
hcAss = histcounts(M488.A_tot(3,IdxAss),-20:10:600);
hcNAss = histcounts(M488.A_tot(3,IdxNAss),-20:10:600);

hcAssNorm = hcAss/max(hcAss);
hcNAssNorm = hcNAss/max(hcNAss);

figure,
bar(hcAssNorm)
hold on
bar(hcNAssNorm)
hold off
xticklabels({-20:100:600})
title('G10-VHH association to VSV EboV particles')
ylim([0 1.1])
ylabel('Normalized Counts (#)')
set(gcf,'Color','w')
xlabel('VHH-647 fluorescence (a.u.)')
leg = legend({'AP2-pit associated','Not AP2-associated'});
set(leg,'Box','Off')

%% G10: 488 vs. 647, with (BLUE) or without AP2 (RED)

% Idx488
Idx488gfp = find(M488.A_tot(1,ExpIdx488)>200 & M488.A_tot(1,ExpIdx488)<6000); % 488gfp
Idx488rfp = find(M488.A_tot(2,ExpIdx488)>10 & M488.A_tot(2,ExpIdx488)<100); % 488rfp

I4GRa = intersect(Idx488gfp, Idx488rfp);
I4G = setdiff(Idx488gfp, Idx488rfp);

IdxAss = ExpIdx488(I4GRa);
IdxNAss = ExpIdx488(I4G);

% total values
figure,
scatter(M488.A_tot(1,IdxAss),M488.A_tot(3,IdxAss),'r')
hold on
scatter(M488.A_tot(1,IdxNAss),M488.A_tot(3,IdxNAss),'b')
hold off
xlabel('VSV EboV fluorescence (a.u.)')
ylabel('VHH fluorescence (a.u.)')
set(gcf,'Color','w')
ylim([0 1800])

% relative values
figure,
scatter(M488.A_tot(1,IdxAss),M488.A_tot(3,IdxAss)./M488.A_tot(1,IdxAss),'r')
hold on
scatter(M488.A_tot(1,IdxNAss),M488.A_tot(3,IdxNAss)./M488.A_tot(1,IdxNAss),'b')
hold off
xlabel('VSV EboV fluorescence (a.u.)')
ylabel('Ratio (VHH/VSV)')
set(gcf,'Color','w')
ylim([0 2.5])

%% G84 Overview of Cells

rtdir = 'G:\LLSM_VSV_VHH\TEMPORARYDATA\CS_G84';
cd(rtdir)

load('M488.mat') % [488,560,642]
load('M560.mat') % [560,488,642], needed to add from DataProc_211102
Exps = unique(M488.ExpNum);

% Plot all 5 560-Cells with 488-GFP-VHH attached
Idx488gfp = find(M488.A_tot(1,:)>200 & M488.A_tot(1,:)<6000); % 488gfp

% FIRST CHECK THE CUTOFF VALUE FOR VHH-VALUES
figure, histogram(M488.A_tot(3,Idx488gfp),-500:10:2000)
set(gcf,'Color','w'),title('G10 VHH-values associated with VSV EboV')
% use VHH cut-off at 80

figure,
for g = 1:numel(Exps)
    subplot(1,numel(Exps),g)
    % Idx Exp
    ExpIdx = find(M488.ExpNum==Exps(g));
    % 488
    Idx488gfp = find(M488.A_tot(1,ExpIdx)>200 & M488.A_tot(1,ExpIdx)<6000); % 488gfp
    I488Exp = ExpIdx(Idx488gfp);
    % Modify VHH-values (everything below 80 == 0)
    VHHvals = M488.A_tot(3,I488Exp);
    VHHvals(VHHvals<80) = 0.1;
    % 560
    Idx560rfp = find(M560.A_tot(1,ExpIdx)>10 & M560.A_tot(1,ExpIdx)<100); % 560rfp
    I560Exp = ExpIdx(Idx560rfp);
    % plot them
    scatter3(M560.x_tot(1,I560Exp)',M560.y_tot(1,I560Exp)',M560.z_tot(1,I560Exp)',0.5*M560.A_tot(1,I560Exp)','r','filled')
    hold on
    scatter3(M488.x_tot(1,I488Exp)',M488.y_tot(1,I488Exp)',M488.z_tot(1,I488Exp)',0.01*M488.A_tot(1,I488Exp)','g','filled')
    scatter3(M488.x_tot(3,I488Exp)',M488.y_tot(3,I488Exp)',M488.z_tot(3,I488Exp)',0.1*VHHvals','b','filled')
    hold off
    % OBS: preventing negative VHH-values
    set(gca,'Zdir','reverse')
    view(98.0882,-64.9092)
    grid off
    title(['Ex' num2str(Exps(g))]) 
end
set(gcf,'Color','w')

%% G84: Histogram and pie-diagram 

% using only Ex5 or Ex9 as example, as there is very little precipitation outside cell

ExpIdx488 = find(M488.ExpNum==Exps(1)); % ind = 1, Exp5, 
ExpIdx560 = find(M560.ExpNum==Exps(1)); % ind = 5, Exp9,

% Pie Diagram
Idx488gfp = find(M488.A_tot(1,ExpIdx488)>200 & M488.A_tot(1,ExpIdx488)<6000); % 488gfp
Idx488rfp = find(M488.A_tot(2,ExpIdx488)>10 & M488.A_tot(2,ExpIdx488)<100); % 488rfp

Idx560rfp = find(M560.A_tot(1,ExpIdx560)>10 & M560.A_tot(1,ExpIdx560)<100); % 560rfp
Idx560gfp = find(M560.A_tot(2,ExpIdx560)>200 & M560.A_tot(2,ExpIdx560)<6000); % 560gfp

I4GRa = intersect(Idx488gfp, Idx488rfp);
I4G = setdiff(Idx488gfp, Idx488rfp);

I4RGa = intersect(Idx560rfp, Idx560gfp);
I4R = setdiff(Idx560rfp, Idx560gfp);

% Quantify all, base on M488 w/RFP-ass. or base on M560 w/GFP-ass.
figure,
pie([numel(setdiff(Idx488gfp,Idx488rfp)),numel(intersect(Idx488gfp,Idx488rfp))],{'VSV EboV alone','VSV in AP2-pits'})
ax = gca;
newColors = [0, 1, 0; 0, 1, 1];
ax.Colormap = newColors; 
set(gcf,'Color','w')
title('VSV G-EboV G84-VHH added')

% Histograms
IdxAss = ExpIdx488(I4GRa);
IdxNAss = ExpIdx488(I4G);

% Normalize Histogram -100:600
hcAss = histcounts(M488.A_tot(3,IdxAss),-20:10:600);
hcNAss = histcounts(M488.A_tot(3,IdxNAss),-20:10:600);

hcAssNorm = hcAss/max(hcAss);
hcNAssNorm = hcNAss/max(hcNAss);

figure,
bar(hcAssNorm)
hold on
bar(hcNAssNorm)
hold off
xticklabels({-20:100:600})
title('G84-VHH association to VSV EboV particles')
ylim([0 1.1])
ylabel('Normalized Counts (#)')
set(gcf,'Color','w')
xlabel('VHH-647 fluorescence (a.u.)')
leg = legend({'AP2-pit associated','Not AP2-associated'});
set(leg,'Box','Off')
