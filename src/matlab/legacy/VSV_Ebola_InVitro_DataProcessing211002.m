%% VSV MeGFP Ebola SDCM

% Data Processing 2021, starting over


%% 1: Check VSV MeGFP Ebola-G in vitro (SDCM)

% C:\Users\Ras\Dropbox\TK lab\PROJECTS\VSV Ebola (Jason)\SDCM\190324 VSVEboV with VHHs 44nM\VSVEboV_alone
% Ex24:Ex36

% Check Homogeneity of VSV MeGFP population, subplots

rtdir = 'C:\Users\Ras\Dropbox\TK lab\PROJECTS\VSV Ebola (Jason)\SDCM\190324 VSVEboV with VHHs 44nM\VSVEboV_alone';

FNs = dir(rtdir);

n = 0;
figure,
for g = 3:numel(FNs)
    n = n+1;
    load(fullfile(rtdir,FNs(g).name,'Ch0','Detection','detection_v2.mat'))
    VSVs = frameInfo.A(1,:);
    subplot(numel(FNs)-2,1,n)
    histogram(VSVs,1:1:200)
    xlim([0 80])
end
set(gcf,'Color','w')
suptitle('VSV MeGFP G-Ebola integrity')
xlabel('Fluorescence (a.u.)')

%% Pool MeGFP and AF647 fluorescence from individual experiments

% from 190324 VSVEboV with VHHs 44nM
rtdir = 'C:\Users\Ras\Dropbox\TK lab\PROJECTS\VSV Ebola (Jason)\SDCM\190324 VSVEboV with VHHs 44nM';

FldNames = dir(rtdir);
FldNames = FldNames([3 5 7 10]);

for g = 1:numel(FldNames)
    ExpNs = dir(fullfile(FldNames(g).folder,FldNames(g).name));
    ExpNs = ExpNs(3:end);
    Fluor = [];
    for h = 1:numel(ExpNs)
        load(fullfile(FldNames(g).folder,FldNames(g).name,ExpNs(h).name,'Ch0','Detection','detection_v2.mat'))
        F_temp = frameInfo.A;
        Fluor = cat(2,Fluor,F_temp);
    end
    VSV_binding.(FldNames(g).name).A = Fluor;
end

% from 190322 VSVEboV with VHHs 440nM
rtdir = 'C:\Users\Ras\Dropbox\TK lab\PROJECTS\VSV Ebola (Jason)\SDCM\190322 VSVEboV with VHHs 440nM';

FldNames = dir(rtdir);
FldNames = FldNames(3:6);

for g = 1:numel(FldNames)
    ExpNs = dir(fullfile(FldNames(g).folder,FldNames(g).name));
    ExpNs = ExpNs(3:end);
    Fluor = [];
    for h = 1:numel(ExpNs)
        load(fullfile(FldNames(g).folder,FldNames(g).name,ExpNs(h).name,'Ch0','Detection','detection_v2.mat'))
        F_temp = frameInfo.A;
        Fluor = cat(2,Fluor,F_temp);
    end
    VSV_binding.(FldNames(g).name).A = Fluor;
end

save(fullfile('C:\Users\Ras\Dropbox\TK lab\PROJECTS\VSV Ebola (Jason)\SDCM','VSV_binding.mat'),'VSV_binding')


%% Histograms of MeGFP sizes, and AF647/MeGFP-ratios

% OBS!! Compensate VHH-crossbleeding effect of AF647 in VSV fluorescence

figure,
% VSV fluorescence with 44nM
subplot(5,4,1)
histogram(VSV_binding.VSVEboV_alone.A(1,:),1:1:200,'FaceColor',[0 1 0])
title('VSV alone')
subplot(5,4,2)
histogram(VSV_binding.G10_44nM.A(1,:),1:1:200,'FaceColor',[0 1 0])
title('G10 44nM')
subplot(5,4,3)
histogram(VSV_binding.G84_44nM.A(1,:),1:1:200,'FaceColor',[0 1 0])
title('G84 44nM')
subplot(5,4,4)
title('G68flu 44nM')
% VHH-binding ratio with 44nM
subplot(5,4,5)
histogram(VSV_binding.VSVEboV_alone.A(2,:)./VSV_binding.VSVEboV_alone.A(1,:),'FaceColor',[0 0 1])
title('VSV alone (VHH)')
subplot(5,4,6)
histogram(VSV_binding.G10_44nM.A(2,:)./VSV_binding.G10_44nM.A(1,:),'FaceColor',[0 0 1])
title('G10 44nM ')
subplot(5,4,7)
histogram(VSV_binding.G84_44nM.A(2,:)./VSV_binding.G84_44nM.A(1,:),'FaceColor',[0 0 1])
title('G84 44nM')
subplot(5,4,8)
title('G68flu 44nM')
% VSV fluorescence with 440nM
subplot(5,4,9)
title('VSV alone')
subplot(5,4,10)
histogram(VSV_binding.G10_440nM.A(1,:),'FaceColor',[0 1 0])
title('G10 440nM')
subplot(5,4,11)
histogram(VSV_binding.G84_440nM.A(1,:),'FaceColor',[0 1 0])
title('G84 440nM')
subplot(5,4,12)
histogram(VSV_binding.G68flu_440nM.A(1,:),'FaceColor',[0 1 0])
title('G68flu 440nM')
% VHH-binding ratio with 440nM
subplot(5,4,13)
title('VSV alone (VHH)')
subplot(5,4,14)
histogram(VSV_binding.G10_440nM.A(2,:)./VSV_binding.G10_440nM.A(1,:),'FaceColor',[0 0 1])
title('G10 440nM (VHH)')
subplot(5,4,15)
histogram(VSV_binding.G84_440nM.A(2,:)./VSV_binding.G84_440nM.A(1,:),'FaceColor',[0 0 1])
title('G84 440nM (VHH)')
subplot(5,4,16)
histogram(VSV_binding.G68flu_440nM.A(2,:)./VSV_binding.G68flu_440nM.A(1,:),'FaceColor',[0 0 1])
title('G68flu 440nM (VHH)')
% VSV fluorescence with 440nM increased exposure
subplot(5,4,17)
title('VSV alone')
subplot(5,4,18)
histogram(VSV_binding.G10_440nM_incr_exp.A(1,:),'FaceColor',[0 1 0])
title('G10 440nM Incr. Exp')
subplot(5,4,19)
histogram(VSV_binding.G84_440nM_incr_Exp.A(1,:),'FaceColor',[0 1 0])
title('G84 440nM incr. Exp')
subplot(5,4,20)
title('G68flu 440nM')
set(gcf,'Color','w')

%% CHECK DATA FROM G84 440nM

FNs = dir('C:\Users\Ras\Dropbox\TK lab\PROJECTS\VSV Ebola (Jason)\SDCM\190322 VSVEboV with VHHs 440nM\G84_440nM');
figure('Position',[544 78 1229 1032])
for f = 3:12
    load(fullfile(FNs(f).folder,FNs(f).name,'Ch0','Detection','detection_v2'))
    % try with indices for MeGFP > 12000
    idx = find(frameInfo.A(1,:)>12000);
    scatter(frameInfo.x(1,idx),frameInfo.y(1,idx),0.01*frameInfo.A(1,idx),'g')
    hold on
    SlaveSz = frameInfo.A(2,idx);
    SlaveSz(SlaveSz<=0) = 1;
    scatter(frameInfo.x(2,idx),frameInfo.y(2,idx),0.01*SlaveSz,'b')
    title(FNs(f).name)
    set(gcf,'Color','w')
    hold off
    pause
end

XY = [234.58 247.271];
XYs(:,1) = frameInfo.x(1,:)';
XYs(:,2) = frameInfo.y(1,:)';

Dist = sqrt(sum(bsxfun(@minus, XYs, XY).^2,2));
[~,I] = min(Dist);
MeGFP = frameInfo.A(1,I);
AF647 = frameInfo.A(2,I);

% then do distance-analysis between AF_647s and MeGFP_488s

% Distance and fluorescence levels

clearvars -except FNs
XYgfp = [];
XY647 = [];
GFP = [];
AF647 = [];
for f = 3:12
    clear XYgfp_t XY647_t
    load(fullfile(FNs(f).folder,FNs(f).name,'Ch0','Detection','detection_v2'))
    % try with indices for MeGFP > 12000
    idx = find(frameInfo.A(1,:)>12000);
    % and indices with AF647 < 20000
    idxAF = find(frameInfo.A(2,:)<20000);
    idxtot = intersect(idx,idxAF);
    % XYs and Fluor
    XYgfp_t(:,1) = frameInfo.x(1,idxtot)';
    XYgfp_t(:,2) = frameInfo.y(1,idxtot)';
    XY647_t(:,1) = frameInfo.x(2,idxtot)';
    XY647_t(:,2) = frameInfo.y(2,idxtot)';
    % concatenate
    XYgfp = cat(1,XYgfp,XYgfp_t);
    XY647 = cat(1,XY647,XY647_t);
    GFP = cat(1,GFP,frameInfo.A(1,idxtot)');
    AF647 = cat(1,AF647,frameInfo.A(2,idxtot)');
end

Dist = sqrt(sum(bsxfun(@minus, XYgfp, XY647).^2,2));

% plot it all (limits

figure,
scatter(GFP,AF647)
xlabel('GFP'), ylabel('AF647'), title('G84 440nM')
set(gcf,'Color','w')

figure,
histogram(GFP/15000,200,'FaceColor',[0 1 0])
xlabel('GFP (a.u.)')
set(gcf,'Color','w')
title('G84 440nM - 15000au/SglVir')

% Idx into SglVir ranges
figure
for g = 1:6
    Idx = find(GFP>(g-0.5)*15000 & GFP<=(g+0.5)*15000);
    SglVir(g).Number = g;
    SglVir(g).AF647 = AF647(Idx);
    SglVir(g).AF647_mean = mean(AF647(Idx));
    SglVir(g).AF647_std = std(AF647(Idx));
end

% Boxplot
figure,
boxplot([SglVir(1).AF647; SglVir(2).AF647; SglVir(3).AF647; SglVir(4).AF647;...
    SglVir(5).AF647; SglVir(6).AF647],[repmat({'One Virion'},numel(SglVir(1).AF647),1);...
    repmat({'Two'},numel(SglVir(2).AF647),1);repmat({'Three'},numel(SglVir(3).AF647),1);...
    repmat({'Four'},numel(SglVir(4).AF647),1);repmat({'Five'},numel(SglVir(5).AF647),1);...
    repmat({'Six'},numel(SglVir(6).AF647),1)])
title('G84 440nM - AF647 binding')
set(gcf,'Color','w')


%% Bleaching analysis (?)

clearvars -except VSV_binding

cd('C:\Users\Ras\Dropbox\TK lab\PROJECTS\VSV Ebola (Jason)\SDCM\190324 VSVEboV with VHHs 44nM\G10_44nM_bleach\Ex23\Ch0\Detection')

figure,
for h = 1:numel(frameInfo)
    subplot(2,1,1)
    plot(frameInfo(h).A(1,:))
    subplot(2,1,2)
    plot(frameInfo(h).A(2,:))
    pause
end
% Detection, not trackin

% Check detection 3 in frame 1
figure,
scatter(frameInfo(1).x(1,:),frameInfo(1).y(1,:),'r')
set(gcf,'Color','w')
title('FOV 1, Detection 3')
hold on
scatter(frameInfo(1).x(1,3),frameInfo(1).y(1,3),'k')
hold off

figure,
for n = 1:numel(frameInfo(1).A(1,:))
    % Follow with X and Y
    XY(1) = frameInfo(1).x(1,n);
    XY(2) = frameInfo(1).y(1,n);
    % find index for XY in each frame, extract PSF-amplitude
    for g = 1:numel(frameInfo)
        clear XYs
        XYs(:,1) = frameInfo(g).x(1,:)';
        XYs(:,2) = frameInfo(g).y(1,:)';
        % find closest point in each frame
        Dist = sqrt(sum(bsxfun(@minus, XYs, XY).^2,2));
        [~,I] = min(Dist);
        MeGFP(g) = frameInfo(g).A(1,I);
    end
    plot(MeGFP)
    title(['Particle #' num2str(n)])
    pause
end

%% Single Molecules and Single Viruses

% if SM is really ~45au:
SglMol = VSV_binding.G84_44nM.A(1,:)/45;

figure, histogram(SglMol)
title('G84_44nM')

% Maybe 3160au/SglVir

%% Cross-bleeding of fluorescence from 647 to eGFP channel

figure,
scatter(VSV_binding.G84_440nM.A(2,:), VSV_binding.G84_440nM.A(1,:))
xlabel('AF647')
ylabel('MeGFP')
title('Cross-Bleeding G84 440nM')
ylim([0 40000])

%% 44nM alone

figure,
scatter(VSV_binding.G10_44nM.A(1,:),VSV_binding.G10_44nM.A(2,:))
title('G10 44nM'), xlabel('GFP'),ylabel('AF647'), set(gcf,'Color','w')

figure,
scatter(VSV_binding.G84_44nM.A(1,:),VSV_binding.G84_44nM.A(2,:))
title('G84 44nM'), xlabel('GFP'),ylabel('AF647'), set(gcf,'Color','w')

%% Detections of 647nm

FNs = fieldnames(Detections);
for t = 1:5
    Mean_647(t) = mean(Detections.(FNs{t}).Ch1);
    Mean_647_corr(t) = mean(Detections.(FNs{t}).Ch1_corr);
    Med_647(t) = median(Detections.(FNs{t}).Ch1);
    Med_647_corr(t) = median(Detections.(FNs{t}).Ch1_corr);
end
x = [5 7.5 10 15 20];

figure,
subplot(2,2,1)
scatter(x,Mean_647)
title('Mean 647')
subplot(2,2,2)
scatter(x,Mean_647_corr)
title('Mean 647 corr')
subplot(2,2,3)
scatter(x,Med_647)
title('Median 647')
xlabel('time (sec)')
subplot(2,2,4)
scatter(x,Med_647_corr)
title('Median 647 corr')
xlabel('time (sec)')
set(gcf,'Color','w')

% continue with median 647
Ratio = Med_647./x; % 130au/(SM*sec), at lp-80

Ratio = (Med_647(2:3)-Med_647(1:2))./2.5; % 120au/(SM*sec), at lp-80

% 100ms Exp = 12au/SM

%% FINAL FIGURES FOR 44nM

% 1) Make histograms, 2) Find representative images (search artificially)

% 1) Make histograms and SM-calibration

orange      = [1 0.5 0];
blue        = [0 0.5 1];
green       = [0 0.6 0.3];
red         = [1 0.2 0.2];

f = polyfit(x(1:3),Med_647(1:3),1);
Fxy = polyval(f,[4 11]);

% Single Molecule Calibration (In Vitro)
figure('Position',[230 632 585 441])
scatter(x(1),Med_647(1),'MarkerFaceColor',[0.25 0.05 0.05],'MarkerEdgeColor',[0.75 0.75 0.75])
hold on
scatter(x(2),Med_647(2),'MarkerFaceColor',[0.67 0.13 0.13],'MarkerEdgeColor',[0.75 0.75 0.75])
scatter(x(3),Med_647(3),'MarkerFaceColor',[1 0.2 0.2],'MarkerEdgeColor',[0.75 0.75 0.75])
xlim([4 11])
xlabel('Exposure Time (sec)')
ylabel('Fluorescence (a.u.)')
set(gcf,'Color','w')
hold on
plot([4 11],Fxy,'k--')
title('Single Molecule Calibration')

% histograms for SM-cal
figure('Position',[681 927 610 172])
histogram(Detections.G84_10s.Ch1,800,'FaceColor',[1 0.2 0.2])
hold on
histogram(Detections.G84_07p5.Ch1,800,'FaceColor',[0.67 0.13 0.13])
histogram(Detections.G84_05s.Ch1,800,'FaceColor',[0.25 0.05 0.05])
xlim([200 2000])
set(gca,'Xdir','reverse')
set(gcf,'Color','w')
xticks([])
yticks([])

% Single Virions Alone
figure('Position',[230 632 585 441])
% idx = find(VSV_binding.VSVEboV_alone.A
histogram(VSV_binding.VSVEboV_alone.A(1,:),0:1:80,'FaceColor',[0 0.6 0.3],'EdgeColor',[0.75 0.75 0.75])

[N,E] = histcounts(VSV_binding.VSVEboV_alone.A(1,:),0:1:80);

f = fit(E(2:end)',N','gauss1');

data = VSV_binding.VSVEboV_alone.A(1,:);

idx = find(data<80);

data = data(idx);

mu = mean(data);
sd = std(data);
ndfcn = @(mu,sd,x) exp(-(x-mu).^2 ./ (2*sd^2)) /(sd*sqrt(2*pi));    % Standard Normal Distribution
[hc,edges] = histcounts(data, 25);                                  % Histogram
ctrs = edges(1:length(edges)-1) + mean(diff(edges))/2;              % Calculate Centres
ctrsx = linspace(min(ctrs), max(ctrs));                             % High-Resolution Vector
sdnd = ndfcn(mu,sd,ctrsx);                                          % Calculate Standard Normal Distribution

figure('Position',[230 632 585 441])
bar(ctrs, hc,'FaceColor',[0 0.6 0.3])                                                       % Plot Histogram
hold on
plot(ctrsx, sdnd*max(hc)/max(sdnd), 'k--', 'LineWidth',1)            % Plot Scaled Standard Normal Distribution
hold off
xlim([0 80])
set(gcf,'Color','w')
title('Single Virion Calibration')
xlabel('Fluorescence (a.u.)')
ylabel('Counts')

% Histograms of Single Virions and VHH-numbers

% G10
data = VSV_binding.G10_44nM.A(1,:);
idx = find(data<90);
data = data(idx)/70;

mu = mean(data);
sd = std(data);
ndfcn = @(mu,sd,x) exp(-(x-mu).^2 ./ (2*sd^2)) /(sd*sqrt(2*pi));    % Standard Normal Distribution
[hc,edges] = histcounts(data, 25);                                  % Histogram
ctrs = edges(1:length(edges)-1) + mean(diff(edges))/2;              % Calculate Centres
ctrsx = linspace(min(ctrs), max(ctrs));                             % High-Resolution Vector
sdnd = ndfcn(mu,sd,ctrsx);      

figure('Position',[230 632 585 441])
histogram(VSV_binding.G10_44nM.A(1,:)/70,0:0.05:2,'FaceColor',[0 0.6 0.3],'EdgeColor',[0.25 0.25 0.25])
hold on
plot(ctrsx, sdnd*max(hc)/max(sdnd)*2, 'k--', 'LineWidth',1)           % Plot Scaled Standard Normal Distribution
hold off
set(gcf,'Color','w')
title('VSV MeGFP G-EboV w/G10')
xlabel('Single Virions')
ylabel('Counts')

% G84
data = VSV_binding.G84_44nM.A(1,:);
idx = find(data<70);
data = data(idx)/50;

mu = mean(data);
sd = std(data);
ndfcn = @(mu,sd,x) exp(-(x-mu).^2 ./ (2*sd^2)) /(sd*sqrt(2*pi));    % Standard Normal Distribution
[hc,edges] = histcounts(data, 25);                                  % Histogram
ctrs = edges(1:length(edges)-1) + mean(diff(edges))/2;              % Calculate Centres
ctrsx = linspace(min(ctrs), max(ctrs));                             % High-Resolution Vector
sdnd = ndfcn(mu,sd,ctrsx);   

figure('Position',[230 632 585 441])
histogram(VSV_binding.G84_44nM.A(1,:)/50,0:0.05:2,'FaceColor',[0 0.6 0.3],'EdgeColor',[0.25 0.25 0.25])
hold on
plot(ctrsx, sdnd*max(hc)/max(sdnd)*1.5, 'k--', 'LineWidth',1)           % Plot Scaled Standard Normal Distribution
hold off
set(gcf,'Color','w')
title('VSV MeGFP G-EboV w/G84')
xlabel('Single Virions')
ylabel('Counts')

% G68
data = VSV_binding.G68flu_440nM.A(1,:);
idx = find(data<16);
data = data(idx)/11;

mu = mean(data);
sd = std(data);
ndfcn = @(mu,sd,x) exp(-(x-mu).^2 ./ (2*sd^2)) /(sd*sqrt(2*pi));    % Standard Normal Distribution
[hc,edges] = histcounts(data, 25);                                  % Histogram
ctrs = edges(1:length(edges)-1) + mean(diff(edges))/2;              % Calculate Centres
ctrsx = linspace(min(ctrs), max(ctrs));                             % High-Resolution Vector
sdnd = ndfcn(mu,sd,ctrsx);   

figure('Position',[230 632 585 441])
histogram(VSV_binding.G68flu_440nM.A(1,:)/11,0:0.05:2,'FaceColor',[0 0.6 0.3],'EdgeColor',[0.25 0.25 0.25])
hold on
plot(ctrsx, sdnd*max(hc)/max(sdnd)*1.3, 'k--', 'LineWidth',1)           % Plot Scaled Standard Normal Distribution
hold off
set(gcf,'Color','w')
title('VSV MeGFP G-EboV w/G68')
xlabel('Single Virions')
ylabel('Counts')

%% VHH plots

% G10 AF647 : VHH/SV (12au/VHH)

data = VSV_binding.G10_44nM.A(1,:);
idx = find(data<140);

VHH_G10 = VSV_binding.G10_44nM.A(2,idx);
VHH_G10 = VHH_G10(VHH_G10>11);
VHH_G10 = VHH_G10/12;

figure('Position',[230 632 585 441])
histogram(VHH_G10,0:1:175,'FaceColor',[0.67 0.13 0.13])
ylim([0 15])
xlabel('VHHs per Single Virion')
ylabel('Counts')
title('G10 VHH per Virion')
set(gcf,'Color','w')

% G84 AF647 : VHH/SV (12au/VHH)

data = VSV_binding.G84_44nM.A(1,:);
idx = find(data<100);

VHH_G84 = VSV_binding.G84_44nM.A(2,idx);
VHH_G84 = VHH_G84(VHH_G84>11);
VHH_G84 = VHH_G84/12;

figure('Position',[230 632 585 441])
histogram(VHH_G84,0:1:175,'FaceColor',[0.67 0.13 0.13])
ylim([0 25])
xlabel('VHHs per Single Virion')
ylabel('Counts')
title('G84 VHH per Virion')
set(gcf,'Color','w') 

% G84 AF647 : VHH/SV (12au/VHH)

data = VSV_binding.G68flu_440nM.A(1,:);
idx = find(data<22);

VHH_G68 = VSV_binding.G68flu_440nM.A(2,idx);
VHH_G68 = VHH_G68(VHH_G68>0);
VHH_G68 = VHH_G68/12;

figure('Position',[230 632 585 441])
histogram(VHH_G68,0:1:175,'FaceColor',[0.67 0.13 0.13])
% ylim([0 25])
xlabel('VHHs per Single Virion')
ylabel('Counts')
title('G68 VHH per Virion')
set(gcf,'Color','w') 

%% Locate positions for Representative Images (!)

% Plot only VSVs within a reasonable range, and VHHs of reasonable number

% VSVs < 2

% VHHs 

%% G10 - VSVs < 140, VHH ~50, ~12*50 = 600 a.u.

RtDir = 'C:\Users\Ras\Dropbox\TK lab\PROJECTS\VSV Ebola (Jason)\SDCM\190324 VSVEboV with VHHs 44nM\G10_44nM';
FldNMs = dir(RtDir);

for g = 3:numel(FldNMs)
    load(fullfile(FldNMs(g).folder,FldNMs(g).name,'Ch0','Detection','detection_v2.mat'))
    % find criteria
    Idx1 = find(frameInfo.A(1,:)<140);
    Idx2 = find(frameInfo.A(2,:)>400 & frameInfo.A(2,:)<800);
    Idx3 = intersect(Idx1,Idx2);
    % plot them
    figure('Position',[833 464 620 592])
    scatter(frameInfo.x(1,Idx3),frameInfo.y(1,Idx3),frameInfo.A(1,Idx3)*0.5,'g')
    hold on
    scatter(frameInfo.x(2,Idx3),frameInfo.y(2,Idx3),frameInfo.A(2,Idx3)*0.25,'b')
    % Navigation
    % [~,I] = maxk(frameInfo.A(1,:),140);
    % scatter(frameInfo.x(1,I),frameInfo.y(1,I),frameInfo.A(1,I)*0.05,'MarkerFaceColor',[1 0 0])
    hold off
    xlim([1 512]), ylim([1 512])
    title(FldNMs(g).name)
    pause
end








