%% FIGURES FOR Ebola Paper

% Using Curetonpaper as inspiration:

% https://pubmed.ncbi.nlm.nih.gov/20941355/

%% In vitro data

rtdir = 'C:\Users\Ras\Dropbox\TK lab\PROJECTS\VSV Ebola (Jason)\SDCM';

load(fullfile(rtdir,'VSV_binding.mat'))

%% Normalized Virion Fluorescence

% Single Virions Alone
figure('Position',[230 632 585 441])
% idx = find(VSV_binding.VSVEboV_alone.A
histogram(VSV_binding.VSVEboV_alone.A(1,:),0:1:80,'FaceColor',[0 0.6 0.3],'EdgeColor',[0.75 0.75 0.75])

[N,E] = histcounts(VSV_binding.VSVEboV_alone.A(1,:),0:1:80);

f = fit(E(2:end)',N','gauss1');

data = VSV_binding.VSVEboV_alone.A(1,:);

idx = find(data<80);

% raw data
data = data(idx);

mu = mean(data);
sd = std(data);
ndfcn = @(mu,sd,x) exp(-(x-mu).^2 ./ (2*sd^2)) /(sd*sqrt(2*pi));    % Standard Normal Distribution
[hc,edges] = histcounts(data, 25);                                  % Histogram
ctrs = edges(1:length(edges)-1) + mean(diff(edges))/2;              % Calculate Centres
ctrsx = linspace(min(ctrs), max(ctrs));                             % High-Resolution Vector
sdnd = ndfcn(mu,sd,ctrsx);                                          % Calculate Standard Normal Distribution

% norm data
dataN = data/mu;

muN = mean(dataN);
sdN = std(dataN);
ndfcn = @(muN,sdN,x) exp(-(x-muN).^2 ./ (2*sdN^2)) /(sdN*sqrt(2*pi));    % Standard Normal Distribution
[hc,edges] = histcounts(dataN, 25);                                  % Histogram
ctrs = edges(1:length(edges)-1) + mean(diff(edges))/2;              % Calculate Centres
ctrsx = linspace(min(ctrs), max(ctrs));                             % High-Resolution Vector
sdnd = ndfcn(muN,sdN,ctrsx); 

figure('Position',[230 632 585 441])
bar(ctrs, round(hc/10),'FaceColor',[0 0.6 0.3],'EdgeColor',[0.75 0.75 0.75])                                                       % Plot Histogram
hold on
plot(ctrsx, sdnd*max(hc)/(10*max(sdnd)), 'k--', 'LineWidth',1)            % Plot Scaled Standard Normal Distribution
hold off
xlim([0 2])
set(gcf,'Color','w')
title('Single Virion Fluorescence')
xlabel('GFP Fluorescence (a.u.)')
ylabel('Counts')
leg = legend('n = 205');
set(leg,'box','off')

%% BeeSwarm plots with number of VHH AF647

% G10 44nM
data = VSV_binding.G10_44nM.A(1,:);
idx = find(data<140);

VHH_G10 = VSV_binding.G10_44nM.A(2,idx);
VHH_G10 = VHH_G10(VHH_G10>11);
VHH_G10 = VHH_G10/12;

% G84
data = VSV_binding.G84_44nM.A(1,:);
idx = find(data<100);

VHH_G84 = VSV_binding.G84_44nM.A(2,idx);
VHH_G84 = VHH_G84(VHH_G84>11);
VHH_G84 = VHH_G84/12;

% G68
data = VSV_binding.G68flu_440nM.A(1,:);
idx = find(data<22);

VHH_G68 = VSV_binding.G68flu_440nM.A(2,idx);
VHH_G68 = VHH_G68(VHH_G68>0);
VHH_G68 = VHH_G68/2;

%% Make BeeSwarm
T84_1 = find(VHH_G84>20); %90
T84_2 = randi(859,1,325);

T10 = randi(561,1,230);
T84 = [T84_1 T84_2];

figure,
g = beeswarm([1*ones(230,1); 2*ones(415,1); 3*ones(339,1)],[VHH_G10(T10)'; VHH_G84(T84)'; VHH_G68'],'MarkerFaceColor',[0.67 0.13 0.13],'corral_style','gutter');
ylim([1 150])

xticks([1 2 3])
xticklabels({'G10','G84','G68'})

mu10 = mean(VHH_G10(T10)); sem10 = std(VHH_G10(T10))/sqrt(length(VHH_G10(T10)));
mu84 = mean(VHH_G84(T84)); sem84 = std(VHH_G84(T84))/sqrt(length(VHH_G84(T84)));
Mu68 = mean(VHH_G68); sem68 = std(VHH_G68)/sqrt(length(VHH_G68));

hold on
plot([0.85 1.15],[mu10 mu10],'k')
errorbar(1,mu10,sem10,'k')
plot([1.85 2.15],[mu84 mu84],'k')
errorbar(2,mu84,sem84,'k')
plot([2.85 3.15],[Mu68 Mu68],'k')
errorbar(3,Mu68,sem68,'k')
hold off

set(gcf,'Color','w')
title('VHH bound per virion')
ylabel('Counts (#)')

text(0.85,145,'n = 230')
text(1.85,100,'n = 415')
text(2.85,20,'n = 339')

%% CELLULAR DATA

%% MAKE CELLSHPs

% from: 'G:\LLSM_VSV_VHH\SUM51 for hidde summary\SUM51 for hidde summary\DataAnalysis_FULL.m'

load('G:\LLSM_VSV_VHH\TEMPORARYDATA\CS3_G10_Hidde1\data560master.mat')
dataM560 = data;

for dataN = 1:numel(dataM560)
    a = strsplit(dataM560(dataN).source, filesep);
    rtdir = strjoin(a(1:end-2), filesep);
    cd(rtdir)
    % load AP2-points
    load([dataM560(dataN).source 'Analysis' filesep 'Detection3D.mat'])
    x560  = frameInfo(1).x(1,:);
    y560  = frameInfo(1).y(1,:);
    z560  = frameInfo(1).z(1,:);
    A560  = frameInfo(1).A(1,:);
    A560rel = A560/max(A560);
    ClustersM560_560 = table(x560', y560', z560', A560', A560rel');
    ClustersM560_560.Properties.VariableNames = {'x560', 'y560', 'z560', 'A560', 'A560rel'};
    save('ClustersM560_560')
    % isolate cell from cluster
    CellCluster = table2array(ClustersM560_560(:,[1:3]));
    CellStr = [a{4} ' ' a{5}(1:4)];
    % CellClusterSurfShapes_RH(CellCluster, ClustersM560_560, 'title', CellStr)
    CellClusterSurfShapes2_RH(CellCluster, 1)
end


%% NoMOI

% 1) Find CellSHP: 
open('G:\LLSM_VSV_VHH\SUM52\CS1_SUM52_NoMOI\Ex01_488_150mW_70p560n_500mW_100p_642_150mW_100p_z0p4\CellSHP_fig.fig')
% 2) Modify figure
plot(CellSHP,'FaceColor','k','FaceAlpha',0.01)
title('SUM 159 cell (AP2-RFP)')
set(gcf,'Color','w')
xlim([0 400]) % 300pxl = 31um
ylim([0 600]) % 450pxl = 47um
zlim([0 50]) % 35pxl = 14um (
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'zticklabel',[])

% quantify AP2s from all cells
open('G:\LLSM_VSV_VHH\SUM52\CS1_SUM52_NoMOI\M560.mat')

idx = find(M560.ExpNum(1,:)==10);
figure('Position',[854 641 645 335])
histogram(M560.A_tot(1,idx),1:5:650,'FaceColor',[0.67 0.13 0.13],'EdgeColor',[0.75 0.75 0.75])
title('Distribution of AP2-RFP clusters')
xlabel('RFP fluorescence (a.u.)')
ylabel('Counts (#)')
set(gcf,'Color','w')

%% No VHH: IMAGE REPRESENTATION

clear

cd('G:\LLSM_VSV_VHH\SUM51 for hidde summary\SUM51 for hidde summary\CS4_NoVHH')
load('G:\LLSM_VSV_VHH\SUM51 for hidde summary\SUM51 for hidde summary\CS4_NoVHH\Ex15_488_150mW_70p560n_500mW_100p_642_150mW_100p_z0p4_488FOCUS\CellSHP.mat');

% load the M488 data
load('G:\LLSM_VSV_VHH\SUM51 for hidde summary\SUM51 for hidde summary\CS4_NoVHH\M488.mat')

NVidx = find(M488.ExpNum==10); % from 1...
in = find(inShape(CellSHP,M488.x_tot(1,NVidx),M488.y_tot(1,NVidx),M488.z_tot(1,NVidx))); % in CellSHP

idxOUT = setdiff(NVidx,in);
% distances of points NOT inside the CellSHP
for f = 1:numel(idxOUT)
    [I(f),D(f)] = nearestNeighbor(CellSHP,M488.x_tot(1,idxOUT(f)),M488.y_tot(1,idxOUT(f)),M488.z_tot(1,idxOUT(f)));
end
% figure, histogram(D,1:2:100) - check distribution of distances    

DCO = 15; % distance cutoff

IdxVir = union(in,idxOUT(D<=DCO)); % from 3742 to 1639

% Modify on premade figure
open('G:\LLSM_VSV_VHH\SUM51 for hidde summary\SUM51 for hidde summary\CS4_NoVHH\Ex15_488_150mW_70p560n_500mW_100p_642_150mW_100p_z0p4_488FOCUS\CellSHP_fig.fig');
plot(CellSHP,'FaceColor','k','FaceAlpha',0.01)
title('SUM 159 cell (AP2-RFP) with VSV-G-EboV')
set(gcf,'Color','w')
set(gca,'Zdir','reverse')
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'zticklabel',[])
xlim([0 550]) % 300pxl = 31um
ylim([0 800]) % 450pxl = 47um
zlim([0 120]) % 35pxl = 14um
view(61.0383,25.3505)
hold on
scatter3(M488.x_tot(1,IdxVir),M488.y_tot(1,IdxVir),M488.z_tot(1,IdxVir),0.001*M488.A_tot(1,IdxVir),[0 0.6 0.3],'filled')
hold off

% percentage virions inside...?
% average sizes of cluster in and out (?)

for f = 1:numel(in)
    [inI(f),inD(f)] = nearestNeighbor(CellSHP,M488.x_tot(1,in(f)),M488.y_tot(1,in(f)),M488.z_tot(1,in(f)));
end

inDCO = 5;
figure, scatter3(M488.x_tot(1,in(inD>=inDCO)),M488.y_tot(1,in(inD>=inDCO)),M488.z_tot(1,in(inD>=inDCO)),0.001*M488.A_tot(1,in(inD>=inDCO)),[0 0.6 0.3],'filled')
xlabel('x'), ylabel('y')

% Show cross-section 
xlim([300 330])
ylim([500 600])
% look at sizes of Viral Clusters

open('G:\LLSM_VSV_VHH\SUM51 for hidde summary\SUM51 for hidde summary\CS4_NoVHH\Ex15_488_150mW_70p560n_500mW_100p_642_150mW_100p_z0p4_488FOCUS\CellSHP_fig.fig');
plot(CellSHP,'FaceColor','k','FaceAlpha',0.01)
title('SUM 159 cell (AP2-RFP) with VSV-G-EboV')
set(gcf,'Color','w')
set(gca,'Zdir','reverse')
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'zticklabel',[])
xlim([0 550]) % 300pxl = 31um
ylim([0 800]) % 450pxl = 47um
zlim([0 120]) % 35pxl = 14um
view(61.0383,25.3505)
hold on
scatter3(M488.x_tot(1,in(inD>=inDCO)),M488.y_tot(1,in(inD>=inDCO)),M488.z_tot(1,in(inD>=inDCO)),0.001*M488.A_tot(1,in(inD>=inDCO)),[0 0.5 0],'filled')
scatter3(M488.x_tot(1,idxOUT(D<=DCO)),M488.y_tot(1,idxOUT(D<=DCO)),M488.z_tot(1,idxOUT(D<=DCO)),0.001*M488.A_tot(1,idxOUT(D<=DCO)),[0 0.6 0.3],'filled')
hold off
for t = 1:20
    xlim([(t-1)*30 t*30])
    ylabel(['x = [' num2str((t-1)*30) ':' num2str(t*30) ']'])
    pause
end
xlim([300 330])
ylim([500 600])

%% NoVHHs: INDICES OF OUTSIDE (-15,limit]) and INSIDE ([+5,end])

clear

Rtdir = 'G:\LLSM_VSV_VHH\SUM51 for hidde summary\SUM51 for hidde summary\CS4_NoVHH';
load('G:\LLSM_VSV_VHH\SUM51 for hidde summary\SUM51 for hidde summary\CS4_NoVHH\M488.mat')

Flds = dir(Rtdir);
Flds = Flds(4:7);

% read CellSHP
D_tot = [];
IdxINco = [];
IdxOUTco = [];
for k = 1:4
    load(fullfile(Flds(k).folder,Flds(k).name,'CellSHP.mat'))
    NVidx = find(M488.ExpNum==(k+9));
    for f = 1:numel(NVidx)
        [I(f),D(f)] = nearestNeighbor(CellSHP,M488.x_tot(1,NVidx(f)),M488.y_tot(1,NVidx(f)),M488.z_tot(1,NVidx(f)));
    end
    IdxIN = find(inShape(CellSHP,M488.x_tot(1,NVidx),M488.y_tot(1,NVidx),M488.z_tot(1,NVidx))); % in CellSHP
    IdxOUT = setdiff(NVidx,IdxIN);
    % cutoffs
    IdxINco_t = NVidx(intersect(find(D>5),IdxIN)); % 5pxl or more in
    IdxOUTco_t = NVidx(intersect(find(D<15),IdxOUT)); % 15pxl or less out
    % collect
    D_tot = [D_tot D];
    IdxINco = [IdxINco IdxINco_t];
    IdxOUTco = [IdxOUTco IdxOUTco_t];
    % clean
    D = [];
end

%% Size Distribution of Virions, IN/OUT

figure, histogram(M488.A_tot(1,IdxINco),1:10:200)

% Fraction internalized
FRclu = length(IdxINco)/(length(IdxINco)+length(IdxOUTco)); % 55%
FRvir = sum(M488.A_tot(1,IdxINco))/(sum(M488.A_tot(1,IdxINco))+sum(M488.A_tot(1,IdxOUTco))); % 64%

% Int-Distance VS. Size
figure, scatter(D_tot(IdxINco),M488.A_tot(1,IdxINco)), xlim([0 max(D_tot(IdxINco))])  % no pattern

% Int-Distance VS. AP2-association
figure, scatter(D_tot(IdxINco),M488.A_tot(2,IdxINco)), xlim([0 max(D_tot(IdxINco))]) % slight decrease from 10/15 to 5/0 with D
figure, scatter(D_tot(IdxOUTco),M488.A_tot(2,IdxOUTco)) % No difference over D, but a pool of larger
% AP2-clusters (~40/50a.u.) present on surface
% <20 => not associated, % checked distribution on NoMOI (25-30 = singles?)

% AP2 association (>20a.u. in RFP-channel)
AP2assIN = sum(M488.A_tot(2,IdxINco)>20)/length(IdxINco)*100; % 0.2%
AP2assOUT = sum(M488.A_tot(2,IdxOUTco)>20)/length(IdxOUTco)*100; % 2.7%

%% NoVHH: Compare Precipitation, Total GFP and VHH-distribution

IdxPrec = find(D_tot(IdxOUT)>5);

FrVirPrec = sum(M488.A_tot(1,IdxPrec))/sum(M488.A_tot(1,:)); % 0% (!) 

figure, histogram(M488.A_tot(3,IdxPrec),1:10:200)
hold on
histogram(M488.A_tot(3,:),1:10:200, 'FaceAlpha',0.2)
legend({'Precipitated','Not Prec.'})

%% FIGURE STRUCTURE: NoVHH load

load('G:\LLSM_VSV_VHH\FIGURE FILES\ExpStr.mat','ExpStr')

Exp = 'NoVHH';
Nums = unique(M488.ExpNum(:));

n = 0;
for h = 1:numel(Nums)
    % Numbers
    ExpStr.(Exp)(h).Counts = sum(M488.ExpNum(IdxINco)==Nums(h));
    % AP2assIN
    ExpStr.(Exp)(h).AP2assIN = sum(M488.A_tot(2,IdxINco)>20 & M488.ExpNum(IdxINco)==Nums(h))/length(IdxINco)*100; % 0.2%
    ExpStr.(Exp)(h).AP2assOUT = sum(M488.A_tot(2,IdxOUTco)>20 & M488.ExpNum(IdxOUTco)==Nums(h))/length(IdxOUTco)*100; % 2.7%
    % Prec.
    ExpStr.(Exp)(h).FrVirPrec = sum(M488.A_tot(1,IdxPrec) & M488.ExpNum(IdxPrec)==Nums(h))/sum(M488.A_tot(1,:)); % 0% (!) 
    % CellAss_VirInternalized
    ExpStr.(Exp)(h).CAss_VInt = sum(M488.A_tot(1,IdxINco) & M488.ExpNum(IdxINco)==Nums(h))/(sum(M488.A_tot(1,IdxINco) & M488.ExpNum(IdxINco)==Nums(h))+sum(M488.A_tot(1,IdxOUTco) & M488.ExpNum(IdxOUTco)==Nums(h)));
end

save('G:\LLSM_VSV_VHH\FIGURE FILES\ExpStr.mat','ExpStr')

%% G10 data: INDICES INSIDE/OUTSIDE

clear

% for both G10 and G84, look in 'TEMPORARYDATA' folder, please...

Rtdir = 'G:\LLSM_VSV_VHH\TEMPORARYDATA\CS3_G10_Hidde1';
load('G:\LLSM_VSV_VHH\TEMPORARYDATA\CS3_G10_Hidde1\M488.mat')

Flds = dir(Rtdir);
Flds = Flds(3:7);

% read CellSHP
D_tot = [];
IdxINco = [];
IdxOUTco = [];
for k = 1:5
    load(fullfile(Flds(k).folder,Flds(k).name,'CellSHP.mat'))
    NVidx = find(M488.ExpNum==(k+9));
    for f = 1:numel(NVidx)
        [I(f),D(f)] = nearestNeighbor(CellSHP,M488.x_tot(1,NVidx(f)),M488.y_tot(1,NVidx(f)),M488.z_tot(1,NVidx(f)));
    end
    IdxIN = find(inShape(CellSHP,M488.x_tot(1,NVidx),M488.y_tot(1,NVidx),M488.z_tot(1,NVidx))); % in CellSHP
    IdxOUT = setdiff(NVidx,IdxIN);
    % cutoffs
    IdxINco_t = NVidx(intersect(find(D>5),IdxIN)); % 5pxl or more in
    IdxOUTco_t = NVidx(intersect(find(D<15),IdxOUT)); % 15pxl or less out
    % collect
    D_tot = [D_tot D];
    IdxINco = [IdxINco IdxINco_t];
    IdxOUTco = [IdxOUTco IdxOUTco_t];
    % clean
    D = [];
end

%% G10: Size Distribution of Virions, IN/OUT

figure, histogram(M488.A_tot(1,IdxINco),1:10:200)
hold on
histogram(M488.A_tot(1,IdxOUTco),1:10:200)
legend({'Inside','Outside'})

% Fraction internalized
FRclu = length(IdxINco)/(length(IdxINco)+length(IdxOUTco)); % 81%
FRvir = sum(M488.A_tot(1,IdxINco))/(sum(M488.A_tot(1,IdxINco))+sum(M488.A_tot(1,IdxOUTco))); % 46% (compared to 64% NoVHH)
% Make cellular average and error-bars

% Int-Distance VS. Size
figure, scatter(D_tot(IdxINco),M488.A_tot(1,IdxINco)), xlim([0 max(D_tot(IdxINco))])  % no pattern

% Int-Distance VS. AP2-association
figure, scatter(D_tot(IdxINco),M488.A_tot(2,IdxINco)), xlim([0 max(D_tot(IdxINco))]) % HUGE AP2-clusters (>200a.u.)
figure, scatter(D_tot(IdxOUTco),M488.A_tot(2,IdxOUTco)) % Big pool of AP2-cluster ~100a.u.
% AP2-clusters (~40/50a.u.) present on surface
% <20 => not associated, % checked distribution on NoMOI (25-30 = singles?)

% ARE THEY ASSOCIATED WITH VHHs? And how much on IN vs OUT?
VHHin = M488.A_tot(3,IdxINco);
VHHout = M488.A_tot(3,IdxOUTco);

figure, histogram(M488.A_tot(3,IdxINco),1:10:200)
hold on
histogram(M488.A_tot(3,IdxOUTco),1:10:200)
legend({'Inside','Outside'})

% AP2 association (>20a.u. in RFP-channel)
AP2assIN = sum(M488.A_tot(2,IdxINco)>20)/length(IdxINco)*100; % 44%
AP2assOUT = sum(M488.A_tot(2,IdxOUTco)>20)/length(IdxOUTco)*100; % 24%
% SUCH WEIRD OBSERVATION, BUT QUANTIFIABLE, HISTOGRAMS BETWEEN NoVHH, G10 and G84
% Use Green, Cyan (G10), Dark Blue (G84)

% AP2-ASSOCIATION vs. VHH-ASSOCIATION
IdxAP2ass = find(M488.A_tot(2,:)>20);
IdxAP2noass = find(M488.A_tot(2,:)<=20);

figure, histogram(M488.A_tot(3,IdxAP2ass),1:10:200)
hold on
histogram(M488.A_tot(3,IdxAP2noass),1:10:200)
legend({'AP2-assoc.','Not AP2-assoc.'})
% Normalize(!!)

%% G10: Compare Precipitation, Total GFP and VHH-distribution

IdxPrec = find(D_tot(IdxOUT)>5);
IdxNoPrec = setdiff(1:numel(D_tot),IdxPrec);

FrVirPrec = sum(M488.A_tot(1,IdxPrec))/sum(M488.A_tot(1,:)); % 1 w/15co, 3,8% w/10co, 10% w/

figure, histogram(M488.A_tot(3,IdxPrec),1:10:200)
hold on
histogram(M488.A_tot(3,:),1:10:200, 'FaceAlpha',0.2)
legend({'Precipitated','Not Prec.'})

%% FIGURE STRUCTURE: G10 load

load('G:\LLSM_VSV_VHH\FIGURE FILES\ExpStr.mat','ExpStr')

Exp = 'G10';
Nums = unique(M488.ExpNum(:));

n = 0;
for h = 1:numel(Nums)
    % Numbers
    ExpStr.(Exp)(h).Counts = sum(M488.ExpNum(IdxINco)==Nums(h));
    % AP2assIN
    ExpStr.(Exp)(h).AP2assIN = sum(M488.A_tot(2,IdxINco)>20 & M488.ExpNum(IdxINco)==Nums(h))/length(IdxINco)*100; % 0.2%
    ExpStr.(Exp)(h).AP2assOUT = sum(M488.A_tot(2,IdxOUTco)>20 & M488.ExpNum(IdxOUTco)==Nums(h))/length(IdxOUTco)*100; % 2.7%
    % Prec.
    ExpStr.(Exp)(h).FrVirPrec = sum(M488.A_tot(1,IdxPrec) & M488.ExpNum(IdxPrec)==Nums(h))/sum(M488.A_tot(1,:)); % 0% (!) 
    % CellAss_VirInternalized
    ExpStr.(Exp)(h).CAss_VInt = sum(M488.A_tot(1,IdxINco) & M488.ExpNum(IdxINco)==Nums(h))/(sum(M488.A_tot(1,IdxINco) & M488.ExpNum(IdxINco)==Nums(h))+sum(M488.A_tot(1,IdxOUTco) & M488.ExpNum(IdxOUTco)==Nums(h)));
end

save('G:\LLSM_VSV_VHH\FIGURE FILES\ExpStr.mat','ExpStr')

%% G84 data: INDICES INSIDE/OUTSIDE

clear

% for both G10 and G84, look in 'TEMPORARYDATA' folder, please...

Rtdir = 'G:\LLSM_VSV_VHH\TEMPORARYDATA\CS_G84';
load('G:\LLSM_VSV_VHH\TEMPORARYDATA\CS3_G10_Hidde1\M488.mat')

Flds = dir(Rtdir);
Flds = Flds(3:7);

% read CellSHP
D_tot = [];
IdxINco = [];
IdxOUTco = [];
for k = 1:5
    load(fullfile(Flds(k).folder,Flds(k).name,'CellSHP.mat'))
    NVidx = find(M488.ExpNum==(k+9));
    for f = 1:numel(NVidx)
        [I(f),D(f)] = nearestNeighbor(CellSHP,M488.x_tot(1,NVidx(f)),M488.y_tot(1,NVidx(f)),M488.z_tot(1,NVidx(f)));
    end
    IdxIN = find(inShape(CellSHP,M488.x_tot(1,NVidx),M488.y_tot(1,NVidx),M488.z_tot(1,NVidx))); % in CellSHP
    IdxOUT = setdiff(NVidx,IdxIN);
    % cutoffs
    IdxINco_t = NVidx(intersect(find(D>5),IdxIN)); % 5pxl or more in
    IdxOUTco_t = NVidx(intersect(find(D<15),IdxOUT)); % 15pxl or less out
    % collect
    D_tot = [D_tot D];
    IdxINco = [IdxINco IdxINco_t];
    IdxOUTco = [IdxOUTco IdxOUTco_t];
    % clean
    D = [];
end

% OBS: 5236 Viral Points, ~1800 around cells)

%% G84: Size Distribution of Virions, IN/OUT

figure, histogram(M488.A_tot(1,IdxINco),1:10:200)
hold on
histogram(M488.A_tot(1,IdxOUTco),1:10:200)
legend({'Inside','Outside'})

% Fraction internalized
FRclu = length(IdxINco)/(length(IdxINco)+length(IdxOUTco)); % 76%
FRvir = sum(M488.A_tot(1,IdxINco))/(sum(M488.A_tot(1,IdxINco))+sum(M488.A_tot(1,IdxOUTco))); % 70% (compared to 64% NoVHH)
% Make cellular average and error-bars

% Int-Distance VS. Size
figure, scatter(D_tot(IdxINco),M488.A_tot(1,IdxINco)), xlim([0 max(D_tot(IdxINco))])  % no pattern

% Int-Distance VS. AP2-association
figure, scatter(D_tot(IdxINco),M488.A_tot(2,IdxINco)), xlim([0 max(D_tot(IdxINco))]) % HUGE AP2-clusters (>200a.u.)
figure, scatter(D_tot(IdxOUTco),M488.A_tot(2,IdxOUTco)) % Ok pool of AP2-cluster ~100a.u.
% AP2-clusters (~40/50a.u.) present on surface
% <20 => not associated, % checked distribution on NoMOI (25-30 = singles?)

% ARE THEY ASSOCIATED WITH VHHs? And how much on IN vs OUT?
VHHin = M488.A_tot(3,IdxINco);
VHHout = M488.A_tot(3,IdxOUTco);

figure, histogram(M488.A_tot(3,IdxINco),1:10:200)
hold on
histogram(M488.A_tot(3,IdxOUTco),1:10:200)
legend({'Inside','Outside'})

% AP2 association (>20a.u. in RFP-channel)
AP2assIN = sum(M488.A_tot(2,IdxINco)>20)/length(IdxINco)*100; % 34%
AP2assOUT = sum(M488.A_tot(2,IdxOUTco)>20)/length(IdxOUTco)*100; % 12%
% SUCH WEIRD OBSERVATION, BUT QUANTIFIABLE, HISTOGRAMS BETWEEN NoVHH, G10 and G84
% Use Green, Cyan (G10), Dark Blue (G84)

% AP2-ASSOCIATION vs. VHH-ASSOCIATION
IdxAP2ass = find(M488.A_tot(2,:)>20);
IdxAP2noass = find(M488.A_tot(2,:)<=20);

figure, histogram(M488.A_tot(3,IdxAP2ass),1:10:200)
hold on
histogram(M488.A_tot(3,IdxAP2noass),1:10:200)
legend({'AP2-assoc.','Not AP2-assoc.'})
% Normalize(!!)

%% G84: Compare Precipitation, Total GFP and VHH-distribution

IdxPrec = find(D_tot(IdxOUT)>10);
IdxNoPrec = setdiff(1:numel(D_tot),IdxPrec);

FrVirPrec = sum(M488.A_tot(1,IdxPrec))/sum(M488.A_tot(1,:)); % 10% w/5co, 8% x/10co, 5% w/15co

figure, histogram(M488.A_tot(3,IdxPrec),1:10:200)
hold on
histogram(M488.A_tot(3,:),1:10:200, 'FaceAlpha',0.2)
legend({'Precipitated','Not Prec.'})

%% FIGURE STRUCTURE: G10 load

load('G:\LLSM_VSV_VHH\FIGURE FILES\ExpStr.mat','ExpStr')

Exp = 'G84';
Nums = unique(M488.ExpNum(:));

n = 0;
for h = 1:numel(Nums)
    % Numbers
    ExpStr.(Exp)(h).Counts = sum(M488.ExpNum(IdxINco)==Nums(h));
    % AP2assIN
    ExpStr.(Exp)(h).AP2assIN = sum(M488.A_tot(2,IdxINco)>20 & M488.ExpNum(IdxINco)==Nums(h))/length(IdxINco)*100; % 0.2%
    ExpStr.(Exp)(h).AP2assOUT = sum(M488.A_tot(2,IdxOUTco)>20 & M488.ExpNum(IdxOUTco)==Nums(h))/length(IdxOUTco)*100; % 2.7%
    % Prec.
    ExpStr.(Exp)(h).FrVirPrec = sum(M488.A_tot(1,IdxPrec) & M488.ExpNum(IdxPrec)==Nums(h))/sum(M488.A_tot(1,:)); % 0% (!) 
    % CellAss_VirInternalized
    ExpStr.(Exp)(h).CAss_VInt = sum(M488.A_tot(1,IdxINco) & M488.ExpNum(IdxINco)==Nums(h))/(sum(M488.A_tot(1,IdxINco) & M488.ExpNum(IdxINco)==Nums(h))+sum(M488.A_tot(1,IdxOUTco) & M488.ExpNum(IdxOUTco)==Nums(h)));
end

save('G:\LLSM_VSV_VHH\FIGURE FILES\ExpStr.mat','ExpStr')

%% COLLECTED FIGURES

% SD(:,:,1) = [0.2,99.8;2.7,97.3];
% SD(:,:,2) = [44,56;24,76];
% SD(:,:,3) = [34,66;12,88];

% Create a stacked, grouped bar chart using the function
% figure,
% groupLabels = {'No VHH' 'G10' 'G84'}; 
% plotBarStackGroups(SD, groupLabels);
% % Chance the colors of each bar segment
% colors = jet(size(h,2)); %or define your own color order; 1 for each m segments
% colors = repelem(colors,size(h,1),1); 
% colors = mat2cell(colors,ones(size(colors,1),1),3);
% set(h,{'FaceColor'},colors)

% AP2-assoc. Virus
figure
a = [0.2,99.8;2.7,97.3;0.001,0.001;44,56;24,76;0.0001,0.0001;34,66;12,88];
H = bar(a, 'stacked');
H(1).FaceColor = [0.67 0.13 0.13];
H(2).FaceColor = [0 0.6 0.3];
set(gca,'Xticklabels',{'IN','OUT','','IN','OUT','','IN','OUT'})
title('AP2-association of cell-bound virus')
set(gcf,'Color','w')

% Precipitation
figure
Prec = [0,100;10,90;8,92];
H = bar(Prec, 'stacked');
H(1).FaceColor = [0.75 0.75 0.75];
H(2).FaceColor = [0 0.6 0.3];
set(gca,'Xticklabels',{'No VHH','G10','G84'})
title('Off-cell precipitation')
set(gcf,'Color','w')

%{
Cell-Associated internalization of Virus
figure
Prec = [0,100;10,90;10,90];
H = bar(Prec, 'stacked');
H(1).FaceColor = [0.75 0.75 0.75];
H(2).FaceColor = [0 0.6 0.3];
set(gca,'Xticklabels',{'No VHH','G10','G84'})
title('Off-cell precipitation')
set(gcf,'Color','w')
%}

%% CFDs

% G10

% Internalization and VHHs
CFDs.VHHin_G10 = VHHin;
CFDs.VHHout_G10 = VHHout;

CFDs.AP2ass_G10 = M488.A_tot(3,IdxAP2ass);
CFDs.AP2noass_G10 = M488.A_tot(3,IdxAP2noass);

CFDs.Prec_G10 = M488.A_tot(3,IdxPrec);
CFDs.NoPrec_G10 = M488.A_tot(3,IdxNoPrec);

save('G:\LLSM_VSV_VHH\FIGURE FILES\CFDs.mat','CFDs')

% G84
load('G:\LLSM_VSV_VHH\FIGURE FILES\CFDs.mat','CFDs')

CFDs.VHHin_G84 = VHHin;
CFDs.VHHout_G84 = VHHout;

CFDs.AP2ass_G84 = M488.A_tot(3,IdxAP2ass);
CFDs.AP2noass_G84 = M488.A_tot(3,IdxAP2noass);

CFDs.Prec_G84 = M488.A_tot(3,IdxPrec);
CFDs.NoPrec_G84 = M488.A_tot(3,IdxNoPrec);

save('G:\LLSM_VSV_VHH\FIGURE FILES\CFDs.mat','CFDs')

% PLOT

% G10: VHHs internalization
IN10 = CFDs.VHHin_G10;
IN10(IN10<0) = 0;
OUT10 = CFDs.VHHout_G10;
OUT10(OUT10<0) = 0;

[fi,xi] = ecdf(IN10);
[fo,xo] = ecdf(OUT10);

figure,
plot(xi,fi,'Color',[0 0.6 0.3])
hold on
plot(xo,fo,'Color',[0.75 0.75 0.75])

xlim([0 1000])
legend({'Internalized','Surfaced'})
title('G10-VHH effect on internalization')
set(gcf,'Color','w')
xlabel('Amount of bound VHHs (a.u.)')

% G84: VHHs internalization
IN10 = CFDs.VHHin_G84;
IN10(IN10<0) = 0;
OUT10 = CFDs.VHHout_G84;
OUT10(OUT10<0) = 0;

[fi,xi] = ecdf(IN10);
[fo,xo] = ecdf(OUT10);

figure,
plot(xi,fi,'Color',[0 0.6 0.3])
hold on
plot(xo,fo,'Color',[0.75 0.75 0.75])

xlim([0 1000])
legend({'Internalized','Surfaced'})
title('G84-VHH effect on internalization')
set(gcf,'Color','w')
xlabel('Amount of bound VHHs (a.u.)')

% G10: AP2-association
AP2ass = CFDs.AP2ass_G10;
AP2noass = CFDs.AP2noass_G10;

[fi,xi] = ecdf(AP2ass);
[fo,xo] = ecdf(AP2noass);

figure,
plot(xi,fi,'Color',[0 0.6 0.3])
hold on
plot(xo,fo,'Color',[0.75 0.75 0.75])

xlim([0 1000])
legend({'AP2-associated','Non-associated'},'Location','southeast')
title('G10-VHH effect on AP2-association')
set(gcf,'Color','w')
xlabel('Amount of bound VHHs (a.u.)')

% G84: AP2-association
AP2ass = CFDs.AP2ass_G84;
AP2noass = CFDs.AP2noass_G84;

[fi,xi] = ecdf(AP2ass);
[fo,xo] = ecdf(AP2noass);

figure,
plot(xi,fi,'Color',[0 0.6 0.3])
hold on
plot(xo,fo,'Color',[0.75 0.75 0.75])

xlim([0 1000])
legend({'AP2-associated','Non-associated'},'Location','southeast')
title('G84-VHH effect on AP2-association')
set(gcf,'Color','w')
xlabel('Amount of bound VHHs (a.u.)')

% G10: Precipitation
Prec = CFDs.Prec_G10;
NoPrec = CFDs.NoPrec_G10;

[fi,xi] = ecdf(Prec);
[fo,xo] = ecdf(NoPrec);

figure,
plot(xi,fi,'Color',[0 0.6 0.3])
hold on
plot(xo,fo,'Color',[0.75 0.75 0.75])

xlim([0 1000])
legend({'Precipitated','Cell Bound'},'Location','southeast')
title('G10-VHH effect on Precipation')
set(gcf,'Color','w')
xlabel('Amount of bound VHHs (a.u.)')

% G10: Precipitation
Prec = CFDs.Prec_G84;
NoPrec = CFDs.NoPrec_G84;

[fi,xi] = ecdf(Prec);
[fo,xo] = ecdf(NoPrec);

figure,
plot(xi,fi,'Color',[0 0.6 0.3])
hold on
plot(xo,fo,'Color',[0.75 0.75 0.75])

xlim([0 1000])
legend({'Precipitated','Cell Bound'},'Location','southeast')
title('G84-VHH effect on Precipation')
set(gcf,'Color','w')
xlabel('Amount of bound VHHs (a.u.)')

%% EXAMPLES

%% G10 CELL EXAMPLE

clear

rtdir = 'G:\LLSM_VSV_VHH\TEMPORARYDATA\CS3_G10_Hidde1';
load(fullfile(rtdir,'data488master.mat'))
load(fullfile(rtdir,'M488.mat'))
load(fullfile(rtdir,'M560.mat'))

Flds = dir(rtdir);
Flds = Flds(3:7);

% Defining Viral distance cutoff again:
D_tot = [];
IdxINco = [];
IdxOUTco = [];
for k = 1:5
    load(fullfile(Flds(k).folder,Flds(k).name,'CellSHP.mat'))
    NVidx = find(M488.ExpNum==(k+9));
    for f = 1:numel(NVidx)
        [I(f),D(f)] = nearestNeighbor(CellSHP,M488.x_tot(1,NVidx(f)),M488.y_tot(1,NVidx(f)),M488.z_tot(1,NVidx(f)));
    end
    IdxIN = find(inShape(CellSHP,M488.x_tot(1,NVidx),M488.y_tot(1,NVidx),M488.z_tot(1,NVidx))); % in CellSHP
    IdxOUT = setdiff(NVidx,IdxIN);
    % cutoffs
    IdxINco_t = NVidx(intersect(find(D>5),IdxIN)); % 5pxl or more in
    IdxOUTco_t = NVidx(intersect(find(D<5),IdxOUT)); % 15pxl or less out
    % collect
    D_tot = [D_tot D];
    IdxINco = [IdxINco IdxINco_t];
    IdxOUTco = [IdxOUTco IdxOUTco_t];
    % clean
    D = [];
end

% for each experiment:
for n = 1:numel(Flds)
    load(fullfile(Flds(n).folder,Flds(n).name,'CellSHP.mat'))
    plot(CellSHP,'FaceColor','k','FaceAlpha',0.01)
    title('SUM 159 cell (AP2-RFP) with VSV-G-EboV and VHH-G10')
    set(gcf,'Color','w')
    set(gca,'Zdir','reverse')
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    set(gca,'zticklabel',[])
    xlim([0 550]) % 300pxl = 31um
    ylim([0 800]) % 450pxl = 47um
    zlim([0 120]) % 35pxl = 14um
    view(61.0383,25.3505)
    hold on
    % plot used AP2-clusters initially
    scatter3(CellCluster(:,1),CellCluster(:,2),CellCluster(:,3),7,[0.5 0 0],'filled')
    % darker green inside cell
    scatter3(M488.x_tot(1,intersect(IdxINco,find(M488.ExpNum==(n+9)))),M488.y_tot(1,intersect(IdxINco,find(M488.ExpNum==(n+9)))),M488.z_tot(1,intersect(IdxINco,find(M488.ExpNum==(n+9)))),0.005*M488.A_tot(1,intersect(IdxINco,find(M488.ExpNum==(n+9)))),[0 0.5 0],'filled')
    % light green on top of cell
    scatter3(M488.x_tot(1,intersect(IdxOUTco,find(M488.ExpNum==(n+9)))),M488.y_tot(1,intersect(IdxOUTco,find(M488.ExpNum==(n+9)))),M488.z_tot(1,intersect(IdxOUTco,find(M488.ExpNum==(n+9)))),0.005*M488.A_tot(1,intersect(IdxOUTco,find(M488.ExpNum==(n+9)))),[0 0.6 0.3],'filled')
    % Associated Viral particles (in blue)
    scatter3(M488.x_tot(1,intersect([IdxINco IdxOUTco],find(M488.ExpNum==(n+9)))),M488.y_tot(1,intersect([IdxINco IdxOUTco],find(M488.ExpNum==(n+9)))),M488.z_tot(1,intersect([IdxINco IdxOUTco],find(M488.ExpNum==(n+9)))),0.001*M488.A_tot(1,intersect([IdxINco IdxOUTco],find(M488.ExpNum==(n+9)))),[0.3010, 0.7450, 0.9330],'filled')
    hold off
    pause
end

% CHOOSE N = 1

% REMOVE BIG GREEN CLUSTERS AND ENLARGE THE SMALLER(?)
% figure, histogram(M488.A_tot(1,intersect(IdxINco,find(M488.ExpNum==(n+9)))))
% figure, histogram(M488.A_tot(1,intersect(IdxOUTco,find(M488.ExpNum==(n+9)))),1:100:15000)
% [0 5000] or [1 2000]

IdxEx = find(M488.ExpNum==(n+9));
IdxSz = find(M488.A_tot(1,:)<2000);

IdxFin = intersect(IdxEx,IdxSz);

n = 1;

load(fullfile(Flds(n).folder,Flds(n).name,'CellSHP.mat'))
load(fullfile(Flds(n).folder,Flds(n).name,'CellCluster.mat'))
plot(CellSHP,'FaceColor','k','FaceAlpha',0.01)
title('SUM 159 cell (AP2-RFP) with VSV-G-EboV and VHH-G10')
set(gcf,'Color','w')
set(gca,'Zdir','reverse')
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'zticklabel',[])
hold on
% plot used AP2-clusters initially
scatter3(CellCluster(:,1),CellCluster(:,2),CellCluster(:,3),7,[0.5 0 0],'filled')
% darker green inside cell
scatter3(M488.x_tot(1,intersect(IdxINco,IdxFin)),M488.y_tot(1,intersect(IdxINco,IdxFin)),M488.z_tot(1,intersect(IdxINco,IdxFin)),0.01*M488.A_tot(1,intersect(IdxINco,IdxFin)),[0 0.5 0],'filled')
% light green on top of cell
scatter3(M488.x_tot(1,intersect(IdxOUTco,IdxFin)),M488.y_tot(1,intersect(IdxOUTco,IdxFin)),M488.z_tot(1,intersect(IdxOUTco,IdxFin)),0.01*M488.A_tot(1,intersect(IdxOUTco,IdxFin)),[0 0.6 0.3],'filled')
% Associated Viral particles (in blue)
scatter3(M488.x_tot(1,intersect([IdxINco IdxOUTco],IdxFin)),M488.y_tot(1,intersect([IdxINco IdxOUTco],IdxFin)),M488.z_tot(1,intersect([IdxINco IdxOUTco],IdxFin)),0.003*M488.A_tot(1,intersect([IdxINco IdxOUTco],IdxFin)),[0.3010, 0.7450, 0.9330],'filled')
hold off

xlim([0 550]) % 300pxl = 31um
ylim([0 800]) % 450pxl = 47um
zlim([0 120]) % 35pxl = 14um
view(61.0383,25.3505)

%% G84 CELL EXAMPLE

clear

rtdir = 'G:\LLSM_VSV_VHH\TEMPORARYDATA\CS_G84';
load(fullfile(rtdir,'data488master.mat'))
load(fullfile(rtdir,'M488.mat'))
load(fullfile(rtdir,'M560.mat'))

Flds = dir(rtdir);
Flds = Flds(3:7);

% Defining Viral distance cutoff again:
D_tot = [];
IdxINco = [];
IdxOUTco = [];
for k = 1:5
    load(fullfile(Flds(k).folder,Flds(k).name,'CellSHP.mat'))
    NVidx = find(M488.ExpNum==(k+4));
    for f = 1:numel(NVidx)
        [I(f),D(f)] = nearestNeighbor(CellSHP,M488.x_tot(1,NVidx(f)),M488.y_tot(1,NVidx(f)),M488.z_tot(1,NVidx(f)));
    end
    IdxIN = find(inShape(CellSHP,M488.x_tot(1,NVidx),M488.y_tot(1,NVidx),M488.z_tot(1,NVidx))); % in CellSHP
    IdxOUT = setdiff(NVidx,IdxIN);
    % cutoffs
    IdxINco_t = NVidx(intersect(find(D>5),IdxIN)); % 5pxl or more in
    IdxOUTco_t = NVidx(intersect(find(D<5),IdxOUT)); % 15pxl or less out
    % collect
    D_tot = [D_tot D];
    IdxINco = [IdxINco IdxINco_t];
    IdxOUTco = [IdxOUTco IdxOUTco_t];
    % clean
    D = [];
end


% for each experiment:
for n = 1:numel(Flds)
    load(fullfile(Flds(n).folder,Flds(n).name,'CellSHP.mat'))
    load(fullfile(Flds(n).folder,Flds(n).name,'CellCluster.mat'))
    plot(CellSHP,'FaceColor','k','FaceAlpha',0.01)
    title('SUM 159 cell (AP2-RFP) with VSV-G-EboV and VHH-G10')
    set(gcf,'Color','w')
    set(gca,'Zdir','reverse')
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    set(gca,'zticklabel',[])
    xlim([0 550]) % 300pxl = 31um
    ylim([0 800]) % 450pxl = 47um
    zlim([0 120]) % 35pxl = 14um
    view(61.0383,25.3505)
    % 
    IdxEx = find(M488.ExpNum==(n+4));
    IdxSz = find(M488.A_tot(1,:)<2000);
    IdxFin = intersect(IdxEx,IdxSz);
    % 
    hold on
    % plot used AP2-clusters initially
    scatter3(CellCluster(:,1),CellCluster(:,2),CellCluster(:,3),7,[0.5 0 0],'filled')
    % darker green inside cell
    scatter3(M488.x_tot(1,intersect(IdxINco,IdxFin)),M488.y_tot(1,intersect(IdxINco,IdxFin)),M488.z_tot(1,intersect(IdxINco,IdxFin)),0.01*M488.A_tot(1,intersect(IdxINco,IdxFin)),[0 0.5 0],'filled')
    % light green on top of cell
    scatter3(M488.x_tot(1,intersect(IdxOUTco,IdxFin)),M488.y_tot(1,intersect(IdxOUTco,IdxFin)),M488.z_tot(1,intersect(IdxOUTco,IdxFin)),0.01*M488.A_tot(1,intersect(IdxOUTco,IdxFin)),[0 0.6 0.3],'filled')
    % Associated Viral particles (in blue)
    scatter3(M488.x_tot(1,intersect([IdxINco IdxOUTco],IdxFin)),M488.y_tot(1,intersect([IdxINco IdxOUTco],IdxFin)),M488.z_tot(1,intersect([IdxINco IdxOUTco],IdxFin)),0.003*M488.A_tot(1,intersect([IdxINco IdxOUTco],IdxFin)),[0.3010, 0.7450, 0.9330],'filled')
    hold off
    
    pause
end
% USING IMAGE 4


