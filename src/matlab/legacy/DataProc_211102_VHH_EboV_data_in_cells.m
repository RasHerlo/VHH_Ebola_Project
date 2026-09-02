%% Data Processing for VHH-647, EboV-GFP, RFP-AP2 cells

% F:\LLSM_VSV_VHH\2021 Processing In Vivo

%{
OBS: G10 data has bad analysis file in:

F:\LLSM_VSV_VHH\SUM51 for hidde summary\SUM51 for hidde summary\Cs3_G10\

Checking
"Ex10_488_150mW_70p560n_500mW_100p_642_150mW_100p_z0p4\ch488nmCamA\DS"
- Stacks look ok after DS (see images from PPT

We need an analysis file that will give "Detection3D" as output
Detection3D: "frameInfo-structure w/26 fields

Check GU_runDetection3D in "E:\MATLAB\Matlab
Repositories\GU_Repository\GokulScripts"

Input: data-structure

Use GU_loadConditionData3D

fx. "data = loadConditionData('Parameters', [1.45 160 1e-7],'MovieSelector', 'Image');"
...but in updated 3D-version (selectDS-data to load)...

PATH: "C:\Users\Ras\Desktop\Code Repositories\VHH_CODES\GokulScripts"

%}

%% Only two images in each

%% F:\LLSM_VSV_VHH\SUM51 for hidde summary\SUM51 for hidde summary\Cs3_G10\

[data] = GU_loadConditionData3D;

fs = strfind(data(1).source,filesep);
rtdir = data(1).source(1:fs(end-2));
save(fullfile(rtdir,'data488master.mat'),'data')

runDetection3D(data,'Sigma',[1.26, 1.32; 1.41, 1.38;  1.58, 1.608;],'Overwrite', true)

% OBS: FILES DAMAGED (lacking single stacks)
% Using 'F:\LLSM_VSV_VHH\TEMPORARYDATA\CS3_G10_Hidde1' instead

% Concatenate the data
x_tot = [];
y_tot = [];
z_tot = [];
A_tot = [];
ExpNum_all = [];
for g = 1:numel(data)
    load(fullfile(data(g).source,'Analysis','Detection3D.mat'))
    x_temp = frameInfo(1).x;
    y_temp = frameInfo(1).y;
    z_temp = frameInfo(1).z;
    A_temp = frameInfo(1).A;
    ExpNum = (g+9)*ones(1,size(A_temp,2));
    x_tot = [x_tot x_temp];
    y_tot = [y_tot y_temp];
    z_tot = [z_tot z_temp];
    A_tot = [A_tot A_temp];
    ExpNum_all = [ExpNum_all ExpNum];
end
M488.x_tot = x_tot;
M488.y_tot = y_tot;
M488.z_tot = z_tot;
M488.A_tot = A_tot;
M488.ExpNum = ExpNum_all;
% For each virus particle, check out the association with RFP as a function
% of VHH-binding

figure, histogram(M488.A_tot(1,:))

for h = 10:14
    idx = find(M488.ExpNum == h);
    figure, scatter3(M488.x_tot(1,idx),M488.y_tot(1,idx),M488.z_tot(1,idx),M488.A_tot(1,idx)/250,'g','filled')
end

% smaller limit can be found at:
x = 297.1130;
y = 339.3214;
z = 25.1932;
ExpNum = 4; idx = find(M488.ExpNum == 13);

[Mi,Xi] = mink(abs(M488.x_tot(1,:)-x),10);

M488.A_tot(1,3819) % A = 1347.8

save('M488.mat','M488')

% Limit set at 1000
figure, scatter(M488.A_tot(3,M488.A_tot(1,:)>1000),M488.A_tot(2,M488.A_tot(1,:)>1000))
xlabel('VHH-Alexa647-values')
ylabel('RFP-AP2-values')
title('EboV-GFP-particles > 1000a.u.')
set(gcf,'Color','w')

% Consider three things:
% A) What is the expected size of a viral particle, and normalize VHHs to
% number of particles
% B) What is the normal RFP-value for a pit (?)
% - run 560-master of data, and do "single-pit" fit
% C) ...Check fractions of particles associating to "at least one pit" (binary)

%% Making 560-master and Data-file of G10-data

[data] = GU_loadConditionData3D;

fs = strfind(data(1).source,filesep);
rtdir = data(1).source(1:fs(end-2));
save(fullfile(rtdir,'data560master.mat'),'data')

runDetection3D(data,'Sigma',[1.41, 1.38;  1.26, 1.32; 1.58, 1.608;],'Overwrite', true)

% Concatenate the data
x_tot = [];
y_tot = [];
z_tot = [];
A_tot = [];
ExpNum_all = [];
for g = 1:numel(data)
    load(fullfile(data(g).source,'Analysis','Detection3D.mat'))
    x_temp = frameInfo(1).x;
    y_temp = frameInfo(1).y;
    z_temp = frameInfo(1).z;
    A_temp = frameInfo(1).A;
    ExpNum = (g+9)*ones(1,size(A_temp,2));
    x_tot = [x_tot x_temp];
    y_tot = [y_tot y_temp];
    z_tot = [z_tot z_temp];
    A_tot = [A_tot A_temp];
    ExpNum_all = [ExpNum_all ExpNum];
end
M560.x_tot = x_tot;
M560.y_tot = y_tot;
M560.z_tot = z_tot;
M560.A_tot = A_tot;
M560.ExpNum = ExpNum_all;

save('M560.mat','M560')

%% Single Virions in LLSM

rtdir = 'F:\LLSM_VSV_VHH\20190415_p5_p55_sCMOS_RH_SM_MeGFP_AF647\CS1_VSV_MeGFP_EboV_VHHG84_AF647_44nM';

cd(rtdir)

% [data] = GU_loadConditionData3D; Error

FN = dir(rtdir); FN = FN(3:9);

% save structure
save(fullfile(rtdir,'M488.mat'),'M488')

figure, histogram(A_tot) % Cannot find the data files ??

% USE ANYTHING FOR Single Virions, fx. 1000

%% 

A_tot = M488.A_tot;

VHH = M488.A_tot(3,:);
VHH(VHH<0) = 0;
VHH_norm = VHH./M488.A_tot(1,:);

figure, scatter3(VHH_norm(M488.A_tot(1,:)>500),M488.A_tot(2,M488.A_tot(1,:)>500),M488.A_tot(1,M488.A_tot(1,:)>500))
xlabel('VHH-Alexa647-values')
ylabel('RFP-AP2-values')
zlabel('VSV-EboV-eGFP-values')
title('EboV-GFP-particles > 500a.u.')
set(gcf,'Color','w')

% Try using only values GFP <= 6000
Idx = find(M488.A_tot(1,:)<6000 & M488.A_tot(1,:)>200);

figure, scatter3(VHH_norm(Idx),M488.A_tot(2,Idx),M488.A_tot(1,Idx))
xlabel('VHH-Alexa647-values')
ylabel('RFP-AP2-values')
zlabel('VSV-EboV-eGFP-values')
title('EboV-GFP-particles [200:6000]')
set(gcf,'Color','w')

%% G84 

rtdir = 'F:\LLSM_VSV_VHH\SUM52\CS2_G84';
rtdir = 'G:\LLSM_VSV_VHH\SUM52\CS2_G84';


[data] = GU_loadConditionData3D;

save('dataM488.mat','data')

% Concatenate the data
x_tot = [];
y_tot = [];
z_tot = [];
A_tot = [];
ExpNum_all = [];
for g = 1:numel(data)
    load(fullfile(data(g).source,'Analysis','Detection3D.mat'))
    x_temp = frameInfo(1).x;
    y_temp = frameInfo(1).y;
    z_temp = frameInfo(1).z;
    A_temp = frameInfo(1).A;
    ExpNum = (g+9)*ones(1,size(A_temp,2));
    x_tot = [x_tot x_temp];
    y_tot = [y_tot y_temp];
    z_tot = [z_tot z_temp];
    A_tot = [A_tot A_temp];
    ExpNum_all = [ExpNum_all ExpNum];
end
M488.x_tot = x_tot;
M488.y_tot = y_tot;
M488.z_tot = z_tot;
M488.A_tot = A_tot;
M488.ExpNum = ExpNum_all;

save('M488.mat','M488')

clearvars -except M488

VHH = M488.A_tot(3,:);
VHH(VHH<0) = 0;
VHH_norm = VHH./M488.A_tot(1,:);

% Try using only values GFP <= 6000
Idx = find(M488.A_tot(1,:)<6000 & M488.A_tot(1,:)>200);

figure, scatter3(VHH_norm(Idx),M488.A_tot(2,Idx),M488.A_tot(1,Idx))
xlabel('VHH-Alexa647-values')
ylabel('RFP-AP2-values')
zlabel('VSV-EboV-eGFP-values')
title('EboV-GFP-particles [200:6000]')
set(gcf,'Color','w')

%% Making 560-master and Data-file of G84-data

% OBS: Files lacking from DS-files, again, so continue with only two first

% FROM SCRATCH AGAIN WITH TEMPORARY FILES.... including only two first
% images in each experiment

rtdir = 'G:\LLSM_VSV_VHH\TEMPORARYDATA\CS_G84';
cd(rtdir)

% data488master
[data] = GU_loadConditionData3D; % [488; 560; 647]

fs = strfind(data(1).source,filesep);
rtdir = data(1).source(1:fs(end-2));
save(fullfile(rtdir,'data488master.mat'),'data')

runDetection3D(data,'Sigma',[1.26, 1.32; 1.41, 1.38; 1.58, 1.608;],'Overwrite', true)

% Concatenate the data
x_tot = [];
y_tot = [];
z_tot = [];
A_tot = [];
ExpNum_all = [];
for g = 1:numel(data)
    load(fullfile(data(g).source,'Analysis','Detection3D.mat'))
    x_temp = frameInfo(1).x;
    y_temp = frameInfo(1).y;
    z_temp = frameInfo(1).z;
    A_temp = frameInfo(1).A;
    ExpNum = (g+4)*ones(1,size(A_temp,2));
    x_tot = [x_tot x_temp];
    y_tot = [y_tot y_temp];
    z_tot = [z_tot z_temp];
    A_tot = [A_tot A_temp];
    ExpNum_all = [ExpNum_all ExpNum];
end
M488.x_tot = x_tot;
M488.y_tot = y_tot;
M488.z_tot = z_tot;
M488.A_tot = A_tot;
M488.ExpNum = ExpNum_all;

save('M488.mat','M488')

% data560master
[data] = GU_loadConditionData3D; % [560; 488; 647]

fs = strfind(data(1).source,filesep);
rtdir = data(1).source(1:fs(end-2));
save(fullfile(rtdir,'data560master.mat'),'data')

runDetection3D(data,'Sigma',[1.41, 1.38;  1.26, 1.32; 1.58, 1.608;],'Overwrite', true)

% Concatenate the data
x_tot = [];
y_tot = [];
z_tot = [];
A_tot = [];
ExpNum_all = [];
for g = 1:numel(data)
    load(fullfile(data(g).source,'Analysis','Detection3D.mat'))
    x_temp = frameInfo(1).x;
    y_temp = frameInfo(1).y;
    z_temp = frameInfo(1).z;
    A_temp = frameInfo(1).A;
    ExpNum = (g+4)*ones(1,size(A_temp,2));
    x_tot = [x_tot x_temp];
    y_tot = [y_tot y_temp];
    z_tot = [z_tot z_temp];
    A_tot = [A_tot A_temp];
    ExpNum_all = [ExpNum_all ExpNum];
end
M560.x_tot = x_tot;
M560.y_tot = y_tot;
M560.z_tot = z_tot;
M560.A_tot = A_tot;
M560.ExpNum = ExpNum_all;

save('M560.mat','M560')






%% G10 second dataset is the same as first G10 data set (!!)

%% NoVHH

rtdir = 'F:\LLSM_VSV_VHH\SUM51 for hidde summary\SUM51 for hidde summary\CS4_NoVHH';

[data] = GU_loadConditionData3D; % Not working for some reason

FN = dir(rtdir); FN = FN([4 5 7]);

% OBS: Ex17 lacking Detection3D-file

% Concatenate the data
x_tot = [];
y_tot = [];
z_tot = [];
A_tot = [];
ExpNum_all = [];
for g = 1:numel(FN)
    load(fullfile(FN(g).folder,FN(g).name,'ch488nmCamA','Analysis','Detection3D.mat'))
    x_temp = frameInfo(1).x;
    y_temp = frameInfo(1).y;
    z_temp = frameInfo(1).z;
    A_temp = frameInfo(1).A;
    ExpNum = (g+9)*ones(1,size(A_temp,2));
    x_tot = [x_tot x_temp];
    y_tot = [y_tot y_temp];
    z_tot = [z_tot z_temp];
    A_tot = [A_tot A_temp];
    ExpNum_all = [ExpNum_all ExpNum];
end
M488.x_tot = x_tot;
M488.y_tot = y_tot;
M488.z_tot = z_tot;
M488.A_tot = A_tot;
M488.ExpNum = ExpNum_all;

save('M488.mat','M488')

clearvars -except M488

VHH = M488.A_tot(3,:);
VHH(VHH<0) = 0;
VHH_norm = VHH./M488.A_tot(1,:);

% Try using only values GFP <= 6000
Idx = find(M488.A_tot(1,:)<6000 & M488.A_tot(1,:)>200);

figure, scatter3(VHH_norm(Idx),M488.A_tot(2,Idx),M488.A_tot(1,Idx))
xlabel('VHH-Alexa647-values')
ylabel('RFP-AP2-values')
zlabel('VSV-EboV-eGFP-values')
title('EboV-GFP-particles [200:6000]')
set(gcf,'Color','w')

%% No VHH - Checking 560nm master population

rtdir = 'F:\LLSM_VSV_VHH\SUM51 for hidde summary\SUM51 for hidde summary\CS4_NoVHH';

FN = dir(rtdir); FN = FN([4 5 7]);

% Concatenate the data
x_tot = [];
y_tot = [];
z_tot = [];
A_tot = [];
ExpNum_all = [];
for g = 1:numel(FN)
    load(fullfile(FN(g).folder,FN(g).name,'ch560nmCamB','Analysis','Detection3D.mat'))
    x_temp = frameInfo(1).x;
    y_temp = frameInfo(1).y;
    z_temp = frameInfo(1).z;
    A_temp = frameInfo(1).A;
    ExpNum = (g+9)*ones(1,size(A_temp,2));
    x_tot = [x_tot x_temp];
    y_tot = [y_tot y_temp];
    z_tot = [z_tot z_temp];
    A_tot = [A_tot A_temp];
    ExpNum_all = [ExpNum_all ExpNum];
end
M560.x_tot = x_tot;
M560.y_tot = y_tot;
M560.z_tot = z_tot;
M560.A_tot = A_tot;
M560.ExpNum = ExpNum_all;

save('M560.mat','M560')

clearvars -except M560

figure, histogram(M560.A_tot(1,:))
title('RFP values master560'), set(gcf,'Color','w')

%% NoMOI - Check distribution of RFP-clusters naturally

% Find minimum sizes, if exists

rtdir = 'F:\LLSM_VSV_VHH\SUM52\CS1_SUM52_NoMOI';
rtdir = 'G:\LLSM_VSV_VHH\SUM52\CS1_SUM52_NoMOI';

[data] = GU_loadConditionData3D;
% M560, 488, 647

save('dataM560.mat','data')

% Concatenate the data
x_tot = [];
y_tot = [];
z_tot = [];
A_tot = [];
ExpNum_all = [];
for g = 1:numel(data)
    load(fullfile(data(g).source,'Analysis','Detection3D.mat'))
    x_temp = frameInfo(1).x;
    y_temp = frameInfo(1).y;
    z_temp = frameInfo(1).z;
    A_temp = frameInfo(1).A;
    ExpNum = (g+9)*ones(1,size(A_temp,2));
    x_tot = [x_tot x_temp];
    y_tot = [y_tot y_temp];
    z_tot = [z_tot z_temp];
    A_tot = [A_tot A_temp];
    ExpNum_all = [ExpNum_all ExpNum];
end
M560.x_tot = x_tot;
M560.y_tot = y_tot;
M560.z_tot = z_tot;
M560.A_tot = A_tot;
M560.ExpNum = ExpNum_all;

save('M560.mat','M560')

clearvars -except M560

figure, histogram(M560.A_tot(1,:))
title('RFP measures in NoMOI samples')
set(gcf, 'Color','w')

%% COMPARE FRACTION OF VSVs with RFP

% G10:
G10dir = 'F:\LLSM_VSV_VHH\TEMPORARYDATA\CS3_G10_Hidde1';

% G84: 
G84dir = 'F:\LLSM_VSV_VHH\SUM52\CS2_G84';

% NoVHH:
NoVHHdir = 'F:\LLSM_VSV_VHH\SUM51 for hidde summary\SUM51 for hidde summary\CS4_NoVHH';

FN = fieldnames(M488all);

% make piediagrams of association
figure,
for f = 1:numel(FN)
    subplot(1,numel(FN),f)
    Idx488 = find(M488all.(FN{f}).A_tot(1,:)>200 & M488all.(FN{f}).A_tot(1,:)<6000);
    A = size(M488all.(FN{f}).A_tot(:,Idx488),2);
    B = length(find(M488all.(FN{f}).A_tot(2,Idx488)>20));
    pie([A-B B])
    title(FN{f})
end
set(gcf,'Color','w')

% VHH distributions
figure,
for f = 1:numel(FN)
    subplot(2,numel(FN),f)
    Idx488 = find(M488all.(FN{f}).A_tot(1,:)>200 & M488all.(FN{f}).A_tot(1,:)<6000);
    Idx560 = find(M488all.(FN{f}).A_tot(2,:)>20);
    IdxVHH = intersect(Idx488,Idx560);
    IdxnoVHH = setdiff(1:size(M488all.(FN{f}).A_tot,2),IdxVHH);
    histogram(M488all.(FN{f}).A_tot(3,IdxVHH));
    title([FN{f} ' EboV in AP2'])
    subplot(2,numel(FN),f+numel(FN))
    histogram(M488all.(FN{f}).A_tot(3,IdxnoVHH));
    title([FN{f} ' EboV outside AP2'])
end
set(gcf,'Color','w')




