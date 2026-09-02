%% Single Molecule Calibration of VHH AF647 on VSV MeGFP EboV GP

%{
Strategy:
- make FlatField corrected images by BG-extracted Flatfields
- use 40% mask from BG-FlatFields

    - to Detect SM @
        - 05s
        - 07p5s
        - 10s
        - 15s
        - 20s
    - to Detect VHH bound to EboV VHHs
        - VHH G10 440nM
        - VHH G10 44nM
        - VHH G84 440nM
        - VHH G84 44nM
        - VHH 68flu 440nM

- Use VSV MeGFP EboV GP alone to estimate the Virus population
    - use the eGFP SM_calibration to estimarte # eGFP/Virion
    - use this information to estimate the intensity of a single virion

- altogether this gives the "Binding-distribution" of VHH on VSV EboV(!)

%}

%% Sort files and extract overview of different VSV samples

SortFiles_RH.m

% data loading
data_G10_440nM = loadConditionData('Parameters', [1.4 94.5 16e-6], 'MovieSelector', 'Ex');
data_G10_440nM_incrExp = loadConditionData('Parameters', [1.4 94.5 16e-6], 'MovieSelector', 'Ex');
data_G10_44nM = loadConditionData('Parameters', [1.4 94.5 16e-6], 'MovieSelector', 'Ex');

data_G84_440nM = loadConditionData('Parameters', [1.4 94.5 16e-6], 'MovieSelector', 'Ex');
data_G84_440nM_incrExp = loadConditionData('Parameters', [1.4 94.5 16e-6], 'MovieSelector', 'Ex');
data_G84_44nM = loadConditionData('Parameters', [1.4 94.5 16e-6], 'MovieSelector', 'Ex');

data_G68flu_440nM = loadConditionData('Parameters', [1.4 94.5 16e-6], 'MovieSelector', 'Ex');

data_VSVEboV_alone = loadConditionData('Parameters', [1.4 94.5 16e-6], 'MovieSelector', 'Ex');

save('datafiles')

%% Use common BG-derived FlatField to correct every imaging 

% PS: BG FF can only be derived at very large exposures (i.e. at SM)
% needs to be done for each channel though

save('IllumProf_median.mat','IllumProf_median')
save('Masks.mat','Masks')

IP_med_Norm = IllumProf_median/max(IllumProf_median(:));
save('IP_med_Norm.mat','IP_med_Norm')

data_FldNames = fieldnames(datafiles);
for i = 1:numel(datafiles)
    for j = 1:numel(datafiles.(data_FldNames{i}))
        rtdir = datafiles.(data_FldNames{i})(j).source;
        pos = find(rtdir==filesep);
        Chdir = rtdir(1:pos(end-1));
        cd(Chdir)
        for k = 1:numel(datafiles.(data_FldNames{i})(j).channels)
            mkdir(['Ch' num2str(k-1) 'FFcorr'])
            Imdir = datafiles.(data_FldNames{i})(j).channels{k};
            ImName = dir(fullfile(Imdir, '*.tif'));
            A = imread([Imdir ImName.name]);
            A_FFcorr = im2double(A)./IP_med_Norm; % completely overcompensates!!
            
        end
    end
end
dir(fullfile(data(1).source, '*.tif'));

%% Using the defined 40% mask to run Detections

% loading masks and data
Mask20 = Masks.Msk5;

sigma488 = getGaussianPSFsigma(1.4, 94.5, 16e-6, 'gfp');
sigma560 = getGaussianPSFsigma(1.4, 94.5, 16e-6, 'rfp');
sigma647 = getGaussianPSFsigma(1.4, 94.5, 16e-6, 'alexa647');

% run Detection in 40% mask on all images

data_FldNames = fieldnames(datafiles);
for i = 1:length(fieldnames(datafiles))
    tic
    DataSpec = data_FldNames{i};
    data = datafiles.(DataSpec);
    if numel(data(1).channels)==1
        runDetection(data,'Sigma', sigma488, 'CellMask', Mask20, 'Overwrite', true);
    elseif numel(data(1).channels)==2
        runDetection(data,'Sigma', [sigma488, sigma647], 'CellMask', Mask20, 'Overwrite', true);
    end
    toc
end

% NB: compare Detections with Mask20 and full detections

% inspection
figure,
A = dir(fullfile(data(2).channels{1}, '*.tif'));
Im488 = imread([data(2).source, A.name]);
imagesc(Im488)
axis image
hold on
scatter(frameInfo.x(1,:), frameInfo.y(1,:), 'o', 'r')

% PS: x_init and y_init marks the initial peaks (for signal deletion)

% collection of intensity data and coordinates from 2chan experiments
data_FldNames = fieldnames(datafiles);
for i = [1, 3, 4, 5, 7, 8]
    DataName = data_FldNames(i);
    GFP = [];
    AF647 = [];
    Dpfr = []; % detections per frame (for tracking positions)
    x = [];
    y = [];
    data = datafiles.(DataName{1}); % temp structure
    for j = 1:numel(data)
        cd([data(j).channels{1} 'Detection'])
        load('detection_v2.mat')
        GFP = [GFP frameInfo.A(1,:)];
        AF647 = [AF647 frameInfo.A(2,:)];
        x = [x frameInfo.x(1,:)];
        y = [y frameInfo.y(1,:)];
        Dpfr = [Dpfr length(frameInfo.x)];
    end
    DataStrct.(DataName{1}).GFP = GFP;
    DataStrct.(DataName{1}).AF647 = AF647;
    DataStrct.(DataName{1}).x = x;
    DataStrct.(DataName{1}).y = y;
    DataStrct.(DataName{1}).Dpfr = Dpfr;
end

% save data structure
pos = find(data(1).source==filesep);
rtdir = data(1).source(1:pos(end-5));
cd(rtdir)

save('DataStrct_2chan.mat', 'DataStrct')

%% DarkCurrent, BG, Noise and SNR

data_DC = loadConditionData('Parameters', [1.4 94.5 16e-6], 'MovieSelector', 'Ex');

% concatenate all 10 DC images
for i = 1:numel(data_DC)
    ImFileName = dir(fullfile([data(i).source '*.tif']));
    SD_IDI_DCim_temp = readtiff([data(i).source ImFileName.name]);
    SD_IDI_DCim(:,:,i) = SD_IDI_DCim_temp;
end

% averaging over 10 FOVs with s.d.
SD_IDI_avgBK = mean(SD_IDI_DCim(:));
tmpstd = std(single(SD_IDI_DCim),0,3);
SD_IDI_avgSig = mean(tmpstd(:));
SD_IDI_DC_FWHM = 2.355 * SD_IDI_avgSig;

% show in figures, surface plot
figure,
SCm = surfc(SD_IDI_DCim(:,:,1));
title('SD IDI Dark Current')
ax1 = gca;
zticks(ax1, [(SD_IDI_avgBK-SD_IDI_DC_FWHM) SD_IDI_avgBK (SD_IDI_avgBK+SD_IDI_DC_FWHM)])
zticklabels(ax1, {'BG - FWHM(noise)', 'BG', 'BG + FWHM(noise)'})
zlim([(SD_IDI_avgBK-4*SD_IDI_DC_FWHM) (SD_IDI_avgBK+4*SD_IDI_DC_FWHM)])
Smcontour = SCm(2);
Smcontour.ContourZLevel = (SD_IDI_avgBK-4*SD_IDI_DC_FWHM);

% show in figures, histogram
figure,
histogram(SD_IDI_DCim(:))
xlim([(SD_IDI_avgBK-4*SD_IDI_DC_FWHM) (SD_IDI_avgBK+4*SD_IDI_DC_FWHM)])
hold on
vline(SD_IDI_avgBK, 'k', 'avr BG')
vline((SD_IDI_avgBK-SD_IDI_DC_FWHM), 'r', 'BG - s.d.')
vline((SD_IDI_avgBK+SD_IDI_DC_FWHM), 'r', 'BG + s.d.')
title('SD IDI Dark Current (+- FWHM)')
xlabel('Intensity (a.u.)')


%% Assess Virus Distributions and potential "Bleed-Through"

histogram(DataStrct.data_VSVEboV_alone.GFP(:))
% majority of particles below 60a.u.
GFP_all = DataStrct.data_VSVEboV_alone.GFP(:);
GFP_below60 = GFP_all(GFP_all<60);
GFP_above60 = GFP_all(GFP_all>=60);

AF647_all = DataStrct.data_VSVEboV_alone.AF647(:);
AF647_GFPbelow60 = AF647(GFP_all<60);
AF647_GFPabove60 = AF647(GFP_all>=60);

figure,
scatter(GFP_above60, AF647_GFPabove60)
hold on
hline(0, 'k', 'background')
hline(44.136, 'r', 'BG + FWHM (noise)')
hline(-44.136, 'r', 'BG - FWHM (noise')
title('Bleed-through of large GFP clusters')
xlabel('GFP intensity (a.u.)')
ylabel('AF647 intensity (a.u.)')

figure,
hist(GFP_above60, 100);
title('Large GFP clusters above SNR = 3')
xlabel('GFP intensities (a.u.)')
      
figure,
hist(GFP_below60, 100);
title('Small GFP clusters, SNR < 3')
xlabel('GFP intensities (a.u.)')  

% Single viruses = 110 (?)

%% Checking Binding to viruses

cut0ff = 80; % GFP cutoff

FIGdir = uigetdir; % FIGURES folder

data_FldNames = fieldnames(DataStrct);
for i = 1:length(fieldnames(DataStrct))
    cd(FIGdir)
    GFP = DataStrct.(data_FldNames{i}).GFP(:);
    VSVEboV = GFP/110;
    AF647 = DataStrct.(data_FldNames{i}).AF647(:);
    VHH647 = AF647/15;
    idx = find(GFP>=80);
    % GFP histogram
    figVSVEboV = figure;
    hist(VSVEboV, 100)
    title([data_FldNames{i} 'VSVEboV number'])
    xlabel('# VSV EboV virions (a.u.)')
    xlim([0.5 10])
    mkdir(data_FldNames{i})
    cd([FIGdir filesep data_FldNames{i}])
    savefig([data_FldNames{i} 'VSVEboV'])
    % AF647 histogram
    figVHH647 = figure;
    hist(VHH647, 100)
    title([data_FldNames{i} 'VHHAF647'])
    xlabel('VHH-AF647-numbers (a.u.)')
    savefig([data_FldNames{i} 'VHHAF647'])
    % scatter plot
    figure
    scatter(VSVEboV, VHH647);
    hold on
    hline(0, 'k', 'background')
    hline(1.47, 'r', 'BG + 0.5*FWHM (noise)')
    hline(-1.47, 'r', 'BG - 0.5*FWHM (noise')
    title([data_FldNames{i} 'Scatter Plot'])
    xlabel('# VSV EboV virions (a.u.)')
    ylabel('VHH-AF647-numbers (a.u.)')
    savefig([data_FldNames{i} 'NumericScatterPlot'])
end

%% Convert to single molecule counts

time = [5000, 7500, 10000, 15000, 20000];

mu = [776, 941, 1562, 2286, 2976];

f=fittype('a*x');
ft=fit(time',mu',f);
figure,
plot(time,mu,'.','Markersize',10), hold on, plot(ft)
set(gca,'Fontsize',20)
xlabel('Exposure Time (msec)')
ylabel('Intensity (mu/a.u.)')
a = ft.a;

xlim([0 25000])
ylim([0 3250])
legend(['a = ' num2str(a)])
title('Single Molecule Calibration of AF647')

% Expected Intensity per VHH = 15



% are there to many aggregates??

% is the binding specific??




