%% setup
% -- preset 1: original Lck. 2: Jen Lck. 5: GluSnFR
% 
%
% Read Me:
% 'p0' is the folder containing tifs you want to deal with.
% Suggest sort the files in order, so that you can set the parameters 
% conviniently. AQuA/cfg/parameters_for_batch is the parameters excel.
% The script will read the parameters from that excel to deal with data.
% How many files you have, how many parameter settings should be in that excel.

close all;
clc;
clearvars
startup;  % initialize

pIn = 'X:\Methoxy\20211029_5mon_CCKCre5xFAD_14dpiu_DIOhM3DmCh_GFAPGCAmpLCK_1dpi_Methoxy\Mxxx1\4_TSereis_Resonant_ACSF\Green\'; %% tif folder
pOut = 'X:\Methoxy\20211029_5mon_CCKCre5xFAD_14dpiu_DIOhM3DmCh_GFAPGCAmpLCK_1dpi_Methoxy\Mxxx1\4_TSereis_Resonant_ACSF\GreenRes\'; %% tif folder
mkdir(pOut);
files = dir(fullfile(pIn,'*.tif'));
for xxx = 1:numel(files)
    p0 = [pIn,files(xxx).name,'\'];
    f1 = files(xxx).name; 
    load('random_Seed.mat');
    rng(s);
    %% Setting for this dataset
    opts = util.parseParam(1);
    opts.singleChannel = true;

    opts.registrateCorrect = 2;
    opts.bleachCorrect = 3;

    opts.noiseEstimation = 2;
    opts.smoXY = 3;
    opts.thrARScl = 3;
    opts.minSize = 100;
    opts.maxSize = 10000;
    opts.minDur = 5;
    opts.circularityThr = 0.1;
    opts.compress = 0.3;

    opts.needTemp = true;
    opts.ratio = 0.5;
    opts.sigThr = 3.5;
    opts.maxDelay = 10;
    opts.needRefine = false;
    opts.needGrow = false;

    opts.needSpa = true;
    opts.cRise = 5;
    opts.cDelay = 2;
    opts.gtwSmo = 1;

    opts.ignoreTau = true;

    [datOrg1,datOrg2,opts] = burst.prep1(pIn,f1,pIn,[],[],opts);
    [H,W,T] = size(datOrg1);

    %% preprocessing
    if(opts.registrateCorrect == 2)
        [datOrg1,~,tforms] = reg.regCrossCorrelation(datOrg1,[]);
    end

    if(opts.bleachCorrect==2)
        [datOrg1] = bleach_correct(datOrg1);
    elseif(opts.bleachCorrect==3)
        [datOrg1] = burst.bleach_correct2(datOrg1,opts);
    end
    opts.sz = size(datOrg1);
    %% active region
    xx = randn(10000,opts.cut); %% bias need to be reconsidered
    xxMA = movmean(xx,opts.movAvgWin,2);
    xxMin = min(xxMA,[],2);
    opts.xBias = nanmean(xxMin(:));
    [dat1,dF1,dFOrg1,opts] = burst.baselineRemoveAndNoiseEstimation(datOrg1,opts);
    opts.maxdF1 = quantile(dF1(:),0.999);
    [arLst1] = burst.acDetect(dF1,opts,true(H,W));

    %% temporal segmentation
    if(opts.needTemp)
         [seLst1,subEvtLst1,seLabel1,majorInfo1,opts,sdLst1,~,~] = se.seDetection(dF1,dFOrg1,arLst1,opts);
    else
        seLst1 = arLst1;
    end

    %% spatial segmentation
    if(opts.needTemp && opts.needSpa)
        [subEvtLst1,seLabel1,majorInfo1] = burst.preSplitSub(subEvtLst1,seLabel1,majorInfo1,size(dF1));
        [riseLst1,datR1,evtLst1,~] = burst.se2evtTop(dF1,seLst1,subEvtLst1,seLabel1,majorInfo1,opts);
    else
        evtLst1 = seLst1;
        datR1 = 255*uint8(ones(size(datOrg1)));
        riseLst1 = [];
    end

    %% global signal detection


    %% feature extraction
    [ftsLstE1, dffMat1, dMat1,dffAlignedMat1] = fea.getFeaturesTop(datOrg1, evtLst1, opts);
    ftsLstE1 = fea.getFeaturesPropTop(dat1,datR1 , evtLst1, ftsLstE1, opts);

    %% save output
    % export 
    evt1 = evtLst1;
    fts1 = ftsLstE1;
    evt2 = []; fts2 = []; dffMat2 = []; dMat2=[]; riseLst2 = []; dF2 = [];kernelLst2 = [];datR2 = ones(size(dF1))*255;dFOrg2 = [];dFOrg1 = [];
    res.maxVal = opts.maxValueDat;
    vSave0 = {...  % basic variables for results analysis
        'opts','ov','datOrg1','datOrg2','evt1','evt2','fts1','fts2','dffMat1','dMat1',...
        'dffMat2','dMat2','riseLst1','riseLst2','dF1','dF2','tforms',...
        'dFOrg1','dFOrg2'};

    ov = containers.Map('UniformValues',0);
    ov('None') = [];
    ovName = 'Events';
    fprintf('Overlay for events...\n')
    ov1 = ui.over.getOv([],evtLst1,opts.sz,datR1,1);
    ov1.name = ovName;
    ov1.colorCodeType = {'Random'};
    ov([ovName,'_Red']) = ov1;
    ov2 = ui.over.getOv([],evt2,opts.sz,datR2,2);
    ov2.name = ovName;
    ov2.colorCodeType = {'Random'};
    ov([ovName,'_Green']) = ov2;

    opts.bitNum = 8;
    opts.maxdF1 = max(dF1(:));
    opts.maxdF2 = opts.maxdF1;

    vSave = vSave0;
    datOrg1 = datOrg1*255;
    datOrg2 = [];
    res = [];
    for ii=1:numel(vSave)
        v0 = vSave{ii};
        res.(v0) = eval(v0);
    end
    bd = containers.Map;
    bd('None') = [];
    res.bd = bd;

    %% save output
    name = [f1(1:end-4)];

    save([pOut,'\',name,'_AQuA.mat'], 'res','-v7.3');   
end
