function res = saveExp(~,~,f,file0,path0,modex)
% saveExp save experiment (and export results)

fts = getappdata(f,'fts1');
if ~exist('modex','var')
    modex = 0;
end
if isempty(fts)
    msgbox('Please save after event detection\n');
    return
end
if ~exist(path0,'file') && ~isempty(path0)
    mkdir(path0);    
end


%% gather results
ff = waitbar(0,'Gathering results ...');

% if do not want to detect again, do not need to save dF
vSave0 = {...  % basic variables for results analysis
    'opts','scl','btSt','ov','bd','datOrg1','datOrg2','evt1','evt2','fts1','fts2','dffMat1','dMat1',...
    'dffMat2','dMat2','riseLst1','riseLst2','featureTable','userFeatures','dF1','dF2','cfuInfo1','cfuInfo2'};
% vSave1 = {...  % extra variables for event detection
%     'arLst','lmLoc','svLst','seLstAll','riseX','riseLstAll','evtLstAll','ftsLstAll',...
%     'dffMatAll','datRAll','evtLstFilterZ','dffMatFilterZ','tBeginFilterZ',...
%     'riseLstFilterZ','evtLstMerge','dF'...
% };
% vSave = [vSave0,vSave1];
vSave = vSave0;

res = [];
for ii=1:numel(vSave)
    v0 = vSave{ii};
    res.(v0) = getappdata(f,v0);
end

% filter features and curves
ov = getappdata(f,'ov');
ov1 = ov('Events_Red');
xSel1 = ov1.sel;
ov2 = ov('Events_Green');
xSel2 = ov2.sel;
btSt = getappdata(f,'btSt');
opts = getappdata(f,'opts');
% xSel = ov0.sel;
xSelFav1 = false(numel(res.evt1),1);
xSelFav1(btSt.evtMngrMsk1) = true;
xSelFav2 = false(numel(res.evt2),1);
xSelFav2(btSt.evtMngrMsk2) = true;

res.ftsFilter1 = util.filterFields(res.fts1,xSel1);
res.ftsFav1 = util.filterFields(res.fts1,xSelFav1);
res.evtFilter1 = res.evt1(xSel1);
res.evtFav1 = res.evt1(xSelFav1);
res.dffMatFilter1 = res.dffMat1(xSel1,:,:);
res.dffMatFav1 = res.dffMat1(xSelFav1,:,:);

res.ftsFilter2 = util.filterFields(res.fts2,xSel2);
res.ftsFav2 = util.filterFields(res.fts2,xSelFav2);
res.evtFilter2 = res.evt2(xSel2);
res.evtFav2 = res.evt2(xSelFav2);
res.dffMatFilter2 = res.dffMat2(xSel2,:,:);
res.dffMatFav2 = res.dffMat2(xSelFav2,:,:);

if ~isempty(res.dMat1)
    res.dMatFilter1 = res.dMat1(xSel1,:,:);
     res.dMatFav1 = res.dMat1(xSelFav1,:,:);
end
if ~isempty(res.dMat2)
    res.dMatFilter2 = res.dMat2(xSel2,:,:);
     res.dMatFav2 = res.dMat2(xSelFav2,:,:);
end

if ~isempty(res.riseLst1)  % rising map is for super events
    res.riseLstFilter1 = res.riseLst1(xSel1);
    res.riseLstFav1 = res.riseLst1(xSelFav1);
end
if ~isempty(res.riseLst2)  % rising map is for super events
    res.riseLstFilter2 = res.riseLst2(xSel2);
    res.riseLstFav2 = res.riseLst2(xSelFav2);
end

res.evtSelectedList1 = find(xSel1>0);
res.evtFavList1 = find(xSelFav1>0);
res.evtSelectedList2 = find(xSel2>0);
res.evtFavList2 = find(xSelFav2>0);

% save raw movie with 8 or 16 bits to save space
res.opts.bitNum = 16;
res.maxVal = nanmax(res.datOrg1(:));  
if(~opts.singleChannel)
    res.maxVal = max(res.maxVal,max(res.datOrg2(:)));  
end
res.datOrg1 = res.datOrg1/res.maxVal;
res.datOrg1 = uint16(res.datOrg1*(2^res.opts.bitNum-1));

res.datOrg2 = res.datOrg2/res.maxVal;
res.datOrg2 = uint16(res.datOrg2*(2^res.opts.bitNum-1));


res.stg.post = 1;
res.stg.detect = 1;

if modex>0
    waitbar(1,ff);
    delete(ff);
    return
end

%% export
fh = guidata(f);

% btSt = getappdata(f,'btSt');
favEvtLst1 = btSt.evtMngrMsk1;
favEvtLst2 = btSt.evtMngrMsk2;
fout = [path0,filesep,file0];
[fpath,fname,ext] = fileparts(fout);

if fh.expEvt.Value==1
    waitbar(0.25,ff,'Saving res file...');
    if ~strcmp(fout(end-3:end),'.mat')
        fout = [fout,'.mat'];
    end
    save(fout,'res','-v7.3');
end

if fh.expEvt2.Value==1
    waitbar(0.25,ff,'Saving res file...');
    if ~strcmp(fout(end-3:end),'.mat')
        fout = [fout,'.mat'];
    end
    field = {'datOrg1','datOrg2','dF1','dF2','ov'};
    res = rmfield(res,field);
    save(fout,'res','-v7.3');
end



% export movie
if fh.expMov.Value==1
    waitbar(0.5,ff,'Writing movie ...');
    ov1 = zeros(opts.sz(1),opts.sz(2),3,opts.sz(3));
    for tt=1:opts.sz(3)
        if mod(tt,100)==0
            fprintf('Frame %d\n',tt); 
        end
        ov1(:,:,:,tt) = ui.movStep(f,tt,1);
    end
    ui.movStep(f);
    fmov = [fpath,filesep,fname,'.tif'];
    io.writeTiffSeq(fmov,ov1,8);
end

if fh.expFt.Value==1
    waitbar(0.75,ff,'Writing feature table ...');
    % export feature table
    ftTb = getappdata(f,'featureTable1');
    if isempty(ftTb)
        ui.detect.getFeatureTable(f);
        ftTb = getappdata(f,'featureTable1');
    end
    cc = ftTb{:,1};

    % all selected events
    cc1 = cc(:,res.evtSelectedList1);
    ftTb1 = table(cc1,'RowNames',ftTb.Row);
    ftb = [fpath,filesep,fname,'_Ch1.csv'];
    writetable(ftTb1,ftb,'WriteVariableNames',0,'WriteRowNames',1);
    
    % for favorite events
    if ~isempty(favEvtLst1)
        cc00 = cc(:,favEvtLst1);
        ftTb00 = table(cc00,'RowNames',ftTb.Row);
        ftb00 = [fpath,filesep,fname,'_Ch1_favorite.xlsx'];
        writetable(ftTb00,ftb00,'WriteVariableNames',0,'WriteRowNames',1);
    end

    if(~opts.singleChannel)
        ftTb = getappdata(f,'featureTable2');
        if isempty(ftTb)
            ui.detect.getFeatureTable(f);
            ftTb = getappdata(f,'featureTable2');
        end
        cc = ftTb{:,1};

        % all selected events
        cc1 = cc(:,res.evtSelectedList2);
        ftTb1 = table(cc1,'RowNames',ftTb.Row);
        ftb = [fpath,filesep,fname,'_Ch2.csv'];
        writetable(ftTb1,ftb,'WriteVariableNames',0,'WriteRowNames',1);
        
            % for favorite events
        if ~isempty(favEvtLst2)
            cc00 = cc(:,favEvtLst2);
            ftTb00 = table(cc00,'RowNames',ftTb.Row);
            ftb00 = [fpath,filesep,fname,'_Ch2_favorite.xlsx'];
            writetable(ftTb00,ftb00,'WriteVariableNames',0,'WriteRowNames',1);
        end
    end
    
    %% curves
    dffAlignedMat1 = getappdata(f,'dffAlignedMat1');
    rowName = cell(size(dffAlignedMat1,1),1);
    for k = 1:numel(rowName)
       rowName{k}  = ['Event ',num2str(k)];
    end
    curves1 = table(num2cell(dffAlignedMat1),'RowNames',rowName);
    writetable(curves1,[fpath,filesep,fname,'_Ch1_curves.xlsx'],'WriteVariableNames',0,'WriteRowNames',1);
    if(~opts.singleChannel)
        dffAlignedMat2 = getappdata(f,'dffAlignedMat2');
        rowName = cell(size(dffAlignedMat2,1),1);
        for k = 1:numel(rowName)
           rowName{k}  = ['Event ',num2str(k)];
        end
        curves2 = table(num2cell(dffAlignedMat2),'RowNames',rowName);
        writetable(curves2,[fpath,filesep,fname,'_Ch2_curves.xlsx'],'WriteVariableNames',0,'WriteRowNames',1);
    end
    
    bd = getappdata(f,'bd');
    
    % for each region
    if ~isempty(fts.region) && isfield(fts.region.cell,'memberIdx') && ~isempty(fts.region.cell.memberIdx)
        bdcell = bd('cell');
        fpathRegion = [fpath,'\Regions'];
        if ~exist(fpathRegion,'file') && ~isempty(fpathRegion)
            mkdir(fpathRegion);    
        end

        memSel = fts.region.cell.memberIdx(res.evtSelectedList1,:);
        for ii=1:size(memSel,2)
            mem00 = memSel(:,ii);
            Name = 'None';
            if numel(bdcell{ii})>=4
                Name = bdcell{ii}{4};
            end
            if strcmp(Name,'None')
               Name = num2str(ii); 
            end
            if(sum(mem00>0)==0)
                continue;
            end
            cc00 = cc(:,mem00>0);
            ftTb00 = table(cc00,'RowNames',ftTb.Row);
            ftb00 = [fpathRegion,filesep,fname,'ch1_region_',Name,'.xlsx'];
            writetable(ftTb00,ftb00,'WriteVariableNames',0,'WriteRowNames',1);
        end
        
        if(~opts.singleChannel)
            fts = getappdata(f,'fts2');
            memSel = fts.region.cell.memberIdx(res.evtSelectedList2,:);
            for ii=1:size(memSel,2)
                mem00 = memSel(:,ii);
                Name = 'None';
                if numel(bdcell{ii})>=4
                    Name = bdcell{ii}{4};
                end
                if strcmp(Name,'None')
                   Name = num2str(ii); 
                end
                if(sum(mem00>0)==0)
                    continue;
                end
                cc00 = cc(:,mem00>0);
                ftTb00 = table(cc00,'RowNames',ftTb.Row);
                ftb00 = [fpathRegion,filesep,fname,'ch2_region_',Name,'.xlsx'];
                writetable(ftTb00,ftb00,'WriteVariableNames',0,'WriteRowNames',1);
            end
        end
    end

    
    % region and landmark map
    f00 = figure('Visible','off');
    dat = getappdata(f,'datOrg1');
    dat = mean(dat,3);
    dat = dat/max(dat(:));
    Low_High = stretchlim(dat,0.001);
    dat = imadjust(dat,Low_High);
    axNow = axes(f00);
    image(axNow,'CData',flipud(dat),'CDataMapping','scaled');
    axNow.XTick = [];
    axNow.YTick = [];
    axNow.XLim = [0.5,size(dat,2)+0.5];
    axNow.YLim = [0.5,size(dat,1)+0.5];
    axNow.DataAspectRatio = [1 1 1];
    colormap gray
    ui.mov.addPatchLineText(f,axNow,0,1)
    % saveas(f00,[fpath,filesep,fname,'_landmark.fig']);
    saveas(f00,[fpath,filesep,fname,'_landmark.png'],'png');
    delete(f00);

    % rising maps
%     riseLst = getappdata(f,'riseLst');
%     if ~isempty(favEvtLst)
%         f00 = figure('Visible','off');
%         axNow = axes(f00);
%         fpathRising = [fpath,filesep,'risingMaps'];
%         if ~exist(fpathRising,'file')
%             mkdir(fpathRising);
%         end
%         for ii=1:numel(favEvtLst)
%             rr = riseLst{favEvtLst(ii)};
%             imagesc(axNow,rr.dlyMap);
%             colorbar(axNow);
%             xx = axNow.XTickLabel;
%             for jj=1:numel(xx)
%                 xx{jj} = num2str(str2double(xx{jj})+min(rr.rgw)-1);
%             end
%             axNow.XTickLabel = xx;
%             xx = axNow.YTickLabel;
%             for jj=1:numel(xx)
%                 xx{jj} = num2str(str2double(xx{jj})+min(rr.rgw)-1);
%             end
%             axNow.YTickLabel = xx;
%             axNow.DataAspectRatio = [1 1 1];
%             saveas(f00,[fpathRising,filesep,num2str(favEvtLst(ii)),'.png'],'png');
%         end
%     end
end

delete(ff);

end







