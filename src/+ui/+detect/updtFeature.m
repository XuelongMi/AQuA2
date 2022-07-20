function updtFeature(~, ~, f, stg)
    % updtFeature update network features after user draw regions
    % regions are all in x,y coordinate, where y need to be flipped for matrix manipulation

    fprintf('Updating basic, network, region and landmark features\n')

    % read data
    ov = getappdata(f, 'ov');
    opts = getappdata(f, 'opts');
    evtLst1 = getappdata(f, 'evt1');
    evtLst2 = getappdata(f, 'evt2');
    
    fh = guidata(f);
    fh.nEvtName.String = 'nEvt';
    if(~opts.singleChannel)
        fh.nEvt.String = [num2str(numel(evtLst1)),' | ',num2str(numel(evtLst2))];
    else
        fh.nEvt.String = [num2str(numel(evtLst1))];
    end
    
    gg = waitbar(0, 'Updating features ...');
    sz = opts.sz;

    % gather data
    fprintf('Gathering data ...\n')
    ov0 = ov('Events_Red');
    datR1 = zeros(sz, 'uint8');

    for tt = 1:sz(3)
        ov00 = ov0.frame{tt};
        dRecon00 = zeros(sz(1), sz(2));
        if isempty(ov00)
            continue
        end 
        for ii = 1:numel(ov00.idx)
            pix00 = ov00.pix{ii};
            val00 = ov00.val{ii};
            dRecon00(pix00) = uint8(val00 * 255);
        end 
        datR1(:, :, tt) = dRecon00;
    end 
    
    if(~opts.singleChannel)
        ov0 = ov('Events_Green');
        datR2 = zeros(sz, 'uint8');

        for tt = 1:sz(3)
            ov00 = ov0.frame{tt};
            dRecon00 = zeros(sz(1), sz(2));
            if isempty(ov00)
                continue
            end 
            for ii = 1:numel(ov00.idx)
                pix00 = ov00.pix{ii};
                val00 = ov00.val{ii};
                dRecon00(pix00) = uint8(val00 * 255);
            end 
            datR2(:, :, tt) = dRecon00;
        end 
    end

    % basic features
    waitbar(0.2, gg);
    dat1 = getappdata(f, 'dat1');
    datOrg1 = getappdata(f, 'datOrg1');
%     dF1 = getappdata(f, 'dF1');

    if stg == 0
        fprintf('Updating basic features ...\n')
        [ftsLstE1, dffMat1, dMat1,dffAlignedMat1] = fea.getFeaturesTop(datOrg1, evtLst1, opts,gg);
        setappdata(f, 'dffMat1', dffMat1);
        setappdata(f, 'dMat1', dMat1);
        setappdata(f, 'dffAlignedMat1', dffAlignedMat1);
        % filter table init
        setappdata(f,'fts1',ftsLstE1);
    else 
        ftsLstE1 = getappdata(f, 'fts1');
    end 
    
    if(~opts.singleChannel)
        dat2 = getappdata(f, 'dat2');
        datOrg2 = getappdata(f, 'datOrg2');
%         dF2 = getappdata(f, 'dF2');

        if stg == 0
            fprintf('Updating basic features ...\n')
            [ftsLstE2, dffMat2, dMat2,dffAlignedMat2] = fea.getFeaturesTop(datOrg2, evtLst2, opts,gg);
            setappdata(f, 'dffMat2', dffMat2);
            setappdata(f, 'dMat2', dMat2);
            setappdata(f, 'dffAlignedMat2', dffAlignedMat2);
            % filter table init
            setappdata(f,'fts2',ftsLstE2);
        else 
            ftsLstE2 = getappdata(f, 'fts2');
        end 
    end
    ui.detect.filterInit([],[],f);
    
    waitbar(0.6, gg);
    
%     % ROI info % connect event with ROI
%     bd = getappdata(f,'bd');
%     if bd.isKey('roi')
%         % ROI map
%         ROIinfo = bd('roi');
%         H = sz(1);
%         W = sz(2);
%         ROImap = zeros(H,W);
%         for i = 1:numel(ROIinfo)
%            ROImap(ROIinfo{i}.pix)  = i;
%         end
%         ROIlabel = cell(1,numel(ftsLstE.basic.center));
%         for i = 1:numel(ftsLstE.basic.center)
%             center = ftsLstE.basic.center{i};
%             numLabel = ROImap(center(1),center(2));
% %             if(numLabel>0)
% %                 ROIlabel{i} = ROIinfo{numLabel}.name;
% %             else
% %                 ROIlabel{i} = 'NaN';
% %             end
%             ROIlabel{i} = num2str(numLabel);
%         end     
%         ftsLstE.ROIlabel = ROIlabel;
%         setappdata(f,'fts',ftsLstE);
%     end

    % propagation features
    ftsLstE1 = fea.getFeaturesPropTop(dat1, datR1, evtLst1, ftsLstE1, opts);
%     ftsLstE1.propDir = getFeaturesPropDirection(dat1, evtLst1, opts);
    ftsLstE1.channel = 1;
    setappdata(f,'fts1',ftsLstE1);
    if(~opts.singleChannel)
        ftsLstE2 = fea.getFeaturesPropTop(dat2, datR2, evtLst2, ftsLstE2, opts);
        ftsLstE2.channel = 2;
%         ftsLstE2.propDir = getFeaturesPropDirection(dat2, evtLst2, opts);
    else
        ftsLstE2 = [];
    end
    setappdata(f,'fts2',ftsLstE2);
    
     waitbar(0.8, gg);
    % region, landmark, network and save results
    ui.detect.updtFeatureRegionLandmarkNetworkShow(f, datR1, evtLst1, ftsLstE1, gg,1);
    if(~opts.singleChannel)
        ui.detect.updtFeatureRegionLandmarkNetworkShow(f, datR2, evtLst2, ftsLstE2, gg,2);
    end
    ui.over.updtEvtOvShowLst([],[],f);
    waitbar(1, gg);
    % feature table
    ui.detect.getFeatureTable(f);
    fprintf('Done.\n')
    delete(gg)

end 



