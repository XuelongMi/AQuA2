function drawReg(~,~,f,op,lbl)
    % updtFeature update network features after user draw regions
    
    fh = guidata(f);
    bd = getappdata(f,'bd');
    btSt = getappdata(f,'btSt');
    opts = getappdata(f,'opts');
    
    if bd.isKey(lbl)
        bd0 = bd(lbl);
    else
        bd0 = [];
    end
    
    ax = fh.mov;
    if btSt.sbs==0
        ax = fh.mov;
    end
    if btSt.sbs==1
        ax = fh.movL;
    end
    
  
    if strcmp(op,'add')
        tmp = [];
        hh = impoly(ax);
        if ~isempty(hh)
            nPts = size(hh.getPosition,1);
            if nPts>2
                msk = flipud(hh.createMask);
                tmp{1} = bwboundaries(msk);
                tmp{2} = find(msk>0);
                tmp{3} = 'manual';
                tmp{4} = 'None';
                bd0{end+1} = tmp;
                delete(hh)
            end
        end
        
    end
    
    if strcmp(op,'extract') 
        if bd.isKey('roi')
            hh = impoly(ax);
            if ~isempty(hh)
                nPts = size(hh.getPosition,1);
                if nPts>2
                    msk = flipud(hh.createMask);
                    delete(hh);
                end
            end
            ROIinfo = bd('roi');
            ROImap = zeros(size(msk));
            for i = 1:numel(ROIinfo)
               ROImap(ROIinfo{i}.pix)  = i;
            end
            validROI = setdiff(ROImap(msk>0),0);
%             selectROIname = cell(numel(validROI),1);
%             for i = 1:numel(validROI)
%                 selectROIname{i} = ROIinfo{i}.name;
%             end
            selectROIname = cell(numel(validROI),1);
            for i = 1:numel(validROI)
                selectROIname{i} = num2str(validROI(i));
            end
            selectROITable = cell2table(selectROIname);
            setappdata(f,'featureTable',selectROITable);
            selpath = uigetdir(opts.filePath,'Choose output folder for exporting ROIs');
            selpath = [selpath,filesep,'selectingROIs.csv']
            writetable(selectROITable,selpath,'WriteVariableNames',0,'WriteRowNames',1);
        end
    end
    
    if strcmp(op,'arrow')
        opts = getappdata(f,'opts');
        hh = imline(ax);
        if ~isempty(hh)
            bd0 = hh.getPosition;
            opts.northx = bd0(2,1)-bd0(1,1);
            opts.northy = bd0(2,2)-bd0(1,2);
            setappdata(f,'opts',opts);
            delete(hh)
        end
    end
    
    bd(lbl) = bd0;
    setappdata(f,'bd',bd);
    f.Pointer = 'arrow';
    ui.movStep(f,[],[],1);
    
end







