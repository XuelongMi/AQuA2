function datxCol = movStep(f,n,ovOnly,updtAll)
    % use btSt.sbs, btSt.leftView and btSt.rightView to determine what to show
    opts = getappdata(f,'opts');
    fh = guidata(f);
    scl = getappdata(f,'scl');
    btSt = getappdata(f,'btSt');
    
    if ~exist('ovOnly','var') || isempty(ovOnly)
        ovOnly = 0;
    end

    if ~exist('updtAll','var') || isempty(updtAll)
        updtAll = 0;
    end
        
    if ~exist('n','var') || isempty(n)
        n = round(fh.sldMov.Value);
    end
    
    %% channel 1
    if ~isfield(btSt,'GaussFilter') ||(btSt.GaussFilter==0) 
        dat1 = getappdata(f,'datOrg1');
    else
        dat1 = getappdata(f,'dat1');  
    end
    if isempty(dat1) dat1 = getappdata(f,'dat1'); end
    dat1 = dat1(:,:,n);
    
    if (scl.map==1) dat1 = dat1.^2;  end
    dat1 = (dat1-scl.min)/max(scl.max-scl.min,0.01)*scl.bri1;
    
    %% channel 2
    if(~opts.singleChannel)
        if ~isfield(btSt,'GaussFilter') ||(btSt.GaussFilter==0) 
            dat2 = getappdata(f,'datOrg2');
        else
            dat2 = getappdata(f,'dat2');   
        end
        if isempty(dat2) dat2 = getappdata(f,'dat2'); end
        dat2 = dat2(:,:,n);
        
        if (scl.map==1) dat2 = dat2.^2;  end
        dat2 = (dat2-scl.min)/max(scl.max-scl.min,0.01)*scl.bri2;
        datx = cat(3,dat2,dat1,zeros(size(dat1)));
    else
        datx = cat(3,dat1,dat1,dat1);
    end

    datxCol = ui.over.addOv(f,datx,n);
    
    if ovOnly>0
        return
    end
    
    %% overlay
    if btSt.sbs==0
        fh.ims.im1.CData = flipud(datxCol);
        ui.mov.addPatchLineText(f,fh.mov,n,updtAll);
    end
    if btSt.sbs==1
        datxL = datx*scl.briL;
        datxR = datx*scl.briR;
        viewName = {'leftView','rightView'};
        imName = {'im2a','im2b'};
        axLst = {fh.movL,fh.movR};
        for ii=1:2
            curType = btSt.(viewName{ii});
            axNow = axLst{ii};
            
            % clean all patches
            if updtAll>0
                types = {'quiver','line','patch','text'};
                for jj=1:numel(types)
                    h00 = findobj(axNow,'Type',types{jj});
                    if ~isempty(h00)
                        delete(h00);
                    end
                end
            end
            switch curType
                case 'Raw'
                    if ii==1
                        fh.ims.(imName{ii}).CData = flipud(datxL);
                    else
                        fh.ims.(imName{ii}).CData = flipud(datxR);
                    end
%                     ui.mov.addPatchLineText(f,axNow,n,updtAll);
                case 'dF'
                        dF1 = getappdata(f,'dF1');
                        if (~isempty(dF1)) dF1 = dF1(:,:,n);  else dF1 = dat1; end
                        dF1 = dF1/opts.maxdF1;  % re-scale
                        dF1 = (dF1-scl.min)/max(scl.max-scl.min,0.01)*scl.bri1;
                        if(~opts.singleChannel)
                            dF2 = getappdata(f,'dF2');
                            if (~isempty(dF2)) dF2 = dF2(:,:,n); else dF2 = dat2; end
                            dF2 = dF2/opts.maxdF2;
                            dF2 = (dF2-scl.min)/max(scl.max-scl.min,0.01)*scl.bri2;
                            dFx = cat(3,dF2,dF1,zeros(size(dF1)));
                        else
                            dFx = cat(3,dF1,dF1,dF1);
                        end
                        if ii==1
                            dFxL = dFx*scl.briL;
                            fh.ims.(imName{ii}).CData = flipud(dFxL);
                        else
                            dFxR = dFx*scl.briR;
                            fh.ims.(imName{ii}).CData = flipud(dFxR);
                        end
                case 'Raw + overlay'
                    if ii==1
                        datxColL = ui.over.addOv(f,datxL,n);
                        fh.ims.(imName{ii}).CData = flipud(datxColL);
                    else
                        datxColR = ui.over.addOv(f,datxR,n);
                        fh.ims.(imName{ii}).CData = flipud(datxColR);
                    end
                    ui.mov.addPatchLineText(f,axNow,n,updtAll);
                case 'Rising map'
                    ui.mov.showRisingMap(f,imName{ii},n);
                case 'Maximum Projection'
                    if scl.map==1
                        datM = fh.maxPro.^2;
                    end
                    if ii==1
                        datML = (datM-scl.min)/max(scl.max-scl.min,0.01)*scl.briL;
                        datMxL = cat(3,datML,datML,datML);
                        fh.ims.(imName{ii}).CData = flipud(datMxL);
                    else
                       datMR = (datM-scl.min)/max(scl.max-scl.min,0.01)*scl.briR;
                        datMxR = cat(3,datMR,datMR,datMR);
                        fh.ims.(imName{ii}).CData = flipud(datMxR);
                    end
                    ui.mov.addPatchLineText(f,axNow,n,updtAll);
                case 'Average Projection'
                    if scl.map==1
                        datM = fh.averPro.^2;
                    end
                    if ii==1
                        datML = (datM-scl.min)/max(scl.max-scl.min,0.01)*scl.briL;
                        datMxL = cat(3,datML,datML,datML);
                        fh.ims.(imName{ii}).CData = flipud(datMxL);
                    else
                       datMR = (datM-scl.min)/max(scl.max-scl.min,0.01)*scl.briR;
                        datMxR = cat(3,datMR,datMR,datMR);
                        fh.ims.(imName{ii}).CData = flipud(datMxR);
                    end
                    ui.mov.addPatchLineText(f,axNow,n,updtAll);    
                case 'Rising map'
                    ui.mov.showRisingMap(f,imName{ii},n);
            end
        end
    end
    
    %% finish
    % adjust area to show
    if btSt.sbs==0
        fh.mov.XLim = scl.wrg;
        fh.mov.YLim = scl.hrg;
    else
        fh.movL.XLim = scl.wrg;
        fh.movL.YLim = scl.hrg;
        fh.movR.XLim = scl.wrg;
        fh.movR.YLim = scl.hrg;
    end
    
    % frame number
    ui.mov.updtMovInfo(f,n,opts.sz(3));
    % f.Pointer = 'arrow';
    
end



