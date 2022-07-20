function datxCol = addOv(f,datx,n)

datxCol = datx;
[H,W,~] = size(datx);

scl = getappdata(f,'scl');
btSt = getappdata(f,'btSt');
sclOv = scl.briOv;
if ~strcmp(btSt.overlayDatSel,'None')
    % show movie with overlay
    rPlane = zeros(H,W);
    gPlane = rPlane;
    bPlane = rPlane;
    
    ov = getappdata(f,'ov');
    if(btSt.channelSelect~=3)
        ov1 = ov([btSt.overlayDatSel,'_Red']);
        x1 = ov1.frame{n};
        c1 = ov1.col;
        if isfield(ov1,'colVal') && strcmp(btSt.overlayColorSel,'Random')==0
            v1 = ov1.colVal;
            c1 = ui.over.reMapCol(btSt.mapNow,v1,scl);
        end
        if ~isempty(x1)
            for ii=1:numel(x1.idx)
                idx1 = x1.idx(ii);
                if ov1.sel(idx1)>0
                    pix0 = x1.pix{ii};
                    val0 = x1.val{ii};
                    col0 = c1(idx1,:);
                    rPlane(pix0) = rPlane(pix0) + val0*col0(1);
                    gPlane(pix0) = gPlane(pix0) + val0*col0(2);
                    bPlane(pix0) = bPlane(pix0) + val0*col0(3);
                end
            end
        end
    end
    if(btSt.channelSelect~=2)
        ov2 = ov([btSt.overlayDatSel,'_Green']);
        x2 = ov2.frame{n};
        c2 = ov2.col;

        % remap color
        if isfield(ov2,'colVal') && strcmp(btSt.overlayColorSel,'Random')==0
            v2 = ov2.colVal;
            c2 = ui.over.reMapCol(btSt.mapNow,v2,scl);
        end

        if ~isempty(x2)
            for ii=1:numel(x2.idx)
                idx2 = x2.idx(ii);
                if ov2.sel(idx2)>0
                    pix0 = x2.pix{ii};
                    val0 = x2.val{ii};
                    col0 = c2(idx2,:);
                    rPlane(pix0) = rPlane(pix0) + val0*col0(1);
                    gPlane(pix0) = gPlane(pix0) + val0*col0(2);
                    bPlane(pix0) = bPlane(pix0) + val0*col0(3);
                end
            end
        end
    end
    datxCol(:,:,1) = rPlane*sclOv + datxCol(:,:,1);
    datxCol(:,:,2) = gPlane*sclOv + datxCol(:,:,2);
    datxCol(:,:,3) = bPlane*sclOv + datxCol(:,:,3);
end

end