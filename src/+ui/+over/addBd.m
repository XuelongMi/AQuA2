function addBd(f,axNow,lst,bdCol,n,nCh)

btSt = getappdata(f,'btSt');
opts = getappdata(f,'opts');
flexLst = getappdata(f,'flexLst');

if ~isfield(opts,'minShowEvtGUI')
    opts.minShowEvtGUI = 0.5;
end

H = opts.sz(1);
name00 = btSt.overlayDatSel;
if(nCh==1)
   ovName = 'Events_Red';
   name00 = [name00,'_Red'];
   ftsName = 'fts1';
else
   ovName = 'Events_Green';
   name00 = [name00,'_Green'];
   ftsName = 'fts2';
end
if ~isempty(lst) && strcmp(name00,ovName) && n>0
    ov = getappdata(f,'ov');
    ov0 = ov(ovName);
    x0 = ov0.frame{n};
    if ~isempty(x0)
        idx = x0.idx;
        fts = getappdata(f,ftsName);
        bds = fts.bds;
        loc2D = fts.loc.x2D;
        for ii=1:numel(idx)
            if sum(idx(ii)==lst)>0
                % only draw when area is large enough
                nPixTot = numel(loc2D{idx(ii)});
                nPixNow = numel(x0.pix{ii});
                if nPixNow/nPixTot>opts.minShowEvtGUI  % FIXME: do not make duration too long
                    xyC = bds{idx(ii)};
                    for jj=1:numel(xyC)
                        xy = xyC{jj};
                        flexLst{end+1} = patch(axNow,'XData',xy(:,2),'YData',H-xy(:,1)+1,...
                            'FaceColor','none','EdgeColor',bdCol,'Tag','flex','LineWidth',1); %#ok<AGROW>
                        if jj==1
                            flexLst{end+1} = text(axNow,xy(1,2)+1,H-xy(1,1),num2str(idx(ii)),...
                                'Color',bdCol,'FontSize',18,'Tag','flex'); %#ok<AGROW>
                        end
                    end
                end
            end
        end
    end
end

setappdata(f,'flexLst',flexLst);

end
