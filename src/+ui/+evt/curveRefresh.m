function curveRefresh(~,~,f,evtIdxVec1,evtIdxVec2)
% curveRefresh draw single or multiple dff curves

fh = guidata(f);
opts = getappdata(f,'opts');
sz = opts.sz;
fh.showcurves = [];
ofstGap = 0.3;
ofst0 = 0;

ax = fh.curve;
% delete existing curves
hh = findobj(ax,'Type','line');
delete(hh);
hh = findobj(ax,'Type','text');
delete(hh);
ax.XLim = [0,sz(3)+1];

xxMax = -inf;
xxMin = inf;

if(~isempty(evtIdxVec1))
    dffMat = getappdata(f,'dffMat1');
    fts = getappdata(f,'fts1');
    evtIdxVec = evtIdxVec1;
    % save the status
    fh.showcurves = [evtIdxVec,1];
    guidata(f,fh);
    n = round(fh.sldMov.Value);

    xx = double(reshape((dffMat(evtIdxVec,:,1)),numel(evtIdxVec),[]));
    xxMin = min(xxMin,min(xx(:)));
    xxMax = max(xxMax,max(xx(:)));

    % draw new curves
    for ii=1:numel(evtIdxVec)
        evtIdx = evtIdxVec(ii);
%         t0 = fts.curve.rgt1(evtIdx,1);
%         t1 = fts.curve.rgt1(evtIdx,2);
        t0 = fts.curve.tBegin(evtIdx);
        t1 = fts.curve.tEnd(evtIdx);
        x = xx(ii,:);

        if sum(numel(evtIdxVec1) + numel(evtIdxVec2))==1
            xAll = dffMat(evtIdx,:,1);

    %         line(ax,1:sz(3),xAll,'Color',[0.5 0.5 0.5],'LineWidth',1);
            line(ax,1:sz(3),x*0,'Color','k','LineStyle','--');
    %         line(ax,1:sz(3),x,'Color','b');
    %         line(ax,t0:t1,x(t0:t1),'Color','r');
            line(ax,1:sz(3),xAll,'Color','b','LineWidth',1);
            line(ax,t0:t1,xAll(t0:t1),'Color','r','LineWidth',1);
        else
            ofst0 = ofst0 + ofstGap;
            x = x + ofst0;
            col = rand(1,3); col = col/sum(col);
            line(ax,1:sz(3),x,'Color',col);
            line(ax,t0:t1,x(t0:t1),'Color',col,'LineWidth',1.5);
        end
        xSel = x(t0:t1);
        [xm,ixm] = max(xSel);
        txt0 = [num2str(evtIdx),' dff:',num2str(fts.curve.dffMax(evtIdx))];
        text(ax,ixm+t0-1,xm,txt0);
    end
end

%% channel 2
if(~isempty(evtIdxVec2))
    dffMat = getappdata(f,'dffMat2');
    fts = getappdata(f,'fts2');
    evtIdxVec = evtIdxVec2;
    % save the status
    fh.showcurves = [fh.showcurves;evtIdxVec,2];
    guidata(f,fh);
    n = round(fh.sldMov.Value);

    ofstGap = 0.3;

    xx = double(reshape((dffMat(evtIdxVec,:,1)),numel(evtIdxVec),[]));
    xxMin = min(xxMin,min(xx(:)));
    xxMax = max(xxMax,max(xx(:)));
  
    % draw new curves
    for ii=1:numel(evtIdxVec)
        evtIdx = evtIdxVec(ii);
        % t0 = fts.curve.rgt1(evtIdx,1);
        % t1 = fts.curve.rgt1(evtIdx,2);
        t0 = fts.curve.tBegin(evtIdx);
        t1 = fts.curve.tEnd(evtIdx);
        x = xx(ii,:);

        if sum(numel(evtIdxVec1) + numel(evtIdxVec2))==1
            xAll = dffMat(evtIdx,:,1);

    %         line(ax,1:sz(3),xAll,'Color',[0.5 0.5 0.5],'LineWidth',1);
            line(ax,1:sz(3),x*0,'Color','k','LineStyle','--');
    %         line(ax,1:sz(3),x,'Color','b');
    %         line(ax,t0:t1,x(t0:t1),'Color','r');
            line(ax,1:sz(3),xAll,'Color','b','LineWidth',1);
            line(ax,t0:t1,xAll(t0:t1),'Color','r','LineWidth',1);
        else
            ofst0 = ofst0 + ofstGap;
            x = x + ofst0;
            col = rand(1,3); col = col/sum(col);
            line(ax,1:sz(3),x,'Color',col);
            line(ax,t0:t1,x(t0:t1),'Color',col,'LineWidth',1.5);
        end
        xSel = x(t0:t1);
        [xm,ixm] = max(xSel);
        txt0 = [num2str(evtIdx),' dff:',num2str(fts.curve.dffMax(evtIdx))];
        text(ax,ixm+t0-1,xm,txt0);
    end 
end
xxMax = xxMax + ofst0;
xxRg = xxMax-xxMin;
ax.YLim = [xxMin-xxRg*0.1,xxMax+xxRg*0.2];
line(ax,[n,n],ax.YLim,'Color','k','LineStyle','--');
end



