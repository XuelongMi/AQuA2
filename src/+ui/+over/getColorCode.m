function [col0,cMap] = getColorCode(f,nEvt,cType,cVal,nCh)
% getColorCode color codee each event
%
% nEvt: number of events
% cType: color code name
% cVal: strength of each event for the given feature
% cMap: obsolete
if(~isempty(f))
    opts = getappdata(f,'opts');
else
    opts.singleChannel = false; 
end
cMap = [];

if ~exist('cType','var') || isempty(cType)
    cType = 'Random';
end

col0 = zeros(nEvt,3);

if(opts.singleChannel)
    if strcmp(cType,'Random')
        for nn=1:nEvt
            x = rand(1,3);
            while (x(1)>0.8 && x(2)>0.8 && x(3)>0.8) || sum(x)<1
                x = rand(1,3);
            end
            x = x/max(x);
            col0(nn,:) = x;
        end
        return
    end
else
   if(nCh==1) 
       if strcmp(cType,'Random')
            for nn=1:nEvt
                x = rand(1,3);
                while x(2)<0.8 || x(1)>0.5 || sum(x)<1
                    x = rand(1,3);
                end
                x = x/max(x);
                col0(nn,:) = x;
            end
            return
       end
   else
       if strcmp(cType,'Random')
            for nn=1:nEvt
                x = rand(1,3);
                while x(1)<0.8 || x(2)>0.5 || sum(x)<1
                    x = rand(1,3);
                end
                x = x/max(x);
                col0(nn,:) = x;
            end
            return
       end
   end
end

% other color schemes
btSt = getappdata(f,'btSt');
colName = btSt.colNames;
colMapLst = btSt.colVals;

idxGood = strcmp(cType,colName);
idx = find(idxGood,1);
cMap = colMapLst{idx};

sclx = [];
sclx.minOv = min(cVal);
sclx.maxOv = max(cVal);

col0 = ui.over.reMapCol(cMap,cVal,sclx);

end







