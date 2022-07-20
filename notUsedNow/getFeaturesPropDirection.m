function propDir = getFeaturesPropDirection(dat, evtLst, opts)
    [H,W,T] = size(dat);
%     evtMap = zeros(H,W,T);
%     for i = 1:numel(evtLst)
%         evtMap(evtLst{i}) = i;
%     end
    
%     datSmo = imgaussfilt3(dat,[0.001,0.001,2]);
%     datVec = reshape(datSmo,[],T);
% %     evtMap = reshape(evtMap,[],T);
%     thrVec = [0.4,0.5,0.6];
%     % 50% rising time
%     for i = 1:numel(evtLst)
%         [ih,iw,it] = ind2sub([H,W,T],evtLst{i});
%         t0 = min(it);
%         t1 = max(it);
%         ihw = unique(sub2ind([H,W],ih,iw));
%         dFCurves = datVec(ihw,t0:t1);
%         F0 = min(dFCurves,[],2);
%         dFCurves = dFCurves - repmat(F0,1,t1-t0+1);
%         [maxV] = max(dFCurves,[],2);
%         dFCurves = dFCurves./repmat(maxV,1,t1-t0+1);
%         
%         risingTime = T*ones(numel(ihw),3);
%         for k = 1:numel(thrVec)
%             select = (dFCurves>thrVec(k));
%             for t = 1:t1-t0+1
%                 risingTime(select(:,t),k) = min(risingTime(select(:,t),k),t);
%             end
%         end
%         risingTime = mean(risingTime,2);
%         propMap = nan(H,W);
%         propMap(ihw) = risingTime;
%         figure;imagesc(propMap);colorbar;
%         
%         t0 = ceil(min(risingTime));
%         t1 = max(risingTime);
%         curProp = zeros(1,4);
%         initialR = 1;
%         for t = t0:t1
%             ihw0 = ihw(risingTime<=t);
%             map = false(H,W);
%             map(ihw0) = true;
%             map = bwareaopen(map,20,4);  
%             CC = bwconncomp(map);
%             CC = CC.PixelIdxList;
%             sz = cellfun(@numel,CC);
%             [~,id] = max(sz);
%             ihw0 = CC{id};
% %             ihw0 = find(map);
%             [ih0,iw0] = ind2sub([H,W],ihw0);
%             if(numel(ihw0)>numel(ihw)/10)
%                 if(initialR)
%                    upB = min(ih0);
%                    downB = max(ih0);
%                    westB = min(iw0);
%                    eastB = max(iw0); 
%                    initialR = 0;
%                    ts = t;
%                 else
%                    curProp = max(curProp,[upB-min(ih0),max(ih0)-downB,westB-min(iw0),max(iw0)-eastB]);
%                 end
%             end
%         end
%         propDirV(i,:) = curProp;
%     end
    
    
    
    
    
    
    
    %% propagation direction
    propDirV = zeros(numel(evtLst),4);
    for i = 1:numel(evtLst)
        [ih,iw,it] = ind2sub([H,W,T],evtLst{i});
        ihw = unique(sub2ind([H,W],ih,iw));
        upB = min(ih);
        downB = max(ih);
        westB = min(iw);
        eastB = max(iw); 
        for t = min(it):max(it)
            select = it==t;
            ih0 = ih(select);
            iw0 = iw(select);
            ihw0 = sub2ind([H,W],ih0,iw0);
            map = false(H,W);
            map(ihw0) = true;
            map = bwareaopen(map,20,4);  
            CC = bwconncomp(map);
            CC = CC.PixelIdxList;
            sz = cellfun(@numel,CC);
            [maxSz,id] = max(sz);
            minSz = min(sz);
            if(maxSz/minSz>20)
                ihw0 = CC{id};
            else
                ihw0 = find(map);
            end
            [ih0,iw0] = ind2sub([H,W],ihw0);
            if(numel(ih0)>max(numel(ihw)/10,opts.minSize))
               upB = min(ih0);
               downB = max(ih0);
               westB = min(iw0);
               eastB = max(iw0);               
               break;
            end
        end
        propDirV(i,:) = [upB-min(ih),max(ih)-downB,westB-min(iw),max(iw)-eastB];
    end
    
    minDistThr = 50;
    %% classification
    propDir = cell(numel(evtLst),1);
    for i = 1:numel(evtLst)
        if(max(propDirV(i,:))<minDistThr)
            propDir{i} = 'Stationary Event';
        else
            upDown = max(propDirV(i,1:2));
            leftRight = max(propDirV(i,3:4));
            if(upDown>=leftRight/3)
                propDir{i} = 'Has vertical propagation';
            else
                west = propDirV(i,3);
                east = propDirV(i,4);
                if(west/east<1.5 && east/west<1.5)
                    propDir{i} = 'Horizontal: no obvious direction';
                else
                    if(west>east)
                        propDir{i} = 'Towards left';
                    else
                        propDir{i} = 'Towards right';
                    end
                end
            end
        end
    end
end