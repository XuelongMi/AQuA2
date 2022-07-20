function [evt2Lst,major] = splitEvt_Ac(dF,dFOrg,evtLst,major,dFSmoVec,opts,Ttestopen)

    [H,W,T] = size(dFOrg);
    opts.sigThr = opts.sigThr-0.5;
    opts.ratio = 0.05;
    %% get Majority
    N = numel(major);    
    trivial = false(N,1);
    gap = opts.smoT*2+1;
    %% Split
    Map = zeros(H,W,T,'uint16');
    for i = 1:numel(evtLst)
       Map(evtLst{i})  = i;
    end

    dFVec = reshape(dFOrg,[],T);
    clear dFOrg;
%     dFSmoVec = reshape(imgaussfilt(dFOrg,4),[],T);
    Map = reshape(Map,[],T);
    evtVecLst = label2idx(Map);
    nEvt = N;
    minDur = opts.minDur;

    sourceEvt = zeros(N,1);
    for i = 1:N
       pix = evtVecLst{i} ;
       [~,it] = ind2sub([H*W,T],pix);
       ihw = major{i}.ihw;
       curve = mean(dFVec(major{i}.ihw,:),1);
       s0 = sqrt(median((curve(gap+1:end)-curve(1:end-gap)).^2)/0.9133);
       curve = curve/s0;    % curve = curve*sqrt(numel(ihw));
       curveSmo = imgaussfilt(curve,2);
       major{i}.TW = se.getMajorityTem(curveSmo,major{i}.TW,pix,[H,W,T]);
       t00 = min(it);
       t11 = max(it);
       curve00 = mean(dFVec(major{i}.ihw,:),1); curve00 = curve00(major{i}.TW);
%         curve00 = curve00 - ([0:numel(curve00)-1]/numel(curve00)*(curve00(end)-curve00(1))+curve00(1));
        curve00 = curve00 - min(curve00);
        curve00 = curve00/max(curve00);
        major{i}.curve = curve00;
       
       % initialization
       thrs = floor(max(curveSmo(t00:t11))):-opts.ratio:floor(min(curveSmo(t00:t11)));
       labels = true(T,1);
       labels(t00:t11) = false;
       TWcandidate = cell(0);
       cnt = 0;
       %% check peak
       for k = 1:numel(thrs)
           thr = thrs(k);
          
           curTWs = bwconncomp(curveSmo(t00:t11)>thr);
           curTWs = curTWs.PixelIdxList;
           sz = cellfun(@numel,curTWs);
           curTWs = curTWs(sz>=minDur);
            
           for j = 1:numel(curTWs)
               TW = curTWs{j} + t00 - 1;
               t0 = TW(1);
               t1 = TW(end);
               if(t0==t00 && t00>1 && curveSmo(t00-1)<thr) % left side
                   continue;
               end
               if(t1==t11 && t11<T && curveSmo(t11+1)<thr) % right side
                   continue;
               end
               
               % whether detected
               if(~isempty(find(labels(t0:t1),1)))
                  continue; 
               end
               fg = curve(TW);
               [bg1,bg2,nv1,nv2] = se.findNeighbor(curve,t0,t1,T,thr,curveSmo);
               [Tscore1,Tscore2] = se.calTScore(fg,bg1,bg2);
               if(Ttestopen && (Tscore1<=opts.sigThr || Tscore2<=opts.sigThr)) % T-test check
                  continue; 
               end
               [score1,score2] = se.calOrderScore(fg,bg1,bg2,nv1,nv2);           
%                    figure;plot(curve);hold on;plot(TW,curve(TW),'r')
                if(min(score1,score2)>opts.sigThr)    % Update
                    cnt = cnt + 1;
                    TWcandidate{cnt} = t0:t1;
                    labels(TW) = true;
                end
           end
       end

       if(numel(TWcandidate)==0)
          continue; 
       end
       
       maxSignal = max(curveSmo)-min(curveSmo);
       
       if(numel(TWcandidate)==1)
           major{i}.TW = se.getMajorityTem(curveSmo,TWcandidate{1},pix,[H,W,T]);
           curve00 = mean(dFVec(major{i}.ihw,:),1); curve00 = curve00(major{i}.TW);
%            curve00 = curve00 - ([0:numel(curve00)-1]/numel(curve00)*(curve00(end)-curve00(1))+curve00(1));
           curve00 = curve00 - min(curve00);
           curve00 = curve00/max(curve00);
           major{i}.curve = curve00;
           continue;
       end
       
        % find splitting margin
        tStart = cellfun(@(x)x(1),TWcandidate);
        [~,id] = sort(tStart);
        TWcandidate = TWcandidate(id);
        GapPoint = [];
        for k = 1:(numel(TWcandidate)-1)
            t0 = TWcandidate{k}(end) + 1;
            t1 = TWcandidate{k+1}(1) - 1;
            [~,tGap] = min(curveSmo(t0:t1));
            tGap = tGap + t0 - 1;
            GapPoint = [GapPoint;tGap];
        end

       % Extend
       seedTW = major{i}.TW;
       GapPoint = [t00;GapPoint;t11];
       validGap = true(numel(GapPoint),1);
       validPeak = true(numel(GapPoint)-1,1);
       % check
       for k = 1:(numel(GapPoint)-1)
            TW = (GapPoint(k)+1):GapPoint(k+1);
            if(k==1)
                TW = [t00,TW];
            end
            %% check whether significant
            t0 = min(TW);
            t1 = max(TW);
            [maxV,tPeak] = max(curveSmo(TW));
            tPeak = t0 + tPeak - 1;
            [minVLeft] = min(curveSmo(t0:tPeak));
            [minVRight] = min(curveSmo(tPeak:t1));
            minV = max(minVLeft,minVRight);
            if(maxV-minV<max(3,0.05*maxSignal))
                validPeak(k) = false;
                if(k+1 == numel(GapPoint))
                    t = k;
                    while(~validGap(t))
                        t = t-1;
                    end
                    validGap(t) = false;
                else
                    validGap(k+1) = false;
                end
            end
       end
       TWcandidate = TWcandidate(validPeak);
       GapPoint = GapPoint(validGap);
       
       %
       if(numel(TWcandidate)==0)
          continue; 
       elseif(numel(TWcandidate)==1)
            curTW = TWcandidate{1};
            t0 = curTW(1);   % update time window
            t1 = curTW(end);
            major{i}.TW = se.getMajorityTem(curveSmo,t0:t1,pix,[H,W,T]);
            curve00 = mean(dFVec(major{i}.ihw,:),1); curve00 = curve00(major{i}.TW);
%             curve00 = curve00 - ([0:numel(curve00)-1]/numel(curve00)*(curve00(end)-curve00(1))+curve00(1));
            curve00 = curve00 - min(curve00);
            curve00 = curve00/max(curve00(:));
            major{i}.curve = curve00;
            continue;
       end

%        szT = zeros(numel(TWcandidate),1);
%        for k = 1:numel(TWcandidate)
%             szT(k) = numel(intersect(TWcandidate{k},seedTW));
%        end
       szT = cellfun(@(x)numel(intersect(x,seedTW)),TWcandidate);
       [~,id] = max(szT);
       
       for k = 1:numel(TWcandidate)            
           select = it>=GapPoint(k)&it<=GapPoint(k+1);
           pix2 = pix(select);
           t0 = TWcandidate{k}(1);   % update time window
           t1 = TWcandidate{k}(end);
           %% majorPeak
           if(k==id)
                major{i}.TW = se.getMajorityTem(curveSmo,t0:t1,pix2,[H,W,T]);
                Map(pix2) = i;
                major{i}.needUpdatePeak = true;
                curve00 = mean(dFVec(major{i}.ihw,:),1); curve00 = curve00(major{i}.TW);
%                 curve00 = curve00 - ([0:numel(curve00)-1]/numel(curve00)*(curve00(end)-curve00(1))+curve00(1));
                curve00 = curve00 - min(curve00);
                curve00 = curve00/max(curve00(:));
                major{i}.curve = curve00;
           else
               [mIhw,TW,delays] = se.getRefineSpaMajority_Ac(ihw,dFSmoVec,[H,W,T],pix2,t0:t1,opts,false);
               % check trivial
                n0 = numel(intersect(ihw,mIhw));
                n1 = numel(ihw);
                n2 = numel(mIhw); 
                nEvt = nEvt+1;  
                trivial(nEvt) = ~(n0/n1>=opts.overlap & n0/n2>=opts.overlap);
                sourceEvt(nEvt) = i;
                major{nEvt}.TW = TW;
                major{nEvt}.ihw = mIhw;
                major{nEvt}.delays = delays;
                major{nEvt}.needUpdatePeak = true;
%                 major{nEvt}.orgSeed = major{i}.orgSeed;
                curve00 = mean(dFVec(mIhw,:),1); curve00 = curve00(TW);
%                 curve00 = curve00 - ([0:numel(curve00)-1]/numel(curve00)*(curve00(end)-curve00(1))+curve00(1));
                curve00 = curve00 - min(curve00);
                curve00 = curve00/max(curve00(:));
                major{nEvt}.curve = curve00;
                Map(pix2) = nEvt;
           end

       end 
    end
    
    Map = reshape(Map,[H,W,T]);
    evt2Lst = label2idx(Map);
    dh = [-1 0 1 -1 1 -1 0 1];
    dw = [-1 -1 -1 0 0 1 1 1];
    % trivial should be merged before merging step
    % to avoid wrong forbidden pair
    % And from principle, it should have large overlap with its source
    for i = (N+1):nEvt
       if(trivial(i))    % is trivial event
            pix = evt2Lst{i};
            [ih,iw,it] = ind2sub([H,W,T],pix);
            curTW = major{i}.TW;
            neib0 = [];
            % Find neighbors
            for ii=1:numel(dh)
                ih1 = min(max(ih + dh(ii),1),H);
                iw1 = min(max(iw + dw(ii),1),W);
                vox1 = sub2ind([H,W,T],ih1,iw1,it);
                idxSel = setdiff(Map(vox1),[0,i]);
                neib0 = union(neib0,idxSel);
            end
%             neib0 = setdiff(neib0,[0,i]);
            neib0 = neib0(neib0<=N);
            
            mergedLabel = [];
            if(~isempty(neib0))
                tOs = zeros(numel(neib0),1); 
                for k = 1:numel(neib0)
                    nLabel = neib0(k);
                    nTW = major{nLabel}.TW;
                    tOs(k) = numel(intersect(curTW,nTW))/numel(union(curTW,nTW));
                end
                [~,id] = max(tOs);
                if(tOs(id)>0.8)
                    mergedLabel = neib0(id);
                end
            end
            
            if(~isempty(mergedLabel))
%                 label = mergedLabel;
                Map(pix) = mergedLabel;
                evt2Lst{i} = [];
                evt2Lst{mergedLabel} = [evt2Lst{mergedLabel};pix];
                % Update
                [mIhw,TW,delays] = se.getRefineSpaMajority_Ac(major{mergedLabel}.ihw,dFSmoVec,[H,W,T],evt2Lst{mergedLabel},major{mergedLabel}.TW,opts,true,major{mergedLabel}.delays);
                major{mergedLabel}.TW = TW;
                major{mergedLabel}.ihw = mIhw;
                major{mergedLabel}.delays = delays;
                curve00 = mean(dFVec(mIhw,:),1); curve00 = curve00(TW);
%                 curve00 = curve00 - ([0:numel(curve00)-1]/numel(curve00)*(curve00(end)-curve00(1))+curve00(1));
                curve00 = curve00 - min(curve00);
                curve00 = curve00/max(curve00(:));
                major{mergedLabel}.curve = curve00;
            else    % no neighbor
                trivial(i) = false;
                if(numel(major{i}.ihw) < opts.minSize)
                    major{i}.ihw = se.getMajoritySpa(pix,curTW,[H,W,T],opts);
                    if(numel(major{i}.ihw) < opts.minSize)
                        trivial(i) = true;  % Still small
                        OrgLabel = sourceEvt(i);
                        evt2Lst{OrgLabel} = [evt2Lst{OrgLabel};pix];
                    end
                    major{i}.delays = zeros(numel(major{i}.ihw),1);
                    curve00 = mean(dFVec(mIhw,:),1); curve00 = curve00(major{i}.TW);
%                     curve00 = curve00 - ([0:numel(curve00)-1]/numel(curve00)*(curve00(end)-curve00(1))+curve00(1));
                    curve00 = curve00 - min(curve00);
                    curve00 = curve00/max(curve00(:));
                    major{i}.curve = curve00;
                end
            end
       end
    end
    
    
    major = major(~trivial);
    evt2Lst = evt2Lst(~trivial);
end