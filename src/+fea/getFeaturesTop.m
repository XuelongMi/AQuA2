function [ftsLst,dffMat,dMat,dffAlignedMat] = getFeaturesTop(dat,evtLst,opts,ff)
    % getFeaturesTop extract curve related features, basic features and propagation
    % dat: single (0 to 1)
    % evtMap: single ( integer)
    
    [H,W,T] = size(dat);
    evtMap = zeros(size(dat),'uint16');
    for ii=1:numel(evtLst)
        evtMap(evtLst{ii}) = ii;
    end
    
    if opts.usePG
        dat = dat.^2;
    end
    
    secondPerFrame = opts.frameRate;
    muPerPix = opts.spatialRes;
    
    % impute events
    fprintf('Imputing ...\n')
    datx = dat;
    datx(evtMap>0) = nan;
    datx = img.imputeMov(datx);
    
    if ~isfield(opts,'maxValueDat')
        opts.maxValueDat = 1;
    end
    if ~isfield(opts,'correctTrend')
        opts.correctTrend = 0;
    end
    if ~isfield(opts,'bgFluo')
        opts.bgFluo = 0;
    end
    
    Tww = min(opts.movAvgWin,T/4);
    opts.major = 0.8;
    
    % bias in moving average minimum
    %xx = randn(1000,T);
    %xxMovAvg = movmean(xx,Tww,2);
    %bbm = mean(min(xxMovAvg,[],2));
    bbm = 0;
    
    % bb = zeros(1,100);
    % for ii=1:100
    %     xx = randn(T,1);
    %     xxMovAvg = movmean(xx,Tww);
    %     bb(ii) = min(xxMovAvg);
    % end
    % bbm = mean(bb);
    
    ftsLst = [];
    ftsLst.basic = [];
    
    waitbar(0.3, ff);

    
    ftsLst.propagation = [];
    
    foptions = fitoptions('exp1');
    foptions.MaxIter = 100;
    
    dMat = zeros(numel(evtLst),T,2,'single');
    dffMat = zeros(numel(evtLst),T,2,'single');
    alignedrgT = -10:100;
    dffAlignedMat = nan(numel(evtLst),110,'single');
    datVec = reshape(dat,[],T);

    for ii=1:numel(evtLst)
        if mod(ii,100)==0
            fprintf('%d/%d\n',ii,numel(evtLst))
            waitbar(0.3 + 0.3*ii/numel(evtLst), ff);
        end
        pix0 = evtLst{ii};
        if isempty(pix0)
            continue
        end
        [ih,iw,it] = ind2sub([H,W,T],pix0);
        ihw = unique(sub2ind([H,W],ih,iw));
        rgH = max(min(ih)-1,1):min(max(ih)+1,H);
        rgW = max(min(iw)-1,1):min(max(iw)+1,W);
        rgT = max(min(it)-1,1):min(max(it)+1,T);
        
        if numel(rgT)==1
            continue
        end
        
        % dff
        
        voxi1 = evtMap(rgH,rgW,:);
        voxi1 = reshape(voxi1,[],T);               
        
        sigxy = sum(voxi1==ii,2);
        sigz = sum(voxi1,1)>0;
        if sum(sigz)>T/2
            sigz = sum(voxi1==ii,1)>0;
        end
        
        voxd1 = dat(rgH,rgW,:);
        voxd1 = reshape(voxd1,[],T);
        charxIn1 = mean(voxd1(sigxy>0,:),1);
        idx = sub2ind([numel(rgH),numel(rgW),T],ih-min(rgH)+1,iw-min(rgW)+1,it);
        evtData = voxd1(idx);
        clear voxd1;
%         charx1 = fea.curvePolyDeTrend(charxIn1,sigz,opts.correctTrend);
        charx1 = charxIn1;
        sigma1 = sqrt(median((charx1(2:end)-charx1(1:end-1)).^2)/0.9113);
        
%         charxBg1 = min(movmean(charx1,Tww));
        
        % correct baseline method
        charxBg1 = getBaseline(charx1,Tww,opts.cut);
        charxBg1 = charxBg1 - opts.xBias*sigma1;
%         charxBg1 = charxBg1 - bbm*sigma1 - opts.bgFluo^2;
        dff1 = (charx1-charxBg1)./charxBg1;
        sigma1dff = sqrt(median((dff1(2:end)-dff1(1:end-1)).^2)/0.9113);
        

        dff1Sel = dff1(rgT);
        dffMax1= max(dff1Sel);
        
        % dff without other events
        voxd2 = reshape(datx(rgH,rgW,:),[],T);
        voxd2(idx) = evtData;  % bring current event back
        charxIn2 = nanmean(voxd2(sigxy>0,:),1);
        clear voxd2;
%         charx2 = fea.curvePolyDeTrend(charxIn2,sigz,opts.correctTrend);
        charx2 = charxIn2;
%         charxBg2 = min(movmean(charx2,Tww));
        charxBg2 = getBaseline(charx2,Tww,opts.cut);
        charxBg2 = charxBg2 - opts.xBias*sigma1;
%         charxBg2 = charxBg2 - bbm*sigma1 - opts.bgFluo^2;
        dff2 = (charx2-charxBg2)./charxBg2;
        
        if 1  % for p values
            dff2Sel = dff2(rgT);
            [dffMax2,tMax] = max(dff2Sel);
            xMinPre = max(min(dff2Sel(1:tMax)),sigma1dff);
            xMinPost = max(min(dff2Sel(tMax:end)),sigma1dff);
            dffMaxZ = max((dffMax2-xMinPre+dffMax2-xMinPost)/sigma1dff/2,0);
            dffMaxPval = 1-normcdf(dffMaxZ);
        end
        
        % extend event window in the curve
        voxi1(voxi1==ii) = 0;
        sigxOthers = sum(voxi1>0,1)>0;
        [dff2e,rgT1] = fea.extendEventTimeRangeByCurve(dff2,sigxOthers,it);
        
        % curve features
        [ rise19,fall91,width55,width11,decayTau,pp,riseTime] = fea.getCurveStat( ...
            dff2e, secondPerFrame, foptions, opts.ignoreTau );
        riseTime = riseTime + min(rgT1) - 1;

        dffMat(ii,:,1) = single(dff1);
        dffMat(ii,:,2) = single(dff2);
        dMat(ii,:,1) = single(charx1*opts.maxValueDat);
        dMat(ii,:,2) = single(charx2*opts.maxValueDat);

%         %% modification
%         curve = mean(datVec(ihw,:),1);
%         s0 = sqrt(mean((curve(2:end)-curve(1:end-1)).^2)/2);
%         
%         sz2D = zeros(max(it)-min(it)+1,1);
%         for t = min(it):max(it)
%             sz2D(t-min(it)+1) = sum(it==t);
%         end
% %         [~,tPeak] = max(sz2D);
%         curve0 = curve(min(it):max(it));
%         curve0(sz2D<0.5*max(sz2D)) = 0;
% 
%         [~,tPeak] = max(curve0);
%         tPeak = tPeak + min(it) - 1;
%         ts = tPeak;minV = curve(ts);
%         for t = tPeak:-1:min(it)
%             if(curve(t)<minV)
%                 minV = curve(t);
%                 ts = t;
%             else
%                 if(curve(t)-minV>=3*s0)
%                     break;
%                 end
%             end
%         end
% 
%         % end time
%         te = tPeak;minV = curve(te);
%         for t = tPeak:max(it)
%             if(curve(t)<minV)
%                 minV = curve(t);
%                 te = t;
%             else
%                 if(curve(t)-minV>=3*s0)
%                     break;
%                 end
%             end
%         end

        
        ftsLst.loc.t0(ii) = min(it);
        ftsLst.loc.t1(ii) = max(it);
        ftsLst.loc.x3D{ii} = pix0;
        ftsLst.loc.x2D{ii} = ihw;
        ftsLst.curve.rgt1(ii,:) = [min(rgT1),max(rgT1)];
        ftsLst.curve.dffMax(ii) = dffMax1;
        ftsLst.curve.dffMax2(ii) = dffMax2;
        ftsLst.curve.dffMaxFrame(ii) = (tMax+min(rgT)-1);
        ftsLst.curve.dffMaxZ(ii) = dffMaxZ;
        ftsLst.curve.dffMaxPval(ii) = dffMaxPval;
        ftsLst.curve.tBegin(ii) = min(it);%ts;
        ftsLst.curve.tEnd(ii) = max(it);%te;
        ftsLst.curve.duration(ii) = (max(it)-min(it)+1)*secondPerFrame;
        ftsLst.curve.rise19(ii) = rise19;
        ftsLst.curve.fall91(ii) = fall91;
        ftsLst.curve.width55(ii) = width55;
        ftsLst.curve.width11(ii) = width11;
        ftsLst.curve.riseTime(ii) = riseTime;
        ftsLst.curve.dff1Begin(ii) = (pp(1,1)+min(rgT1)-1);
        ftsLst.curve.dff1End(ii) = (pp(1,2)+min(rgT1)-1);
        ftsLst.curve.width11(ii) = width11;
        ftsLst.curve.decayTau(ii) = decayTau;
        
        %%
        % extend leftwards to calculate 
        t0 = min(rgT1) + pp(1,1) - 1;
        minV = dff1(t0);
        t = t0;
        while(t>max(1,min(it)-10))
            if(dff1(t)<minV)
                t0 = t;
                minV = dff1(t);
            else
                if(dff1(t)>minV+2*sigma1dff)
                   break;
                end
            end
            t = t-1;            
        end
        range = alignedrgT+t0;
        valid = range>0 & range<T;
        dffAlignedMat(ii,11+alignedrgT(valid)) = dff1(range(valid));
        
%         smodff1 = imgaussfilt(dff1,5);
%         d1dff1 = [0,smodff1(2:end) - smodff1(1:end-1)];
%         d2dff2 = [0,d1dff1(2:end) - d1dff1(1:end-1)];
%         
%         [~,maxd2] = max(d2dff2(min(rgT1):min(rgT1)+pp(3,1)-1));
%         searchRange = min(rgT1) -1 + [max(1,maxd2-1):min(T,maxd2+1)];
%         [~,minTP] = min(dff1(searchRange));
%         t0 = searchRange(minTP);
%         figure;
%         subplot(3,1,1);
%         plot(dff1);hold on;plot(t0:T,dff1(t0:T));
%         subplot(3,1,2);plot(d1dff1);
%         subplot(3,1,3);plot(d2dff2)
         
        % AUC
        datVec = reshape(dat,[],T);
        datCurve = mean(datVec(ihw,:),1);
        ftsLst.curve.datAUC(ii) = sum(datCurve(min(it):max(it)));
        ftsLst.curve.dffAUC(ii) = sum(dff1(min(it):max(it)));
        
        % basic features
        rgT = min(it):max(it);
        ih1 = ih-min(rgH)+1;
        iw1 = iw-min(rgW)+1;
        it1 = it-min(rgT)+1;
%         voxd = dat(rgH,rgW,rgT);
        voxi = zeros(length(rgH),length(rgW),length(rgT));
        pix1 = sub2ind(size(voxi),ih1,iw1,it1);
        voxi(pix1) = 1;
        ftsLst.basic = fea.getBasicFeatures(voxi,muPerPix,ii,ftsLst.basic);
        ftsLst.basic.center{ii} = [round(mean(ih)), round(mean(iw))];
        ftsLst.basic.spaBoundaryEvt(ii) = min(ih)<3 | max(ih)>H-2 | min(iw)<3 | max(iw)>W-2;
        ftsLst.basic.tempBoundaryEvt(ii) = min(it)<3 | max(it)>T-2;
        % p values
        %[p0,z0] = fea.getPval(voxd,voxi,1,0,0,4,4,sqrt(opts.varEst));
        %ftsLst.basic.p0(ii) = p0;
        %ftsLst.basic.z0(ii) = z0;
    end
    
    ftsLst.bds = img.getEventBorder(evtLst,[H,W,T]);
    
end

function F0 = getBaseline(x0,window,cut)
    datMA = movmean(x0,window);
    T = numel(datMA);
    step = cut;
    nSegment = ceil(T/step);

    F0 = zeros(size(x0));
    for k = 1:nSegment
        t0 = 1 + (k-1)*step;
        t1 = min(T,t0+cut);
        
        [curMinV,curMinT] = min(datMA(t0:t1));
        curMinT = curMinT + t0 - 1;
        if(k==1)
            F0(1:curMinT) = curMinV;
        else
            F0(preMinT:curMinT) = preMinV + (curMinV-preMinV)/(curMinT-preMinT)*[0:curMinT-preMinT]; 
        end      
        if(k==nSegment)
            F0(curMinT:end) = curMinV;
        end
        preMinT = curMinT;
        preMinV = curMinV;
    end

end









