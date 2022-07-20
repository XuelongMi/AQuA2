function [spLst,cx,dlyMap,distMat,rgt00,isFail,evtMemC,evtMemCMap] = spgtw_MuYuProject(...
    dF,seMap,seSel,smoBase,maxStp,cDelay,spSz,spT,superVoxels,opts)
% spgtw super pixel GTW 
% make one burst to super pixels and run gtw

if ~isfield(opts,'gtwGapSeedMin') || ~isfield(opts,'gtwGapSeedRatio')
    opts.gtwGapSeedRatio = 4;
    opts.gtwGapSeedMin = 5;
end

%% setting
[H,W,T] = size(dF);
isFail = 0;
maxStp = max(min(maxStp,round(T/2)),1);
s00 = 1; % already normalized

%% super pixels
m0Msk = seMap==seSel;
validMap = sum(m0Msk,3)>0;

dFip = dF;
dFip(seMap~=seSel) = nan;
dFip = gtw.imputeMov_Fast(dFip,validMap);

dFAvg = max(dFip,[],3);
dFAvg(dFAvg<0) = 0;
dFAvg = medfilt2(dFAvg);  

%% super pixel
spSzMinLimitation = spSz/4;
nSpEstimate = max(1,round(H*W/spSz));
[L] = superpixels(dFAvg,nSpEstimate);
L(dFAvg==0) = 0;
spLst = label2idx(L);

%% get super pixel
if sum(validMap(:))< spSzMinLimitation || numel(spLst)<2
    spLst = {find(validMap>0)};
    cx = [];
    dlyMap = [];
    distMat = [];
    evtMemC = [];
    evtMemCMap = [];
    rgt00 = 1:T;
    isFail = 1;
    return
end





%% connectivity
nSp = numel(spLst);
nAdd = nSp + 1;
for i = 1:nSp
    pix = spLst{i};
   if(~isempty(pix)) 
       cMap = false(H,W);
       cMap(pix) = true;
       cc = bwconncomp(cMap);
       if(cc.NumObjects>1)
           comps = cc.PixelIdxList;
           for j = 1:numel(comps)
               if(j==1)
                   spLst{i} = comps{j};
               else
                   spLst{nAdd} = comps{j};
                   nAdd = nAdd+1;
               end
           end
       end
   end
end

dh = [-1,-1,-1,0,0,1,1,1];
dw = [-1,0,1,-1,1,-1,0,1];
isolatedRegion = [];
% merge small pixels
while(1)
    for i = 1:numel(spLst)
        pix = spLst{i};
        curSz = numel(spLst{i});
        if(curSz>0 && curSz<spSzMinLimitation)
            [ih,iw] = ind2sub([H,W],pix);
           for k = 1:numel(dw)
               ih1 = min(H,max(1,ih+dh(k)));
               iw1 = min(W,max(1,iw+dw(k)));
               pixShift = sub2ind([H,W],ih1,iw1);
               neiL = setdiff(L(pixShift),[0,i]);
               if(~isempty(neiL))
                  break; 
               end
           end
           if(isempty(neiL))
              isolatedRegion = [isolatedRegion,i];
              continue;
           end
           neiL = neiL(1);
           % update
           spLst{neiL} = [spLst{neiL};pix];
           L(pix) = neiL;
           spLst{i} = [];
        end
    end
    sz = cellfun(@numel,spLst);
    smallRegions = find(sz>0 & sz<spSz/4);
    smallRegions = setdiff(smallRegions,isolatedRegion);
    if(isempty(smallRegions))
       break; 
    end
end
sz = cellfun(@numel,spLst);
spLst = spLst(sz>0);
[~,seedIn] = max(dFAvg);

% signal part
idx0 = find(m0Msk>0);
[~,~,it0] = ind2sub(size(m0Msk),idx0);
rgt00 = max(min(it0)-5,1):min(max(it0)+5,T);
dat = double(dFip(:,:,rgt00));

fprintf('Node %d\n',numel(spLst));

%% alignment
% graph
[ih0,iw0] = find(validMap>0);
gapSeed = max(ceil(max(max(ih0)-min(ih0),max(iw0)-min(iw0))/opts.gtwGapSeedRatio),opts.gtwGapSeedMin);
[ref,tst,refBase,s,t,idxGood] = gtw.sp2graph(dat,validMap,spLst,seedIn,gapSeed);

% gtw
spLst = spLst(idxGood);
nSp = numel(spLst);
sz = cellfun(@numel,spLst);
s2 = 1./sqrt(sz);
for i = 1:nSp
    spSeedVec(i) = spLst{i}(1);
end
if numel(spLst)>3 && numel(refBase)>5
%     tic
    [ ss,ee,gInfo ] = gtw.buildGTWGraph( ref, tst, s, t, smoBase, maxStp, s2);
    [~, labels1] = aoIBFS.graphCutMex(ss,ee);
    path0 = gtw.label2path4Aosokin( labels1, ee, ss, gInfo );
    t00 = toc;
    if numel(spLst)>1000
        fprintf('Time %fs\n',t00);
    end
else
    [nPix,nTps] = size(tst);
    path0 = cell(nPix,1);
    rg = (0:nTps)';
    p0 = [rg,rg,rg+1,rg+1];
    for ii=1:nPix
        path0{ii} = p0;
    end
end

% warped curves
pathCell = cell(H,W);
vMap1 = zeros(H,W);
vMap1(spSeedVec) = 1:numel(spSeedVec);
for ii=1:numel(spLst)
    [ih0,iw0] = ind2sub([H,W],spSeedVec(ii));
    pathCell{ih0,iw0} = path0{ii};
end
datWarp = gtw.warpRef2Tst(pathCell,refBase/max(refBase(:)),vMap1,[H,W,numel(refBase)]);
dVec = reshape(datWarp,[],numel(refBase));
cx = dVec(spSeedVec,:);

% time to achieve different levels for each seed
nSp = numel(spLst);
thrVec = 0.5:0.05:0.95;
tAch = nan(nSp,numel(thrVec));
for nn=1:nSp
    x = cx(nn,:);
    [~,t0] = max(x);
    x = x(1:t0);
    for ii=1:numel(thrVec)
        t1 = find(x>=thrVec(ii),1);
        if isempty(t1)
            t1 = t0;
        end
        tAch(nn,ii) = t1;
    end
end
tDly = mean(tAch,2);

% direction for each pair
nPair = numel(s);
distMat = nan(nSp,nSp);
% cDelay = 1e8;  % !!
for nn=1:nPair
    s0 = s(nn);
    t0 = t(nn);
    d0 = tAch(s0,:)-tAch(t0,:);  % negative is earlier
    d0 = sum(d0)/numel(thrVec);
    %d0a = abs(d0)-cDelay;
    %d0a(d0a<0) = 0;
    %if sum(d0a)==0
    distMat(s0,t0) = d0;
    distMat(t0,s0) = -d0;
    %end
end

% delay map
dlyMap = inf(H,W);
for nn=1:numel(spLst)
    dlyMap(spLst{nn}) = tDly(nn);
end

% partition by continuity
A = Inf(nSp,nSp);
[ia,ib] = find(~isnan(distMat));
for ii=1:numel(ia)
    ia0 = ia(ii);
    ib0 = ib(ii);
    A(ia0,ib0) = min(abs(distMat(ia0,ib0)),A(ia0,ib0));
end
A(A>cDelay) = Inf;

B = A;
B(A<inf) = 1;
B(isinf(A)) = 0;
B(eye(nSp)>0) = 1;
B = max(B,B');

G = graph(B);
evtMemC = zeros(nSp,1);
cc = conncomp(G,'OutputForm','cell');
for ii=1:numel(cc)
    cc0 = cc{ii};
    evtMemC(cc0) = ii;
end

evtMemCMap = zeros(H,W);
for ii=1:numel(evtMemC)
    evtMemCMap(spLst{ii}) = evtMemC(ii);
end

end







