function [ ss,ee,mapMatrix ] = buildGTWGraphMyVersion( ref, tst, s, t, smoBases, maxShift, stds, path0)
%buildGraph4Aosokin Build the dual graph of dynamic time warping problem
% For Aosokin's Matlab wrappers, like https://github.com/aosokin/graphCutMex_IBFS
%
% Consider graph constraint information and windowing
% Use coordinate to specify nodes in template pair, but use integer for spatial
% Support arbitrary graph
%
% ref: 1xT or NxT reference curve
% tst: NxT curves for pixels
% s,t: graph edge pair, undirected, if has pair (i,j), do not include (j,i), do not use (i,i) 
%
% smoBase: smoothness between edges
% winSize: window size, (winSize-1) lines above and below diagonal
% s2: noise variance for each test curve, or use single s2 for all test curves
%
% Basic constraint:
% No skipping/turn left/go down: otherwise has infinite cost (capacity)
% Start NEAR (1,1) and stop NEAR (T,T)
%
% The graphs has x axis as REF and y axis as TST
%
% ALL INPUTS SHOULD BE DOUBLE
%
% Xuelong Mi, CBIL@VT
% mixl18@vt.edu
% Dec.30, 2020

[nNode,T] = size(tst);
cap = 1e10;


if(~exist('path0','var') || isempty(path0))
    path0 = cell(nNode,1);
    for i = 1:nNode
       path0{i}  = repmat([1:T]',1,2);
    end
end

if(numel(smoBases)==1)
   smoBases = ones(nNode,1) *smoBases;
end

if(size(ref,1)==1)
    ref = repmat(ref,nNode,1);
end

%%  Full DTW dual graph template
% return the edges in dual graph, and start point of each edge in primal
% graph for mapping the weight
[TemplateEdge,source,sink] = getFullPairTemplate(T);
% set each DTW graph has one substitute source and sink.
% first, set each DTW graph and its dual graph
% still, use coordinates

% the matrix which record the mapping (Or the status whether the node is valid)
nNodeInDualDTW = (T-1)*(T-1)*2+2;
mapMatrix = zeros(nNode,nNodeInDualDTW);

% tic;
%% the dual node coordinate relative in primal graph
[x,y,z] = ind2sub([T-1,T-1,2],(1:(T-1)^2*2)');
y_location = y - 0.5 + (z-1.5)/2;
nEdgeInGraph = size(TemplateEdge,1);
ss = zeros(nNode*nNodeInDualDTW,2);
ee = zeros(nNode*nEdgeInGraph,4);
nValidEdges = 0;
for k = 1:nNode
   [lB,uB] = dilatePath(path0{k},maxShift);
   lBcurve = zeros(T-1,1);
   uBcurve = zeros(T-1,1);
   mapToSrcSink = 1:nNodeInDualDTW;
   %% set validDualGraph
   i = 1; j =1;
   for ii = 1:T-1
       % lower bound
       while(lB(i,2)<=ii)
           i = i + 1;
       end
       lBcurve(ii) = (lB(i-1,1) + lB(i,1))/2 - 1; % the midPoint of current interval
       % upper bound
       while(uB(j,2)<=ii)
           j = j+1;
       end
       uBcurve(ii) = (uB(j-1,1) + uB(j,1))/2 - 1; % the midPoint of current interval
   end
   mapToSrcSink(y_location<lBcurve(x)) = source;
   mapToSrcSink(y_location>uBcurve(x)) = sink;
   
   %% %%  Plan for faster construction
   % 1. for each DTW dual graph, construct the full graph, that is one template
   % 2. set non-valid dual nodes to source or sink
   % 3. remove the edges which both start and end are source or sink
   %% replace the non-valid nodes to source or sink.
   curEdge = TemplateEdge;
   curEdge(:,1) = mapToSrcSink(curEdge(:,1));
   curEdge(:,2) = mapToSrcSink(curEdge(:,2));
   %% remove non-valid edges (both points are source or sink)
   validEdges = curEdge(:,2)-curEdge(:,1)~=0;
   curEdge = curEdge(validEdges,:);
   %% weights
   d0 = (gtw.getDistMat(ref(k,:),tst(k,:))/stds(k))';
   weights = d0(curEdge(:,3));
%    drawDualGraph(curEdge,T,weights);
   %% Update
   % link substitute source and sink to primal source and sink
   ss(source + (k-1)*nNodeInDualDTW,:) = [cap,0];
   ss(sink + (k-1)*nNodeInDualDTW,:) = [0,cap];
   % update node label
   curEdge(:,1:2) = curEdge(:,1:2) + (k-1)*nNodeInDualDTW;
   % update ee
   nCurValidEdge = size(curEdge,1);
   ee(nValidEdges + 1:nCurValidEdge + nValidEdges,:) = [curEdge(:,1:2),weights,ones(numel(weights),1)*cap];
   nValidEdges = nValidEdges + nCurValidEdge;
   % update others
   mapMatrix(k,:) = mapToSrcSink;
end
ee = ee(1:nValidEdges,:);
% toc;
% tic;
if(sum(smoBases)>0)
    %% neighbor dual graph Linking.
    % for each neighbor pair, set bidirected edges with smoBase weight
    eeSpa = zeros(nNodeInDualDTW*numel(s),4);
    validEdges = false(nNodeInDualDTW*numel(s),1);
    for k = 1:numel(s)
        id1 = s(k);
        id2 = t(k);
        smoBase = max(smoBases(id1),smoBases(id2));
        curEE = [mapMatrix(id1,:)',mapMatrix(id2,:)',ones(nNodeInDualDTW,1)*smoBase,ones(nNodeInDualDTW,1)*smoBase];
        indexRange = (k-1)*nNodeInDualDTW+1 : k*nNodeInDualDTW;
        validEdges(indexRange) = curEE(:,1)<source | curEE(:,2)<source | curEE(:,1)~=curEE(:,2);
        curEE(:,1) = curEE(:,1) + (id1-1)*nNodeInDualDTW;
        curEE(:,2) = curEE(:,2) + (id2-1)*nNodeInDualDTW;
        eeSpa(indexRange,:) = curEE;
    end
    eeSpa = eeSpa(validEdges,:);
    ee = [ee;eeSpa];
end
% toc;
end