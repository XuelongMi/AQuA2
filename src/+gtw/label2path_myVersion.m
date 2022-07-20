function resPath = label2path_myVersion( labels, mapMatrix, T)
% label2path: Convert label of src and sink to path in primal graph
% For the graph representation of Aosokin's codes
% Src and sink edges do not use explicit node names, but other edges do

    % cut to mapping pattern in primal graph
    [nNode,nNodeInDual] = size(mapMatrix);
    resPath = cell(nNode,1);
    src = (T-1)^2*2+1;
    for k=1:nNode
        nodeInCurDual = (k-1)*nNodeInDual+[1:nNodeInDual-2];
        labelInCurDual = labels(nodeInCurDual);
        map2srcsink = mapMatrix(k,1:src-1);
        validNodes = find(map2srcsink<src);
        validLabels = labelInCurDual(validNodes);
        
%         drawNodesDualGraph([1:src-1]',labelInCurDual,T);
        
        validDualGraph = zeros(T-1,T-1,2,'uint8');
        validDualGraph(map2srcsink==src) = 1;
        validDualGraph(validNodes) = 1-validLabels;
        pathMap = sum(validDualGraph,3); % 0|2: horizontal 1: inclined
        i = 1;
        curPath = [];
        for ii = 1:T-1
            while(i<T && pathMap(ii,i)>0)
                i = i+1;
            end
            if(i==1 || pathMap(ii,i-1)==2)
                curPath = [curPath;i,ii;i,ii+1];
            else
                curPath = [curPath;i-1,ii;i,ii+1];
            end                
        end
        % add missing points
        curPath = [(1:curPath(1,1)-1)',ones(curPath(1,1)-1,1);curPath];
        curPath = [curPath;(curPath(end,1)+1:T)', T*ones(T-curPath(end,1),1)];
        % remove repeat point
        pathIdx = sub2ind([T,T],curPath(:,1),curPath(:,2));
        [pathY,pathX] = ind2sub([T,T],unique(pathIdx));
        curPath = [pathY,pathX];
        resPath{k} = curPath;
    end
end



