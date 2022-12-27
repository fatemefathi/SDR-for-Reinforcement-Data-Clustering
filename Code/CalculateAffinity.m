function [affinity] = CalculateAffinity(data)
% set the parameters
sigma = 1;
for i=1:size(data,1)    
    for j=1:size(data,1)
        dist = sqrt((data(i,1) - data(j,1))^2 + (data(i,2) - data(j,2))^2); 
        affinity(i,j) = exp(-dist/(2*sigma^2));
    end
end

% sigma=1;
% knn    = ceil(0.03*size(data,1));
% m=size(data,1);
% %for knn=1:5
% dist= squareform(pdist(data));
% [srtdDist,srtdIdx] = sort(dist,'ascend');
% dist = srtdDist(1:knn+1,:);
% nidx= srtdIdx(1:knn+1,:);
% tempW  = exp(-dist.^2/sigma); 
% i = repmat(1:m,knn+1,1);
% affinity = sparse(i(:),double(nidx(:)),tempW(:),m,m); 
% affinity= max(affinity,affinity'); % for undirected graph.
% affinity=full(affinity);