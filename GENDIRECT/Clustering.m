function Dat = Clustering(Dat)

% Check if: fevals >= 100n + 1 and (local search was not performed or
% DIRECT sufficiently improved Fmin)
if Dat.VAL.I >= Dat.VAL.GLC && (Dat.VAL.GLCcount == 0 || (abs(Dat.VAL.Fmin - Dat.VAL.Fmin_old) > 0.01*abs(Dat.VAL.Fmin)))

    % Use hierarchical clustering using standart Ward's method to cluster
    % sampled points
    Dat.VAL.GLCcount = Dat.VAL.GLCcount + 1;
    CC = Dat.MSS.CC(:, 1:Dat.VAL.I)';
    X = pdist2( CC, CC, "euclidean" );
    Z = linkage( X , 'ward' );
    maxClusters = 10*Dat.Problem.n;
    silhouetteValues = zeros(maxClusters, 1);
    for k = 2:maxClusters
        clusters = cluster(Z, 'maxclust', k);
        silhouetteValues(k) = mean(silhouette(X, clusters));
    end

    % Find the optimal number of clusters as the one with the highest silhouette value
    [~, optimalK] = max(silhouetteValues);
    
    % Find starting points from each cluster
    Clusters = cluster( Z, 'maxclust', optimalK );
    POHLocal = zeros(1, optimalK);
    for i = 1:optimalK
        indexs = find(Clusters == i);
        POHLocal(i) = indexs(find(Dat.MSS.FF(indexs) == min(Dat.MSS.FF(indexs)), 1));
    end
    
    % Run local search from each cluster
    [xminlocs,fminlocs,~,outputs]=arrayfun(@(idx)fmincon(@(x)Dat.VAL.TPD(x),...
    Dat.MSS.CC(:,Dat.POHLocal(idx)),[],[],[],[],zeros(Dat.Problem.n,1), ...
    ones(Dat.Problem.n,1),[],Dat.Hybridization.options),1:length(POHLocal),'UniformOutput',false);

    % Extract data from local solver
    Dat.VAL.nLocSearch = Dat.VAL.nLocSearch + length(POHLocal);
    Dat.VAL.fevalLocal = Dat.VAL.fevalLocal + sum(cellfun(@(out) out.funcCount, outputs));

    %  Update Fmin if it has improved
    [mins, idxs] = min(cell2mat(fminlocs));
    if mins < Dat.VAL.Fmin
        Dat.VAL.Fmin = mins;
        Dat.VAL.Xmin = (xminlocs{idxs} - Dat.Problem.xl)./Dat.VAL.delta;
    end
end

end
