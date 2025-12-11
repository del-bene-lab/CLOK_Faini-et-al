%%% THIS Script show an exampel of data treatment to perfom the early and late born neuron clusteirng used in Fgiure 6a and Supp. Figure 20.

%IT  first load the x y position and fluo itensnityy of cell ROI extracted from
%%% the two dye images (obtained with 561 and 642 nm excitation laser )
%%% (contained in fluoANDpositionWHALOCAMP.mat file). Then matches the
%%% coordinate in the two channels and define the dye ration of each resulting neurons. 
%%% Then, based on k-means or GMM clustering, accoridgn to ratio, the full
%%% popualtion is divided int otwo group , earry and late born neurons 


load fluoANDpositionWHALOCAMP

%% -------- PARAMETERS --------
pair_thr  = 6;        % max distance (px) to consider ROIs a match
eps_denom = 1e-12;    % avoid divide-by-zero
METRIC    = 'ratio';  % 'ratio' or 'ratioVector'

% Required inputs (must already be in workspace):
% x561, y561, x642, y642

% -------- 1) Pair ROIs by nearest neighbor --------
n1 = numel(x561);
n2 = numel(x642);

% Distance matrix between ROIs in the two channels
DX   = x561(:) - x642(:)';  
DY   = y561(:) - y642(:)';  
Dist = hypot(DX, DY);       % n1 x n2

% Greedy matching: smallest distance first, no duplicates
pairs = [];
used1 = false(n1,1);
used2 = false(n2,1);

while true
    DistMasked = Dist;
    DistMasked(used1, :) = inf;
    DistMasked(:, used2) = inf;

    [mn, idx] = min(DistMasked(:));
    if ~isfinite(mn) || mn > pair_thr
        break;
    end

    [i1, i2] = ind2sub([n1 n2], idx);
    pairs(end+1,:) = [i1 i2]; 
    used1(i1) = true;
    used2(i2) = true;
end

if isempty(pairs)
    warning('No ROI pairs found within threshold %.1f px.', pair_thr);
end

% Indices of matched ROIs
i561 = pairs(:,1);
i642 = pairs(:,2);

%% -------- 3) Compute metric & normalization --------
switch lower(METRIC)
    case 'ratio'
        % Xraw in [0,1]: 561 / (561 + 642)
        Xraw   = fluo561 ./ max(fluo561 + fluo642, eps_denom);
        Xnorm  = min(max(Xraw,0),1);     % for colormap
        cbLabel = 'RATIO = 561 / (561 + 642)';

    case 'ratiovector'
        % Xraw in [0, pi/2]: atan(561/642)
        Xraw   = atan2(fluo561, max(fluo642, eps_denom));
        Xnorm  = Xraw / (pi/2);          % normalize to [0,1]
        cbLabel = 'atan(561/642), normalized';

    otherwise
        error('Unknown METRIC "%s". Use ''ratio'' or ''ratioVector''.', METRIC);
end

% -------- Purple → cyan colormap --------
nC  = 256;
c1  = [0.5 0.00 0.5];  % purple
c2  = [0.0 1.00 1.0];  % cyan
cmap = [linspace(c1(1),c2(1),nC)', ...
        linspace(c1(2),c2(2),nC)', ...
        linspace(c1(3),c2(3),nC)'];

ci = max(1, min(nC, round(1 + Xnorm*(nC-1))));  % index into colormap

% --------  Plot: colored ROIs over 561 image --------
figure('Name','ROI Metric','Position',[100 100 700 600]);

%imshow(img561, [], 'InitialMagnification','fit'); 
%hold on;

scatter(x561(i561), y561(i561), 45, cmap(ci,:), ...
    'filled', 'MarkerEdgeColor','k', 'LineWidth',0.5);

title(sprintf('ROIs colored by %s', METRIC));
colormap(gca, cmap);
caxis([0 1]);                     % normalized color scale
cb = colorbar;
cb.Label.String = cbLabel;

hold off;


%% 
% Clustering on METRIC (Xraw): K-means or GMM + colored ROI map + per-cluster rasters
% ---------------------------

% --- Choose clustering method and K
CLUST_METHOD = 'kmeans';    % 'kmeans' or 'gmm'
K_user = [];                % [] → auto-select (silhouette for kmeans, BIC for GMM)
Kmin = 2; 
Kmax = 4;                   % search range for auto

X = Xraw(:);                % metric vector for matched ROIs
N = numel(X);

% --- Compute clusters
switch lower(CLUST_METHOD)
    case 'kmeans'
        if isempty(K_user)
            bestK = []; 
            bestSil = -inf; 
            bestIdx = [];
            opts = statset('MaxIter', 1000, 'Display', 'off');
            for K = Kmin:min(Kmax, N)
                try
                    idx0 = kmeans(X, K, 'Replicates', 10, 'Options', opts, 'Distance', 'sqeuclidean');
                    s = silhouette(X, idx0);
                    sil = mean(s, 'omitnan');
                    if sil > bestSil
                        bestSil = sil; 
                        bestK = K; 
                        bestIdx = idx0;
                    end
                catch
                end
            end
            if isempty(bestIdx)
                bestK = min(max(2, Kmin), N);
                if bestK <= 1
                    warning('Not enough data to cluster.'); 
                    return;
                end
                bestIdx = kmeans(X, bestK, 'Replicates', 5, ...
                    'Options', statset('MaxIter',500,'Display','off'));
            end
        else
            bestK = K_user;
            bestIdx = kmeans(X, bestK, 'Replicates', 10, ...
                'Options', statset('MaxIter',1000,'Display','off'));
        end

        % order clusters by ascending cluster mean
        centers = accumarray(bestIdx, X, [bestK 1], @mean, NaN);
        [centers, ord] = sort(centers, 'ascend');
        map = zeros(size(bestIdx));
        for k = 1:bestK
            map(bestIdx==ord(k)) = k;
        end
        clustIdx = map;
        K = bestK;

    case 'gmm'
        opts = statset('MaxIter', 2000, 'Display', 'off');
        if isempty(K_user)
            gmBest = []; 
            bestBIC = inf; 
            bestK = NaN;
            for K = Kmin:min(Kmax, max(2,N))
                try
                    gm = fitgmdist(X, K, 'RegularizationValue', 1e-6, ...
                        'Options', opts, 'Replicates', 10);
                    if gm.BIC < bestBIC
                        bestBIC = gm.BIC; 
                        gmBest = gm; 
                        bestK = K;
                    end
                catch
                end
            end
            if isempty(gmBest)
                warning('GMM failed to fit.'); 
                return;
            end
            [idx, ~, ~] = cluster(gmBest, X);
            [muSorted, ord] = sort(gmBest.mu, 'ascend');
            map = zeros(size(idx));
            for k = 1:bestK
                map(idx==ord(k)) = k;
            end
            clustIdx = map;
            K = bestK;
        else
            gm = fitgmdist(X, K_user, 'RegularizationValue', 1e-6, ...
                'Options', opts, 'Replicates', 10);
            [idx, ~, ~] = cluster(gm, X);
            [muSorted, ord] = sort(gm.mu, 'ascend');
            map = zeros(size(idx));
            for k = 1:K_user
                map(idx==ord(k)) = k;
            end
            clustIdx = map;
            K = K_user;
        end

    otherwise
        error('Unknown CLUST_METHOD: use ''kmeans'' or ''gmm''.');
end

% --- Colors for up to 6 clusters
baseC = lines(max(6,K)); 
C = baseC(1:K,:);

% --- Summary print
fprintf('Clustering (%s): K = %d\n', lower(CLUST_METHOD), K);
for k = 1:K
    fprintf('  Cluster %d: n = %d, mean metric = %.3f\n', ...
        k, sum(clustIdx==k), mean(X(clustIdx==k), 'omitnan'));
end
% ---------------------------
% Simplified clustering summary (2 clusters only)
% ---------------------------

% Force to 2 clusters maximum
K = min(2, numel(unique(clustIdx)));
X = Xraw(:);

% limit to first two clusters (ignore extra)
validIdx = clustIdx <= K;
clustIdx = clustIdx(validIdx);
X = X(validIdx);
i561 = i561(validIdx);
i642 = i642(validIdx);

% Define colors specifically for the FIRST PANEL:
% Cluster 1 → magenta, Cluster 2 → cyan
C_firstPanel = [1 0 1;   % magenta
                0 0.7 0.7];  % cyan

% Also keep generic colors for other plots if you want (optional)
C_generic = lines(max(K,2));

% Define x-axis for metric
switch lower(METRIC)
    case 'ratio'
        xLabel = 'RATIO = 561 / (561 + 642)';
        xMin = 0; xMax = 1; xs = linspace(xMin, xMax, 300)';
    case 'ratiovector'
        xLabel = 'RATIOVECTOR = atan(561/642) [rad]';
        xMin = 0; xMax = pi/2; xs = linspace(xMin, xMax, 300)';
    otherwise
        xLabel = 'Metric value';
        xMin = min(X); xMax = max(X); xs = linspace(xMin, xMax, 300)';
end

% --- Gaussian fits (approximate)
mus = zeros(1,K); sigs = zeros(1,K); weights = zeros(1,K);
for k = 1:K
    Xk = X(clustIdx==k);
    mus(k) = mean(Xk,'omitnan');
    sigs(k) = std(Xk,[],'omitnan');
    if ~isfinite(sigs(k)) || sigs(k)==0, sigs(k) = (xMax-xMin)/100; end
    weights(k) = numel(Xk)/numel(X);
end

% --- Figure layout
figure;

% (1) ROI MAP — 561 image (USE MAGENTA/CYAN HERE)
subplot(2,3,[1 2]);
imshow(img561, [], 'InitialMagnification','fit'); hold on;
for k = 1:K
    ii = i561(clustIdx==k);
    % choose magenta for cluster 1, cyan for cluster 2 (fallback to gray if K=1)
    col = (K>=2) * C_firstPanel(k,:) + (K==1) * [0.5 0.5 0.5];
    scatter(x561(ii), y561(ii), 25, 'o', 'MarkerEdgeColor', col, ...
        'LineWidth', 1.5, 'DisplayName', sprintf('Cluster %d', k));
end
title(sprintf('ROI map (561 image) — %s clustering', upper(CLUST_METHOD)));
legend('Location','southeastoutside'); legend boxoff;
colormap(gray); hold off;

% (2) METRIC DISTRIBUTION + GAUSSIANS (keep generic colors)
subplot(2,3,3); hold on;
nb = max(10, min(40, round(numel(X)/5)));
edges = linspace(xMin, xMax, nb+1);
for k = 1:K
    histogram(X(clustIdx==k), edges, 'Normalization','pdf', ...
        'FaceColor', C_generic(k,:), 'EdgeColor','none', 'FaceAlpha', 0.4);
    plot(xs, weights(k)*normpdf(xs, mus(k), sigs(k)), '-', ...
        'Color', C_generic(k,:), 'LineWidth', 2);
end
histogram(X, edges, 'Normalization','pdf', 'EdgeColor','k', 'FaceColor','none', 'LineWidth',1);
xlabel(xLabel); ylabel('PDF'); grid on;
title(sprintf('%s distribution + Gaussian fits (K=%d)', METRIC, K));
xlim([xMin xMax]); hold off;
