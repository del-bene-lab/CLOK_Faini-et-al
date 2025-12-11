% this is a exemplary script to plot dye ratio distribution and dispotion in the ioptic tectum 
% (graph of Supp.Figure )  (download alos the file "Exemple_fish1_ratio_position"


%% read values of dye ratio and cell position from csv file 

filename='Exemple_fish1_ratio_position.csv'
currentFolder = pwd;
fullpath = fullfile(currentFolder, filename);

% Create import options based on the file
opts = detectImportOptions(fullpath);

% Make sure we use the header names and select only the needed columns
opts.SelectedVariableNames = {'X', 'Y', 'ratio'};

% Read the table
dataTbl = readtable(fullpath, opts);

% Extract columns into separate variables
X     = dataTbl.X;
Y     = dataTbl.Y;
ratio = dataTbl.ratio;

%%  plot scatter dot plot and circultr fit to extract angualr position of each cell 
[xc, yc, ~] = fitCircleXY([X,Y], true, ratio);


%% Plot distribution of dye ratio values 
figure;

    [f, xi] = ksdensity(ratio);
    plot(xi, f, 'LineWidth', 2, 'Color', colors(i,:));

xlabel('Dye Ratio');
ylabel('Probability Density');
title('Distribution of Ratios');
grid on;
hold off;



function [xc, yc, R, rmse] = fitCircleXY(XY, showPlot, ratio, axisLim, savePdf)
%FITCIRCLEXY  Best-fit circle to 2D points, with optional color-coded visualization.

    % --- Default arguments ---
    if nargin < 2 || isempty(showPlot)
        showPlot = false;
    end
    if nargin < 3
        ratio = [];
    end
    if nargin < 4
        axisLim = [];
    end
    if nargin < 5
        savePdf = [];
    end

    % --- Input checks ---
    if size(XY,2) ~= 2 || size(XY,1) < 3
        error('XY must be N×2 with N ≥ 3.');
    end
    if ~isempty(ratio) && numel(ratio) ~= size(XY,1)
        error('ratio must have same number of elements as rows in XY.');
    end

    x = XY(:,1);
    y = XY(:,2);

    % ---------- (1) Algebraic initial estimate ----------
    xmean = mean(x);
    ymean = mean(y);
    u = x - xmean;
    v = y - ymean;

    z = u.^2 + v.^2;
    D = [u, v, ones(size(u))];
    abc = D \ z;
    A = abc(1); B = abc(2); C = abc(3);

    uc0 = A/2; vc0 = B/2;
    R0  = sqrt(max(0, C + uc0^2 + vc0^2));

    xc0 = uc0 + xmean;
    yc0 = vc0 + ymean;

    % ---------- (2) Nonlinear refinement ----------
    obj = @(p) sum((hypot(x - p(1), y - p(2)) - ...
                   mean(hypot(x - p(1), y - p(2)))).^2);
    opts = optimset('Display','off','TolX',1e-10,'TolFun',1e-10);
    pOpt = fminsearch(obj, [xc0, yc0], opts);
    xc = pOpt(1); yc = pOpt(2);

    r_i = hypot(x - xc, y - yc);
    R   = mean(r_i);
    rmse = sqrt(mean((r_i - R).^2));

    % ---------- (3) Optional visualization ----------
    if showPlot
        figure; hold on; axis equal; box on
        title(sprintf('Best-Fit Circle: center=(%.3f, %.3f), R=%.3f', xc, yc, R))
        xlabel('X'); ylabel('Y');

        if ~isempty(ratio)
            % --- Continuous colormap (magenta → cyan) ---
            N = 256;
            cmap = [linspace(1,0,N)', linspace(0,1,N)', ones(N,1)]; 
            colormap(cmap);

            % Semi-transparent dots
            scatter(x, y, 60, ratio, 'filled', 'MarkerFaceAlpha', 0.5);

            cb = colorbar;
            ylabel(cb, 'Ratio');
        else
            plot(x, y, 'ko', 'MarkerFaceColor',[0.2 0.6 0.9], 'MarkerFaceAlpha', 0.5);
        end

        % Plot fitted circle
        theta = linspace(0, 2*pi, 500);
        x_fit = xc + R * cos(theta);
        y_fit = yc + R * sin(theta);
        plot(x_fit, y_fit, 'r-', 'LineWidth', 1.5);

        % Plot center
        plot(xc, yc, 'rx', 'MarkerSize', 10, 'LineWidth', 2);

        % ------- Axis limits control -------
        if ~isempty(axisLim)
            if isscalar(axisLim)
                L = axisLim;
                xlim([xc-L, xc+L]);
                ylim([yc-L, yc+L]);
            elseif numel(axisLim)==4
                axis(axisLim(:)');
            else
                warning('axisLim must be scalar or [xmin xmax ymin ymax]. Ignored.');
            end
            axis equal;
        end

    end
end
