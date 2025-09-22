% Batch run test script
clear; clc; close all;

%gamma_list = [logspace(log10(0.001), log10(1), 8)]; % n = 10
%kr_list = linspace(0, 20, 20);
%tgo_list = linspace(5, 20, 10);
gamma_list = 1.0;
kr_list = 6;
tgo_list = 9.70249498737309;

[Xg, Xk, Xt] = ndgrid(gamma_list, kr_list, tgo_list); % Grid for all combinations of IC's
X0_all = [Xg(:), Xk(:), Xt(:)];

valid = X0_all(:,2) >= 2*(X0_all(:,1)+2); % Use only valid gamma/kr combos
X0_all = X0_all(valid,:);

numTests = size(X0_all,1);
Results(numTests,1) = struct('ok', false, 'err', "", 'S', []);

cfg = struct(); % Leave blank if no special config of problem values

% Parallel Setup
%Parallel setup dropped runtime from 477 to 78 seconds, defintiely use this
%later during montecarlo sims
if contains(struct2array(ver), 'Parallel Computing Toolbox') && numTests >= 100
    usingParallel = true;
    if isempty(gcp("nocreate"))
        parpool %Create parallel pool if none
    end
else
    usingParallel = false; % Default to not using parallel processing
end
%usingParallel = false;
% Batch Run
tStart = tic;
if usingParallel
    parfor i = 1:numTests
        x0 = X0_all(i,:)';
        try
            S = runLEMMassOptOpenLoop(x0);
            %S = runLEMMassOptClosedLoop(x0);
            Results(i).ok = true;
            Results(i).S = S;
        catch ME
            Results(i).ok = false;
            Results(i).err = string(ME.message);
        end
    end
else
    for i = 1:numTests
        x0 = X0_all(i,:)';
        try
            S = runLEMMassOptOpenLoop(x0);
            %S = runLEMMassOptClosedLoop(x0);
            Results(i).ok = true;
            Results(i).S = S;
        catch ME
            Results(i).ok = false;
            Results(i).err = string(ME.message);
        end
    end
end
elapsed = toc(tStart);
fprintf("Finished %d sims in %.1f s\n", numTests, elapsed);

% Build a summary table
idx_ok = find([Results.ok]);
n_ok = numel(idx_ok);

fprintf("%d tests were OK out of the total %d tests\n", n_ok, numTests);

%Pre allocate
gamma   = nan(numTests,1);
kr      = nan(numTests,1);
tgo     = nan(numTests,1);
tgoVirt = nan(numTests,1);
cost    = nan(numTests,1);
prop_kg = nan(numTests,1);
errMsg  = strings(numTests,1);

for i = 1:numTests
    if Results(i).ok %if ok, but so far all tests have been succesful, if it continues to do so then this if else can be removed
        S = Results(i).S;
        gamma(i)   = S.opt.gamma;
        kr(i)      = S.opt.kr;
        tgo(i)     = S.opt.tgo;
        tgoVirt(i) = S.opt.tgoVirtual;
        cost(i)    = S.opt.costEval;

        % propellant used in kg
        M_ref = S.refs.M_ref;            % equals massInitDim in current setup
        massInitDim = S.masses.massInitDim;
        mf_kg = S.stateTraj(end,7) * M_ref;
        prop_kg(i) = massInitDim - mf_kg;
    else
        errMsg(i) = Results(i).err;
    end
end

Summary = table(gamma, kr, tgo, tgoVirt, cost, prop_kg, errMsg); %Table helps for easy mass inspection
Summary.x0_gamma = X0_all(:,1);
Summary.x0_kr    = X0_all(:,2);
Summary.x0_tgo   = X0_all(:,3);

% Keep only successful rows for ranking, probably not necessary
SummaryOK = Summary( [Results.ok].', : );
SummaryOK = sortrows(SummaryOK, 'prop_kg');   % or 'cost'

% 7) Save everything
% ts = datestr(now, 'yyyymmdd_HHMMSS');
% md = matlabdrive;
% folder = "/Thesis Research/Batch Runs";
% fname = ['LEM_batch_' ts '.mat'];
% fullPath = fullfile(md, folder, fname);
% save(fullPath, 'Summary', 'SummaryOK', 'Results', 'X0_all', 'cfg', 'elapsed');
% fprintf('Saved batch to %s\n', fname);

%% 8) Plotting
% Run this section after loading desired results

% Find index of min propellant in the OK table
[~, minRow] = min(SummaryOK.prop_kg);
best = SummaryOK(minRow,:);

% Find this row in the original Summary, same idx as Results
% Probably unnecessary step if all tests are OK
bestIdx = find(Summary.gamma == best.gamma & ...
               Summary.kr    == best.kr & ...
               Summary.tgo   == best.tgo, 1, 'first');

if ~isempty(bestIdx) && Results(bestIdx).ok % Found a result and it actually ran
    plotLEMMassOpt(Results(bestIdx).S);
    fprintf('Plotted best run: gamma=%.4f, kr=%.4f, tgo=%.4f, prop=%.1f kg\n', ...
        best.gamma, best.kr, best.tgo, best.prop_kg);
else
    warning('Best run not found in Results.');
end


[minGamma, maxGamma] = bounds(SummaryOK.gamma);
rangeGamma = maxGamma - minGamma;
[minkr, maxkr] = bounds(SummaryOK.kr);
rangekr = maxkr - minkr;

% Figure 8: Gamma vs Kr spread for batch runs
figure(); hold on;
scatter(SummaryOK.gamma, SummaryOK.kr,10)
xscale('log');
title('Gamma vs Kr');
subtitle(sprintf('Gamma Spread: %.3e \n Kr Spread: %.3e', rangeGamma, rangekr));
xlabel('Gamma');
ylabel('Kr');

%Figure 9: Case Number vs Dimensional Tgo
figure();
plot(1:numTests, Summary.tgo*S.refs.T_ref,'.','MarkerSize',30);
title('Case Number vs Tgo (Dimensional)');
xlabel('Case Number');
ylabel('Tgo (seconds)');
grid;


% waitforbuttonpress;
% for i = 1:7
%     fig = figure(i);
%     waitforbuttonpress;
%     savefig(sprintf('C:/Users/casey/MATLAB Drive/Thesis Research/Media/8_15_25_Figures/Gamma 0.5/Figure_%d',i));
% end