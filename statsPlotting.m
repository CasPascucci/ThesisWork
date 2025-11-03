function statsPlotting(Results)
ok   = [Results.exit_ok]' == 1;
gamma= [Results.gamma]';
kr   = [Results.kr]';
tgo  = [Results.tgo]';        % seconds
fuel = [Results.fuel_opt]';
coeff1   = [Results.coeff1]';
coeff2   = [Results.coeff2]';
coeff3   = [Results.coeff3]';
coeff4   = [Results.coeff4]';

nOK = sum(ok); N = numel(Results);
fprintf('Success: %d / %d  (%.1f%%)\n', nOK, N, 100*nOK/max(N,1));

vars = {'kr','tgo','fuel_opt','coeff1','coeff2','coeff3','coeff4'};
Xall = [kr,tgo,fuel,coeff1,coeff2,coeff3,coeff4];
%% Histograms
Xs = Xall(ok,:);
figure('Name','Histograms (success only)');
tiledlayout(2,4,'Padding','compact','TileSpacing','compact');
for i = 1:numel(vars)
    nexttile; histogram(Xs(:,i), 50); grid on;
    xlabel(vars{i}); ylabel('count'); title(['Hist ' vars{i}]);
end
%% Boxplots
figure('Name','Boxplots by exit_ok');
tiledlayout(3,3,'Padding','compact','TileSpacing','compact');

vars = {'kr','tgo','fuel_opt','coeff1','coeff2','coeff3','coeff4'};
okcat  = categorical([Results.exit_ok]', [false true], {'fail','success'});

for i = 1:numel(vars)
    v = [Results.(vars{i})]';
    nexttile;
    boxplot(v, okcat);                 % <-- key change
    grid on; title(['Box ' vars{i}]); ylabel(vars{i});
end
end