clear all; clc;
% All values are dimensional
PDINom = struct;
PDINom.altitude_km        = 13.36;
PDINom.lonInitDeg         = 41.85;
PDINom.latInitDeg         = -71.59;
PDINom.inertialVelocity   = 1693.8;
PDINom.flightPathAngleDeg = 0;
PDINom.azimuth            = pi;

planetaryParams = struct;
planetaryParams.rPlanet = 1736.01428 * 1000; % m
planetaryParams.gPlanet = 1.622; % m/s^2
planetaryParams.gEarth = 9.81; % m/s^2

vehicleNom = struct;
vehicleNom.massInit = 15103.0; % kg
vehicleNom.dryMass = vehicleNom.massInit - 8248; % kg
vehicleNom.isp = 311; % s
vehicleNom.maxThrust = 45000; % N
vehicleNom.minThrust = 4500; % N

targetState = struct;
targetState.landingLonDeg      = 41.85;
targetState.landingLatDeg      = -90.0;
targetState.rfLanding = [0;0;0];
targetState.vfLanding = [0;0;-1];
targetState.afLanding = [0;0;2*planetaryParams.gPlanet];
%targetState.afLanding = [0;0;vehicleParams.maxThrust/(6710)]; % This
                                                               %afStar seems to be too high for a solvable path
targetState.delta_t   = 5; % seconds dim, for btt

optimizationParams = struct;
optimizationParams.nodeCount = 997; %Count must be odd for Simpson
optimizationParams.glideSlopeFinalTheta = 45; %deg
optimizationParams.glideSlopeEnabled = true;
optimizationParams.pointingEnabled = false;
optimizationParams.maxTiltAccel = 2; % deg/s^2
optimizationParams.maxTiltRate = 5; %deg/s

beta = 0.65;
doPlotting = false; % disable this to not plot results
verboseOutput = false;


seedDir = fileparts(mfilename("fullpath"));
seedFilenames = struct( ...
    'alt', fullfile(seedDir,'Seeds/alt_seeds.dat'), ...
    'lon', fullfile(seedDir,'Seeds/lon_seeds.dat'), ...
    'lat', fullfile(seedDir,'Seeds/lat_seeds.dat'), ...
    'v', fullfile(seedDir,'Seeds/V_seeds.dat'), ...
    'fpa', fullfile(seedDir,'Seeds/fpa_seeds.dat'), ...
    'azmth', fullfile(seedDir,'Seeds/azmth_seeds.dat'), ...
    'mass', fullfile(seedDir,'Seeds/mass_seeds.dat'));

seeds = struct(...
    'alt_seeds', load(seedFilenames.alt),...
    'lon_seeds', load(seedFilenames.lon),...
    'lat_seeds', load(seedFilenames.lat),...
    'v_seeds', load(seedFilenames.v),...
    'fpa_seeds', load(seedFilenames.fpa),...
    'azmth_seeds', load(seedFilenames.azmth),...
    'mass_seeds', load(seedFilenames.mass));

% 3-sigma ranges
r_disp = 200; % m
lon_disp = 0.25; % deg
lat_disp = 0.25; % deg
v_disp = 3; % m/s
fpa_disp = 0.1; % deg
azmth_disp = 0.2; % deg
mass_disp = 0.0066; % fraction of nominal mass

caseCount = numel(seeds.alt_seeds);

Results = struct('k',{},'gamma',{},'kr',{},'tgo',{},'fuel_opt',{},...
    'coeff1',{},'coeff2',{},'coeff3',{},'coeff4',{},'exit_ok',{},'msg',{}, 'exitflag',{});

L_ref = 10000; A_ref = planetaryParams.gPlanet; T_ref = sqrt(L_ref/A_ref); V_ref = L_ref/T_ref; M_ref = 15103.0;
dispTime = tic;
parfor idx = 1:caseCount
    try
        dalt = seeds.alt_seeds(idx) * (r_disp / 3);
        dlon = seeds.lon_seeds(idx) * (lon_disp / 3);
        dlat = seeds.lat_seeds(idx) * (lat_disp / 3);
        dv = seeds.v_seeds(idx) * (v_disp / 3);
        dfpa = seeds.fpa_seeds(idx) * (fpa_disp / 3);
        dazmth = seeds.azmth_seeds(idx) * (azmth_disp / 3);
        mass_mult = seeds.mass_seeds(idx) * (mass_disp / 3);

        PDI = PDINom;
        PDI.altitude_km = PDINom.altitude_km + dalt/1000;
        PDI.lonInitDeg = PDINom.lonInitDeg + dlon;
        PDI.latInitDeg = PDINom.latInitDeg + dlat;
        PDI.inertialVelocity = PDINom.inertialVelocity + dv;
        PDI.flightPathAngleDeg = PDINom.flightPathAngleDeg + dfpa;
        PDI.azimuth = PDINom.azimuth + dazmth*pi/180;
        
        vehicle = vehicleNom;
        vehicle.massInit = vehicleNom.massInit * (1+ mass_mult);
        vehicle.dryMass = vehicle.massInit - 8248; % this way makes the dry mass variable but fuel amount constant
        

        [gammaOpt, krOpt, tgoOpt, aTOpt, exitflag, optFuel] = getParams(PDI, planetaryParams, targetState, vehicle, optimizationParams, beta, false, false, true);
        
        Results(idx).k = idx;
        Results(idx).gamma = gammaOpt;
        Results(idx).kr = krOpt;
        Results(idx).tgo = tgoOpt;
        Results(idx).fuel_opt = optFuel;
        Results(idx).coeff1 = gammaOpt*(krOpt/(2*gammaOpt +4) -1);
        Results(idx).coeff2 = (gammaOpt*krOpt/(2*gammaOpt+4)-gammaOpt-1);
        Results(idx).coeff3 = ((gammaOpt+1)/tgoOpt)*(1-krOpt/(gammaOpt+2));
        Results(idx).coeff4 = (krOpt/tgoOpt^2);
        Results(idx).exit_ok = true;
        Results(idx).msg = '';
        Results(idx).exitflag = exitflag;
    catch ME
        Results(idx).k = idx;
        Results(idx).exit_ok = false;
        Results(idx).msg = ME.message;
    end
end
elapsed = toc(dispTime)
if ~exist('Dispersion','dir')
    mkdir('Dispersion');
end

timeRun = datestr(now,'yyyymmdd_HHMMSS');
runDir = fullfile('DispersionOnlyThrust', timeRun);
mkdir(runDir);
matfile = fullfile(runDir, ['results_dispersion_' timeRun '.mat']);
csvfile = fullfile(runDir, ['results_dispersion_' timeRun '.csv']);

save(matfile,'Results');

T = table( (1:numel(Results)).', [Results.exit_ok].',[Results.gamma].', [Results.kr].', [Results.tgo].', [Results.fuel_opt].', ...
           [Results.coeff1].', [Results.coeff2].', [Results.coeff3].', [Results.coeff4].', ...
           'VariableNames', {'k','exit_ok','gamma','kr','tgo','fuel_opt','coeff1','coeff2','coeff3','coeff4'} );
writetable(T, csvfile);


statsPlotting(Results);