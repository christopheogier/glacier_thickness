% Biafo Hewitt GPR points exported as 13 individual points
% Fields stored as column vectors with one entry per point:
%   Hgpr    : ice thickness (m)
%   Hgprbot : Hgpr - 50 m
%   Hgprtop : Hgpr + 50 m
%   xgpr    : longitude coordinate for each point
%   ygpr    : latitude coordinate for each point
%   zgpr    : bed elevation coordinate (NaN, unknown)

clear;

% Concatenate all Hewitt points across the three profiles
Hgpr = [450 450 450 450 1400 1300 1200 300 350 700 650 500 400]';

% Uncertainty bounds
Hgprbot = Hgpr - 50;
Hgprtop = Hgpr + 50;

% Coordinates derived from the shapefile line geometries and transformed to
% WGS84 / UTM zone 43N (EPSG:32643).
xgpr = [
    555051.806;
    554968.464;
    554890.677;
    552630.948;
    553376.730;
    553656.398;
    553936.066;
    565760.005;
    566029.346;
    566328.613;
    566777.514;
    567196.488;
    567256.341
];

ygpr = [
    3983544.184;
    3984289.265;
    3984984.674;
    3980805.674;
    3981094.363;
    3981202.622;
    3981310.881;
    3967032.778;
    3967393.076;
    3967793.409;
    3968393.908;
    3968954.375;
    3969034.442
];

zgpr = nan(size(Hgpr));

outdir = fullfile(pwd, 'output');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

save(fullfile(outdir, 'biafo_GPR_data.mat'), 'Hgpr', 'Hgprbot', 'Hgprtop', 'xgpr', 'ygpr', 'zgpr');

fprintf('Saved %s\n', fullfile(outdir, 'biafo_GPR_data.mat'));
