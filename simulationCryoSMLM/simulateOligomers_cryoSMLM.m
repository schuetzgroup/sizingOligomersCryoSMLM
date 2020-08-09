function [sim_result,par] = simulateOligomers_cryoSMLM(par)

%% Set general parameters
if isfield(par,'roiSize')
    if isempty(par.roiSize)
        par.roiSize = round(sqrt(par.numberOligomers*1e5));
    end
else
    par.roiSize = round(sqrt(par.numberOligomers*1e5));
end

%% Create shape (regular polygon)
switch par.sizeSpecification
    case 'radius'
        r = par.radius;
    case 'sidelength'
        r = par.sidelength / (2*sin( deg2rad(180/par.oligoDegree) ));
end
allPoints = fft(eye(par.oligoDegree)*-r*1i);
shape = allPoints(:,2);
par.oligomerShape = [real(shape),imag(shape)];


%% Create molecule positions

oligomers_position = par.roiSize.*rand( par.numberOligomers, 2);

% Rotation matrices (R=[cos(t),-sin(t); sin(t),cos(t)])
theta = rand(size(oligomers_position,1),1)*2*pi; % angle between 0 and 360 degrees;
R_all = NaN(size(oligomers_position,1)*2,2);
R_all(1:2:end,1)=cos(theta);
R_all(1:2:end,2)=-sin(theta);
R_all(2:2:end,1)=sin(theta);
R_all(2:2:end,2)=cos(theta);

tmp_rotated = R_all*par.oligomerShape';
mols_oligo_position = reshape(tmp_rotated,2,[])'; % reshape, molecule order: 1st molecule of all oligomers, 2nd protein of all oligomers etc.
mols_oligo_oligomerID = repmat((1:par.numberOligomers)',par.oligoDegree,1); % oligomer ID (ID of corresponding oligomer)
mols_oligo_position = mols_oligo_position + repmat(oligomers_position,par.oligoDegree,1); % shift oligomers to their respective position

% Combined list of molecules in oligomers and molecules present as monomers
molecules_position = mols_oligo_position;
molecules_oligomerID = mols_oligo_oligomerID;

% Create table
molecules = table;
molecules.position = molecules_position;
molecules.oligomerID = molecules_oligomerID;
molecules.moleculeID = (1:size(molecules_position,1))';


%% Label molecules

labels = table;
numMols = size(molecules.position,1);
numLabeledMols = floor(par.labelEff*numMols);
idx_labeled = randsample((1:numMols)',numLabeledMols);

labels.position = molecules.position(idx_labeled,:);
labels.oligomerID = molecules.oligomerID(idx_labeled,:);
labels.labelID = (1:size(labels.position,1))';
labels.moleculeID = molecules.moleculeID(idx_labeled,:);

% Assign random dipole to labels
numLabels = size(labels,1);

% Uniform distribution on sphere
dipole_azimuth = 2*pi*rand(numLabels,1); % azimuth angle in [0,2*pi)
rvals = rand(numLabels,1);
dipole_inclination = asin(rvals); % elevation

% Radial distance (brightness)
dipole_distance = repmat(par.maxIntensity,[numLabels,1]);

% Transform to cartesian coordinates
labels.dipole = NaN(numLabels,3);
labels.dipole = [dipole_azimuth,dipole_inclination,dipole_distance];


%% Create detections

num_label = size(labels.position,1);
blink_dist = par.blink_dist{1};
blinkInfo_numBlinks = ceil(random(blink_dist.num,num_label,1));

locs_tmp_position = repelem(labels.position, blinkInfo_numBlinks, 1); % repeat rows of label_position according to blinkInfo_numBlinks times
locs_tmp_labelID = repelem((1:num_label)', blinkInfo_numBlinks, 1); % repeat ID (index) of label according to blinkInfo_numBlinks times
locs_tmp_moleculeID = repelem(labels.moleculeID, blinkInfo_numBlinks, 1); % repeat corresponding molecule ID of label according to blinkInfo_numBlinks times

% add oligomer ID
locs_tmp_oligomerID = repelem(labels.oligomerID, blinkInfo_numBlinks, 1); % repeat corresponding molecule ID of label according to blinkInfo_numBlinks times

% add dipole orientation
locs_tmp_dipole = repelem(labels.dipole, blinkInfo_numBlinks, 1); % repeat dipole

% create table
locs_tmp = table;
locs_tmp.position = locs_tmp_position;
locs_tmp.oligomerID = locs_tmp_oligomerID;
locs_tmp.labelID = locs_tmp_labelID;
locs_tmp.moleculeID = locs_tmp_moleculeID;

locs_tmp.dipole = locs_tmp_dipole;

locs = locs_tmp;
blink_info = table;
blink_info.numBlinks = blinkInfo_numBlinks;


%% Alternating excitation with differently polarized excitation light

% First polarization direction (Nx)
Nx = locs.dipole(:,3) .* cos(locs.dipole(:,2)).^2 .* cos(locs.dipole(:,1)).^2;
Nx = poissrnd(Nx); % include Poisson shot noise

% Second polarization direction (Ny)
Ny = locs.dipole(:,3) .* cos(locs.dipole(:,2)).^2 .* sin(locs.dipole(:,1)).^2;
Ny = poissrnd(Ny); % include Poisson shot noise

brightness = [Nx,Ny];

% Calculate estimated intensities

% Background parameter
tau = (2.*pi.*par.background.*(par.sigma^2+par.pxSize^2/12)) ./ (brightness.*(par.pxSize^2));

% Determine deltaN
deltaN = sqrt( brightness.*(1+4.*tau+sqrt(tau./(14.*(1+2.*tau)))) );

% estimate intensity
brightness_error = normrnd(0,deltaN);
brightness = brightness + brightness_error;
brightness(brightness<0) = 0;

locs.brightness = brightness;


%% Displace localization according to localization precision

N = sum(locs.brightness,2); % sum along 2nd dimension

par.background = sqrt(2)*par.background; % combined noise for both polarizations (for uncorrelated, Gaussian noise)
tau = (2.*pi.*par.background.*(par.sigma^2+par.pxSize^2/12)) ./ (N.*(par.pxSize^2));

deltaX = sqrt( (par.sigma^2+par.pxSize^2/12)./N .* (1+4.*tau + sqrt(2.*tau./(1+4.*tau))) );
loc_error = normrnd(0,repmat(deltaX,1,2),size(locs.position,1),2);
locs_tmp_locPrec = deltaX;

% Displace each localisation with given localisation precision (pa)
locs.position = locs.position+loc_error;

% Store localisation precision
locs.locPrec = locs_tmp_locPrec;

% Filter data based on localization precision
locs = locs( locs.locPrec<=par.locPrecThreshold,: );


%% Results

% Adjust table columns
molecules = splitvars(molecules,'position','NewVariableNames',{'x','y'});
labels = splitvars(labels,'position','NewVariableNames',{'x','y'});
locs = splitvars(locs,'position','NewVariableNames',{'x','y'});
locs = removevars(locs,'dipole');
locs = splitvars(locs,'brightness','NewVariableNames',{'brightness1','brightness2'});

% save simulation results and parameters
sim_result.molecules = molecules; % molecules
sim_result.labels = labels;       % labels
sim_result.locs = locs;           % localizations
sim_result.parameters = par;

end