%% Sizing oligomeric structures based on cryo-SMLM data
% Script for running algorithm on a test dataset

% Corresponding publication:
% 'A workflow for sizing oligomeric biomolecules based on cryo single
% molecule localization microscopy'
% Magdalena C. Schneider, Roger Telschow, Gwenael Mercier, Montserrat
% Lopez-Martinez, Otmar Scherzer, and Gerhard J. Schütz
% Affiliations:
% TU Wien, Vienna, Austria
% University of Vienna, Vienna, Austria

% Input data: Localization map of oligomer particles given as mat- or csv-file.
%             Localizations need to be assigned to oligomers already.
%             The file must contain the following columns:
%             
%             x           ... x-coordinates of localizations
%             y           ... y-coordinates of localizations
%             brightness1 ... signal brightness for first polarization
%                             direction
%             brightness2 ... signal brightness for second polarization
%                             direction
%             locPrec     ... localization precision
%             oligomerID  ... ID of corresponding oligomer
%             
%             The oligomeric structure is supposed to be a regular polygon.
%             The oligomerization degree needs to be specified as input
%             parameter.
%
% Output data: Estimated oligomer radius and side length.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add folder of this script and subfolders to path
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));


%% Set input parameters

oligomerDegree = 4;

% Specify input file (.mat or .csv)
fileOligomerData = 'testdata/simulatedOligomers_1000particles_1e5photons.mat';
[~,~,fileExtension] = fileparts(fileOligomerData);

switch fileExtension
     case '.mat'
         load(fileOligomerData);
         locs = sim_result.locs;
     case '.csv'
         locs = readtable(fileOligomerData);
end
columnNames = locs.Properties.VariableNames;

% Plotting options
plotInput = true;

%% Plot input data
if plotInput
    figure
    plot(locs.x,locs.y,'.')
    axis equal
    xlabel('x')
    ylabel('y')
    title('Localization map of oligomer particles')
end



%% Run analysis
disp('Running analysis...')

[sidelength,radius] = estimateOligomerSidelength( locs,oligomerDegree );


%% Display results
fprintf('\n')
fprintf('Estimated oligomer radius: %.4f\n',radius)
fprintf('Estimated oligomer side length: %.4f\n',sidelength)
fprintf('(Units corresponding to input data.)\n')


