function [sidelength,radius] = estimateOligomerSidelength(locs,oligomerDegree,analysisMethod)
%  estimateOligomerSidelength
%  INPUT:
%         - locs: Table of localizations from oligomer particles.
%                 Localizations need to be assigned to oligomers already.
%                 The tatble must contain the following columns:
%                 x           ... x-coordinates of localizations
%                 y           ... y-coordinates of localizations
%                 brightness1 ... signal brightness for first polarization
%                                 direction
%                 brightness2 ... signal brightness for second polarization
%                                 direction
%                 locPrec     ... localization precision
%                 oligomerID  ... ID of corresponding oligomer
%
%         - oligomerDegree: oligomeriation degree (i.e. number of protomers
%                           constituting an oligomer)
%
%         - analysisMethod: method one or method two as explained in the
%                           corresponding publication
%
%  OUTPUT:
%          - sidelength:    estimated sidelength of oligomer structure
%          - radius:        estimated radius of oligomer structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get column indices
columnNames = locs.Properties.VariableNames;
[columnXpos,columnYpos,columnNx,columnNy,columnLocPrec,columnOligomerID] = getColumnIndices( columnNames );
locs = table2array(locs);

% Set parameters
% Threshold for outlier removal
rmax = 100;
% Threshold for clustering of (Nx,Ny)
Nmax = max(locs(:,[columnNx,columnNy]),[],'all');
clusterThresholdN = 300 + Nmax/100;


% Sort localizations
sortlocs = sortrows(locs,columnOligomerID);
ind_oligos = [0; find(diff(sortlocs(:,columnOligomerID))); size(sortlocs,1)];

numberOligomers = size( unique(sortlocs(:,columnOligomerID)),1 );

fprintf(['Number of input oligomers: ' num2str(numberOligomers) '\n']);

% Initialize
isEligible = zeros(numberOligomers,1);
radii = NaN(numberOligomers,1);   % radii calculated for each oligomer
sum_var = NaN(numberOligomers,1); % sum of variances

for i=1:numberOligomers
    locs_cluster = sortlocs((ind_oligos(i)+1):ind_oligos(i+1),:);
    oligoID = locs_cluster(1,columnOligomerID);
    
    [clusterIDs,numberClusters] = cluster_intensity(locs_cluster(:,[columnNx,columnNy]),clusterThresholdN);
    
    if(numberClusters == oligomerDegree)
        isEligible(oligoID) = 1;
        
        corners = NaN(oligomerDegree,2); % coordinates of 'corners' (centers of mass of protomer clusters)
        varianceCorners = NaN(oligomerDegree,1); % variance of protomer cluster
        
        for j=1:numberClusters
            isInCorner = (clusterIDs == j);        % get points belonging to corner j
            numberLocsInCorner = sum(isInCorner);  % number of blinks corresponding to corner j
            
            corners(j,:)  = sum(locs_cluster(isInCorner,[columnXpos,columnYpos]),1)/sum(isInCorner);
            varianceCorners(j)  = sum(locs_cluster(isInCorner,columnLocPrec).^2)/(numberLocsInCorner^2);
        end
        
        centerInitial = sum(corners,1)/oligomerDegree; % initial guess for center (= center of mass)
        fittedCircle = LM(corners,[centerInitial(1),centerInitial(2),1]); % Levenberg-Marquardt
        
        sum_var(oligoID) = (1/(2*oligomerDegree)+1/(oligomerDegree^2))*sum(varianceCorners);
        
        radii(oligoID) = fittedCircle(3);
    end
end

% Outlier removal
radii(radii>rmax) = NaN;


%% Calculate corrected radii

switch analysisMethod
    case 'one'
        rhalf = radii./2;
        hasRealSolution = (rhalf.^2 > sum_var); % check that solution is not imaginary
        radii = rhalf + sqrt(rhalf.^2 - sum_var); % subtracting computed bias from estimation
        radii(~hasRealSolution) = NaN;
    case 'two'
        overallRadius_uncorr = nanmedian(radii);
        radii = radii - sum_var./overallRadius_uncorr;
end


%% Calculate overall estimation for radius and side length

fprintf(['Number of eligible oligomers used: ' num2str(sum(isEligible)) '\n']);

radius = nanmedian(radii);
sidelength = 2* radius * sin(pi/oligomerDegree);

end

