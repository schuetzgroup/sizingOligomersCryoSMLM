function [clboth, ncl] = cluster_spatial(B,d)
%B is a 2 times l matrix with point coordinates and d is the distance we use to do the clustering.

l = size(B,1);

%Clustering in x
[clusterx, ~] = cluster_1d(B(:,1), d);

% Clustering in y
[clustery,~] = cluster_1d(B(:,2), d);
%fprintf(['there are ' num2str(nclx) ' clusters in the x direction and ' num2str(ncly) ' in the y direction \n']);

clustercon = [clusterx , clustery]; % Concatenate the two columns obtaining by clustering in x and y
[clussort,idx] = sortrows(clustercon, [1 2]); % Sort this matrix with respect to the first column, then with respect to the second one.
diff = clussort(2:l,:) - clussort(1:l-1,:);
match = ismember( diff, [0,0], 'rows') .* (min(clussort(2:l,:),[],2) > 0); %decide that two points are in the same cluster if they both share their cluster number in x and y (and also, are not alone in their cluster).
%fprintf(['there are ' num2str(sum(match)) ' matches \n']);
clboth = zeros(l,1);
ncl = 0; incl = 0;
for i=1:l-1
    if match(i) 
        if (1-incl)
            ncl = ncl+1;
            clboth(idx(i)) = ncl;
            clboth(idx(i+1)) = ncl;
            incl = 1;
        else
            clboth(idx(i+1))= ncl;
        end
    else
        incl = 0;
    end
end