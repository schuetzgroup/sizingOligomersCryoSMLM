function [clustera,ncla] = cluster_intensity(locs_brightness,D)
% Cluster localizations based on brightness
[clustera,ncla]=cluster_spatial(locs_brightness,D);
end

