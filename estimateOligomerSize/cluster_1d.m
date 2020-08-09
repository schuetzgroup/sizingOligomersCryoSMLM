function [clusterx,ncl] = cluster_1d(X,d)
% Function that returns a column where the ith element is the cluster number corresponding to the point X(i) 

l = size(X,1);
[Xsort, ind] = sort(X);
diff = Xsort(2:l) - Xsort(1:l-1);
match = (diff < d);
clusterx = -ones(l,1);
ncl = 0; incl = 0;
for i=1:l-1
    if match(i) 
        if (1-incl)
            ncl = ncl+1;
            clusterx(ind(i)) = ncl;
            clusterx(ind(i+1)) = ncl;
            incl = 1;
        else
            clusterx(ind(i+1))= ncl;
        end
    else
        incl = 0;
    end
end
