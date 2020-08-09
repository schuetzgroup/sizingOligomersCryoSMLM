function [columnXpos,columnYpos,columnNx,columnNy,columnLocPrec,columnOligomerID] = getColumnIndices( columnNames )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getColumnIndices

% getColumnIndices checks if columns with certain headers exist and returns
% respective column indices

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(columnNames)
    if strcmp(columnNames{i},'x')
        columnXpos = i;
    end
    if strcmp(columnNames{i},'y')
        columnYpos = i;
    end
    if strcmp(columnNames{i},'brightness1')
        columnNx = i;
    end
    if strcmp(columnNames{i},'brightness2')
        columnNy = i;
    end
    if strcmp(columnNames{i},'locPrec')
        columnLocPrec = i;
    end
    if strcmp(columnNames{i},'oligomerID')
        columnOligomerID = i;
    end
end

if ~exist('columnXpos','var')
    error('Input data does not contain x-position of localizations!')
end
if ~exist('columnYpos','var')
    error('Input data does not contain y-position of localizations!')
end
if ~exist('columnNx','var')
    error('Input data does not contain brightness values for first polarization direction!')
end
if ~exist('columnNy','var')
    error('Input data does not contain brightness values for second polarization direction!')
end
if ~exist('columnLocPrec','var')
    error('Input data does not contain localization precision!')
end
if ~exist('columnOligomerID','var')
    error('Input data does not contain assignment to oligomers!')
end

end

