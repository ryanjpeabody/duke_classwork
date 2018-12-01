% ve_dineof_north_atlantic_chl
%--------------------------------------------------------------------------
% Uses Variable EOFs Data INterpolating Orthogonal Functions (VE-DINEOF) to
% interpolate missing North Atlantic satellite chlorophyll a data. Methods
% from Beckers and Rixen (2003) and Ping et al. (2016)
%
%--------------------------------------------------------------------------
% AUTHOR INFORMATION
%   Ryan Peabody
%   Duke University
%   Nicholas School of the Environment
%
%--------------------------------------------------------------------------
% FILE LOCATION
%   /space/ryan/Projects/Natl_OceanColor/North_Atlantic_Variability/Matlab
%
%--------------------------------------------------------------------------
clearvars
close all

% USER INPUTS
%
% PLOT FIGURES?
  PlotFigures = 'Yes';

% SAVE OUTPUTS?
  SaveOutputs  ='Yes';
%
% LON BOUNDARIES
  LonBound = [-82 0];
%
% LAT BOUNDARIES
  LatBound = [15 65];

% ADDITIONAL EXCLUSIONS
  LatExclude = 55;
  LonExclude = -50;
%
%--------------------------------------------------------------------------
% HOUSEKEEPING
FilePath.home = ['/space/ryan/Projects/Natl_OceanColor/',...
    'North_Atlantic_Variability/Matlab'];
FilePath.chl = '/space/ryan/Data/Chlorophyll/OCCCI_V3/Processed';
FilePath.out = ['/space/ryan/Projects/Natl_OceanColor/',...
    'North_Atlantic_Variability/Outputs'];
FilePath.tool = '/space/ryan/Toolbox';

addpath(genpath(FilePath.tool))
    load NorthAtlanticTopography

cd(FilePath.chl)
    load Chl8DayDataGrid1.mat
    
cd(FilePath.home)

%--------------------------------------------------------------------------
% We need to reshape the data from a [x by y by t] sized array to a [x*y by
% t] sized array.
%
% But first, we need to differentiate data that are missing because they
% are actually land, from data that are missing because cloud cover, etc.
% The latter are what we are trying to interpolate.

chl = Chl8DayDataGrid1.chl;
lon = Chl8DayDataGrid1.lon;
    lonOld = lon;
lat = Chl8DayDataGrid1.lat;
    latOld = lat;
mDate = Chl8DayDataGrid1.mDate;
[latDim, lonDim, tDim] = size(chl);

[x, y] = meshgrid(NorthAtlanticTopography.lon(:,1),...
    NorthAtlanticTopography.lat(1,:));
[xQ, yQ] = meshgrid(Chl8DayDataGrid1.lon, Chl8DayDataGrid1.lat);
zInterp = interp2(x, y, NorthAtlanticTopography.z', xQ, yQ);

% Select those points with an elevation greater than 0 and flag chlorophyll
% values
isShallow = (zInterp >= -250);
isShallow = repmat(isShallow, [1 1 tDim]);
chlFlag = chl;
chlFlag(isShallow) = 1e36;

% And remove points north of N degrees and west of W degrees
isInside = (yQ >= LatBound(1) & yQ <= LatBound(2)) &...
    (xQ >= LonBound(1) & xQ <= LonBound(2));
isInside = repmat(isInside, [1 1 tDim]);
isOutside = ~isInside;
chlFlag(isOutside) = 1e36;

isOutside = (yQ > LatExclude) & (xQ <= LonExclude);
chlFlag(isOutside) = 1e36;


clear Chl8DayDataGrid1 isLand x y xQ yQ zInterp

% Now reshape the matrix
[lonGrid, latGrid] = meshgrid(lon, lat);
lon3d(:, :, 1:tDim) = repmat(lonGrid, [1 1 tDim]);
lat3d(:, :, 1:tDim) = repmat(latGrid, [1 1 tDim]);
t3d = repmat(permute(mDate, [3 2 1]), [latDim lonDim 1]);

chlX = reshape(chlFlag, [lonDim*latDim, tDim]);
lonX = reshape(lon3d, [lonDim*latDim, tDim]);
latX = reshape(lat3d, [lonDim*latDim, tDim]);
tX = reshape(t3d, [lonDim*latDim, tDim]);

% And remove land-flagged values, as well as values over N north 
isRemove = (chlX(:,1) == 1e36);
chlX(isRemove, :) = [];
lonX(isRemove, :) = [];
latX(isRemove, :) = [];
tX(isRemove, :) = [];

clear isNorth isRemove 

pctMissing = sum(double(isnan(chlX)), 2)/tDim*100;

if strcmpi(PlotFigures, 'yes')
    
    figure(1)
        clf
        hold on
    m_proj('Miller', 'lon', LonBound, 'lat', LatBound)
    m_grid
    m_scatter(lonX(:,1), latX(:,1), 25, pctMissing, 'Filled')
        colorbar
    m_coast('patch', 'k');
    title('Percent of data missing')
end


% Randomly sample 3% of the valid data

isPresent = ~isnan(chlX);
chlXNonNan = chlX(isPresent);
lonXNonNan = lonX(isPresent);
latXNonNan = latX(isPresent);
tXNonNan = tX(isPresent);
idxFull = find(~isnan(chlX));

numData = length(chlXNonNan);
[~, idx] = datasample(chlXNonNan, ceil(0.03*numData), 'Replace', false);
lonSamp = lonXNonNan(idx);
latSamp = latXNonNan(idx);
chlSamp = chlXNonNan(idx);
tSamp = tX(idx);
idxSamp = idxFull(idx);

clear isPresent numData idxFull

% Replace sampled data with nans
chlX(idxSamp) = nan;

% Subtract spatio-temporal mean from X
chlMean = nanmean(nanmean(chlX));
chlXcorr = chlX - chlMean;

% And replace all nan values with 0s
isMissing = isnan(chlXcorr);
chlXcorr(isMissing) = 0;

% Compute singular value decomposition

ii = 0;
rmseTotal = 1e36;
XrIn = chlXcorr;
while ii < 500 
    ii = ii+1;
    
    % Compute SVD of reconstructed matrix
    [U, S, V] = svd(XrIn);
    
    % Determine optimal number of EOFs to use for this reconstruction
    jj=0;
    rmseOld = 1e36;
    while jj <= 300
        jj = jj+1;    
        Un = U(:, 1:jj);
        Sn = S(1:jj, 1:jj);
        Vn = V(:, 1:jj);
        Xr = Un*Sn*Vn';
        chlSampNew = Xr(idxSamp);
    
        rmseNew = sqrt(mean((chlSamp - chlSampNew).^2));
        rmse(ii,jj) = rmseNew;
        eofCount(ii,jj) = jj;
        if rmseNew < rmseOld
            rmseOld = rmseNew;
        elseif (rmseNew > rmseOld) && (jj > 10)
            break
        end
    end
    rmse(rmse==0) = nan;
    numEof = find(rmse(ii,:)==min(rmse(ii,:)));
    eofUsed(ii) = numEof;
    
    % Compute Xr using optimal number of EOFs
    Un = U(:, 1:numEof);
    Sn = S(1:numEof, 1:numEof);
    Vn = V(:, 1:numEof);
    Xr = Un*Sn*Vn';
    
    XrIn(isMissing) = Xr(isMissing);
    
    
    % And calculate rmse for sample points
    chlSampNew = XrIn(idxSamp);
    rmseTotalNew = sqrt(mean((chlSamp - chlSampNew).^2));
    if (abs(rmseTotal-rmseTotalNew) < 1e-4) %&& (ii > 10)
        rmseTotal = rmseTotalNew;
        XrFinal = Xr;
        break
    end
    rmseTotal = rmseTotalNew;
    
    disp(['Iteration ',num2str(ii),' : ',num2str(numEof), ' EOFs',...
        ' : rmse = ', num2str(rmseTotal)])

end
rmse(rmse==0) = nan;
eofCount(eofCount==0) = nan;

clear Xr
lon = LonBound(1) : LonBound(2);
lat = LatBound(1) : LatBound(2);
Xr = nan(length(lat), length(lon), tDim);
missing = nan(length(lat), length(lon));
for iLat = 1:length(lat)
    disp(lat(iLat))
    for jLon = 1:length(lon)
        try
            I = (latX == lat(iLat)) & (lonX == lon(jLon));
                Xr(iLat, jLon, 1:tDim) = XrIn(I);
        end
        try
            I = (latX(:,1) == lat(iLat)) & (lonX(:,1) == lon(jLon));
                missing(iLat, jLon) = pctMissing(I);
        end
    end
end

chlR = Xr+chlMean;
chlR(chlR<0) = 0;

if strcmpi(SaveOutputs, 'yes')
    ChlVeDineof.lon = lon;
    ChlVeDineof.lat = lat;
    ChlVeDineof.chl = chlR;
    ChlVeDineof.mDate = mDate;
    ChlVeDineof.pctMissing = missing;
    EofCounts.rmse = rmse;
    EofCounts.num = eofCount;
    
    cd(FilePath.out)
    save ChlVeDineof ChlVeDineof
    save EofCounts EofCounts
    cd(FilePath.home)
end


        
        




    
    