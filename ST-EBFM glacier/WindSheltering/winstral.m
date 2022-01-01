function [SX,SB,SI,SO,DD] = winstral(Z,cx,xy,windDir,sDist,rDist,rDistMax) 
% Derives the terrain parameter Sx following Winstral 2002. High values
% means sheltering effects, negative indicates topographically exposed.
% -------------------------------------------------------------------------
% SYNTAX:
%
%   [SX,SB,SI,SO,DD] = winstral(X,Y,Z,windDir,sDist,rDist,rDistMax) 
%
% INPUT:
% Z         Elevation data as matrix.
% cx,cy     Cell size of elevation matrix. 
% windDir   Direction of prevailing wind given in degrees. 
% sDist     Sheltering distance. The value determines the upwind distance 
%           to search for sheltering terrain features.
% rDist     Range to separation between inner and outer zone that
%           determines up-wind slope break.
% rDistMax  Maximum range of the outer zone for up-wind slope break 
%           determination. 
%
% OUTPUIT:
% SX        Sheltering index. Determined as the maximum angle between the 
%           current cell and the upwind cells up to dDist distance. 
%           Positive value indicate a sheltering terrain feature and a
%           negative value indicates topographically exposed cell. 
%           Magnitude is given in degrees. Maximum angle between the cell
%           and the the upwind cells up to
% SI        Maximum slope in inner zone that is used to determine the slope
%           break. Given in degrees.
% SO        Maximum slope in outer zone that is used to determine the slope
%           break. Given in degrees.
% SB        SB gives the upwind slope break. It is defined as the 
%           difference between maximum slope in inner zone and outer
%           zone. Negative value indicates an acute angle and positive
%           indicate obtuse angle. Value indicate severness of the 
%           angle in degrees.
%
%	Rickard Pettersson, Uppsala University
% -------------------------------------------------------------------------

% TESTING

% nwaves = 1;
% t = [0:99*nwaves];
% x = sin(2.*pi().*1/99*t);
% x = [zeros(1,200) x zeros(1,200)]; 
% Z = repmat(x',1,500)*500;
% cx = 10;
% 
% % Define parameters
% windDir = 180;
% 
% % Distances
% sDist = 200;         % Largest distance (dmax)
% rDist = 100;         % Separation distance
% rDistMax = 1000;     % Upper limit

% ------------------------------------------------------------------------
% Paramter for number of rays within prevaling wind direction

rayInc = 7;     % Number of rays including limits of up-wind sector  
binSize = 30;   % Width in degrees of the up-wind sector

% ------------------------------------------------------------------------

% Size of DEM
[nrows, ncols] = size(Z);
cellsize = cx;

% Find mximum kernel size in number of pixels
dMax = max(sDist, rDist + rDistMax);
kMax = ceil((dMax*2)/cellsize);
kMax = kMax + ~bitand(kMax, 1); % Odd kernel size

% Construct distance matrix for kernel
kx = [-(kMax-1)/2:(kMax-1)/2] * cellsize;
ky = flipud(kx');
[kX, kY] = meshgrid(kx,ky);
kD = sqrt(kX.^2 + kY.^2);

% convert kD double to kD int16 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kD = int16(kD);

% Index of center cell in kernel
indxck = sub2ind([kMax kMax],kMax/2 + 0.5,kMax/2 + 0.5);

% Calculate the ray directions
llim = windDir - 0.5 * binSize;
ulim = windDir + 0.5 * binSize;
RayAngles =llim:(rayInc - 2):ulim;

% Calculate rays for the kernel
Rays = calculateRays(RayAngles, [kMax kMax]);

% Construct kernel templates for different distances 
d = [sDist, rDist, rDist + rDistMax];

% Limit each ray depending on distance
[distMax, bSepDistIn, bSepDistOut] = distLimitRays(Rays, kD, d);

% Get number of rays
Raylen = size(Rays,3);

for n = 1:Raylen
    % Get distance limited ray 
    A = distMax(:,:,n);       % Range limit
    B = bSepDistIn(:,:,n);    % Inner range for slope breaks    
    C = bSepDistOut(:,:,n);   % Outer range for slope breaks   
    
    % Find indices of ray within range limit
    indxSx{n} = find(A);
    
    % Find indices to within separation distance
    indxSxi{n} = find(B);
   
    % Find indices on the ray outside the separation distance 
    % but below the maximum limit
    indxSxo{n} = find(C); 
    
    % Find fist pixel along ray outside the separation disntance 
    [~,i] = min(abs(indxSxo{n} - indxck));
    indxSxof(n) = indxSxo{n}(i);
end

% Pad the DEM to allow for the kernel at borders 
Zp = padGrid(Z, (kMax-1)/2);

% Initialize output terrain parameter grids
SX = nan(nrows, ncols);
DD = nan(nrows, ncols);
SB = nan(nrows, ncols);
SI = nan(nrows, ncols);
SO = nan(nrows, ncols);

% Loop over all pixels in DEM, starting and ending
% at first and last pixel of the un-padded DEM
fprintf(['Start time: ' datestr(now()) '\n']);
RevStr = '';
for r = ((kMax-1)/2 + 1):((kMax-1)/2 + nrows)
    % Verbose
    fprintf(RevStr);
    msg = ['Processing row: ' num2str(r-(kMax-1)/2) ' of ' num2str(nrows)];
    fprintf(msg);
    RevStr = repmat('\b',1,length(msg));
    
    for c = ((kMax-1)/2 + 1):((kMax-1)/2 + ncols)
        % Get indices in the kernel
        ir = [-(kMax-1)/2:(kMax-1)/2] + r;
        ic = [-(kMax-1)/2:(kMax-1)/2] + c;

        % Extract the elevation and elevation difference
        % for the kernel
        kZ = Zp(ir,ic);
        kdZ = kZ - Zp(r,c);
              
        % Prepare variables
        Sx = nan(Raylen,1);
        Si = nan(Raylen,1);
        So = nan(Raylen,1);
        
        for n = 1:Raylen              
            % Indices to within separation distance for this ray
            iSi = indxSxi{n};
            
            if ~isempty(iSi)
                % Indices on the ray outside the separation distance 
                % but below the maximum limit
                iSo = indxSxo{n};
                
                % Calculate maximum slope outside the separation distance
                if ~isempty(iSo)  
                    x = (kZ(iSo) - kZ(indxSxof(n)))./...
                        (kD(iSo) - kD(indxSxof(n)));            
                    So(n) = max(atand(x));   
                end
            else
                continue
            end
            
            % Calculate maximum slope for points
            % within separation distance 
            Si(n) = max(atand(kdZ(iSi) ./ kD(iSi)));
                        
            % Indices of ray within range limit
            iSx = indxSx{n};

            % Nothing to do
            if isempty(iSx)
                continue
            end
            
            % Calculate maximum slope for points
            % within range limit (eq 1. Winstral & Marks, 2002)
            Sx(n) = max(atand(kdZ(iSx) ./ kD(iSx)));         
        end
        
        % Average all rays
        Sx_Avg = mean(Sx);
        Sb_Avg = mean(Si-So);
        Si_Avg = mean(Si);
        So_Avg = mean(So);
        
        % Drift delineator
        %Sb_Avg > deg2rad(7) && SxO_Avg < deg2rad(5) 
        if Sb_Avg > 7 && So_Avg < 5 
            bDrift = 1;
        elseif isnan(Sb_Avg) || isnan(So_Avg)
            bDrift = NaN;
        else
            bDrift = 0;
        end
       
        % Output grid
        ix = (r - (kMax-1)/2);
        iy = (c - (kMax-1)/2);
        SX(ix,iy) = Sx_Avg;         % Sheltering
        SB(ix,iy) = Sb_Avg;         % Slope breaks
        SI(ix,iy) = Si_Avg;       % Slope inner range
        SO(ix,iy) = So_Avg;       % Slope outer range
        DD(ix,iy) = bDrift;         % Drift delineator
    end
end

% Verbose
fprintf(RevStr);
msg = ['Processing row: ' num2str(ix) ' of ' num2str(nrows)];
fprintf([msg '\n']);
fprintf(['Stop time: ' datestr(now()) '\n']);



%axis equal
