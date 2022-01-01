%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INIT_alpha_deltaH %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determination of runnout zones following Lied & Bakkehoi (1980)

% needs a digital elevation model as input. 
% output is:
%   - a .mat file with the runout angle alpha, quatifying exposure to
%   overhead terrain and avalanche hazard. 
%   - a .mat file with the height delta_H, used in quantifying the potential energy
%   of graviationally transported snowmass
% for further details, see Terleth MSc. thesis 
% 
% Script is partially based on numerical approach by Rickard Petterson (winstral
% script). Adjusted for current purposes in 2021 by Yoram Terleth.  

%% ----------------- USER INPUT ------------------------------------------- 

% INPUT DEM of glacier and surroundings
load tarfala.mat 
Z = tarfala.elev ;
E = tarfala.UTMeast ;
N = tarfala.UTMnorth ; 

% Optional: put boundaries on the area of interest (e.g. just the
% glacierized catchment) - saves processing time. 
spatial_boundaries = 'yes' ; % restrict in space?

% boundaries (can be left empty if previous is 'no'). 
bound_elev = 1000 ;                 % [m a.s.l.] lowest elevation to consider
bound_easting = [396150 399700] ;   % [m] limits in easting 
bound_northing = [7533150 7535550] ;% [m] limits in northing

% cell size of the DEM in metres
cx = 10 ; 

% print figure of output?
print_figs = 'yes' ; 

%% ------------------------------------------------------------------------

% alpha: runnout angle (Lied & Bakkehoi 1980)
%alpha = 25 ; 

% 
% Define the four cardinal directions
windDirs = [0.1, 90.1, 180.1, 270.1];

% Distances to range over
 sDist = 500;         % Largest distance (dmax)
 rDist = 100;         % Separation distance
 rDistMax = 1000;     % Upper limit
 
% Implement spatial boundaries 
if strcmp(spatial_boundaries,'yes')
    Z(Z < bound_elev)=nan; 
    Z(tarfala.UTMeast < bound_easting(1)) = nan ;
    Z(tarfala.UTMeast > bound_easting(2)) = nan ; 
    Z(tarfala.UTMnorth < bound_northing(1)) = nan ;
    Z(tarfala.UTMnorth > bound_northing(2)) = nan ; 
end 


%% Loop over cardinal directions: 
% calculate delta_H and alpha for each direction

 for j = 1:length(windDirs)
% ------------------------------------------------------------------------
% Paramter for number of rays within prevaling wind direction

rayInc = 3;         % Number of rays including limits of up-wind sector  
binSize = 30;       % Width in degrees of the up-wind sector

% ------------------------------------------------------------------------
windDir = windDirs(j); 


% Size of DEM
[nrows, ncols] = size(Z);
cellsize = cx;

% Find mximum kernel size in number of pixels
dMax = max(sDist);
kMax = ceil((dMax*2)/cellsize);
kMax = kMax + ~bitand(kMax, 1); % Odd kernel size

% Construct distance matrix for kernel
kx = [-(kMax-1)/2:(kMax-1)/2] * cellsize;
ky = flipud(kx');
[kX, kY] = meshgrid(kx,ky);
kD = sqrt(kX.^2 + kY.^2);

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
RX = nan(nrows, ncols);
RD = nan(nrows, ncols);
DH = nan(nrows, ncols); 

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
        alphat = nan(Raylen,1); 

        
        for n = 1:Raylen              
            % Indices to within separation distance for this ray
            iSi = indxSxi{n};
                                  
            % Indices of ray within range limit
            iSx = indxSx{n};

            % Nothing to do
            if isempty(iSx)
                continue 
            end
            
            % Calculate maximum slope for points
            % within range limit (eq 1. Winstral & Marks, 2002)
            Sx(n) = max(atand(kdZ(iSx) ./ kD(iSx))); 
            
            % claculate max elev diff 
            deltaH(n) = max(kdZ(iSx)); 
            

         end
                      
        % Output grid
        ix = (r - (kMax-1)/2);
        iy = (c - (kMax-1)/2);
        
        RX(ix,iy) = Sx_Avg;         % alpha angle
        DH(ix,iy) = deltaH_Avg ;    % delta_H
        
    end
end

% Verbose
fprintf(RevStr);
msg = ['Processing row: ' num2str(ix) ' of ' num2str(nrows)];
fprintf([msg '\n']);
fprintf(['Stop time: ' datestr(now()) '\n']);


  save(['runnout_dir_' num2str(windDirs(j)) '.mat'],...
             'RX','RD','DH');
 end
 
%% Combine directions
% we have parameters for each cardinal direction: now we combine them to
% retain only the maximum value. 

disp('Combining all four cardinal directions...')

% intialise variables
RD_total = zeros(size(RD)); 
RX_total = zeros(size(RX(:))); 
DH_total = zeros(size(DH(:))); 

for j = 1:length(windDirs)
    load(['runnout_dir_' num2str(windDirs(j)) '.mat']); 
    RD_total = RD_total + RD ; 
    RX = RX(:);
    for k = 1:length(RX_total)
        if RX(k) > RX_total(k)
           RX_total(k)=RX(k); 
        end
    end 
end 

for j = 1:length(windDirs)
    load(['runnout_dir_' num2str(windDirs(j)) '.mat']); 
    DH = DH(:);
    for k = 1:length(DH_total)
        if DH(k) > DH_total(k)
           DH_total(k)=DH(k); 
        end
    end 
end

% Reshape matrices 
runnout_alpha = reshape(RX_total, size(RD_total));
DH_total = reshape(DH_total, size(RD_total));


%% VISUALISATION

if strcmp(print_figs,'yes')
    figure 
    subplot(1,2,1)
    runnout_alpha(runnout_alpha==0)=nan ; 
    pcolor(E, N, runnout_alpha), shading flat
    hold  on 
    xlabel('UTM easting [m]'), ylabel('UTM northing [m]')
    if strcmp(spatial_boundaries,'yes')
        xlim([bound_easting(1) bound_easting(2)])
        ylim([bound_northing(1),bound_northing(2)])
    end 
    %colormap(viridis(256)); 
    caxis([min(runnout_alpha(:)),max(runnout_alpha(:))]) ;
    axis equal
    c = colorbar ;
    ylabel(c,'\alpha')

    subplot(1,2,2)
    DH_total(DH_total==0)=nan ; 
    pcolor(E,N, DH_total),shading flat
    xlabel('UTM easting [m]'), ylabel('UTM northing [m]')
    if strcmp(spatial_boundaries,'yes')
        xlim([bound_easting(1) bound_easting(2)])
        ylim([bound_northing(1),bound_northing(2)]);
    end 
    caxis([min(DH_total(:)),max(DH_total(:))]); 
    %colormap(viridis(256)); 
    axis equal
    c = colorbar ;
    ylabel(c,'\Delta_{H}')
end 

%% SAVE the variables for usage within ST-EBFM 

% save deltaH the max elevation difference with exposed terrain 
save([pwd '\deltaH.mat'],'DH_total');

% save runout_alpha the cells angle to overhead terrain. 
save([pwd '\runout_alpha.mat'],'runnout_alpha'); 
    



