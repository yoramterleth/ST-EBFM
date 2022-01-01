%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INIT_sheltering_index %%
%%%%%%%%%%%%%%%%%%%%%%%% %% 

% This script should be used to generate the sheltering index files required in
% INIT_wind_drift.m

% The output is one .mat file for each of the considered directions, and
% are stored in a directory named after the considered sheltering distance.
% Apply this script for each sheltering distance that needs to be investigated.  

% Script and winstral.m function initially written by Rickard Pettersson 2003(?). Adapted for usage with 
% ST-EBFM in 2021 by Yoram Terleth

%% -------------- USER INPUT -------------------------------------------%% 

% load digital elevation model of glacier + its surroundings. 
% e.g.: load (DEM)
Z = DEM.elev ; 

% cell sizes of the dem, in metres
cx = 15 ; 
cy = 15 ; 

% Dominant wind direction to calculate. Given in degrees, North = 0.
windDirs = [0.1:10:350.1]; % 10 degree intervals 

% Sheltering distance given in meters. 
sDist = 300; % used to be 100

%% ------- THESE PARAMETERS DO NOT NORMALLY NEED TO BE CHANGED -------------

% Slope break distance in meters
rDist = sDist + 75 ; % used to be + 75;

% Maximium distance for calculating slope breaks given in meters
rDistMax = 1000; % used to be 1000

% fake an R
R = "no R used" ; 

%% -------------------------------------------------------------------------

% Loop over all wind directions
for windDir = 1:length(windDirs)
    
    % Run the Winstral model ( Di = DD ) 
    [Sx,Sb,Si,So,Di] = winstral(Z,cx,cy, windDirs(windDir), ...
                                            sDist, rDist, rDistMax); 
    
    % make folder for current sheltering dist if it doesn't yet exist
    % if it does exist, files will be overwritten. 
    if ~exist([num2str(sDist) 'm'], 'dir')
       mkdir([num2str(sDist) 'm'])
    end
    
    % Save files
    save([pwd '\' num2str(sDist) 'm\sheltering_dir_' num2str(windDirs(windDir)) '.mat'],...
             'Sx','Sb','Si','So','Di','windDir', 'sDist', 'rDist', ...
             'rDistMax','Z','R');
    
end

