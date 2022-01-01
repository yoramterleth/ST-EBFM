%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Windsheletering INIT %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initiation of the windsnhetlering routine based on Winstral et al. 2002

% Some user info is required: wind directions to be considered, which
% should match the sheltering index files generated with INIT_sheltering
% index.
% Provide 2 sheltering index distances along with corresponding files: one 
% large scale (generally between 100m and 1000m, and one as small as the DEM
% resolution allows (e.g. 15m when the DEM res. is 10m)

function [WD] = INIT_wind_drift(io)

if io.WindSheltering
    disp('Reading sheltering files...')

    %% ------------------ USER INPUT -------------------------------------
    WD.thresholdWS = 5 ;            % threshold wind speed for snow drift 
    Ldistance = '750m';             % the distance to which terrain is sheltered (CALIBRATION)
    Mdistance = '15m';              % the small scale distance, min possible with DEM res.

    WD.windDirs = [0.1:22.5:360];   % the wind driections to consider. should match those given 
                                    % in INIT_sheltering_index.m

    % Provide timeseries of wind direction.  
    load([pwd '\Climate\EG.mat'])
    WD.wind_dir.Vindriktning = EG.wd ;  % wind direction vector
    WD.wind_dir.time = EG.time ;        % corresponding time vector
    
    %% --------------------------------------------------------------------

    for n = 1:length(WD.windDirs) 
        WD.S_all(n) = load([pwd '\WindSheltering\' Ldistance '\sheltering_dir_' num2str(WD.windDirs(n)) '.mat']);
        WD.S_10m(n) = load([pwd '\WindSheltering\' Mdistance '\sheltering_dir_' num2str(WD.windDirs(n)) '.mat']);
    end


else
    WD =[];
end

end