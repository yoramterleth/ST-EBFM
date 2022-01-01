%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MODEL PARAMETERS:
%%% User-defined run parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [grid, time, io, phys, H] = INIT_parameters()

clearvars;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Time parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time.ts = '15-Sep-1997 00:00';                                              % Date and time start run
time.te = '20-Sep-1997 00:00';                                              % Date and time end run
time.dt = .25;                                                             % Timestep (days)
time.tn = round((datenum(time.te)- ...                                     % Nr. of time-steps
    datenum(time.ts))/time.dt)+1;
time.dT_UTC = 1;                                                           % time difference relative to UTC (hours)          

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Grid parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grid.utmzone = 34;                                                         % UTM zone
grid.max_subZ = .1;                                                        % Maximum first layer thickness (m)
grid.nl = 5;                                                               % Number of vertical layers
grid.doubledepth = 0;                                                      % Double vertical layer depth at layer grid.split (1=yes, 0=no)
grid.split = [3;5;10;15];                                                  % Vertical layer nr's at which layer depth doubles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Gravitational Transport Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% headwall dem parameters
H.cs = 10 ;                                                                   % cell size of DEM
H.x_bounds = [396150 399700];                                                 % UTM easting limits of zone to compute grav. transport for
H.y_bounds = [7533150 7535550];                                               % UTM northing bounds of zone to compute grav. transport for

% Avalanche parameters beta max and d max 
H.Blim = 35 ;                                                               % [degrees]   % see willibald et al. 2020 for angle of repose
H.Dlim = .05 ;                                                               % [ m w e ]  
                                                                   
H.runout_a = 28 ;                                                            % [degrees]  % runnout angle: this is the lowest angle at which deposits occur, "headwall" only considers terrain above this angle 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Model physics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phys.percolation = 'normal';                                                % Water percolation scheme:
                                                                            %   - 'bucket': tipping-bucket method (all water added at the surface)
                                                                            %   - 'normal': normally distributed deep percolation
                                                                            %   - 'linear': linearly distributed deep percolation
                                                                            %   - 'uniform': uniformly distributed deep percolation
phys.snow_compaction = 'firn+snow';                                         % Snow and firn compaction scheme:
                                                                            %   - 'firn_only': apply Ligtenberg et al. (2011) for all snow and firn layers
                                                                            %   - 'firn+snow': apply Ligtenberg et al. (2011) for firn layers and Kampenhout et al. (2017) for seasonal snow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Input/output parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
io.homedir = pwd;                                                           % Home directory
io.outdir = [io.homedir '/Output/'];                                        % Output directory
io.rebootdir = [io.homedir '/Reboot/'];                                     % Restart file directory
io.example_run = 0;                                                         % Run example case (no user input required)
io.readbootfile = 0;                                                        % REBOOT: read initial conditions from file (1=yes, 0=no)
io.writebootfile = 0;                                                       % REBOOT: write file for rebooting (1=yes, 0=no)  
io.bootfilein = 'boot_final_stakes.mat';                                    % REBOOT: bootfile to be read  
io.bootfileout = 'boot_final_stakes.mat';                                   % REBOOT: bootfile to be written
io.freqout = 1;                                                             % OUTPUT: frequency of storing output (every n-th time-step)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Include/Exclude Acc. Components
io.climate_averaging = 1;                                                   % use mean and sum climatic values for timestep (1), or at precise time (averaged input)(0)
                                                                            % using mean of timestep not ecommende for long timesteps
io.WindSheltering = 1;                                                      % turn wind accumulation on or off                                                                          

io.grav_transport = 1;                                                      % turn grav. transport on or off

io.stakes_only = 0 ;                                                        % run model only at specific stake sites (requires specific grid info: see INIT_grid_read_data.m)
end

