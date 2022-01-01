%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% METEOROLOGICAL FORCING:
%%% - Specify / read meteorological input for current time-step, including
%%%     air temperature, precipitation, relative humidity, cloud cover,
%%%     wind-speed and air pressure
%%% - Determine derived meteorological fields, including snowfall, 
%%%     rainfall, annual accumulation, specific humidity, vapor pressure,
%%%     air density, time since last snowfall and potential temperature

%%% - Read in Wind drift 

%%% - Read in Gravitational transport 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [IN,OUT,H] = TIME_climate_forcing(C,grid,IN,t,time,OUT,io,WD,H)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SPECIFY/READ METEO FORCING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IN = TIME_climate_forcing_read_data(IN,io,C,time,grid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DERIVED METEOROLOGICAL FIELDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Snowfall & rainfall
IN.snow = IN.P .* (IN.T < C.rainsnowT-1);
IN.rain = IN.P .* (IN.T > C.rainsnowT+1);
IN.snow = IN.snow + IN.P .* (C.rainsnowT-IN.T+1)./2 .* ...
    (IN.T < C.rainsnowT+1 & IN.T > C.rainsnowT-1);
IN.rain = IN.rain + IN.P .* (1+IN.T-C.rainsnowT)./2 .* ...
    (IN.T < C.rainsnowT+1 & IN.T > C.rainsnowT-1);

% Annual snow accumulation
if t==1
    OUT.ys =[];
    OUT.ys(1:grid.gpsum,1)= 500 ; 
end 

OUT.ys = (1.0-1.0/(C.yeardays/time.dt)).*OUT.ys + IN.P.*1d3;
logys = log(OUT.ys);
IN.yearsnow = repmat(OUT.ys,[1 grid.nl]);
IN.logyearsnow = repmat(logys,[1 grid.nl]);

% Vapor pressure & specific humidity
VPsat = C.VP0.*exp(C.Lv/C.Rv.*(1.0./273.15-1.0./IN.T)) .* (IN.T>=273.15)...
    + C.VP0.*exp(C.Ls/C.Rv.*(1.0./273.15-1.0./IN.T)) .* (IN.T<273.15);
IN.VP = IN.RH .* VPsat;
IN.q = IN.RH .* (VPsat .* C.eps ./ IN.Pres);
    
% Air density
IN.Dair = IN.Pres./C.Rd./IN.T;

% Time since last snow fall event
OUT.timelastsnow(IN.snow/(time.dt*24*3600)>C.Pthres) = time.TCUR;
if t==1
    OUT.timelastsnow(:) = time.TCUR; 
end
IN.timelastsnow = OUT.timelastsnow ; 

% Potential temperature & lapse rate
IN.Theta = IN.T.*(C.Pref./IN.Pres).^(C.Rd/C.Cp);
IN.f_Pottemp = fit(grid.z,IN.Theta,'poly1');
IN.Theta_lapse = max(IN.f_Pottemp.p1,0.0015);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% READ WIND DRIVEN SNOW DRIFT 
if io.WindSheltering && min(IN.WS(:)) >= WD.thresholdWS
    [IN]=TIME_precipslopebreak(IN,WD,grid,time,C);
else
    IN.SD = zeros(size(IN.P)); 
    IN.vm = 500; 
    IN.SD_2d = zeros(size(IN.P_2d));
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% READ GRAVITATIONAL SNOW TRANSPORT
if io.grav_transport 
    [H,IN] = TIME_headwall(H,IN,C,time,grid); 
else 
    IN.D = zeros(size(IN.P)) ; 
    IN.delta_T = zeros(size(IN.P)); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STORE RELEVANT VARIABLES IN OUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OUT.climT = IN.T;
OUT.climP = IN.P;
OUT.climC = IN.C;
OUT.climRH = IN.RH;
OUT.climWS = IN.WS;
OUT.climPres = IN.Pres;
OUT.climsnow = IN.snow;
OUT.climrain = IN.rain;
%OUT.direct_snow = IN.direct_snow;  
OUT.wind_drift = IN.SD;     
OUT.avy_dep = IN.D ;  
OUT.delta_T = IN.delta_T ; 
OUT.vm = IN.vm ; 
%OUT.total_WD = IN.SD_2d ; 
end