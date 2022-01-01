%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% USER INPUT: Provide grid information:
%%%         - input.x: 2-D array containing UTM easting coordinates (m)
%%%         - input.y: 2-D array containing UTM northing coordinates (m)
%%%         - input.z: 2-D array containing elevation (m)
%%%         - input.mask: 2-D array containing mask (0 = no glacier, 
%%%                       1 = glacier)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function input = INIT_grid_read_data(io)

%%% SPECIFY USER INPUT HERE!
% provide x UTM easting and y UTM northing, z elevation, mask should be
% logical 0 - 1 indicating glacier area. 

 %  example: using the tarfala DEM
    load([io.homedir '\Grid\tarfala_with10mDEM.mat']); 
    input.x = tarfala.UTMeast;
    input.y = tarfala.UTMnorth;
    input.z = tarfala.elev; 
    input.mask = double(tarfala.shading);
    
 % run EBFM model component only at stake locations: stakes.S should be a logical 
 % matrix of same size as DEM, with 1 at stake locations and zeros
 % elsewhere 
    if io.stakes_only
        load([io.homedir '\Grid\stakes.mat']); 
        input.mask = stakes.S; 
    end 


end