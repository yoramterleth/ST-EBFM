%% % INIT headwall %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Initialisation of terrain parameters for gravitational transport %%%%%
%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Topography Initialisation 
function [H] = INIT_headwall(H, grid,io)

disp('Initialising surrounding topography...')

%% Grid 

% load the alpha file: generate prior to model run with INIT_alpha_deltaH.m 
load([io.homedir '\Headwall\runout_alpha.mat']); 

% load the deltaH file: generate prior to model run with INIT_alpha_deltaH.m
load([io.homedir '\Headwall\deltaH.mat']);
H.deltaHr = .75 * DH_total ; 

% select zone of interest in DEM and pad it with one cell:
ZI=(grid.x_2D > H.x_bounds(1)-H.cs) .* (grid.x_2D < H.x_bounds(2)+H.cs) .*...
    (grid.y_2D > H.y_bounds(1)-H.cs).* (grid.y_2D < H.y_bounds(2)+H.cs); 

% apply to all relevant input 
x = ZI .* grid.x_2D ; 
y = ZI .* grid.y_2D ; 
alpha = ZI .* grid.aspect_2D ; 
beta = ZI .* grid.slope_2D ; 

% to gain speed: only consider zones that are above the maximum slope 
% angle of runnout: turn on and off & set threshold angle in 'avalanche'. 
if H.runout_a > 0  
    avy_slopes = zeros(size(beta)); % make logical array:
    avy_slopes(runnout_alpha > H.runout_a)=1; % 1: slope > max angle for runout
    z = ZI .* grid.z_2D .* avy_slopes ; % select only those zones for analysis
    z(z==0)= NaN ;
else 
    z = ZI .* grid.z_2D ; 
    z(z==0)= NaN ;
end 

%% flow field derivation 

% preallocate variables for speed:
Lnb = zeros(length(z(:,1)), length(z(1,:)), 4) ;
dz = zeros(length(z(:,1)), length(z(1,:)), 4) ;
Cnb = zeros(length(z(:,1)), length(z(1,:)), 4) ;
fnb = zeros(length(z(:,1)), length(z(1,:)), 4) ;
Dmax = zeros(length(z(:,1)), length(z(1,:))) ; 

% go over DEM
for j = 2: (length(z(:,1))-1) % x direction
    for i = 2:(length(z(1,:))-1) % y direction 
        
        % only if the cell is in an interesting area (gain speed)
        if ~isnan(z(j,i))
            
        % define flow width 
        Lnb(j,i,1) = cos(deg2rad(alpha(j,i))) * H.cs ; % flow length to NB 1 (above)
        Lnb(j,i,2) = -sin(deg2rad(alpha(j,i))) * H.cs ; % flow lenght to NB 2 (left)
        Lnb(j,i,3) = sin(deg2rad(alpha(j,i))) * H.cs ; % flow length to NB 3 (right)
        Lnb(j,i,4) = -cos(deg2rad(alpha(j,i))) * H.cs ; % flow length to NB 4 (bottom)
        
        % heigth difference with each neighbour
        dz(j,i,1) = z(j,i) - z(j-1,i) ; % elev diff (dz) with NB 1 cell above
        dz(j,i,2) = z(j,i) - z(j,i-1) ; % dz with NB 2 left 
        dz(j,i,3) = z(j,i) - z(j,i+1) ; % dz with NB 3 right
        dz(j,i,4) = z(j,i) - z(j+1,i) ; % dz with NB 4 bottom 
        
        % correct the flow width for uphill cells
        for k = 1:4 
            
            % define H function for corrections: 
            if dz(j,i,k) > 0 && runnout_alpha(j,i)> (H.runout_a+1) 
            H0 = 1;                    % : no further movement for slopes    
            else                       % below runout_a + 1 ; to prevent 
            H0 = 0;                    % mass loss at boundary
            end 
            
            % uphill flows are set to zero flow width: 
            Cnb(j,i,k) = Lnb(j,i,k) * (H0 * dz(j,i,k)) * (H0 * Lnb(j,i,k)) ; 
        end 
                   
        % compute normalised fraction draining into each adjacent cell
        for k2 = 1:4
            fnb(j,i,k2) = Cnb(j,i,k2)/(Cnb(j,i,1)+Cnb(j,i,2)+Cnb(j,i,3)+Cnb(j,i,4));
        end 
              
        % compute the local Dmax: the maximum snowmass the cell can hold
        % before avalanching 
        if beta(j,i) < H.Blim && runnout_alpha(j,i) > H.runout_a 
            Dmax(j,i) = (1-(beta(j,i)/H.Blim))*H.Dlim ; 
        else 
           Dmax(j,i) = 0 ; 
        end
        
    end
   end
end 

%% Set order in which the cells should be considered: from high z to low z

% make x and y matrices 
ii = repmat([1: length(z(1,:))],length(z(:,1)),1) ; % x direction indices
jj = repmat([1: length(z(:,1))]',1,length(z(1,:))) ; % y direction indices

% conctenate and order from high to low along cell elevation z(:)
A = [ii(:),jj(:),z(:)];
A= sortrows(A,3,'descend'); 
A = rmmissing(A) ; % delete the NaNs from A so we have less loop iterations
H.ii = A(:,1);       % to go trhough. Indices are preserved.
H.jj = A(:,2); 
H.A = A ; 
%% Store Important Variables
H.z = z ;
H.beta = beta ;
H.alpha = alpha ; 
H.fnb = fnb ;
H.Dmax = Dmax ; 
H.runnout_alpha = runnout_alpha ; 

%% Initialise remaining variables for speed
H.M = zeros(length(z(:,1)), length(z(1,:)));
H.D = zeros(length(z(:,1)), length(z(1,:)));
H.Fnb = zeros(length(z(:,1)), length(z(1,:)), 4) ;
H.Acc = zeros(length(z(:,1)), length(z(1,:)));

%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end 