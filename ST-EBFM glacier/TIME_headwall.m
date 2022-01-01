%% TIME_Headwall %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% contribution of snow to glacier from surrounding terrain
%%% reads in climate variables, applies winstral, calculate gravitational 
%%% transport based on Gruber, 2007 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H,IN] = TIME_headwall(H,IN,C,time,grid)

%% declare topography  

z = H.z ; 
%beta = H.beta ;
fnb = H.fnb ; 
Dmax = H.Dmax ; 

%% Load variables from previous timestep
Acc= zeros(size(H.Acc)) ; 
M = zeros(length(z(:,1)), length(z(1,:))); 
Mstore = H.M ;

% preallocate matrix size to new variables for speed
D = zeros(length(z(:,1)), length(z(1,:)));
Fnb = zeros(length(z(:,1)), length(z(1,:)), 4) ;
delta_T = zeros(size(z)); 


%% load precipitation: 
I = IN.P_2d .* (IN.T_2d < C.rainsnowT-1);
I = I + IN.P_2d .* (C.rainsnowT-IN.T_2d+1)./2 .* ...
        (IN.T_2d < C.rainsnowT+1 & IN.T_2d > C.rainsnowT-1);
    
% add mass from snow drift
I = I + IN.SD_2d ; 

%% load the order of consideration 
ii = H.ii ;
jj = H.jj ; 

%% If there is new snow:
if nansum(I(:))> 0 
    
    %%% START of LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Snow Tranport...')

    %% gravitational transport     
    for i = 1:length(H.A(:,1))

       if ~isnan(H.A(i,3))       

        %% find location 
        ix = ii(i); % the x coord
        iy = jj(i); % the y coord

        %% calculate the inflow 
        M(iy,ix) = M(iy,ix) + I(iy,ix) ; % add the new precip to the mass 

        %% calculate outflow 

        % is snowmass larger than Dmax?
        if M(iy,ix) < Dmax(iy,ix)
           D(iy,ix) = M(iy,ix) ; 
        else 
           D(iy,ix) = Dmax(iy,ix) ;
        end 
        
        % register the relevant output: new deposits only
        % Acc(iy,ix) = D(iy,ix)-(Mstore(iy,ix)+I(iy,ix)); 
        
        for ik = 1:4 

            if isnan(fnb(iy,ix,ik))
                Fnb(iy,ix,ik) = 0 ; 
            else 
            % calculate the ammount of mass to each surrounding cell
            Fnb(iy,ix,ik) = (M(iy,ix) - D(iy,ix)) * fnb(iy,ix,ik) ; 
            end 
        end

        if iy > min(H.jj(:)) && ix > min(H.ii(:))
        % adjust the mass of each surroudning cell
        M(iy-1,ix) = M(iy-1,ix) + Fnb(iy,ix,1) ; % inflow to cell above NB 1
        M(iy,ix-1) = M(iy,ix-1) + Fnb(iy,ix,2) ; % inflow to cell to the left NB 2
        M(iy,ix+1) = M(iy,ix+1) + Fnb(iy,ix,3) ; % inflow to cell to the right NB 3
        M(iy+1,ix) = M(iy+1,ix) + Fnb(iy,ix,4) ; % inflow to cell to the bottom NB 4

       % Subtract snowfall to keep only avalanche fed accumulation
        Acc(iy,ix) = D(iy,ix)-I(iy,ix); 
        
        % register the relevant output: new deposits only
        Acc(iy,ix) = D(iy,ix)-(Mstore(iy,ix)+I(iy,ix)); 

        else 
            M(iy,ix) = 0 ; % set M to zero at padded boudnaries
            Acc(iy,ix) = 0 ; % set Acc to zero at padded boundaries

        end 
        
        %% Friction driven temperature increase %%
        if Acc(iy,ix) > 0 
            
        % elevation diff to start zone 
        Mr = I(iy,ix)*C.Dwater ; 
        % avergae elevation diff to mass entrainment zone
        Me = (Acc(iy,ix) - I(iy,ix)) * C.Dwater ; 
        % ensure no neg entrainement sum;
        Me(Me<0)= 0 ;    
        
        % from Steinkogler et al. (2015)
        delta_H(iy,ix) = (H.deltaHr(iy,ix) * (Mr + (0.5 * Me)))/(Mr + Me);
        delta_H(delta_H > H.deltaHr(iy,ix)) = 0 ; % ensure stability
        delta_H(delta_H < 0 ) = 0 ; 
        delta_T(iy,ix) = (C.g * delta_H(iy,ix))/C.cp ; 
        
        else 
        delta_T(iy,ix) = 0; 
        end 
        %% end of friction driven temp. increase
       end
    end

else 
    M = M + zeros(size(z));
    D = D + zeros(size(z));
    Acc = zeros(size(z));
    delta_T = zeros(size(z)); 
end

%%% end of loop %%% 

Acc(Acc<0)=0; % delete neg values produced when Dlim = 0 

%% store relevant output %%
IN.D = Acc(grid.mask_2D(:)==1) ; 
IN.delta_T = delta_T(grid.mask_2D(:)==1) ; 
H.M = M ; 
H.D = D ;  
end 
    


