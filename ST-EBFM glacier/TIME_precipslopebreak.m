%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% precip slope break: By Rickard Petterson (modified by Yoram Terleth to fit in EBFM)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% utilises winstral et al.'s sheltering index to recalculate snow on the
% ground 

function [IN] = TIME_precipslopebreak(IN,WD,grid,time,C)

    % load precipitation and determine rain/snow partition
    P = IN.P_2d .* (IN.T_2d < C.rainsnowT-1);
    P = P + IN.P_2d .* (C.rainsnowT-IN.T_2d+1)./2 .* ...
        (IN.T_2d < C.rainsnowT+1 & IN.T_2d > C.rainsnowT-1);
     
    % median wind direction over previous timestep
    vm = nanmedian(WD.wind_dir.Vindriktning((find(WD.wind_dir.time == time.TCUR_DT)-...
       (time.dt*24-1)):find(WD.wind_dir.time == time.TCUR_DT)));  
     
    if isnan(vm)
        disp('replaced NaN in wind direction')
        vm=270; 
    end 
    

    % Determine the dominant wind direction
    [~,j] = min(abs(vm - WD.windDirs));
   
    % Get the specific sheltering index file
    S = WD.S_all(j);
       
     % Rescale the sheltering index
    if verLessThan('matlab','9.3')
        Sx = rescaleSX(S.Sx,abs(min(S.Sx(:))/max(S.Sx(:))).*-1,1);
    else
        Sx = rescale(S.Sx,abs(min(S.Sx(:))/max(S.Sx(:))).*-1,1);
    end
    Sx(isnan(Sx)) = 0;
    S.Di(isnan(S.Di)) = 0;
    
    % Adjust for slope break
    Sx = Sx .* S.Di; 
    
    % Set right side up 
    Sx =flipud(Sx);
    
    % and the short scale sheltering 
    S10m = WD.S_10m(j); 
    Sx10m = rescale(S10m.Sx,abs(min(S10m.Sx(:))/max(S10m.Sx(:))).*-1,1);
    Sx10m(isnan(Sx10m)) = 0 ; 
    Sx10m = flipud(Sx10m); 
    
    % Remove any negative precipitation
    P(P < 0) = 0;
    %Remove NaNs
    P(isnan(P)) = 0 ; 
    
    % Rescale the variability so it can be added to the sheltering index
    if max(P(:)) > 0
        if verLessThan('matlab','9.3')
            Ps = rescaleSX(P,min(P(:))/max(P(:)),1);
        else
            Ps = rescale(P,min(P(:))/max(P(:)),1);
        end
    else
        Ps = zeros(size(P));
    end
    
    %re-remove NaNs
    Ps(isnan(Ps))=0 ;
    
    % Add the rescaled spatial precipitation variability
    % forming an accumulation factor
    Af = Sx + Ps + (Sx10m);
    
    % Negative accumulation factor means erosion, 
    % so set it to zero
    Af(Af < 0) = 0;
    
    % Normalize the accumulation factor to the total 
    % so the sum will equal 1.
    Af = Af./sum(Af(:));

    % Calculate the total precipitation volume
    Pv = sum(P(:));

    % Multiply with precipitation volume to get distributed precipitation
    A = Af .* Pv;
    
    % Now store the snow drift contribution only
    IN.SD_2d = A - IN.P_2d ; 
    
    % more direct precip than snow drift = erosion, = zero wind accum 
    IN.SD_2d(IN.SD_2d<0) = 0 ; 
    
    % set to 1d input vector
    IN.SD = IN.SD_2d(grid.mask_2D(:)==1) ;
    IN.vm = vm .* ones(size(IN.SD)); 

end

%%% END %%%
