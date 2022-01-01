%% stability %%
%%%%%%%% inserted to prevent everlasting loop...
% when testing it is sometimes convenient to run longer timesteps or few
% vertical layers, which occasionally leads to instabilities in the snow
% model. This script keeps values within a realistic range and ensures the
% model doesn't fail. 
% Attention: misestimations can still occur within the realistic ranges. 

OUT.subD(OUT.subD(:)>1000) = C.Dice; 
OUT.subD(OUT.subD(:)<10) = C.Dfreshsnow ; 
OUT.subT(OUT.subT<=0) = 0.1 ; 
OUT.subT(OUT.subT>C.T0+1) = C.T0; 
OUT.subZ(OUT.subZ>4) = 4 ; 
OUT.subZ(OUT.subZ<0) = 0 ; 
OUT.subW(OUT.subW>1000) = 1000 ; 
OUT.subW(OUT.subW<0) = 0 ;
OUT.Tsurf(OUT.Tsurf<0) = 0.1 ; 
OUT.Tsurf(OUT.subT>C.T0+1) = C.T0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%