function varargout = distLimitRays(Rays, Distance, limDistances)
% 
%
% INPUT:
%
% OUTPUT:
%

% Logcial arrays for limited distances
distMax = Distance <= limDistances(1);
distSep = Distance <= limDistances(2);
distSep2 = Distance <= limDistances(3);

% For each ray
for r = 1:size(Rays,3)
    bDistMax(:,:,r) = distMax & Rays(:,:,r);
end

for r = 1:size(Rays,3)
    bSepDist(:,:,r) =  distSep & Rays(:,:,r);
end

for r = 1:size(Rays,3)
    bSepDist2(:,:,r) =  ~distSep & distSep2 & Rays(:,:,r);
end

varargout{1} = bDistMax;
varargout{2} = bSepDist;
varargout{3} = bSepDist2;
