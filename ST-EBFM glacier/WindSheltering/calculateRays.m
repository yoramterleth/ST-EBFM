function [iR] = calculateRays(rayAngles, siz)
% Calcualtes a logical raster containing the pixels that exists in 
% given directions from center pixel. The pixels are determined using 
% Bresenham's line algorithm. 
%
%INPUT:
% rayAngles     Vector or scalar containing angles to calculate. 
%               Given in degrees.     
% siz           Size of raster to project the ray on. Must have odd
%               sizeed rows and columns to define a clear center cell.
%
%OUTPUT:
% iR            Logical raster with cells in the direction of given 
%               direction from the center pixel set to true. The output
%               will be as 2D array with the size [siz(1), siz(2)]. 
%               iR will be a 3D array if rayAngles is given as vector
%               where the third dimesion has the size length(rayAngles)
%

% Prepare output matrix
iR = false(siz(1),siz(2),length(rayAngles));

% Create xoordinate matrix for pixels
[X,Y] = meshgrid((1:siz(1))-(siz(1)/2 + 0.5), (1:siz(2))-(siz(2)/2 + 0.5));

% Construct coordinate matrix for all pixels 
Coords = [X(:),Y(:) zeros(prod(siz),1)].';

% Loop over number of rays
for i = 1:length(rayAngles)  
  
    % Get the slope fractions in x, y, z 
    dc = [cosd(rayAngles(i)), sind((rayAngles(i))), 0].';
    
    % x*cos(a) + y*sin(a) 
    % Determine perpendicular distance to all cell from a perpendicular 
    % line to the line with the ray angle direction (ie the line angle-90), 
    % The distance is signed (negative on lefthand side, positive on 
    % right hand side).
    % All this assumes that the line goes through origo. 
    % Otherwise d = x*cos(a) - (y-m)*sin(a), where m is intercept.
    % This is used to limit the line only towards the given direction
    AcuteAngle = (dc.' * Coords) >= 0;
       
    % Calculate distance for pixels from line
    %Dist = sqrt(sum(cross(Coords , dc * ones(1, prod(siz))).^2));
    
    % Calculate indices of a line with specified angle (going through
    % origo) in a raster
    iLine = bresenham(rayAngles(i), siz);
    
    % Reshape and combine
    img = reshape(AcuteAngle & iLine, siz);
    
    % Rotate so 0 angle is to north
    iR(:,:,i) = rot90(img,1);
end

function ind = bresenham(angle, siz)
% Bresenham line algorithm to pixelize a line at a certain angle. 
%
%INPUT:
%   angle       Direction of line in degrees
%   siz         Size of raster to project the line on
%
%OUTPUT:
%   ind         Indices of line in raster of 'siz' size
%

% Center point of raster
cx = floor(siz(2)/2);
cy = floor(siz(1)/2);

% End points of line to rasterize
x1 = -cx; x2 = cx;
y1 = x1*tand(angle); y2 = x2*tand(angle);

% Delta values
dx = abs(x2-x1);
dy = abs(y2-y1);

steep = abs(dy)>abs(dx);
if steep 
    t = dx;
    dx = dy;
    dy = t; 
end

if dy == 0 
    q = zeros(dx+1,1);
else
    q = [0; diff(mod([floor(dx/2):-dy:-dy*dx+floor(dx/2)]',dx))>=0];
end

if steep
    if y1 <= y2 
        y = [y1:y2]'; 
    else
        y = [y1:-1:y2]';
    end
    if x1 <= x2
        x = x1 + cumsum(q);
    else
        x = x1 - cumsum(q);
    end
else
    if x1 <= x2 
        x = [x1:x2]'; 
    else
        x = [x1:-1:x2]'; 
    end
    if y1 <= y2
        y = y1 + cumsum(q);
    else
        y = y1 - cumsum(q); 
    end
end

% Limit the y within raster
i = (round(y) > cy | round(y) <-cy);
y(i) = [];
x(i) = [];

% Convert into indices into raster of siz size
i = sub2ind(siz,round(y)+cy+1,round(x)+cx+1); 
ind = false(1, prod(siz));
ind(i) = true;

