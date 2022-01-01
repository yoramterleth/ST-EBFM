function [B] = PadGrid(A,n)
% Pads matrix A with n elements around each side.
%
% INPUT:
%   A       Matrix to pad.
%   n       Number of cells to pad.  

% OUTPUT:
%   B       Padded matrix
%

% Pad the matrix
B = [repmat(nan(1,size(A,2)),n,1); A; repmat(nan(1,size(A,2)),n,1)];
B = [repmat(nan(size(B,1),1),1,n), B, repmat(nan(size(B,1),1),1,n)];