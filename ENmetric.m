function [N,RMSE,R2] = calcENmetric(k,X,sX,Y,sY,sU)

% author: C. van der Tol
% date: 11 December 2025
%
% The equivalence ratio of the field reference and satellite data. 
% Source: 
% 1. Personal communication P. de Vis (NPL)
% 2. Section 3.2 of https://qa4eo.org/docs/2_Metrology_Document.pdf). 
% 
% usage: N = ENmetric(k,X,Y,sX,sY,sU)
%
% Input:
%   k: coverage factor (default: 1)
%   X: ground reference data (scalar, vector or matrix)
%   sX: uncertainty of reference data
%   Y: satellite data (scalar, vector or matrix)
%   sY: uncertainty of satellite data
%   sU: uncertainty due to difference in representation (the measurand
%   differs from the validated variable)
% output:
%   N: the equivalence ratio
%   RMSE: the root mean square error
%   R2: the squared correlation coefficient.
%
% Notes: 
% 1. The pair of X and sX may switch place with Y and sY with no
% further consequence for the output.
% 2. The dimensions of X,sX,Y,sY sU are arbitrary, but they should be equal
% for all inputs.

N = abs(X-Y)./(k*sqrt(sX.^2 + sY.^2 + sU.^2));
RMSE = sqrt(mean((X(:)-Y(:)).^2));

R = corrcoef(X(:),Y(:));
R2 = R(2).^2;
