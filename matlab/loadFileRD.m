function [labels_rdpoint rdpoint labels_rdopcurve rdopcurve] = loadFileRD

% plotrd - plots the rate-distortion operational curves
%
% Files:
%   rdpoint.out
%   rdopcurve.out
%
% Output:
%   
%


[labels_rdpoint,x,y] = readColData('rdpoint.out',9,0);

rdpoint = [x y];


[labels_rdopcurve,x,y] = readColData('rdopcurve.out',9,0);

rdopcurve = [x y];
