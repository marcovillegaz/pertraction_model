function [DiffusionCoefficients] = maxwellStefanDiffussion(...
    CompoundLibrary,UNIFACLibrary, xMembrane,Temperature)
%MAXWELLSTEFANDIFFUSSION Summary of this function goes here
%   input:
%       CompoundLibrary (struct): struct that contains the compounds in the
%       system and all their properties. 
%       UNIFACLibrary (struct): structure with the unifac interaaction
%       parameters of the componentes. 
%       xMembrane (array): molar fraction of componentes in the system. The
%       length of this array must be the same as the number of compounds in
%       compoundLibrary. 
%       Temperature: temeprature of the system [K]
%   Output: 
%       DifussionCoefficient (array): n by n matrix that contain all the
%       diffusion coefficientes acoordting to the Maxwell-Stefan Diffusion
%       Theory. [cm2/s]
%
% Reference of Maxwel-Stefan theory:
%   Taylor, R., & Krishna, R. (1993). Multicomponent Mass Transfer (1st ed.).
%       John Wiley & Sons.


end

