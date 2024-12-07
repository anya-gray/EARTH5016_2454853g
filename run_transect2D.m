%***** RUN 2D MODEL FROM IMAGE ***********************************

% clear workspace
clear all; close all; %clc;

% load model setup from image, interpolate to target grid size

% domain width (must correspond to width of image) [m]
W       = 16e3;     

% number of pixels in z direction
Nx      = 200;      

% grid spacing based on image width and target grid size
h       = W/Nx;     

% number of rock units contained in image
n_units = 9;        


% units = value of each pixel
% D = tiff file depth
% Nx = number of pixels in x direction 
[units, D, Nz] = ModelFromImage('section.tiff', n_units, W, Nx);








% material properties for each rock unit (update based on your calibration)
% 1: He1    helmsdale granet phase 1
% 2: Bg     basement gneiss
% 3: He2    helmsdale granet phase 2
% 4: Sa     sand
% 5: Gr     gravel
% 6: Cm     clay, mudstone
% 7: Si     silt
% 8: Ms     mud, silt, sand
% 9: Air/water
matprop = [
% unit  conductivity  density  heat capacity  heat production
%  #        sigma       rho         C           Hr
   1	    3.6788	    2697.6	    600	        4.172               % 1172 other option for C
   2	    1	        2700	    770 	    1
   3	    3.2197	    2703.5	    600	        5.575
   4	    1	        1942.3	    740	        1                   % *used quartz in the paper for C
   5	    1	        2648	    740	        1                   % *
   6	    0.924	    2081.7	    860	        1                   % C from paper, rho from excel
   7	    1	        1916	    910	        1
   8	    0.919	    1909.78	    740	        1                   % *
   9	    1e-6        1	        1000	    0];  % air/water

% get coefficient fields based on spatial distribution of rock units from image
% pay attention if any unit conversion is required!
rho    = reshape(matprop(units,3), Nz, Nx);
Cp     = reshape(matprop(units,4), Nz, Nx);         % need to find values for this ?
sigma  = reshape(matprop(units,2), Nz, Nx);
Hr     = reshape(matprop(units,5), Nz, Nx);

% find diffusivity at each coordinate
k0 = sigma ./ rho ./ Cp;


% remaining model parameters
dTdz_top = 0;                      % flux at top
dTdz_bot = 35/1000;                % flux at bottom
Ttop = 0;       % (C)              % air/water temperature
Tbot = dTdz_bot*5000;              % find temperature at model base


yr    = 3600*24*365;  % seconds per year [s]
tend = 1e7*yr;

CFL   = 1/5;         % Time step limiter
nop   = 5000;          % output figure produced every 'nop' steps


T0 = 5; % surface air temperature

dTdz_boundaries = [0, 35/1000];

% calculate source term
source = Hr ./ rho ./ Cp;



% *****  RUN MODEL
run('transect2D.m');