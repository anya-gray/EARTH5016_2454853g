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
% 1: He1
% 2: Bg
% 3: He2
% 4: Sa
% 5: Gr
% 6: Cm
% 7: Si
% 8: Ms
% 9: Air/water
matprop = [
% unit  conductivity  density  heat capacity  heat production
%  #        sigma       rho         C           Hr
   1	    3.6788	    2697.6	    1000	    4.172
   2	    1	        2000	    1000	    1
   3	    3.2197	    2703.5	    1000	    5.575
   4	    1	        2000	    1000	    1
   5	    1	        2000	    1000	    1
   6	    0.924	    2083.1	    1000	    1
   7	    1	        2000	    1000	    1
   8	    0.919	    1905.9	    1000	    1
   9	    1e-6        1000	    1000	    0];  % air/water

% get coefficient fields based on spatial distribution of rock units from image
% pay attention if any unit conversion is required!
rho    = reshape(matprop(units,3), Nz, Nx);
Cp     = reshape(matprop(units,4), Nz, Nx);
sigma  = reshape(matprop(units,2), Nz, Nx);
Hr     = reshape(matprop(units,5), Nz, Nx);

% find diffusivity at each coordinate
k0 = sigma ./ rho ./ Cp;


% remaining model parameters
dTdz_top = 0;                      % flux at top
dTdz_bot = 35/1000;                % flux at bottom
Ttop = 0;       % (C)              % air/water temperature
Tbot = dTdz_bot*5000;              % find temperature at model base


% should come up with a way to define T at the x boundaries
% T_x0


nop = 50;

tend = 2;
dt = 1;

% *****  RUN MODEL
run('transect2D.m');