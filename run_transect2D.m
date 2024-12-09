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
   1	    3.6788	    2697.6	    600	        4.172               % 1172 other option for C (data_1)
   2	    2.465	    2700	    770 	    2                   % sigma- chatGPT conversion from val in data_3; Hr data_4
   3	    3.2197	    2703.5	    600	        5.575               % ***where are these from??
   4	    0.77	    1942.3	    740	        1                   % C-quartz,data_1; sigma data_2
   5	    0.77	    2648	    740	        1                   % C-quartz,data_1; sigma data_2
   6	    0.924	    2081.7	    860	        1                   % C, data_1, rho, excel, sigma data_2
   7	    1.67	    1916	    910	        1                   % 
   8	    0.919	    1909.78	    740	        1                   % 
   9	    1e-6        1	        1000	    0];  % air/water

% get coefficient fields based on spatial distribution of rock units from image
% pay attention if any unit conversion is required!
rho    = reshape(matprop(units,3), Nz, Nx);
Cp     = reshape(matprop(units,4), Nz, Nx);   
sigma  = reshape(matprop(units,2), Nz, Nx);
Hr     = reshape(matprop(units,5), Nz, Nx)*10^-6;     % test the 1e-6 thing in the excel file


% find diffusivity at each coordinate
k0 = sigma ./ rho ./ Cp*10^3;


% calculate source term
source = Hr ./ rho ./ Cp*10^3;


% remaining model parameters
dTdz_top = 0;                      % flux at top
dTdz_bot = 35/1000;                % flux at bottom


yr    = 3600*24*365;  % seconds per year [s]
tend = 1e7*yr;

CFL   = 1/5;         % Time step limiter
nop   = 5000;          % output figure produced every 'nop' steps


T0 = 5; % surface air temperature


dTdz_boundaries = [0, 35/1000];





% *****  RUN MODEL
run('transect2D.m');