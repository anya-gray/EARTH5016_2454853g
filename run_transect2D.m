%***** RUN 2D MODEL FROM IMAGE ***********************************

% clear workspace
clear all; close all; %clc;


% **** MODEL PARAMTERS ****
yr    = 3600*24*365;    % seconds per year [s]
t_end = 1e6*yr+5000;         % when to stop the simulation (should be 1e6 years)
CFL   = 1/5;           % Time step limiter

plot_freq   = 1000;     % #timesteps between each animation frame


% **** BOOLEANS FOR WHICH MODEL/TEST TO RUN

% kind of temperature gradient to initialise (Gaussian is an unifinished
% test)
linear = true;          
gaussian = false;

% details for 2D Gaussian convergence test
Twidth = 1000;           % initial T peak width
Tpeak = 1000;            % initial T peak amplitude
T0 = 100;                % initial background temperature

% for showing evolution/ results
verification = false;     % energy conservation test? (will not work with Gaussian)
plotAnimation = true;     % show evolution?
plot_result = false;      % plot final result (steatd-state)?




% **** load model setup from image, interpolate to target grid size

W       = 16e3;     % domain width (width of image) [m]
Nx      = 200;      % number of rows 
h       = W/Nx;     % grid spacing based on image width and target grid size
n_units = 9;        % number of rock units contained in image



% ****** CALIBRATE MODEL *******

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
if linear
matprop = [
% unit  conductivity  density  heat capacity  heat production   porosity
%  #        kappa        rho         Cp           Hr (1e-6)       phi (%)
   1	    3.6788	    2697.6	    1172	    4.172             1.03      % C: data_1, phi data_7, others Excel
   2	    2.465	    2700	    979 	    3.7               1.2       % kappa: data_3 chatGPT conversion; Hr data_4; Cp data_1; phi data_8
   3	    3.2197	    2703.5	    1172	    5.575             1.03      % C: data_1, phi data_7, others Excel
   4	    0.77	    1942.3	    740	        0.75              25.33     % C: data_1 (quartz); kappa: data_2; Hr: data_5
   5	    0.77	    2648	    740	        0.95              32.5      % C: data_1 (quartz); kappa: data_2; Hr: data_5
   6	    0.924	    2081.7	    860	        1.43              0.55      % C: data_1; rho: excel avg; kappa: data_2; Hr: data_5 (clay only); phi: data_6 (clay only)
   7	    1.67	    1916	    910	        1.3               17        % kappa,rhoCp: data_1; Hr: data_5 (beach sands)
   8	    0.919	    1909.78	    740	        1.2               21.17     % kappa,rhoCp: data_1; Hr: data_5 (beach sands)
   9	    1e-6        1   	    1000	    0                 100];    
end

% ******* GAUSSIAN TEST CALIBRATION
% use spatially independent material properties for gaussian convergence test
if gaussian
    matprop = [
     % unit  conductivity  density  heat capacity  heat production  porosity
    %  #      kappa (1e3)   rho         Cp           Hr (1e-6)         phi
       1	    1	        2000	    1000	    1                   0
       2	    1	        2000	    1000	    1                   0
       3	    1	        2000	    1000	    1                   0
       4	    1	        2000	    1000	    1                   0
       5	    1	        2000	    1000	    1                   0
       6	    1	        2000	    1000	    1                   0
       7	    1	        2000	    1000	    1                   0
       8	    1	        2000	    1000	    1                   0
       9	    1           2000	    1000	    0                   0];
end

%**************



% get coefficient fields based on spatial distribution of rock units from image
kappa  = reshape(matprop(units,2), Nz, Nx);  
rho    = reshape(matprop(units,3), Nz, Nx);
Cp     = reshape(matprop(units,4), Nz, Nx);   
Hr     = reshape(matprop(units,5), Nz, Nx)*10^-6;     % ensure W/m3
phi    = reshape(matprop(units,6), Nz, Nx)/100;       % percent -> fraction


% ******* account for porosities of the sediments
% ASSUMPTION: all pores are filled with AIR
kappa = (1-phi).*kappa + phi.*1e-6;
rho = (1-phi).*rho + phi.*1000;
Cp = (1-phi).*Cp + phi.*1000;
Hr = (1-phi).*Hr;               % heat rate of air = 0 => phi term vanishes



% find diffusivity for all coordinates
k0 = kappa ./ rho ./ Cp;


% constant term for finding the total thermal energy 
if verification
    energy_factor = sum( rho(:).*Cp(:)*h*h );
end


% calculate source term throughout transect in linear case
if linear
    source = Hr ./ rho ./ Cp;
end



% ***** remaining model parameters
geotherm = 35/1000;                 % 35 degC/km geothermal gradient
base_flux = - k0(end,:)*geotherm;   % for enforcing the flux at the bottom

T_air = 5;                  % surface air temperature to enforce at each t


% *****  RUN MODEL
run('transect2D.m');

