%***** RUN 2D MODEL FROM IMAGE ***********************************

% clear workspace
clear all; close all; %clc;


% **** MODEL PARAMTERS ****
yr    = 3600*24*365;    % seconds per year [s]
t_end = 1e6*yr;         % when to stop the simulation (should be 1e6 years)
CFL   = 1/10;           % Time step limiter
nop   = 3000;           % output figure produced every 'nop' iterations


% number of pixels for the convergence test
% NN = [100, 200, 300];
NN = [200];
% convergenceTest = true;
validation = false;

% show evolution?
plotAnimation = false;

% kind of temperature gradient to initialise
linear = true;        
gaussian = false;


% details for 2D Gaussian convergence test
Twidth = 1000;           % initial T peak width
Tpeak = 1000;            % initial T peak amplitude
T0 = 100;                % initial background temperature


% run model for each specified pixel size
for nn = 1:length(NN)
disp(['Iterating...', nn])



% load model setup from image, interpolate to target grid size

W       = 16e3;     % domain width (width of image) [m]
Nx      = NN(nn);   % number of rows 
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
% unit  conductivity  density  heat capacity  heat production
%  #      sigma (1e3)   rho         Cp           Hr (1e-6)
   1	    3.6788	    2697.6	    600	        4.172               % C: data_1 (1172 other option)
   2	    2.465	    2700	    770 	    3.7                 % sigma: data_3 (conversion from) chatGPT conversion; Hr data_4
   3	    3.2197	    2703.5	    600	        5.575               % ***where are these from??
   4	    0.77	    1942.3	    740	        0.75                % C: data_1 (quartz); sigma: data_2; Hr: data_5
   5	    0.77	    2648	    740	        0.95                % C: data_1 (quartz); sigma: data_2; Hr: data_5
   6	    0.924	    2081.7	    860	        1.43                % C: data_1; rho: excel avg; sigma: data_2; Hr: data_5 (clay only)
   7	    1.67	    1916	    910	        1.3                 % (MISSING REF FOR SIG, RHO, CP); Hr: data_5 (beach sands)
   8	    0.919	    1909.78	    740	        1.2                 % (MISSING REF FOR SIG, RHO, CP); Hr: data_5 (beach sands)
   9	    1e-6        1	        1000	    0];                 % *** should confirm values for simulatoin
end

% use spatially independent material properties for gaussian convergence test
if gaussian
    matprop = [
    % unit  conductivity  density  heat capacity  heat production
    %  #      sigma (1e3)   rho         Cp           Hr (1e-6)
       1	    1	        2000	    1000	    1
       2	    1	        2000	    1000	    1
       3	    1	        2000	    1000	    1
       4	    1	        2000	    1000	    1
       5	    1	        2000	    1000	    1
       6	    1	        2000	    1000	    1
       7	    1	        2000	    1000	    1
       8	    1	        2000	    1000	    1
       9	    1           2000	    1000	    0];     %sigma = 1e-6; rho = 1000 ?
end


% get coefficient fields based on spatial distribution of rock units from image
sigma  = reshape(matprop(units,2), Nz, Nx) *10^3;       % ensure W/m/K
rho    = reshape(matprop(units,3), Nz, Nx);
Cp     = reshape(matprop(units,4), Nz, Nx);   
Hr     = reshape(matprop(units,5), Nz, Nx)*10^-6;       % ensure W/m3


% find diffusivity at each coordinate
k0 = sigma ./ rho ./ Cp;


% calculate source term throughout transect in linear case
if linear
source = Hr ./ rho ./ Cp;
end



% ***** remaining model parameters

dTdz_boundaries = [0, 35/1000];       % heat flux at top and bottom      % this (list) can maybe change later

T_air = 5;                  % surface air temperature to enforce at each t


% *****  RUN MODEL
run('transect2D.m');

Ex(nn) = Errx;
Ez(nn) = Errz;
DH(nn) = h;

end



% an increasing step size should result in an increase in error
% need to run for more years to see if this will work out

if convergenceTest

subplot(1,2,1)
loglog(DH, Ex, 'ro', 'LineWidth',1.5, 'MarkerSize',8); axis tight; box on; hold on
loglog(DH, Ex(1).*[1,1/2,1/4].^1, 'k-','LineWidth', 0.7);
loglog(DH, Ex(1).*[1,1/2,1/4].^2, 'k-','LineWidth', 0.9);
loglog(DH, Ex(1).*[1,1/2,1/4].^3, 'k-','LineWidth', 1.1);
loglog(DH, Ex(1).*[1,1/2,1/4].^4, 'k-','LineWidth', 1.3);
loglog(DH, Ex(1).*[1,1/2,1/4].^5, 'k-','LineWidth', 1.5);
xlabel('Step size', 'FontSize',18)
ylabel('Numerical error', 'FontSize',18)
title('Numerical Convergence in Space x', 'FontSize',20)


subplot(1,2,2)
loglog(DH, Ez, 'ro', 'LineWidth',1.5, 'MarkerSize',8); axis tight; box on; hold on
loglog(DH, Ez(1).*[1,1/2,1/4].^1, 'k-','LineWidth', 0.7);
loglog(DH, Ez(1).*[1,1/2,1/4].^2, 'k-','LineWidth', 0.9);
loglog(DH, Ez(1).*[1,1/2,1/4].^3, 'k-','LineWidth', 1.1);
loglog(DH, Ez(1).*[1,1/2,1/4].^4, 'k-','LineWidth', 1.3);
loglog(DH, Ez(1).*[1,1/2,1/4].^5, 'k-','LineWidth', 1.5);
xlabel('Step size', 'FontSize',18)
ylabel('Numerical error', 'FontSize',18)
title('Numerical Convergence in Space z', 'FontSize',20)

end



% TO DO:
% compare column vals with D site borehole placement

% validation test - closed boundaries:
%   no source term, closed boundaries
%   check T is constant through time
%   check linear gradient consistent through time


% parameter tests:
% use minimum values, medians or maxima from excel file
% see how borehole placement changes
% observe how impactful changing various rock types are
% best or worst case scenario

% spatial distributions of variables
% - where do we see high T and where is high Hr - do they match up

% convergence test:
% run overnight and check results 