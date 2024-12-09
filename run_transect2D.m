%***** RUN 2D MODEL FROM IMAGE ***********************************

% clear workspace
clear all; close all; %clc;

% convergence
NN = [50,100,150];

for nn = 1:length(NN)
disp(['Iterating...', nn])

% load model setup from image, interpolate to target grid size

% domain width (must correspond to width of image) [m]
W       = 16e3;     

% number of pixels in z direction (number of rows) 
Nx = NN(nn);

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
%  #        sigma       rho         Cp           Hr
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
k0 = 10^3 * sigma ./ rho ./ Cp;


% calculate source term
source = Hr ./ rho ./ Cp;


% remaining model parameters
dTdz_top = 0;                      % flux at top
dTdz_bot = 35/1000;                % flux at bottom


yr    = 3600*24*365;  % seconds per year [s]
tend = 1e4*yr;

CFL   = 1/5;         % Time step limiter
nop   = 5000;          % output figure produced every 'nop' steps


T0 = 5; % surface air temperature


dTdz_boundaries = [0, 35/1000];



plotAnimation = false;

% *****  RUN MODEL
run('transect2D.m');

Ex(nn) = Errx;
Ez(nn) = Errz;
DH(nn) = h;

end



% an increasing step size should result in an increase in error
% need to run for more years to see if this will work out
convergenceTest = true;

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


% parameter tests:
% use minimum values, medians or maxima from excel file
% see how borehole placement changes
% observe how impactful changing various rock types are
% best or worst case scenario

% spatial distributions of variables
% - where do we see high T and where is high Hr - do they match up

% convergence test:
% run overnight and check results 