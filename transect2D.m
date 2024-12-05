% model for finding heat gradients through transect defined in section.tiff
% 2D diffusion-source heat transport


%*****  Initialise Model Setup

% coordinate vector for cell centre positions [m]
x_cells = h/2:h:W-h/2;    
z_cells = h/2:h:D-h/2;

% coordinate vector for cell face positions [m]
x_faces = 0:h:W;            
z_faces = 0:h:D;

% create grid of coordinates
[Xc,Zc] = meshgrid(x_cells, z_cells);


% set up index arrays for boundary conditions           (????)

ix = [   1, 1:Nx, Nx  ];          % insulating x boundaries

iz = [   1, 1:Nz, Nz  ];          % isothermal start boundary



% set initial condition for temperature at cell centres

T = T0 + dTdz_boundaries(2)*Zc;  % initialise T array

% Tin = T;                       % store initial condition
% Ta  = T;                       % initialise analytical solution
    

% set time step size
dt = CFL * (h/2)^2 / max(k0, [], "all"); % time step [s]

% indices of air coordinates
air = units == 9;

%*****  Solve Model Equations

t = 0;  % initial time [s]
tau = 0;  % initial time step count

while t <= tend

    T(air) = T0;

    % increment time and step count
    t = t+dt;
    tau = tau+1;
    
    % 4th-order Runge-Kutta time integration scheme
    R1 = diffusion(T, k0, h, ix, iz, dTdz_boundaries(2)) + source;
    R2 = diffusion(T + R1*dt/2, k0, h, ix, iz, dTdz_boundaries(2)) + source;
    R3 = diffusion(T + R2*dt/2, k0, h, ix, iz, dTdz_boundaries(2)) + source;
    R4 = diffusion(T + R3*dt, k0, h, ix, iz, dTdz_boundaries(2)) + source;


    % numerical solution for T
    T = T + (R1 + 2*R2 + 2*R3 + R4)*dt/6 ;


    % get analytical solution at time t


    % plot model progress
    if ~mod(tau,nop)
        makefig(x_cells,z_cells,T);
    end

end



% Function to calculate diffusion rate
function dTdt = diffusion(f, k0, h, ix, iz, base_flux)

% move k0 values to cell centres
kz = ( k0(iz(1:end-1), :) + k0(iz(2:end), :) )/2;
kx = ( k0(:, ix(1:end-1)) + k0(:, ix(2:end)) )/2;


% calculate heat flux by diffusion
qz = - kz .* diff( f(iz,:), 1, 1) /h;
qx = - kx .* diff( f(:,ix), 1, 2) /h;

qz(end,:) = - k0(end,:)*base_flux;

% calculate flux balance for rate of change (2nd derivative)
dTdt = - ( diff(qz,1,1) + diff(qx,1,2) ) /h;

end



% Function to make output figure
function makefig(x, z, T)

clf; 

% plot temperature 

imagesc(x,z,T); axis equal tight; colorbar; hold on

% contour(x,z,T,[100,150,200],'k');

ylabel('z [m]','FontSize',15)
xlabel('x [m]','FontSize',15)
title('Temperature [C]','FontSize',17)

drawnow;

end