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

% set time step size (???)
% dt = CFL * min((h/2)/u0,(h/2)^2/k0); % time step [s]

% set up index arrays for boundary conditions

ix = [   1:Nx   ];          % insulating x boundaries

iz = [    1:Nz   ];          % isothermal start boundary, Nz+1 will be changed based on Nz value



% set initial condition for temperature at cell centres

Tx = Ttop + (Tbot-Ttop) ./ Nx .* Xc;
Tz   = Ttop + (Tbot-Ttop) ./ Nz .* Zc;  % initialise T array on linear gradient

Tin = T;                                         % store initial condition
Ta  = T;                                         % initialise analytical solution


%*****  Solve Model Equations

t = 0;  % initial time [s]
tau = 0;  % initial time step count

while t <= tend

    % increment time and step count
    t = t+dt;
    tau = tau+1;
    
    % 4th-order Runge-Kutta time integration scheme
    % NEED TO MAKE THIS 2D !!!
    R1 = diffusion(Tx, Tz, k0, h, ix, iz);
    R2 = diffusion(T + R1*dt/2, k0, h, ix, iz);
    R3 = diffusion(T + R2*dt/2, k0, h, ix, iz);
    R4 = diffusion(T + R3*dt, k0, h, ix, iz);
    
    % calculate 
    source = Hr ./ rho ./ Cp;
    % 
    % T = T + dt*(R1 + 2*R2 + 2*R3 + R4)/6 + source;
    % 
    % 
    % % get analytical solution at time t
    % 
    % 
    % % plot model progress
    % if ~mod(tau,nop)
    %     makefig(xc,zc,T);
    % end

end



% Function to calculate diffusion rate
function dTdt = diffusion(f, k0, h, ix, iz)

% calculate heat flux by diffusion
qx = - diff( k0(:, ix) ) .* diff( f(:, ix) )/h;
qz = - diff( k0(iz, :) ) .* diff( f(iz, :) )/h;

% disp(size(qx));

% calculate flux balance for rate of change (2nd derivative)
dTdt = - diff(qx + qz)/h;

end



% Function to make output figure
function makefig(x, z, T)

clf; 

% plot temperature 

imagesc(x,z,T); axis equal tight; colorbar; hold on
contour(x,z,T,[100,150,200],'k');

ylabel('z [m]','FontSize',15)
title('Temperature [C]','FontSize',17)


end