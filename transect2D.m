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

T = T0 + dTdz_boundaries(2)*Zc;  % initialise T array (based on what emily wrote)

Tin = T;                                         % store initial condition
Ta  = T;                                         % initialise analytical solution


%*****  Solve Model Equations

t = 0;  % initial time [s]
tau = 0;  % initial time step count

while t <= tend

    % increment time and step count
    t = t+dt;
    tau = tau+1;
    
    R1 = diffusion2D(T, k0, h, dTdz_boundaries(2));
    % 4th-order Runge-Kutta time integration scheme
    % NEED TO MAKE THIS 2D !!!
    % R1 = diffusion(dTdz_boundaries(1), dTdz_boundaries(2), k0);
    % R2 = diffusion(T + R1*dt/2, k0, h, ix, iz);
    % R3 = diffusion(T + R2*dt/2, k0, h, ix, iz);
    % R4 = diffusion(T + R3*dt, k0, h, ix, iz);
    
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
% Dtdzt, dTdzb = top and bottom heat flux
% k0 diffusivity
function dTdt = diffusion(dTdzt, dTdzb, k0)

% calculate heat flux by diffusion
% qz = - diff( k0(iz, :) ) .* diff( f(iz, :) )/h;
% qx = - diff( k0(:, ix) ) .* diff( f(:, ix) )/h;

% use dTdz to 
qz(1,:) = k0(1,:)*dTdzt;
qz(end,:) = k0(1,:)*dTdzb;

size(qz(1,:))
size(qz(end,:))
size(qz)

% qx(end,:) = k0(end,:)*

% disp(size(qx));

% calculate flux balance for rate of change (2nd derivative)
% dTdt = - diff(qx + qz)/h;

end


% SUGGESTION FROM CHATGPT
function dTdt = diffusion2D(f, k0, h, dTdzb)

% Get the size of the temperature grid
[nz, nx] = size(f);

% Preallocate arrays for fluxes
qx = zeros(nz, nx-1); % Flux in the x-direction
qz = zeros(nz-1, nx); % Flux in the z-direction

% Calculate heat flux in the x-direction (horizontal)
for j = 1:nz
    qx(j, :) = -k0(j, 1:end-1) .* diff(f(j, :)) / h;
end

% Calculate heat flux in the z-direction (vertical)
for i = 1:nx
    qz(:, i) = -k0(1:end-1, i) .* diff(f(:, i)) / h;
end

% Apply boundary condition at the bottom (z = depth)
% qz(end, :) corresponds to the flux at the bottom
qz(end+1, :) = k0(end, :) * dTdzb;

% Apply zero-flux boundary condition at the top (z = 0)
% No heat flux means qz(1, :) is effectively 0
qz = [zeros(1, nx); qz];

% Calculate flux divergence
dqx_dx = [qx(:, 1), diff(qx, 1, 2), -qx(:, end)] / h; % Flux change in x-direction
dqz_dz = diff(qz, 1, 1) / h; % Flux change in z-direction

% Rate of temperature change
dTdt = dqx_dx + dqz_dz;

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