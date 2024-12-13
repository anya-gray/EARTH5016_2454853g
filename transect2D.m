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
% if linear
ix = [   1, 1:Nx, Nx  ];          % insulating boundaries

iz = [   1, 1:Nz, Nz  ];          


% set initial condition for temperature at cell centres

% T = T0 + dTdz_boundaries(2)*Zc;  % initialise temperature array

% Tin = T;                       % store initial condition
% Ta  = T;                       % initialise analytical solution (maybe divide by k0 ??)
    

% set time step size
dt = CFL * (h/2)^2 / max(k0(:)); % time step [s]

% indices of air coordinates
air = units == 9;

%*****  Solve Model Equations

t = 0;  % initial time [s]
tau = 0;  % initial time step count


% temperature evolution for:
% - CLOSED, INSULATING BOUNDARIES
% - NO BASAL HEAT FLUX
% - GAUSSIAN INITIALISED T  
if gaussian
    shape = 'gaussian';
    % initialise gaussian temperature array
    T = T0 + Tpeak .* exp( -( ( (Xc - W/2).^2 + (Zc - D/2).^2 ) / (4*Twidth^2) ) );
    k0 = k0(1,1);  % k0 is spatially independent => change to a scalar value
    
    while t <= t_end
        % increment time and step count
        t = t+dt;
        tau = tau+1;
        
        % 4th-order Runge-Kutta time integration scheme
        R1 = diffusion(shape, T, k0, h, ix, iz, dTdz_boundaries(2));
        R2 = diffusion(shape, T + R1*dt/2, k0, h, ix, iz, dTdz_boundaries(2));
        R3 = diffusion(shape, T + R2*dt/2, k0, h, ix, iz, dTdz_boundaries(2));
        R4 = diffusion(shape, T + R3*dt, k0, h, ix, iz, dTdz_boundaries(2));
    
        % numerical solution for T
        T = T + (R1 + 2*R2 + 2*R3 + R4)*dt/6;


        % *** calculate analytical solution at T (with help from ChatGPT)

        % Effective Gaussian width - changes with t due to diffusion
        w_eff = sqrt(Twidth^2 + 4 * k0 * t); 

        % Analytical solution
        Ta = T0 + (Tpeak / (w_eff.^2)) * ...
               exp(-((Xc - W/2).^2 + (Zc - D/2).^2) / w_eff.^2);
    

        % plot model progress
        if plotAnimation
            if ~mod(tau,nop)
                makefig(x_cells, z_cells, T, t, yr);
            end
        end

    end

    % calculate error of final T
    Errz = norm(T-Ta,1)./norm(Ta,1);
    Errx = norm(T-Ta,2)./norm(Ta,2);
    
    disp(' ');
    disp(['Numerical error on x = ', num2str(Errx)]);
    disp(['Numerical error on z = ', num2str(Errz)]);
    disp(' ');


end



% temperature evolution for linear initialisation
if linear
    shape = 'linear';
    
    % initialise linear temperature array
    T = T_air + dTdz_boundaries(2)*Zc;  
    
    % initialise energy to plot
    if validation
        thermal_energy = [];
        temperature = [];
        factor = [];
    end
    
    while t <= t_end
        % disp(t)
        
        if ~validation
            % enfore air temperature at all timesteps
            T(air) = T_air;
        end
        
        % save the total thermal energy at time t
        if validation
            thermal_energy(end+1) = sum( T(:) ) * energy_factor;
            % temperature(end+1) = sum( T(:) );
            % factor(end+1) = sum( rho(:).*Cp(:)*h*h );
        end
    
        % increment time and step count
        t = t+dt;
        tau = tau+1;
        
        % 4th-order Runge-Kutta time integration scheme
        R1 = diffusion(shape, T, k0, h, ix, iz, dTdz_boundaries(2), validation);
        R2 = diffusion(shape, T + R1*dt/2, k0, h, ix, iz, dTdz_boundaries(2), validation);
        R3 = diffusion(shape, T + R2*dt/2, k0, h, ix, iz, dTdz_boundaries(2), validation);
        R4 = diffusion(shape, T + R3*dt, k0, h, ix, iz, dTdz_boundaries(2), validation);
    
        % numerical solution for T
        if validation
            T = T + (R1 + 2*R2 + 2*R3 + R4)*dt/6; % pretend source is 0
        else
            T = T + (R1 + 2*R2 + 2*R3 + R4)*dt/6 + source;
        end
    
        % plot model progress
        if plotAnimation
            if ~mod(tau,nop)
                makefig(x_cells, z_cells, T, t, yr);
            end
        end
    end

    % move out of this loop eventually
    if validation
        clf; 
        % subplot(1,3,1)
        % plot(temperature)
        % xlabel('Time', 'FontSize',18)
        % ylabel('Total Temperature', 'FontSize',18)
        % title('Temperature in the System', 'FontSize',20)
        % 
        % subplot(1,3,2)
        % plot(factor)
        % xlabel('Time', 'FontSize',18)
        % ylabel('rho*Cp*dx*dx', 'FontSize',18)
        % title('That other factor the System', 'FontSize',20)
        % 
        % subplot(1,3,3)
        plot(thermal_energy)
        xlabel('Time', 'FontSize',18)
        ylabel('Total Thermal Energy', 'FontSize',18)
        title('Thermal Energy in the System', 'FontSize',20)
    end
end


% Function to calculate diffusion rate
function dTdt = diffusion(shape, f, k0, h, ix, iz, base_flux, validation)

    switch shape
        case 'linear'
            % move k0 values to cell centres
            kz = ( k0(iz(1:end-1), :) + k0(iz(2:end), :) )/2;
            kx = ( k0(:, ix(1:end-1)) + k0(:, ix(2:end)) )/2;
        
        case 'gaussian'
            % no need to move the k0s in this case as they are homogeneous (?)
            kz = k0;
            kx = k0;
    
        case 'test'
            % move k0 values to cell centres
            kz = ( k0(iz(1:end-1), :) + k0(iz(2:end), :) )/2;
            kx = ( k0(:, ix(1:end-1)) + k0(:, ix(2:end)) )/2;
    end


    % calculate heat flux by diffusion
    qz = - kz .* diff( f(iz,:), 1, 1) /h;
    qx = - kx .* diff( f(:,ix), 1, 2) /h;
    
    
    switch shape
    
        
        case 'linear'
            % ensure basal heat flux in standard linear gradient case
            if ~validation
                qz(end,:) = - k0(end,:)*base_flux;
            end
            
            % close boundaries for validation test
            % if validation
            %     qz(1, :) = 0;
            %     qz(end, :) = 0;
            % 
            %     qx(:, 1) = 0;
            %     qx(:, end) = 0;
            % end
    
        % ensure 0-flux boundaries
        % case 'gaussian'
        %     qz(1, :) = 0;
        %     qz(end, :) = 0;
        %     qx(:, 1) = 0;
        %     qx(:, end) = 0;
        % 
        % case 'test'
        %     qz(1, :) = 0;
        %     qz(end, :) = 0;
        %     qx(:, 1) = 0;
        %     qx(:, end) = 0;
    end

% calculate flux balance for rate of change (2nd derivative)
dTdt = - ( diff(qz, 1, 1) + diff(qx, 1, 2) ) /h;

end



% Function to make output figure
function makefig(x, z, T, t, yr)

clf; 

% plot temperature 

imagesc(x,z,T); axis equal tight; colorbar; hold on

[C,h] = contour(x,z,T,[40, 90, 150],'k');
clabel(C, h, 'Fontsize',12,'Color', 'r')

ylabel('z (m)','FontSize',15, 'FontName','Times New Roman')
ylabel(colorbar, 'Temperature (\circC)', 'FontName','Times New Roman')
xlabel('x (m)','FontSize',15, 'FontName','Times New Roman')
title(['Temperature Distribution (\circC) at ', num2str(floor(t/yr)), ' years.'], 'FontSize',17, 'FontName','Times New Roman')

drawnow;

end