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


% set up index arrays for insulating boundary conditions
ix = [   1, 1:Nx, Nx  ];  
iz = [   1, 1:Nz, Nz  ];          


% set time step size
dt = CFL * (h/2)^2 / max(k0(:)); % time step [s]

% indices of air coordinates
air = units == 9;


%*****  Solve Model Equations

t = 0;  % initial time [s]
tau = 0;  % initial time step count


% temperature evolution for linear initialisation
if linear
    shape = 'linear';
    
    % initialise linear temperature array
    T = T_air + geotherm*Zc;  
    
    % initialise energy to plot
    if verification
        thermal_energy = [];
        time = [];

    end
    
    while t <= t_end

        if ~verification
            % enfore air temperature at all timesteps
            T(air) = T_air;
        end
        
        % save the total thermal energy at time t, and time for plotting
        if verification 
            thermal_energy(end+1) = (sum( T(:) ) * energy_factor);
            time(end+1) = t/yr;
        end
    
        % increment time and step count
        t = t+dt;
        tau = tau+1;
        
        % 4th-order Runge-Kutta time integration scheme
        R1 = diffusion(shape, T, k0, h, ix, iz, base_flux, verification);
        R2 = diffusion(shape, T + R1*dt/2, k0, h, ix, iz, base_flux, verification);
        R3 = diffusion(shape, T + R2*dt/2, k0, h, ix, iz, base_flux, verification);
        R4 = diffusion(shape, T + R3*dt, k0, h, ix, iz, base_flux, verification);
    
        % numerical solution for T, remove source for  verification test
        if verification
            T = T + (R1 + 2*R2 + 2*R3 + R4)*dt/6;     
        else
            T = T + (R1 + 2*R2 + 2*R3 + R4)*dt/6 + source;
        end

        % plot model progress
        if plotAnimation
            if ~mod(tau, plot_freq)
                makefig(x_cells, z_cells, T, t, yr);
            end
        end
    end
   
    % plot final result of simulation
    if plot_result
        makefig(x_cells, z_cells, T, t, yr);
    end
    
    % plot the total thermal energy through time
    if verification
        clf; 
        plot(time, thermal_energy, 'color', [0.5 0 0.5], 'LineWidth',2)
        ylim([2.716e20 2.717e20])           % hard coded for a nice image
        xlabel('Time (years)', 'FontSize', 18, 'FontName', 'Times New Roman')
        ylabel('Energy (J)', 'FontSize',18, 'FontName', 'Times New Roman')
        title('Total Thermal Energy in the System', 'FontSize',20, 'FontName', 'Times New Roman')
    end
end


% Function to calculate diffusion rate dTdt
function dTdt = diffusion(shape, f, k0, h, ix, iz, base_flux, verification)

    % find diffusivity in the linear or gaussian case
    switch shape
        case 'linear'
            % move k0 values to cell centres
            kz = ( k0(iz(1:end-1), :) + k0(iz(2:end), :) )/2;
            kx = ( k0(:, ix(1:end-1)) + k0(:, ix(2:end)) )/2;
        
        case 'gaussian'
            % no need to move the k0s in this case as they are homogeneous
            % (gaussian case below sends in k0 as a scalar)
            kz = k0;
            kx = k0;
    end


    % calculate heat flux by diffusion
    qz = - kz .* diff( f(iz,:), 1, 1) /h;
    qx = - kx .* diff( f(:,ix), 1, 2) /h;
    
    
    % ensure basal heat flux for the normal simulation
    switch shape
        case 'linear'
            if ~verification
                qz(end,:) = base_flux;
            end
    end

% calculate flux balance for rate of change (2nd derivative)
dTdt = - ( diff(qz, 1, 1) + diff(qx, 1, 2) ) /h;

end



% plot animation of temperature dynamics through time
function makefig(x, z, T, t, yr)

clf; 

% Temperature contour plot
imagesc(x,z,T); axis equal tight; colorbar; hold on

% contour lines at 40, 90 and 150 degC
[C,h] = contour(x,z,T,[40, 90, 150],'k');
clabel(C, h, 'Fontsize',15,'Color', 'r')

% axes labels and titles
ylabel('z (m)','FontSize',15, 'FontName','Times New Roman')
ylabel(colorbar, 'Temperature (\circC)', 'FontName','Times New Roman')
xlabel('x (m)','FontSize',15, 'FontName','Times New Roman')
title(['Temperature Distribution (\circC) at ', num2str(floor(t/yr)), ' years.'], 'FontSize',17, 'FontName','Times New Roman')

drawnow;        % animate

end



% ****************** Gaussian case (couldn't get this to fully work)
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
        R1 = diffusion(shape, T, k0, h, ix, iz, base_flux, verification);
        R2 = diffusion(shape, T + R1*dt/2, k0, h, ix, iz, base_flux, verification);
        R3 = diffusion(shape, T + R2*dt/2, k0, h, ix, iz, base_flux, verification);
        R4 = diffusion(shape, T + R3*dt, k0, h, ix, iz, base_flux, verification);
    
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
            if ~mod(tau, plot_freq)
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