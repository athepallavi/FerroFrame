% Parameters
q = 1.60217662e-19;        % Elementary charge in C
epsilon_0 = 8.854187817e-12; % Vacuum permittivity in F/m
epsilon_r = 11.7;           % Relative permittivity for Silicon
epsilon = epsilon_r * epsilon_0; % Permittivity of the semiconductor
k_B = 1.380649e-23;        % Boltzmann constant in J/K
T = 300;                   % Temperature in K
V_T = k_B * T / q;         % Thermal voltage at 300 K
Lx = 1e-6;                 % Length in x-direction (1 micron)
Ly = 1e-6;                 % Length in y-direction (1 micron)
Nx = 100;                  % Number of grid points in x-direction
Ny = 100;                  % Number of grid points in y-direction
dx = Lx / (Nx - 1);        % Grid spacing in x-direction
dy = Ly / (Ny - 1);        % Grid spacing in y-direction
tolerance = 1e-5;          % Convergence tolerance for potential
max_iterations = 1000;     % Maximum number of iterations

% Fermi-Dirac parameters
N_c = 1e19;                % Effective density of states in conduction band (cm^-3)
N_v = 1e19;                % Effective density of states in valence band (cm^-3)
E_c = 0.1;                 % Conduction band energy (eV)
E_v = -0.1;                % Valence band energy (eV)
E_f = 0.05;                % Fermi level (eV)

% Initial potential (e.g., 0 V throughout the domain)
phi = zeros(Nx, Ny);

% Ionized donor and acceptor densities (assuming undoped semiconductor)
N_d_plus = zeros(Nx, Ny);   % Ionized donor concentration
N_a_minus = zeros(Nx, Ny);  % Ionized acceptor concentration

% Boundary conditions (Dirichlet)
phi(:, 1) = 1;            % Left boundary (0 V)
phi(:, end) = 0.4;           % Right boundary (0 V)
phi(1, :) = phi(end,:);             % Bottom boundary (0 V)
phi(end, :) = phi(end-1,:);           % Top boundary (0 V)

% Iterative solver using Gauss-Seidel method
for iter = 1:max_iterations
    phi_old = phi;
    
    % Update potential using finite difference discretization (Gauss-Seidel)
    for i = 2:Nx-1
        for j = 2:Ny-1
            % Calculate electron density using Fermi-Dirac distribution
            n_e = N_c * 1 / (1 + exp((E_c - q * phi(i,j)) / (k_B * T)));
            
            % Calculate hole density using Fermi-Dirac distribution
            n_p = N_v * 1 / (1 + exp((q * phi(i,j) - E_v) / (k_B * T)));
            
            % Calculate charge density (rho)
            rho = q * (n_p - n_e + N_d_plus(i,j) - N_a_minus(i,j));
            
            % Poisson's equation update (finite difference discretization)
            phi(i,j) = (1 / (2 * (1/dx^2 + 1/dy^2))) * ...
                       ((phi(i+1,j) + phi(i-1,j)) / dx^2 + ...
                        (phi(i,j+1) + phi(i,j-1)) / dy^2 - ...
                        rho / epsilon);
        end
    end
    
    % Check for convergence
    if max(max(abs(phi - phi_old))) < tolerance
        disp(['Converged after ', num2str(iter), ' iterations']);
        break;
    end
end

% Plot the resulting potential
figure;
surf(phi);
title('Electric Potential \phi (V)');
xlabel('x');
ylabel('y');
zlabel('Potential (V)');
