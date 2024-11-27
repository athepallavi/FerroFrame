clc;
clear;
epsilon_siO2 = 4*8.854*1e-12; %Permittivity of SiO2 [F/cm]
epsilon_si = 11.68*8.854*1e-12; %Permittivity of Si [F/cm]
epsilon_HfO=24*8.854*1e-12;
kB = 1.380e-23; %Boltzmann constant [J/K]
q = 1.6e-19; %Elementary charge [C]
tox = 10e-9; %The thickness of SiO2 [m]
tfe = 10e-9;    % ferroelectric thickness (m)
Cox = epsilon_siO2/tox; %[F/m2]
E0 = 1; k1 = 1e-10; alpha = -2.5e9; beta = 6e10; k2 = 1e-11; gm = 1.5e11; gamma = 100;
omega = 2*pi; % Define the omega value

%%%%%%%%%% Define the square and grid parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid extension to 3.2 nm
Lx = 3.2e-9;
Ly = 3.2e-9;

% Adjust the number of grid points if necessary
Nx = 160;  % Now each nm has 50 grid points
Ny = 160;

% Compute new grid spacings
hx = 2*Lx / Nx;
hy = 2*Ly / Ny;

% Generate grid coordinates
x = -Lx + (0:Nx) * hx;
y = -Ly + (0:Ny) * hy;
[X, Y] = meshgrid(x, y);
nx=Nx+1;ny=Ny+1;
% Create a permittivity map (assuming permittivity changes at x = 1.6e-9)
epsilon_map = epsilon_HfO * ones(nx, ny);  % Default is HfO
for i = 1:nx
    for j = 1:ny
        if X(i, j) >1.6e-9
            epsilon_map(i, j) = epsilon_siO2;  % Change to SiO2 beyond 1.6 nm
        end
    end
end

%% Diffusion constant and time-step parameters
D=1/gamma;
dt=hx^2/(2*D); %borderline stability of FTCS scheme
alphaa=dt*D/hx^2; %equation parameter
nsteps=1000; %number of time steps

ps=2e-4;
%% Initialize the dipole moment array
%% Set initial conditions for p
p=2*rand(nx,ny)-1;
% p=ps*p;
% mu = 0.002; % Mean of the Gaussian distribution
% sigma = 1; % Standard deviation of the Gaussian distribution
% p = mu + sigma * randn(nx,ny); % Gaussian-distributed values for p
A = 1; % Amplitude of the Gaussian
x0 = 0; % Center of the Gaussian in the x-direction
y0 = 0; % Center of the Gaussian in the y-direction
sigma_x = Lx / 4; % Standard deviation in the x-direction (set to some fraction of the grid size)
sigma_y = Ly / 4; % Standard deviation in the y-direction (set to some fraction of the grid size)

% Initialize the dipole moment array with different conditions for each region
p = zeros(nx, ny);  % Initialize all to zero

% Define ferroelectric region and initialize only in that region
for i = 1:nx
    for j = 1:ny
        if X(i, j) <= 1.6e-9
            % Gaussian initial condition within ferroelectric region
            p(i, j) = ps * exp(-((X(i, j) - x0)^2 / (2 * sigma_x^2) + (Y(i, j) - y0)^2 / (2 * sigma_y^2)));
        end
    end
end
% Initialize the dipole moment array with different conditions for each region
p = zeros(nx, ny);  % Initialize all to zero

% Ferroelectric region boundary
ferro_region_limit = 1.6e-9;

% Define the width of each quarter in the ferroelectric region
quarter_width = ferro_region_limit / 4;

% Define ferroelectric region and initialize only in that region
for i = 1:nx
    for j = 1:ny
        if X(i, j) <= ferro_region_limit
            % Divide the ferroelectric region into four equal parts
            if X(i, j) <= -ferro_region_limit + quarter_width
                % First region (leftmost quarter): Positive polarization
                p(i, j) = ps;
            elseif X(i, j) > -ferro_region_limit + quarter_width && X(i, j) <= -ferro_region_limit + 2 * quarter_width
                % Second region: Negative polarization
                p(i, j) = -ps;
            elseif X(i, j) > -ferro_region_limit + 2 * quarter_width && X(i, j) <= -ferro_region_limit + 3 * quarter_width
                % Third region: Positive polarization
                p(i, j) = ps;
            else
                % Fourth region (rightmost quarter): Negative polarization
                p(i, j) = -ps;
            end
        end
    end
end
% p=ps*p;
% Initial conditions
Phi = zeros(nx, ny); % Initialize potential grid

%% Gauss-Seidel Iteration Parameters
tolerance = 1e-6;  % Convergence tolerance
max_iter = 5000;    % Maximum number of iterations
error = 1;          % Initial error
iter = 0;           % Iteration counter

% Gauss-Seidel Method Loop
% Gauss-Seidel Method with variable permittivity
while error > tolerance && iter < max_iter
    Phi_old = Phi;  % Save old potential

    % Update potentials using Gauss-Seidel method
    for i = 2:nx-1
        for j = 2:ny-1
            % Phi(i,j) = (1/4) * (Phi(i-1,j) + Phi(i+1,j) + Phi(i,j-1) + Phi(i,j+1) ...
            %     - (hx^2/epsilon_map(i,j)) * p(i,j));
            Phi(i,j) = (1/4) * (Phi(i-1,j) + Phi(i+1,j) + Phi(i,j-1) + Phi(i,j+1) - (hx/(2*epsilon_map(i,j))) * (p(i-1,j) - p(i+1,j) -p(i,j-1) + p(i,j+1)));
        end
    end

    % Apply boundary conditions
    Phi(1,:) = Phi(2,:);  % Bottom boundary
    Phi(nx,:) =Phi(nx-1,:); % Top boundary
    Phi(:,1) = 0;%Phi(:,2);  % Left boundary (periodic or symmetric)
    Phi(:,ny) = 0;%Phi(:,ny-1); % Right boundary (periodic or symmetric)
% for j = 2:ny-1
%         Phi(j,1) = (1/4) * (Phi(j,ny-1) + Phi(j,2) + Phi(j-1,1) + Phi(j+1,1) -  (hx/(2*epsilon_map(j,1))) * (p(j,ny-1) - p(j,2) - p(j-1,1) + p(j+1,1))); % Left boundary
%         Phi(j,ny-1) = (1/4) * (Phi(j,ny-2) + Phi(j,1) + Phi(j-1,ny-1) + Phi(j+1,ny) -  (hx/(2*epsilon_map(j,ny-1))) * (p(j,ny-2) - p(j,1) - p(j-1,ny-1) + p(j+1,ny-1))); % Right boundary
% end
    % Compute error and update iteration counter
    error = max(max(abs(Phi - Phi_old)));
    iter = iter + 1;
end

fprintf('Converged after %d iterations with error %e\n', iter, error);

% Plot the dipole moment distribution
figure;
surf(X, Y, p);
shading interp; 
title('Dipole Moment Distribution');
xlabel('x (nm)');
ylabel('y (nm)');
zlabel('Dipole Moment p');

% Plot the electric potential distribution
figure;
surf(X, Y, Phi);
shading interp; 
title('Electric Potential Distribution');
xlabel('x (nm)');
ylabel('y (nm)');
zlabel('\Phi (V)');