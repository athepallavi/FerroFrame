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
Lx=1.6e-9; %
Ly=1.6e-9;
Nx=50; %# of intervals in x  directions
Ny=50;  %# of intervals in y directions
nx=Nx+1; %# of gridpoints in x directions including boundaries
ny=Ny+1; %# of gridpoints in y directions including boundaries
hx=2*Lx/Nx; %grid size in x directions
hy=2*Ly/Ny; %grid size in y directions
x=-Lx + (0:Nx)*hx; %x values on the grid
y=-Ly + (0:Ny)*hy; %y values on the grid
[X,Y]=meshgrid(x,y);

%% boundary_index = [bottom, left, top, right]
boundary_index=[ 1:nx, 1:nx:1+(nx-1)*nx, ...
1+(nx-1)*nx:nx*nx, nx:nx:nx*nx ];

diagonals = [4*ones(nx*ny,1), -ones(nx*ny,4)];
A=spdiags(diagonals,[0 -1 1 -nx nx], nx*ny, nx*ny);
B=spdiags([-1 -1], [-(nx*ny-nx) nx*ny-ny], nx*ny,nx*ny);
A=A+B;

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

% Create the 2D Gaussian function for p
p = A * exp(-((X - x0).^2 / (2 * sigma_x^2) + (Y - y0).^2 / (2 * sigma_y^2)));
% p=randi([0, 1], nx, ny) * 2 - 1;

% Scale by ps if needed
% p = ps * p;
p=ps*p;
% Initial conditions
Phi = zeros(nx, ny); % Initialize potential grid

%% Gauss-Seidel Iteration Parameters
tolerance = 1e-6;  % Convergence tolerance
max_iter = 5000;    % Maximum number of iterations
error = 1;          % Initial error
iter = 0;           % Iteration counter
del_t=1e-12;
% Gauss-Seidel Method Loop
while error > tolerance && iter < max_iter
    Phi_old = Phi; % Save old potential to calculate the error

    for i = 2:nx-1
        for j = 2:ny-1
            % Update each grid point using Gauss-Seidel method
            Phi(i,j) = (1/4) * (Phi(i-1,j) + Phi(i+1,j) + Phi(i,j-1) + Phi(i,j+1) - (hx/(2*epsilon_HfO)) * (p(i-1,j) - p(i+1,j) -p(i,j-1) + p(i,j+1)));
        end
    end

    % Enforce boundary conditions (assuming Dirichlet boundaries)
    % Phi(:,1) = 0;   % Left boundary
    % Phi(:,ny) = 0;  % Right boundary
    Phi(1,:) = 10;   % Bottom boundary
    Phi(nx,:) = Phi(nx-1,:)*exp(-0.002);  % Top boundary
    % Left boundary (i=1) takes values from the right boundary (i=nx)
    for j = 2:ny-1
        Phi(j,1) = (1/4) * (Phi(j,ny) + Phi(j,2) + Phi(j-1,1) + Phi(j+1,1) -  (hx/(2*epsilon_HfO)) * (p(j,ny) - p(j,2) - p(j-1,1) + p(j+1,1))); % Left boundary
        Phi(j,ny) = (1/4) * (Phi(j,ny-1) + Phi(j,1) + Phi(j-1,ny) + Phi(j+1,ny) -  (hx/(2*epsilon_HfO)) * (p(j,ny-1) - p(j,1) - p(j-1,ny) + p(j+1,ny))); % Right boundary
    end
    % Compute the error as the maximum difference between new and old values
    error = max(max(abs(Phi - Phi_old)));
    iter = iter + 1;
end

E_i=Phi/hx;
n=nx; m=ny;

p_i=p;
P_mod=p_i;
        P_old=P_mod;
        p_i1=circshift(p_i,2);
        p_i2=circshift(p_i,2,2);
        P_mod1=p_i1;
        P_mod2=p_i2;

        P_mod(2:end-1,2:end-1)=(-(alpha*p_i(2:end-1,2:end-1)+beta*p_i(2:end-1,2:end-1).^3 +gm*p_i(2:end-1,2:end-1).^5+k1*(2*p_i(2:end-1,2:end-1)-p_i(1:end-2,2:end-1)-p_i(3:end,2:end-1))+k2*(2*p_i(2:end-1,2:end-1)-p_i(2:end-1,1:end-2)-p_i(2:end-1,3:end)))-E_i(2:end-1,2:end-1))*(del_t/gamma)+p_i(2:end-1,2:end-1);
        P_mod1(2:end-1,2:end-1)=(-(alpha*p_i1(2:end-1,2:end-1)+beta*p_i1(2:end-1,2:end-1).^3+gm*p_i1(2:end-1,2:end-1).^5+k1*(2*p_i1(2:end-1,2:end-1)-p_i1(1:end-2,2:end-1)-p_i1(3:end,2:end-1))+k2*(2*p_i1(2:end-1,2:end-1)-p_i1(2:end-1,1:end-2)-p_i1(2:end-1,3:end)))-E_i(2:end-1,2:end-1))*(del_t/gamma)+p_i1(2:end-1,2:end-1);
        P_mod2(2:end-1,2:end-1)=(-(alpha*p_i2(2:end-1,2:end-1)+beta*p_i2(2:end-1,2:end-1).^3+gm*p_i2(2:end-1,2:end-1).^5+k1*(2*p_i2(2:end-1,2:end-1)-p_i2(1:end-2,2:end-1)-p_i2(3:end,2:end-1))+k2*(2*p_i2(2:end-1,2:end-1)-p_i2(2:end-1,1:end-2)-p_i2(2:end-1,3:end)))-E_i(2:end-1,2:end-1))*(del_t/gamma)+p_i2(2:end-1,2:end-1);

        p_i=P_mod;
        p_i(1,:)=P_mod1(3,:);
        p_i(m,:)=P_mod1(2,:);
        p_i(:,1)=P_mod2(:,3);
        p_i(:,n)=P_mod2(:,2);
        p_i(1,1)=(-(alpha*p_i(1,1)+beta*p_i(1,1)^3+gm*p_i(1,1)^5+k1*(2*p_i(1,1)-p_i(m,1)-p_i(2,1))+k2*(2*p_i(1,1)-p_i(1,n)-p_i(1,2)))-E_i(1,1))*(del_t/gamma)+p_i(1,1);
        p_i(1,n)=(-(alpha*p_i(1,n)+beta*p_i(1,n)^3+gm*p_i(1,n)^5+k1*(2*p_i(1,1)-p_i(m,n)-p_i(2,n))+k2*(2*p_i(1,n)-p_i(1,n-1)-p_i(1,1)))-E_i(1,n))*(del_t/gamma)+p_i(1,n);
        p_i(m,1)=(-(alpha*p_i(m,1)+beta*p_i(m,1)^3+gm*p_i(m,1)^5+k1*(2*p_i(m,1)-p_i(m-1,1)-p_i(1,1))+k2*(2*p_i(m,1)-p_i(m,n)-p_i(m,2)))-E_i(m,1))*(del_t/gamma)+p_i(m,1);
        p_i(m,n)=(-(alpha*p_i(m,n)+beta*p_i(m,n)^3+gm*p_i(m,n)^5+k1*(2*p_i(m,n)-p_i(m-1,n)-p_i(1,n))+k2*(2*p_i(m,n)-p_i(m,n-1)-p_i(m,1)))-E_i(m,n))*(del_t/gamma)+p_i(m,n);

fprintf('Converged after %d iterations with error %e\n', iter, error);

% Plot the potential
figure;
surf(X, Y, Phi);
title('Electric Potential \Phi');
xlabel('X [m]');
ylabel('Y [m]');
zlabel('\Phi [V]');
figure;
surf(X, Y, p);
figure;
surf(X, Y, p_i);