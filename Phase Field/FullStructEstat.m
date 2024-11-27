clc;
clear;
epsilon_siO2 = 4*8.854*1e-12; %Permittivity of SiO2 [F/cm]
epsilon_si = 11.68*8.854*1e-12; %Permittivity of Si [F/cm]
epsilon_HfO=24*8.854*1e-12;
k_B = 1.380e-23; %Boltzmann constant [J/K]
q = 1.6e-19; %Elementary charge [C]
tox = 10e-9; %The thickness of SiO2 [m]
tfe = 10e-9;    % ferroelectric thickness (m)
Cox = epsilon_siO2/tox; %[F/m2]
E0 = 1; k1 = 1e-10; alpha = -2.5e9; beta = 6e10; k2 = 1e-11; gm = 1.5e11; gamma = 12;
omega = 2*pi; % Define the omega value
T=300; % room temp in K
V_T = k_B * T / q;         % Thermal voltage at 300 K

% Fermi-Dirac parameters
N_c = 1e25;                % Effective density of states in conduction band (cm^-3)
N_v = 1e25;                % Effective density of states in valence band (cm^-3)
E_c = 0.1;                 % Conduction band energy (eV)
E_v = -0.1;                % Valence band energy (eV)
E_f = 0.05;                % Fermi level (eV)



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
        if X(i, j) >0 && X(i,j)<1.6e-9 
            epsilon_map(i, j) = epsilon_siO2;  % Change to SiO2 beyond 1.6 nm
            nxSi=j+1;
        end
        if X(i, j) >=1.6e-9 
            
            epsilon_map(i, j) = epsilon_si;  % Change to SiO2 beyond 1.6 nm
        end
    end
end

% Ionized donor and acceptor densities (assuming undoped semiconductor)
N_d_plus = zeros(nx, ny);   % Ionized donor concentration
N_a_minus = zeros(nx, ny);  % Ionized acceptor concentration

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


% Ferroelectric region boundary
ferro_region_limit = 0;
Nxf=round((Lx+ferro_region_limit)/hx);
Nyf=Ny;
nxf=Nxf+1;
nyf=Nyf+1;
yf = -Ly + (0:Nyf) * hy;
xf=-Lx+ (0:Nxf) * hx;

[Xf, Yf]=meshgrid(xf,yf);
p = zeros(nxf, nyf);  % Initialize all to zero

% Calculate the quarter width within the ferroelectric region
quarter_width = (Ly) / 4;

% Define the x-coordinates for the division into quarters
quarter_coords = [-Ly, -Ly + 2*quarter_width, -Ly + 2*2 * quarter_width, -Ly + 2*3 * quarter_width, -Ly + ferro_region_limit];

% Iterate through each point in the grid
for i = 1:nxf
    for j = 1:nyf
        if Y(j, i) <= Ly
            % Check which quarter the point belongs to and assign the polarization
            if Y(j, i) <= quarter_coords(2)
                p(i, j) = ps;  % Positive polarization
            elseif Y(j, i) <= quarter_coords(3)
                p(i, j) = -ps;  % Negative polarization
            elseif Y(j, i) <= quarter_coords(4)
                p(i, j) = ps;  % Positive polarization
            else
                p(i, j) = -ps;  % Negative polarization
            end
        end
    end
end
% p(i,:)=p(i-1,:);

p(:,j)=p(:,j-1);
% p=ps*p;
% Initial conditions
Phi = zeros(nx, ny); % Initialize potential grid

%% Gauss-Seidel Iteration Parameters
tolerance = 1e-4;  % Convergence tolerance
max_iter = 1000;    % Maximum number of iterations
error = 1;          % Initial error
iter = 0;           % Iteration counter
del_x=hx;
% Gauss-Seidel Method Loop
% Gauss-Seidel Method with variable permittivity
while error > tolerance && iter < max_iter
    Phi_old = Phi;  % Save old potential

    % Update potentials using Gauss-Seidel method
    for i = 2:nxf-1
        for j = 2:nyf-1
            % Phi(i,j) = (1/4) * (Phi(i-1,j) + Phi(i+1,j) + Phi(i,j-1) + Phi(i,j+1) ...
            %     - (hx^2/epsilon_map(i,j)) * p(i,j));
            Phi(i,j) = (1/4) * (Phi(i-1,j) + Phi(i+1,j) + Phi(i,j-1) + Phi(i,j+1) - (hx/(2*epsilon_map(i,j))) * (p(i-1,j) - p(i+1,j) -p(i,j-1) + p(i,j+1)));
        end
    end
    for i = nxf+1:nxSi-1
            for j = 2:ny-1
                % Phi(i,j) = (1/4) * (Phi(i-1,j) + Phi(i+1,j) + Phi(i,j-1) + Phi(i,j+1) ...
                %     - (hx^2/epsilon_map(i,j)) * p(i,j));
                Phi(i,j) = (1/4) * (Phi(i-1,j) + Phi(i+1,j) + Phi(i,j-1) + Phi(i,j+1) );
            end
    end

    for i = nxSi:nx-1
        for j = 2:ny-1
            % Calculate electron density using Fermi-Dirac distribution
            n_e = N_c * 1 / (1 + exp((E_c - q * Phi(i,j)) / (k_B * T)));
            
            % Calculate hole density using Fermi-Dirac distribution
            n_p = N_v * 1 / (1 + exp((q * Phi(i,j) - E_v) / (k_B * T)));
            
            % Calculate charge density (rho)
            rho = q * (n_p - n_e + N_d_plus(i,j) - N_a_minus(i,j));
            
            % Poisson's equation update (finite difference discretization)
            Phi(i,j) = (1 / (2 * (1/del_x^2 + 1/del_x^2))) * ...
                       ((Phi(i+1,j) + Phi(i-1,j)) / del_x^2 + ...
                        (Phi(i,j+1) + Phi(i,j-1)) / del_x^2 - ...
                        rho / epsilon_si);
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
% Define electric field (triangular waveform)
tspan = linspace(0,1e-8,1e4); % Define the time span
tspan = tspan';
ts = tspan;
t = ts;
tfe=90e-9;
del_t = ts(2) - ts(1);
P_init = repmat(p, 1, 1, length(tspan));

E_tri = zeros(size(t));      % Initialize E_tri vector

t_half = t(end) / 2;         % Midpoint of time vector

VG=1;
for i = 1:length(t)
    if t(i) <= t_half
        E_tri(i) = -VG + 2*VG * (t(i) / t_half);
    else
        E_tri(i) = VG - 2*VG * ((t(i) - t_half) / t_half);
    end
end
V_G=E_tri;

error = 1;

% E1= E_tri*1e9;
del_x=hx;
gamma = 12;
p_in=p;
p_i = p_in;
PT1 = zeros(length(ts), nxf);
Phi_old = Phi;
E_p=zeros(nxf,nyf,length(ts));
for i1 = 1:length(ts)
    error = 1;
    Phi= zeros(nx, ny);  % Save old potential
    iter=0;
    while error > tolerance && iter < max_iter
        Phi_old=Phi;

        % Update potentials using Gauss-Seidel method
        for i = 2:nxf-1
            for j = 2:nyf-1
                % Phi(i,j) = (1/4) * (Phi(i-1,j) + Phi(i+1,j) + Phi(i,j-1) + Phi(i,j+1) ...
                %     - (hx^2/epsilon_map(i,j)) * p(i,j));
                Phi(i,j) = (1/4) * (Phi(i-1,j) + Phi(i+1,j) + Phi(i,j-1) + Phi(i,j+1) - (hx/(2*epsilon_map(i,j))) * (p(i-1,j) - p(i+1,j) -p(i,j-1) + p(i,j+1)));
            end
        end
        for i = nxf:nxSi-1
            for j = 2:ny-1
                % Phi(i,j) = (1/4) * (Phi(i-1,j) + Phi(i+1,j) + Phi(i,j-1) + Phi(i,j+1) ...
                %     - (hx^2/epsilon_map(i,j)) * p(i,j));
                Phi(i,j) = (1/4) * (Phi(i-1,j) + Phi(i+1,j) + Phi(i,j-1) + Phi(i,j+1) );
            end
        end

        for i = nxSi:nx-1
            for j = 2:ny-1
                % Calculate electron density using Fermi-Dirac distribution
                n_e = N_c * 1 / (1 + exp((E_c - q * Phi(i,j)) / (k_B * T)));

                % Calculate hole density using Fermi-Dirac distribution
                n_p = N_v * 1 / (1 + exp((q * Phi(i,j) - E_v) / (k_B * T)));

                % Calculate charge density (rho)
                rho = q * (n_p - n_e + N_d_plus(i,j) - N_a_minus(i,j));

                % Poisson's equation update (finite difference discretization)
                Phi(i,j) = (1 / (2 * (1/del_x^2 + 1/del_x^2))) * ...
                    ((Phi(i+1,j) + Phi(i-1,j)) / del_x^2 + ...
                    (Phi(i,j+1) + Phi(i,j-1)) / del_x^2 - ...
                    rho / epsilon_si);
            end
        end

        % Apply boundary conditions
        Phi(1,:) = 0;%Phi(2,:);  % Bottom boundary
        Phi(nx,:) =V_G(i1);  % Phi(nx-1,:); % Top boundary
        Phi(:,1) = Phi(:,2);  % Left boundary (periodic or symmetric)
        Phi(:,ny) = Phi(:,ny-1);  % Right boundary (periodic or symmetric)
        % for j = 2:ny-1
        %         Phi(j,1) = (1/4) * (Phi(j,ny-1) + Phi(j,2) + Phi(j-1,1) + Phi(j+1,1) -  (hx/(2*epsilon_map(j,1))) * (p(j,ny-1) - p(j,2) - p(j-1,1) + p(j+1,1))); % Left boundary
        %         Phi(j,ny-1) = (1/4) * (Phi(j,ny-2) + Phi(j,1) + Phi(j-1,ny-1) + Phi(j+1,ny) -  (hx/(2*epsilon_map(j,ny-1))) * (p(j,ny-2) - p(j,1) - p(j-1,ny-1) + p(j+1,ny-1))); % Right boundary
        % end
        % Compute error and update iteration counter
        error = max(max(abs(Phi - Phi_old)));
        iter = iter + 1;
    end
    E=Phi(1:nxf,:)*1e9;
    E_p(:,:,i1)=E;
    Psi(:,:,i1)=Phi;
    P_mod=p_i;
    p_i1=circshift(p_i,2);
    p_i2=circshift(p_i,2,2);
    P_mod1=p_i1;
    P_mod2=p_i2;
    P_mod(2:end-1,2:end-1)=(-(alpha*p_i(2:end-1,2:end-1)+beta*p_i(2:end-1,2:end-1).^3 +gm*p_i(2:end-1,2:end-1).^5+del_x*(k1*(2*p_i(2:end-1,2:end-1)-p_i(1:end-2,2:end-1)-p_i(3:end,2:end-1))+k2*(2*p_i(2:end-1,2:end-1)-p_i(2:end-1,1:end-2)-p_i(2:end-1,3:end)))/del_x^2)-E(2:end-1,2:end-1))*(del_t/gamma)+p_i(2:end-1,2:end-1);
    P_mod1(2:end-1,2:end-1)=(-(alpha*p_i1(2:end-1,2:end-1)+beta*p_i1(2:end-1,2:end-1).^3+gm*p_i(2:end-1,2:end-1).^5+del_x*(k1*(2*p_i1(2:end-1,2:end-1)-p_i1(1:end-2,2:end-1)-p_i1(3:end,2:end-1))+k2*(2*p_i1(2:end-1,2:end-1)-p_i1(2:end-1,1:end-2)-p_i1(2:end-1,3:end)))/del_x^2)-E(2:end-1,2:end-1))*(del_t/gamma)+p_i1(2:end-1,2:end-1);
    P_mod2(2:end-1,2:end-1)=(-(alpha*p_i(2:end-1,2:end-1)+beta*p_i2(2:end-1,2:end-1).^3+gm*p_i(2:end-1,2:end-1).^5+del_x*(k1*(2*p_i2(2:end-1,2:end-1)-p_i2(1:end-2,2:end-1)-p_i2(3:end,2:end-1))+k2*(2*p_i2(2:end-1,2:end-1)-p_i2(2:end-1,1:end-2)-p_i2(2:end-1,3:end)))/del_x^2)-E(2:end-1,2:end-1))*(del_t/gamma)+p_i2(2:end-1,2:end-1);
    p_i=P_mod;
    p_i(1,:)=P_mod1(3,:);
    p_i(nxf,:)=P_mod1(2,:);
    p_i(:,1)=P_mod2(:,3);
    p_i(:,nyf)=P_mod2(:,2);
    %  p_i(1,:)=P_mod(2,:);
    % p_i(nxf,:)=P_mod(nxf-1,:);
    % p_i(:,1)=P_mod(:,2);
    % p_i(:,nyf)=P_mod(:,nyf-1);
    p_i(1,1)=(-(alpha*p_i(1,1)+beta*p_i(1,1)^3+gm*p_i(1,1)^5+del_x*(k1*(2*p_i(1,1)-p_i(nxf,1)-p_i(2,1))+k2*(2*p_i(1,1)-p_i(1,nyf)-p_i(1,2)))/del_x^2)-E(1,1))*(del_t/gamma)+p_i(1,1);
    p_i(1,nyf)=(-(alpha*p_i(1,nyf)+beta*p_i(1,nyf)^3+gm*p_i(1,nyf)^5+del_x*(k1*(2*p_i(1,1)-p_i(nxf,nyf)-p_i(2,nyf))+k2*(2*p_i(1,nyf)-p_i(1,nyf-1)-p_i(1,1)))/del_x^2)-E(1,nyf))*(del_t/gamma)+p_i(1,nyf);
    p_i(nxf,1)=(-(alpha*p_i(nxf,1)+beta*p_i(nxf,1)^3+gm*p_i(nxf,1)^5+del_x*(k1*(2*p_i(nxf,1)-p_i(nxf-1,1)-p_i(1,1))+k2*(2*p_i(nxf,1)-p_i(nxf,nyf)-p_i(nxf,2)))/del_x^2)-E(nxf,1))*(del_t/gamma)+p_i(nxf,1);
    p_i(nxf,nyf)=(-(alpha*p_i(nxf,nyf)+beta*p_i(nxf,nyf)^3+gm*p_i(nxf,nyf)^5+del_x*(k1*(2*p_i(nxf,nyf)-p_i(nxf-1,nyf)-p_i(1,nyf))+k2*(2*p_i(nxf,nyf)-p_i(nxf,nyf-1)-p_i(nxf,1)))/del_x^2)-E(nxf,nyf))*(del_t/gamma)+p_i(nxf,nyf);
    p_i(:,nyf)=p_i(:,nyf-1);
    p_in = p_i;
    P_init(:,:,i1)=p_i;
    PT1(i1,:) = sum(p_in');
end
PT=-sum(PT1,2);

% Plot the dipole moment distribution
figure;
surf(Xf, Yf, p');
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
PT=-sum(PT1,2);
% end
% PT=flipud(PT2);
figure()
plot(ts,V_G*1e9,'black',linewidth=2)
xlabel('Time');
ylabel('Electric field');
title('Electric Field vs Time');
figure();
plot(ts, PT,'b',linewidth=2);
xlabel('Time');
ylabel('Dipole Moments');
title('Dipole Moment Evolution Over Time');
figure()
plot(V_G*1e9,PT,'r',linewidth=2)
xlabel('Electric field');
ylabel('Dipole Moments');
title('Dipole Moment vs Electric field');
figure()
surf(Xf,Yf,P_init(:,:,1)');
shading interp;