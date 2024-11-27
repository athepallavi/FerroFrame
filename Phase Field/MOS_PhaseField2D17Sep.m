clc
clear all
%Specify the physical constants
epsilon_siO2 = 4*8.854*1e-12; %Permittivity of SiO2 [F/cm]
epsilon_si = 11.68*8.854*1e-12; %Permittivity of Si [F/cm]
epsilon_HfO=24*8.854*1e-12;
kB = 1.380e-23; %Boltzmann constant [J/K]
q = 1.6e-19; %Elementary charge [C]
%Specify the environmental parameters
T = 300; %Room temperature[K]
ni = 1e16; %Intrinsic semiconductor carrier density [cm-3]
N_A = 1e17*1e6; %Acceptor concentration [cm-3]
tox = 10e-9; %The thickness of SiO2 [m]
V_FB = 0; %The flat-band voltage [V]
tfe = 10e-9;    % ferroelectric thickness (m)
Cox = epsilon_siO2/tox; %[F/m2]
% Ferroelectric Parameters
% E0 = 1; gamma = 1e-8; alpha = -2.5e9; beta = 6e10; k = 25; gm = 1.5e11;
E0 = 1; k1 = 1e-10; alpha = -2.5e9; beta = 6e10; k2 = 1e-11; gm = 1.5e11;gamma =100;
% E0 =1;  gamma =1e-8; alpha = -1; beta = 1; k2= 5;k1=5 ;
omega = 2*pi; % Define the omega value
% n = 50; % Define the number of dipoles
% m=50;
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
% ps=sqrt((-alpha)/beta);
del_tfe=2.8284e-10;
ps=2e-6;
m=Nx; n=Ny;
% Initialize the dipole moment array
% p = zeros(1, n); % Example initialization, modify as needed
% Set initial conditions for p
p=2*rand(m,n)-1;
% p=-1*ones(1,n);
p=ps*p;
% Time span for the simulation
tspan = linspace(0,1e-8,1e4); % Define the time span
tspan=tspan';
ts=tspan;
t=ts;
P_init=repmat(p,1,1, length(tspan));
del_t=ts(2)-ts(1);
Vt= kB*T/q;
Q0 = sqrt(2*epsilon_si*kB*T*ni^2/N_A);
E_tri = zeros(size(t));      % Initialize E_tri vector

t_half = t(end) / 2;         % Midpoint of time vector

for i = 1:length(t)
    if t(i) <= t_half
        E_tri(i) = -10 + 20 * (t(i) / t_half);
    else
        E_tri(i) = 10 - 20 * ((t(i) - t_half) / t_half);
    end
end
% E = 5*sin(omega * ts*1e6);
% E=E_tri*1e9;
p_i= P_init(1,:);
% del_tfe=2.8284e-10;
del_tfe=hx;
% del_tfe=5e-6;
% MOSFET Parameters
L = 1e-6; % channel length (m)
W = 1e-6; % channel width (m)
%Calculate the capacitance
phi_b = Vt*log(N_A/ni); % body-potential (V)
V_G = -0.5:0.01:3; %[V]
A=1e-10;%Area in cm2
V_G =E_tri/4;
%Initialize the vectors
V_length = length(V_G);
Qs = zeros(1,V_length);
psi = zeros(1,V_length);

% Numerical Solver parameters (Newton Method)
iter = 1e5;
tol = 1e-5;
VFE=zeros(m,n,length(V_G));
Qfe=zeros(1,length(V_G));
E=VFE;
p_i= P_init(:,:,1);
VF=VFE(:,:,1);

i=1;
%%%% Initialize phi in FE layer
E(:,:,i)=VFE(:,:,i)/del_tfe;
% VF(1,:)=0.1;
eps=1.e-5;
error=2*eps; ncount=0;
error1=error;
 % VF(1,:)=5;
 VF_mod=VF;
while (error > eps)
    ncount=ncount+1;
    % VF_mod=VF;
    VF_mod(2:end-1,2:end-1)=VF(1:end-2,2:end-1)/2+VF(2:end-1,1:end-2)/2+(p_i(3:end,2:end-1)-p_i(1:end-2,2:end-1)+p_i(2:end-1,3:end)-p_i(2:end-1,1:end-2))*(del_tfe/(4*epsilon_HfO))-Qs(i)*(del_tfe^2/(2*epsilon_HfO*5e-8));
    VF2=circshift(VF,2);
    VF_mod2=VF2;
    VF_mod2(2:end-1,2:end-1)=VF2(1:end-2,2:end-1)/2+VF2(2:end-1,1:end-2)/2+(p_i(3:end,2:end-1)-p_i(1:end-2,2:end-1)+p_i(2:end-1,3:end)-p_i(2:end-1,1:end-2))*(del_tfe/(4*epsilon_HfO))-Qs(i)*(del_tfe^2/(2*epsilon_HfO*5e-8));
    VF1=circshift(VF,2,2);
    VF_mod1=VF1;
    VF_mod1(2:end-1,2:end-1)=VF1(1:end-2,2:end-1)/2+VF1(2:end-1,1:end-2)/2+(p_i(3:end,2:end-1)-p_i(1:end-2,2:end-1)+p_i(2:end-1,3:end)-p_i(2:end-1,1:end-2))*(del_tfe/(4*epsilon_HfO))-Qs(i)*(del_tfe^2/(2*epsilon_HfO*5e-8));
    VF=VF_mod;
    % VF(1,:)=VF_mod2(3,:);
    VF(m,:)=VF_mod2(2,:);
    VF(:,1)=VF_mod1(:,3);
    VF(:,n)=VF_mod1(:,2);
    % VF(1,1)=VF(1,n)/2+VF(1,2)/2+(p_i(1,2)-p_i(1,n)+p_i(2,1)-p_i(n,1))*(del_tfe/(4*epsilon_HfO))-Qs(i)*(del_tfe^2/(2*epsilon_HfO));
    % VF(1,n)=VF(1,n-1)/2+VF(1,1)/2+(p_i(1,1)-p_i(1,n-1)+p_i(2,n)-p_i(n,n))*(del_tfe/(4*epsilon_HfO))-Qs(i)*(del_tfe^2/(2*epsilon_HfO));
    VF(m,1)=VF(m,2)/2+VF(m,n)/2+(p_i(m,2)-p_i(m,n)+p_i(1,1)-p_i(m-1,1))*(del_tfe/(4*epsilon_HfO))-Qs(i)*(del_tfe^2/(2*epsilon_HfO));
    VF(m,n)=VF(m,n-1)/2+VF(m,1)/2+(p_i(m,1)-p_i(m,n-1)+p_i(1,n)-p_i(m-1,n))*(del_tfe/(4*epsilon_HfO))-Qs(i)*(del_tfe^2/(2*epsilon_HfO)) ;
    VFE(:,:,i)=VF;
    %TDGL solution
    E(:,:,i)=VFE(:,:,i)/del_tfe;
    E_i=E(:,:,i);
    error=max(abs(VF_mod(:)-VF(:)))/max(abs(VF_mod(:)));
    VF_mod=VF+0.005*max(abs(VF_mod(:)-VF(:)));
end
VFi=VF;
for i=1:length(V_G)
    %  VF=-0.1*V_G(i);
    VG = V_G(i);
    xp = 10e-5;
    %    VFE(i,:)=0*V_G(i);
    error=2*eps; ncount=0;
    error1=error;
    for m1=1:iter
        xp0 = xp;
        F_psi1 = exp(-xp*q/kB/T) + xp*q/kB/T -1;
        F_psi2 = ni^2/N_A^2*(exp(xp*q/kB/T) - xp*q/kB/T -1);
        F_psi = sqrt(F_psi1 + F_psi2);
        if xp <= 0
            Qs(i) = sqrt(2*epsilon_si*kB*T*N_A)*F_psi;
            %Calculate charge when the surface potential is negative
        else
            Qs(i) = -sqrt(2*epsilon_si*kB*T*N_A)*F_psi;
            %Calculate charge when the surface potential is positive
        end
        %%Ferro effect

        %TDGL solution
        E(:,:,i)=VFE(:,:,i)/del_tfe;
        E_i=E(:,:,i);
        E_i1=circshift(E_i,2);
        E_i2=circshift(E_i,2,2);
        E_mod1=E_i1;
        E_mod2=E_i2;
        VF=VFi;
        VF_mod=VF;
        while (error > eps)
            ncount=ncount+1;
            % VF_mod=VF;
            VF_mod(2:end-1,2:end-1)=VF(1:end-2,2:end-1)/2+VF(2:end-1,1:end-2)/2+(p_i(3:end,2:end-1)-p_i(1:end-2,2:end-1)+p_i(2:end-1,3:end)-p_i(2:end-1,1:end-2))*(del_tfe/(4*epsilon_HfO))-Qs(i)*(del_tfe^2/(2*epsilon_HfO));
            VF2=circshift(VF,2);
            VF_mod2=VF2;
            VF_mod2(2:end-1,2:end-1)=VF2(1:end-2,2:end-1)/2+VF2(2:end-1,1:end-2)/2+(p_i(3:end,2:end-1)-p_i(1:end-2,2:end-1)+p_i(2:end-1,3:end)-p_i(2:end-1,1:end-2))*(del_tfe/(4*epsilon_HfO))-Qs(i)*(del_tfe^2/(2*epsilon_HfO));
            VF1=circshift(VF,2,2);
            VF_mod1=VF1;
            VF_mod1(2:end-1,2:end-1)=VF1(1:end-2,2:end-1)/2+VF1(2:end-1,1:end-2)/2+(p_i(3:end,2:end-1)-p_i(1:end-2,2:end-1)+p_i(2:end-1,3:end)-p_i(2:end-1,1:end-2))*(del_tfe/(4*epsilon_HfO))-Qs(i)*(del_tfe^2/(2*epsilon_HfO));
            VF=VF_mod;
            VF(1,:)=VF_mod2(3,:);
            VF(m,:)=VF_mod2(2,:);
            VF(:,1)=VF_mod1(:,3);
            VF(:,n)=VF_mod1(:,2);
            VF(1,1)=VF(1,n)/2+VF(1,2)/2+(p_i(1,2)-p_i(1,n)+p_i(2,1)-p_i(n,1))*(del_tfe/(4*epsilon_HfO))-Qs(i)*(del_tfe^2/(2*epsilon_HfO));         %  -Qs(i)*(del_tfe^2/(2*epsilon_HfO*5e-8));
            VF(1,n)=VF(1,n-1)/2+VF(1,1)/2+(p_i(1,1)-p_i(1,n-1)+p_i(2,n)-p_i(n,n))*(del_tfe/(4*epsilon_HfO))-Qs(i)*(del_tfe^2/(2*epsilon_HfO));     %-Qs(i)*(del_tfe^2/(2*epsilon_HfO*5e-8));
            VF(m,1)=VF(m,2)/2+VF(m,n)/2+(p_i(m,2)-p_i(m,n)+p_i(1,1)-p_i(m-1,1))*(del_tfe/(4*epsilon_HfO))-Qs(i)*(del_tfe^2/(2*epsilon_HfO));       %-Qs(i)*(del_tfe^2/(2*epsilon_HfO*5e-8));
            VF(m,n)=VF(m,n-1)/2+VF(m,1)/2+(p_i(m,1)-p_i(m,n-1)+p_i(1,n)-p_i(m-1,n))*(del_tfe/(4*epsilon_HfO))-Qs(i)*(del_tfe^2/(2*epsilon_HfO));   %-Qs(i)*(del_tfe^2/(2*epsilon_HfO*5e-8)) ;
            VFE(:,:,i)=VF;
            E(:,:,i)=VFE(:,:,i)/del_tfe;
            E_i=E(:,:,i);
            error=max(abs(VF_mod(:)-VF(:)))/max(abs(VF_mod(:)));
            VF_mod=VF+0.005*max(abs(VF_mod(:)-VF(:)));
        end
        % while (error1 > eps)

        P_mod=p_i;
        P_old=P_mod;
        p_i1=circshift(p_i,2);
        p_i2=circshift(p_i,2,2);
        P_mod1=p_i1;
        P_mod2=p_i2;

        P_mod(2:end-1,2:end-1)=(-(alpha*p_i(2:end-1,2:end-1)+beta*p_i(2:end-1,2:end-1).^3 +gm*p_i(2:end-1,2:end-1).^5+k1*(2*p_i(2:end-1,2:end-1)-p_i(1:end-2,2:end-1)-p_i(3:end,2:end-1))+k2*(2*p_i(2:end-1,2:end-1)-p_i(2:end-1,1:end-2)-p_i(2:end-1,3:end)))-E_i(2:end-1,2:end-1))*(del_t/gamma)+p_i(2:end-1,2:end-1);
        P_mod1(2:end-1,2:end-1)=(-(alpha*p_i1(2:end-1,2:end-1)+beta*p_i1(2:end-1,2:end-1).^3+gm*p_i1(2:end-1,2:end-1).^5+k1*(2*p_i1(2:end-1,2:end-1)-p_i1(1:end-2,2:end-1)-p_i1(3:end,2:end-1))+k2*(2*p_i1(2:end-1,2:end-1)-p_i1(2:end-1,1:end-2)-p_i1(2:end-1,3:end)))-E_i1(2:end-1,2:end-1))*(del_t/gamma)+p_i1(2:end-1,2:end-1);
        P_mod2(2:end-1,2:end-1)=(-(alpha*p_i2(2:end-1,2:end-1)+beta*p_i2(2:end-1,2:end-1).^3+gm*p_i2(2:end-1,2:end-1).^5+k1*(2*p_i2(2:end-1,2:end-1)-p_i2(1:end-2,2:end-1)-p_i2(3:end,2:end-1))+k2*(2*p_i2(2:end-1,2:end-1)-p_i2(2:end-1,1:end-2)-p_i2(2:end-1,3:end)))-E_i2(2:end-1,2:end-1))*(del_t/gamma)+p_i2(2:end-1,2:end-1);

        p_i=P_mod;
        p_i(1,:)=P_mod1(3,:);
        p_i(m,:)=P_mod1(2,:);
        p_i(:,1)=P_mod2(:,3);
        p_i(:,n)=P_mod2(:,2);
        p_i(1,1)=(-(alpha*p_i(1,1)+beta*p_i(1,1)^3+gm*p_i(1,1)^5+k1*(2*p_i(1,1)-p_i(m,1)-p_i(2,1))+k2*(2*p_i(1,1)-p_i(1,n)-p_i(1,2)))-E_i(1,1))*(del_t/gamma)+p_i(1,1);
        p_i(1,n)=(-(alpha*p_i(1,n)+beta*p_i(1,n)^3+gm*p_i(1,n)^5+k1*(2*p_i(1,1)-p_i(m,n)-p_i(2,n))+k2*(2*p_i(1,n)-p_i(1,n-1)-p_i(1,1)))-E_i(1,n))*(del_t/gamma)+p_i(1,n);
        p_i(m,1)=(-(alpha*p_i(m,1)+beta*p_i(m,1)^3+gm*p_i(m,1)^5+k1*(2*p_i(m,1)-p_i(m-1,1)-p_i(1,1))+k2*(2*p_i(m,1)-p_i(m,n)-p_i(m,2)))-E_i(m,1))*(del_t/gamma)+p_i(m,1);
        p_i(m,n)=(-(alpha*p_i(m,n)+beta*p_i(m,n)^3+gm*p_i(m,n)^5+k1*(2*p_i(m,n)-p_i(m-1,n)-p_i(1,n))+k2*(2*p_i(m,n)-p_i(m,n-1)-p_i(m,1)))-E_i(m,n))*(del_t/gamma)+p_i(m,n);

        P_init(:,:,i)=p_i/100;
        PT1(i,:)=sum(P_init(:,:,i));
        % error=max(abs(P_mod(:)-P_old(:)));
        % end

        VFe=mean(VFE(n,:,i));
        % VFe=0;
        % PT(i,:)=sum(P_init(i,:));
        %% Consistency check for surface potential
        if VG<=0
            f = VG+Qs(i)/Cox-VFe;
            del_xp = abs(xp0 - f); % Newton Method
            xp=xp0-del_xp*0.005;
            if abs(del_xp/xp) < tol
                psi(i) = xp;
                break
            end
        end
        if VG>0 && VG<2*phi_b
            f = VG+Qs(i)/Cox-VFe;
            del_xp = abs(xp0 - f); % Newton Method
            xp=(xp0+del_xp*0.005);
            if abs(del_xp/xp) < tol
                psi(i) = xp;
                break
            end
        end
        if  VG>2*phi_b
            f = VG+Qs(i)/Cox-(VFe);
            del_xp = abs(xp0 - f); % Newton Method
            xp=(xp0+del_xp*0.005);
            if abs(del_xp/xp) < tol
                psi(i) = xp;
                break
            end
        end
        %%
    end

    % VFE(i)=(2*(alpha)*abs(Qfe(i))+4*(beta)*abs(Qfe(i))^3+6*(gamma)*abs(Qfe(i))^5)*tfe;
end
% Numerically calculate the derivative of Qs with respect to psi
delta_psi = 1e-5; % Small perturbation for finite difference
dQs_dpsi_numerical = zeros(size(psi));
Csi=zeros(size(psi));
delta_Csi = 1e-5;
delta_Qs=zeros(size(psi));
dCsi_dQs=zeros(size(psi));
Qs=Qs';
for j = 1:length(psi)
    % Perturb the surface potential
    psi_perturbed = psi;
    psi_perturbed(j) = psi_perturbed(j) + delta_psi;

    % Calculate Qs for perturbed psi
    F_psi1_perturbed = exp(-psi_perturbed(j)*q/kB/T) + psi_perturbed(j)*q/kB/T -1;
    F_psi2_perturbed = ni^2/N_A^2*(exp(psi_perturbed(j)*q/kB/T) - psi_perturbed(j)*q/kB/T -1);
    F_psi_perturbed = sqrt(F_psi1_perturbed + F_psi2_perturbed);

    if psi_perturbed(j) <= 0
        Qs_perturbed = sqrt(2*epsilon_si*kB*T*N_A)*F_psi_perturbed;
    else
        Qs_perturbed = -sqrt(2*epsilon_si*kB*T*N_A)*F_psi_perturbed;
    end

    % Numerical derivative using finite difference
    dQs_dpsi_numerical(j) = (Qs_perturbed-Qs(j))/delta_psi;
    delta_Qs(j)=Qs_perturbed-Qs(j);
    Csi(j)=abs(dQs_dpsi_numerical(j)*Cox/(abs(dQs_dpsi_numerical(j))+Cox));

end
for v=1:length(psi)
    Csi_perturbed = Csi;
    Csi_perturbed(v) = Csi_perturbed(v) + delta_Csi;
    dCsi_dQs(v) = (Csi_perturbed(v)-Csi(v))/delta_Qs(v);

end
Qscm=Qs.*1;
semilogy(psi,abs(Qscm)/q,'b','linewidth',3);
hold on
set(gca, 'xlim', [-0.4 1.4], 'ylim', [1e11 1e14]);
set(gca,'fontsize',13);
xlabel('\psi_S (V)');
ylabel('|Qs|/q (cm^-2)');
%Plot x=0
plot(zeros(1,21),logspace(10,16,21),'k--')
%Plot surface potential (psi_S) vs. gate voltage with reference potential 2*psi_B
figure(2)
h1=plot(V_G,psi,'b','linewidth',3);
hold on

xlabel('V_G (V)');
ylabel('\psi_S (V)');
figure(3)
semilogy(V_G,abs(Qscm)/q,'b','linewidth',3);
hold on
set(gca, 'xlim', [-0.4 1.4], 'ylim', [1e11 1e14]);
set(gca,'fontsize',13);
xlabel('VG (V)');
ylabel('|Qs|/q (cm^-2)');
% figure(4)
% plot(V_G,Csi,"b",'LineWidth',3)
% ylabel('Cs');
% xlabel('\psi_S (V)');
% hold on