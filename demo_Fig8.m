%% Script Fig. 8
% The code for computing the bifurcation points in tau-Q plane is provided (Figure 8).
clear; close all; clc;
T      =   2.5e-6;         % tickness ( m )
L      =   160e-6;         % length ( m )
W      =   33e-6;          % width ( m )
rho    =   2.33e3;         % density of the cantilever( kg/m^3 )
f0     =   139.4e3;        % natural frequency in air( Hz )
m0     =   rho*L*W*T*1e9;  % cantilever mass ( mu*g )
w0     =   2*pi*f0*1e-6;   % natural angular frequency in air ( rad/mus )
k      =   w0^2*m0;        % spring constant ( mug/mus^2 ) 
%
% Values of the parameters in air
%
Q0     =   240;            % quality factor
c0     =   w0*m0/Q0;       % damping ( mug/mus ) 
%
% Feedback parameters
%
taul   =   10.7;           % intrinsic delay tau_loop ( mus )
sig    =   1e-1;           % saturation threshold
G      =   0.5;            % Gain (initial value) (mug/mus^2)
%
% Characteristics of the fluid
%
eta_w  =   0.8e-3;         % viscosity of the fluid (change it for different fluid)
rhof   =   998;            % density of the fluid ( kg/m^3 )
%
% We define an average frequency of oscilation to be used in the
% computation of ma and ca which represent the effect of the viscous fluid
% on the oscillating cantilever.
%
fosc   =   0.050;          % average frequency ( muHz)
w      =   2*pi*fosc*1e6;
a1=1.0553; a2=3.7997; b1=3.8018; b2=2.7364;
%
% Computation of ma and ca 
%
ma     =   (pi*rhof*W^2/4)*L*(a1+a2/W*sqrt((2*eta_w)/(rhof*w)))*1e9 ;
ca     =   pi/4*rhof*W^2*L*w*(b1/W*sqrt((2*eta_w)/(rhof*w))+b2/W^2*2*eta_w/(rhof*w))*1e3;
%
% Computation of the overall mass and damping coefficients
%
m      =   m0 + ma; 
c      =   c0 + ca;
wn     =   sqrt( k / m );  % natural frequency in the fluid ( rad/mus )
tau    =   taul*wn;        % initial delay
Q      =   m*wn/c;         % quality factor to be used in (1)

save('data');
% load ddebiftool into path
addpath('dde_biftool/ddebiftool',... 
'dde_biftool/ddebiftool_utilities',...
'dde_biftool/ddebiftool_extra_psol');

funcs = set_funcs('sys_rhs',@(xx,p)sys_general_rhs(xx,p),...
'sys_tau',sys_general_tau,'sys_deri',@sys_general_deri);
           
par     =   [G Q tau]; % bifurcation parameters
ind_G   =   1; % index of G
ind_tau =   3; % index of tau


% 1) Branch F1-F2
load('pitch02.mat')
ind=18; % point near the fold F1
[foldfuncs1,F1_2,suc]=SetupPOfold(funcs,nonsym,ind,'contpar',[3,2],...
    'dir',3,'print_residual_info',1,'step',0.01,'plot_measure',[],...
    'min_bound',[3,0; 2, 0],'max_bound',[3,8; 2,nonsym.point(1).parameter(2)],...
    'max_step',[3,0.1; 2,0.1],'pitchfork',false);
F1_2.method.point.print_residual_info=0;
F1_2.method.continuation.plot=1;% change in 0 to omit the real time computation of branch F1_2
figure(1);clf;
F1_2=br_contn(foldfuncs1,F1_2,500);
xlabel('$\tau$','interpreter','latex', 'Fontsize',15)
ylabel('Q')
storeF1_2=[];
for k=1:length(F1_2.point)
    storeF1_2=[storeF1_2, [F1_2.point(k).parameter(3);F1_2.point(k).parameter(2)]];
end
% 2) branch of the fold points in F3
load('branch_LF061')
figure(2);clf;
ind=279;% point near the fold F3
[foldfuncs,F3,suc]=SetupPOfold(funcs,branch_LF061,ind,'contpar',[3,2],...
    'dir',3,'print_residual_info',1,'step',0.01,'plot_measure',[],...
    'min_bound',[3,0; 2, 0],'max_bound',[3,8; 2,20],...
    'max_step',[3,0.1; 2,0.2],'pitchfork',false);
F3.method.point.print_residual_info=0;
F3.method.continuation.plot=0;
F3=br_contn(foldfuncs,F3,100);
storeF3=[];
for k=1:length(F3.point)
    storeF3=[storeF3, [F3.point(k).parameter(3);F3.point(k).parameter(2)]];
end
% 3) branch of H points
load('brQ28H')
ind_tr=39;% near the secondary Hopf bifurcation
[trfuncs,trbranch,suc]=SetupTorusBifurcation(funcs,brQ28H,ind_tr,...
    'contpar',[3,2],'dir',2,'print_residual_info',1,'step',-0.01,'plot_measure',[],...
    'min_bound',[3,0; 2, 3],'max_bound',[3,8; 2,brQ28H.point(ind_tr).parameter(2)],...
    'max_step',[3,0.1; 2,0.1]);
trbranch.method.point.print_residual_info=0;
trbranch.method.continuation.plot=1;
figure(3);clf;
trbranch=br_contn(trfuncs,trbranch,900);
storeH=[];
for k=1:length(trbranch.point)
    storeH=[storeH, [trbranch.point(k).parameter(3);trbranch.point(k).parameter(2)]];
end

figure(4);clf;
plot(storeH(1,:),storeH(2,:),'k--')
hold on
plot(storeF1_2(1,:),storeF1_2(2,:),'k.-')
hold on
plot(storeF3(1,:),storeF3(2,:),'k-')
xlabel('\tau');ylabel('Q');
ylim([3,11.53])
title('Bifurcations in tau-Q plane');
legend('H','F1-F2','F3')