%% Self-excited microcantilever oscillating in a viscous fluid
%
% From:
%
% "Bifurcation analysis and complex phenomena in self-excited
% microcantilevers", M. Gelli, J. Mouro, P. Paoletti, B.
% Tiribilli and M. Basso
%
% This demo illustrates the bifurcations characterizing the dynamics of a
% self-excited microcantilever subjected to the presence of a variable
% delay in the loop.
%
% Differential equation of the system:
% 
% $$x''(t)+x'(t)/Q+x(t) = -sigma/k*tanh(G/sigma*x(t-tau)) $$ (1)
% 
% where:  
% 
% * x(t) is the deflection of the cantilever and x'(t) and x''(t) are the
% first and the second derivatives respectively
% * Q is the quality factor 
% * k is the spring constant
% * tanh() is used instead of the saturation function sat( ,sig) with sig
% being the saturation threshold. This nonlinear function is necessary to
% limit the oscillation amplitude to few nanometers.
% * G is the feedback gain
% * tau is the variable delay and it includes the constant term tau_loop 
%
%
% Acknowledgements
% We would like to thank Professor Jan Sieber for his precious help in
% computing the symmetry-breaking bifurcation without whom it would not be
% possible.
%%
clear;
%format compact;
close all;clc;

% Dimension of the cantilever and its material
%
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
% We define an average frequency of oscilation (fosc) to be used in the
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

funcs = set_funcs('sys_rhs',@(xx,p)sys_general_rhs(xx,p),...
'sys_tau',sys_general_tau,'sys_deri',@sys_general_deri);
           
par     =   [G Q tau]; % parameters
ind_G   =   1; % index of G
ind_tau =   3; % index of tau
%%
stst.kind='stst';
stst.parameter=par;
stst.x=[0;0];
method=df_mthod(funcs,'stst');
method.stability.minimal_real_part=-1;
stst=p_correc(funcs,stst,[],[],method.point);
stst.stability=p_stabil(funcs,stst,method.stability);
figure(1);clf;
p_splot(stst)
% Preparing for a first continuation of equilibrium points wrt G
branch = df_brnch(funcs,ind_G,'stst');
branch.parameter
branch.parameter.min_bound(2,:) = [ind_G 0];
branch.parameter.max_bound(1,:) = [ind_G 1];
branch.parameter.max_step(1,:) = [ind_G 0.01];
branch.point = stst; % first point of the branch
stst2 = stst;
stst2.parameter(ind_G) = stst2.parameter(ind_G)+0.01;
[stst2,~] = p_correc(funcs,stst2,[],[],method.point);
branch.point(2) = stst2;% second point of the branch
figure(2)
[branch, s, f, r] = br_contn(funcs,branch,300);
branch = br_rvers(branch);
[branch, s, f, r] = br_contn(funcs,branch,300);
xlabel('G', 'Fontsize', 15)
ylabel('First component x(1)','Fontsize', 15)
branch.method.stability.minimal_real_part=-1;
branch = br_stabl(funcs,branch,0,1);
[xm,ym] = df_measr(0,branch);
%%
figure(3);clf;
tiledlayout(2,2)
nexttile([1,2])
br_splot(branch,xm,ym)
xlabel('G', 'Fontsize', 15)
ylabel('x(1)','Fontsize', 15)
nexttile([1,2])
br_splot(branch,xm,ym)
xlabel('G', 'Fontsize', 15)
ylabel('x(1)','Fontsize', 15)
xlim([0 0.1])
sgtitle('Hopf bifurcation points','Fontsize',17)
figure(4);clf;
tiledlayout(2,1)
nexttile
[xm1, ym1]=df_measr(1, ind_G, 'stst');
br_plot(branch,xm1,ym1,'b')
hold on
plot([0 1], [0 0],'k-.')
xlabel('G','Fontsize', 15)
ylabel('$\Re(\lambda)$','interpreter','latex','Fontsize',15)
nexttile
br_plot(branch,[],ym1,'b');
br_plot(branch,[],ym1,'b.');
plot([0 110],[0 0],'k-.');
xlabel('Point number along branch','Fontsize', 15);
ylabel('$\Re(\lambda)$','interpreter','latex','Fontsize',15);
method=df_mthod(funcs, 'hopf',1);
hopf_pL=p_tohopf(funcs, branch.point(100));%first Hopf point (low freq)
[hopf_pL, ~]=p_correc(funcs,hopf_pL, 1, [], method.point);
hopf_pH=p_tohopf(funcs,branch.point(95));% second Hopf point (high freq)
[hopf_pH, ~]=p_correc(funcs,hopf_pH, 1, [], method.point);
% We make a first continuation wrt both G and tau. This tells us
% how the two Hopf points evolves when varying both the gain and
% the delay.
% With branch_hL we define the branch of Hopf points at low
% frequency, instead branch_hH is the one at high frequency.
branch_hL=df_brnch(funcs,[ind_G 3],'hopf');%
branch_hL.parameter
branch_hL.parameter.min_bound(1:2,:)=[[ind_G 0]', [ind_tau tau]']';
branch_hL.parameter.max_bound(1:2,:)=[[ind_G 1]',[ind_tau tau+4]']';
branch_hL.parameter.max_step(1:2,:)=[[ind_G 0.1]',[ind_tau 0.1]']';
branch_hL.point=hopf_pL;% first point of branch_hL
sec_hopfL=hopf_pL;
sec_hopfL.parameter(3)=sec_hopfL.parameter(3)+0.01;
sec_hopfL=p_correc(funcs,sec_hopfL,ind_G,[],method.point);
branch_hL.point(2)=sec_hopfL;% second point of branch_hL
figure(5);clf;
branch_hL=br_contn(funcs,branch_hL,500);
branch_hH=df_brnch(funcs,[ind_G 3],'hopf');
branch_hH.parameter
branch_hH.parameter.min_bound(1:2,:)=[[ind_G 0]', [3 tau]']';
branch_hH.parameter.max_bound(1:2,:)=[[ind_G 1]',[3 tau+4]']';
branch_hH.parameter.max_step(1:2,:)=[[ind_G 0.1]',[3 0.1]']';
branch_hH.point=hopf_pH;
sec_hopfp2=hopf_pH;
sec_hopfp2.parameter(3)=sec_hopfp2.parameter(3)+0.1;
sec_hopfp2=p_correc(funcs,sec_hopfp2,ind_G,[],method.point);
branch_hH.point(2)=sec_hopfp2;
hold on
[branch_hH s f r]=br_contn(funcs,branch_hH,500);
xlabel('G (\mug/\mus^{2})','Fontsize', 15);
ylabel('$\tau$','interpreter','latex','Fontsize', 15);
%

%% High frequency branch
% Here the HF branch of periodic solutions, branch_HF, is
% computed. This branch arises from the second Hopf point hopf_pH. The
% computation is divided into two steps:
%
% 1) Continuation wrt G until G = 1 mug/mus^2 which is the value used in the
% experiments
% 2) Continuation wrt tau in order to understand how the delay causes the
% jumps from the low frequency branch to the high frequency branch 
%
% The computation may take several minutes depending on the values chosen
% for the variables int, degree, max_bound and max_stepsize.

disp('Computing the HF branch')

int=50;
degree=3;
method=df_mthod(funcs,'psol');
method.stability.max_number_of_eigenvalues=5;
[p, step]=p_topsol(funcs,hopf_pH,1e-2,degree,int);
[p, success]=p_correc(funcs,p,ind_G,step,method.point);
deg_p=p_topsol(funcs,hopf_pH,0,degree,int);
method.point.adapt_mesh_after_correct=0;
method.point.newton_max_iterations=7;
method.point.newton_nmon_iterations=2;
branch_1=df_brnch(funcs,ind_G,'psol');
% *1ST STEP*
disp('Continuation wrt G')
G_minH=hopf_pH.parameter(ind_G);
branch_1.parameter.min_bound(1,:)=[ind_G G_minH];
branch_1.parameter.max_bound(1,:)=[ind_G 1.01];
branch_1.parameter.max_step(1,:)=[ind_G 0.1];
branch_1.point=deg_p;
branch_1.point(2)=p;
figure(6);clf;
branch_1=br_contn(funcs,branch_1,500);
xlabel('G','Fontsize',15);
ylabel('Amplitude','Fontsize',15);
% *2ND STEP*
disp('Continuation wrt tau')
branch_HF=df_brnch(funcs,ind_tau,'psol');
branch_HF.parameter.min_bound(1,:)=[ind_tau tau];
branch_HF.parameter.max_bound(1,:)=[ind_tau tau+3];
branch_HF.parameter.max_step(1,:)=[ind_tau 0.05];
branch_HF.point=branch_1.point(end);
punto2=branch_1.point(end);
punto2.parameter(3)=punto2.parameter(3)+0.001;
branch_HF.point(2)=punto2;
figure(7);clf;
[branch_HF, s, f, r]=br_contn(funcs,branch_HF,5000);
xlabel('$\tau$','interpreter','latex','Fontsize',15);
ylabel('Amplitude','Fontsize',15);
disp('Evaluating the stability of HF branch')
% This task may take few minutes
branch_HF=br_stabl(funcs,branch_HF,0,1);
figure(8);clf;
%nunst_2=GetStability(branch_HF,'exclude_trivial',true);
%br_splot2(branch_HF,nunst_2,1.5);
BifDiagram(branch_HF,ind_tau,0);
%% Low frequency branch
% In this section the LF branch, named as branch_pt, is
% computed. This branch arises from the first Hopf point hopf_pL.

disp('Computing the LF branch')
% *1ST STEP*
disp('Continuation wrt G')
int=50;
degree=3;
[psol, step]=p_topsol(funcs,hopf_pL,1e-2,degree,int);
method=df_mthod(funcs,'psol');
method.stability.max_number_of_eigenvalues=5;
psol=p_correc(funcs,psol,ind_G,step,method.point);
method.point.adapt_mesh_after_correct=1;
branch_psol=df_brnch(funcs,ind_G,'psol');
G_min=hopf_pL.parameter(ind_G);
branch_psol.parameter.min_bound(1,:)=[ind_G G_min];
branch_psol.parameter.max_bound(1,:)=[ind_G 1];
branch_psol.parameter.max_step(1,:)=[ind_G 0.1];
deg_psol=p_topsol(funcs,hopf_pL,0,degree,int);%this solution has a zero amplitude
branch_psol.point=deg_psol;
branch_psol.point(2)=psol;
figure(9);clf;
branch_psol=br_contn(funcs,branch_psol,700);
xlabel('G','Fontsize',15);
ylabel('Amplitude (nm)','Fontsize',15);
%
% *2ND STEP*
% Notice that for those values of Q where there are two low frequency
% branches, (Q>=4.5), br_contn is not able to capture the whole solutions.
% In order to do this, one needs to use the SetupBranchSwitch function.
%
disp('Continuation wrt tau')
figure(10);clf;
branch_pt=df_brnch(funcs,ind_tau,'psol');
branch_pt.parameter.min_bound(1,:)=[ind_tau tau];
branch_pt.parameter.max_bound(1,:)=[ind_tau tau+4];
branch_pt.parameter.max_step(1,:)=[ind_tau 0.1];
branch_pt.point=branch_psol.point(end);
punto2=branch_pt.point;
punto2.parameter(3)=punto2.parameter(3)+0.001;
branch_pt.point(2)=punto2;
branch_pt=br_contn(funcs,branch_pt,2000);
xlabel('$\tau$','interpreter','latex','Fontsize',15);
ylabel('Amplitude','Fontsize',15);
% Wait some more minutes for evaluating the stability of the branch
branch_pt=br_stabl(funcs,branch_pt,0,1);
% nunst_1=GetStability(branch_pt,'exclude_trivial',true);
%save('branch_pt_trial','branch_pt')
%%
figure(11);clf;
BifDiagram(branch_pt,ind_tau,0);
hold on
legend off
axis tight 
xlabel('$\tau$','interpreter','latex','Fontsize',15);
ylabel('Amplitude','Fontsize',15);
ylim([7,16]);xlim([3.7,7])
if Q>6.7
    it=295; %try different value for computing more points
    nonsym=pitchfork_branches(branch_pt,it);
    BifDiagram(nonsym,3,0); 
elseif Q>4.6 && Q<6.7
    %closed orbit
    load('bif_eta02.mat')
    it=350;
    nonsym=pitchfork_branches(branch_LF02,it);
    %nonsym=pitchfork_branches(branch_pt,it);
    Qd=Q;% desired value of Q
    br_iso=isolated_lc(nonsym,170,Qd,tau);
    BifDiagram(br_iso,3,0);
end
 
