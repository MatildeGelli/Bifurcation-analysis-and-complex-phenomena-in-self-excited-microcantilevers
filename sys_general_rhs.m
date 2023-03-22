function f=sys_general_rhs(xx,p)

load('data.mat')

% p=[G Q tau]; % bifurcation parameters
%
% Differential equation of the system:
% 
% $$x''(t)+x'(t)/Q+x(t) = -sigma/k*tanh(G/sigma*x(t-tau)) $$ (1)
%
% This form is not suited for the implementation on DDE-Biftool. 
% In fact, it is necessary to rewrite the dynamics in a state-space form. 
% Thus, define the state variable xx which has two rows and two columns.
% The first column contains x(t) and x'(t), instead the second column
% contains a delayed version of those components, i.e: 
%
% * x(t)  =  xx(1,1)
% * x'(t)  =  xx(2,1)
% * x(t-tau)  = xx(1,2)
% * x'(t-tau)  = xx(2,2)
%
% Therefore, equation (1) can be rewritten as:
%
f(1,1) = xx(2,1);
f(2,1) = -sig/k*tanh(p(1)/sig*xx(1,2))-1/p(2)*xx(2,1)-xx(1,1);
end