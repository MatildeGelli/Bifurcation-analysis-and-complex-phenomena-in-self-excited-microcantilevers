%% User-provided partial derivatives of right-hand side |f|
% If a user-provided function is not provided, DDE-Biftool will use
% finite-difference approximation (implemented in |df_deriv.m|).
%%
function J=sys_general_deri(xx,p,nx,np,v)
load('data.mat')
%% Parameter vector
% p=[G Q tau];


%f(1,1)=xx(2,1);
%f(2,1)=-sig/(k1)*tanh(p(1)/sig*xx(1,2))-1/Q*xx(2,1)-xx(1,1);


J=[];
if length(nx)==1 && isempty(np) && isempty(v)
    %% First-order derivatives wrt to state nx+1
    if nx==0 % derivative wrt x(t)
        J(1,1)=0;
        J(1,2)=1;
        J(2,1)=-1;
        J(2,2)=-1/p(2);
    elseif nx==1 % derivative wrt x(t-tau)
        J(1,1)=0;
        J(1,2)=0;
        J(2,1)=-p(1)/k*(1-tanh(p(1)/sig*xx(1,2))^2);
        J(2,2)=0;
    end
elseif isempty(nx) && length(np)==1 && isempty(v)
    %% First order derivatives wrt parameters
    if np==1 % derivative wrt G
        J(1,1)=0;
        J(2,1)=-1/k*(1-tanh(p(1)*xx(1,2)/sig)^2)*xx(1,2);
    elseif np==2 % derivative wrt Q
        J(1,1)=0;
        J(2,1)=1/p(2)^2*xx(2,1);
    elseif np==3 % derivative wrt tau
        J=zeros(2,1);
    end
elseif length(nx)==1 && length(np)==1 && isempty(v)
%Mixed state, parameter derivatives
    if nx==0 % derivative wrt x(t)
        if np==2 % derivative wrt Q
            J(1,1)=0;
            J(1,2)=0;
            J(2,1)=0;
            J(2,2)=1/p(2)^2;
        else %zero for all the others
            J=zeros(2);
        end
    elseif nx==1 % derivative wrt x(t-tau)
        if np==1 % derivative wrt G
            J(1,1)=0;
            J(1,2)=0;
            J(2,1)=-1/k*(1-tanh(p(1)*xx(1,2)/sig)^2)-p(1)/k*(-2*tanh(p(1)*xx(1,2)/sig)*(1-tanh(p(1)*xx(1,2)/sig)^2)*xx(1,2)/sig);
            J(2,2)=0;
        else 
            J=zeros(2);
        end
    end

elseif length(nx)==2 && isempty(np) && ~isempty(v)
    %second order derivatives wrt state variables
    if nx(1)==0 %derivative wrt x(t)
        J=zeros(2);
    elseif nx(1)==1 %derivative wrt x(t-tau)
        if nx(2)==1
            J(1,1)=0;
            J(1,2)=0;
            J(2,1)=(2*p(1)^2)/(k*sig)*tanh(p(1)/sig*xx(1,2))*(1-tanh(p(1)/sig*xx(1,2))^2)*v(1);
            J(2,2)=0;
        else
            J=zeros(2);
        end
    end
end
%% Otherwise raise error
% Raise error if the requested derivative does not exist
if isempty(J)
    error(['SYS_GENERAL_DERI: requested derivative could not be computed!'],nx,np,size(v));
end
end
