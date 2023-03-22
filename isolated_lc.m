function br_iso=isolated_lc(br,it, Q_d, tau)
% isolated_lc computes the closed structure br_iso of the LF branch. This 
% solution exists only for 4.65 < Q < 6.8. This branch cannot be computed 
% using the traditional way (i.e. from a Hopf point), hence, we propose the
% following process. Let us say *Q_d* is the desired quality factor of the
% branch we want to compute and *br* an asymmetric branch born from PF1 with 
% a different quality factor, then:
% 1) Make a continuation wrt Q until Q = Q_d. In this way we get the first
% point of branch br_iso.
% 2) Make a continuation wrt tau. This last continuation allows us to
% compute the solution for different values of the delay.
% 
% INPUTS:
% * br: asymmetric branch computed from the pitchfork in PF1
% * it: number of iterations
% * Q_d: desired value of the quality factor
% * tau: upper bound for the delay
% 
% OUTPUTS: 
% br_iso: closed structure of the LF branch
%%

parnames={'G','Q','tau'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
fun = set_symfuncs(@sym_gelli_rhs,'sys_tau',@()ind.tau);
%continuation wrt Q
br1=df_brnch(fun,2,'psol');
br1.parameter.min_bound(1,:)=[2 Q_d];
br1.parameter.max_bound(1,:)=[2 br.point(110).parameter(2)];
br1.parameter.max_step(1,:)=[2 0.1];
point1=br.point(110);
br1.point=point1;
point2=point1;
point2.parameter(2)=point2.parameter(2)-1e-3;
br1.point(2)=point2;
br1.method.continuation.plot=0;
br1=br_contn(fun,br1,1e3);
%continuation wrt tau
br_iso=df_brnch(fun,3,'psol');
br_iso.parameter.min_bound(1,:)=[3 0];
br_iso.parameter.max_bound(1,:)=[3 tau+4];
br_iso.parameter.max_step(1,:)=[3 0.01];
br_iso.point=br1.point(end);
point2=br1.point(end);
point2.parameter(3)=point2.parameter(3)+0.001;
br_iso.point(2)=point2;
br_iso.method.continuation.plot=0;
br_iso=br_contn(fun,br_iso,it);
br_iso=br_stabl(fun,br_iso);

end