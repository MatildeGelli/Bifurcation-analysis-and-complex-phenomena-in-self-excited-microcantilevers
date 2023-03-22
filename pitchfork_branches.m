function nonsym=pitchfork_branches(branch_pt,it)
%

parnames={'G','Q','tau'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
fun = set_symfuncs(@sym_gelli_rhs,'sys_tau',@()ind.tau);
branch_pt.point=arrayfun(@(p)dde_psol_create('point',p,'stability',p.stability),branch_pt.point);
nunst=GetStability(branch_pt,'exclude_trivial',true);
pf=find(diff(nunst),1);%last stable
mth=df_mthod('psol');
mth.point.matrix='sparse';
mth.stability.matrix='sparse';
branch_pt.method=mth;
[nonsym,suc]=SetupBranchSwitch(fun,branch_pt,pf+(0:1),...
    'print_residual_info',1,'max_step',[ind.tau,0.005]) %[ind.tau,0.01;0,0.1]%#ok<*NOPTS>
nonsym.method.continuation.plot=0;  
nonsym=br_contn(fun,nonsym,it);
nonsym=br_stabl(fun,nonsym,0,0);
end