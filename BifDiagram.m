function BifDiagram(br,par,label)
% BifDiagram plots the bifurcation diagram of branch br wrt par. A
% different color of the points of the branch means a different
% number of stable/unstable multipliers.
%
% INPUT:
% * br is the branch for which we ask the bifurcation diagram
% * par is the bifurcation parameter
% * label decides what to plot, i.e label=0 amplitudes are displayed, else if
% label ~=0 frequency is evaluated.
%
MAX=[];
par_v=[];
branch=br;
for k=1:length(br.point)
    mod=sqrt(br.point(k).profile(1,:).^2.+br.point(k).profile(2,:).^2);
    [mod_M, k_M]=max(mod);
    MAX=[MAX,[mod_M;k_M]];
    par_v=[par_v, br.point(k).parameter(par)];
end
if ~isempty(br.point(1).stability)
    nunst=(GetStability(br,'exclude_trivial',true));
else
    error('compute stability first');
end
[nunst_unique, i_unique]=unique(nunst,'stable');
index_array=cell(length(i_unique),1); % empty array
par_array=cell(length(i_unique),1); % empty array
freq_array=cell(length(i_unique),1);
clrs={[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560],[0.3010 0.7450 0.9330],[.8 .2 .6]};
for k=1:length(i_unique)
    for l=1:length(nunst)
        if nunst(l)==nunst_unique(k)
        index_array{k}=[index_array{k}, l];
        par_array{k}=[par_array{k}, br.point(l).parameter(par)];
        branch.point(l).color=clrs{k};
        freq_array{k}=[freq_array{k},1/(br.point(l).period)];
        end
    end
end


%figure
for k=1:length(i_unique)
    % check if the vector is made of contiguous elements
    if diff(index_array{k}(:))==ones(1,length(index_array{k}(:))-1)
        if label==0
         plot(par_array{k}(:),MAX(1,[index_array{k}(:)]),'-','Color',clrs{nunst_unique(k)+1},'LineWidth',1.5);
        else
            plot(par_array{k}(:),freq_array{k}(:),'-','Color',clrs{nunst_unique(k)+1},'LineWidth',1.5);
        end
            hold on
    else
         ind_c=find(diff(index_array{k}(:))>1);
         v=cell(length(ind_c)+1,1);
         for l=1:length(ind_c)+1
             if l==1
                 i_start=1;
                 i_stop=ind_c(l);
             else
                 i_start=ind_c(l-1)+1;
                 if l<length(ind_c)+1
                     i_stop=ind_c(l);
                 elseif l==length(ind_c)+1
                     i_stop=length(index_array{k});
                 end
             end

          v{l}=[par_array{k,1}([i_start:i_stop])];
              if label==0
              plot(v{l}, MAX(1,[index_array{k,1}([i_start:i_stop])]),'-','Color',clrs{nunst_unique(k)+1},'LineWidth',1.5);
              else 
                  plot(v{l},freq_array{k,1}([i_start:i_stop]),'-','Color',clrs{nunst_unique(k)+1},'LineWidth',1.5);
              hold on
              end
        hold all
        drawnow;
        end
    end

% Add markers at bifurcation points
for k=1:length(branch.point)-1
    if branch.point(k).color~=branch.point(k+1).color
        x=(branch.point(k).parameter(par)+branch.point(k+1).parameter(par))/2;
        if label==0
        y=(MAX(1,k)+MAX(1,k+1))/2;
        else
            y=(1/(br.point(k).period)+1/(br.point(k+1).period))/2;
        end
        plot(x,y,'ko','MarkerFaceColor','k','MarkerSize',4)
    end
end
% for legend
L=zeros(length(nunst_unique),1);
legendText=cell(length(nunst_unique),1);
[n_sorted, i_sorted]=sort(nunst_unique);
for k=1:length(n_sorted)
    L(k)=plot(NaN,NaN,'-','Color',clrs{n_sorted(k)+1},'LineWidth',1.5);
    hold all
    if n_sorted(k)==0
        legendText{k} = sprintf(['stable']);
    else
        legendText{k} = sprintf(['unst=' num2str(n_sorted(k))]);
    end
end
legend(L,legendText,'Location','Best','Fontsize',13);
set(gca,'TickLabelInterpreter','none','Fontsize',13);
end
