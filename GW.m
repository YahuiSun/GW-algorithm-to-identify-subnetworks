clc; clear all;


FinalPrize=zeros(1); % the net-weight of the identified subnetwork
TotalCost=zeros(1); % the net-cost of the identified subnetwork

load D_01_b; % load the drug similarity network
dataname='D_01_b'; % label the saved data

Ti=100; % time increment unit
N=node_num; % the number of vertices in the initial graph
Real_NW=node_weight; % the node weights

a=0;
NW=node_weight;
for i=1:N
    if Terminal(i)==1
        if a==0
            a=1;
            R=i; % choose one terminal to be the root
            NW(i)=1e6;
        else
            NW(i)=1e6; % give other terminals big prizes
        end
    end
end
%

% initialization
cpn_num=N; % numbers of component
surp=zeros(cpn_num,1);% surplus of each component
for i=1:N
    surp(i)=NW(i);
end

active=zeros(cpn_num,1); % 1 means the component is active
for i=1:N
    if surp(i)>0 && i~=R
        active(i)=1;
    end
end

index=zeros(N,1);
for i=1:N
    index(i)=i; % the number is the component number which vertex i belongs to
end

defi=L; % deficit of each edge
paid_set=zeros(N); % included edges in the growth phase

% GW growth
a=0;
while 0<1
    a=a+1;
    % update surplus
    for i=1:cpn_num
        if active(i)==1
            surp(i)=surp(i)-Ti;
        end
    end
    %
    % update deficit
    for i=1:N
        for j=i:N
            if set(i,j)==1 && index(i)~=index(j) % active edge: not joining two vertices in the same component
                choice=active(index(i))+active(index(j));
                if choice==1
                    defi(i,j)=defi(i,j)-Ti;
                elseif choice==2
                    defi(i,j)=defi(i,j)-2*Ti;
                end
            end
        end
    end
    %
    % merge cpn
    for i=1:N
        for j=i:N
            if defi(i,j)<=0 && index(i)~=index(j) % an active edge's deficit is reduced to 0
                [index,surp,cpn_num,active]=Function_merge_cpn(index(i),index(j),N,index,surp,cpn_num,active);% merge component(index(i)) and component(index(j))
                paid_set(i,j)=1; paid_set(j,i)=1;
            end
        end
    end
    % deactive cpn
    for i=1:cpn_num
        if surp(i)<=0
            active(i)=0; % deactived cpn may not be connected to the root cpn, so GW growth may find several trees
        end
    end
    %
    % break condition 1
    finish=1;
    for i=1:N
        if index(i)~=index(R)
            finish=0;
        end
    end
    if finish==1
        break;
    end
    %
    % break condition 2
    finish=1;
    for i=1:N
        if index(i)~=index(R) && active(index(i))==1
            finish=0;
        end
    end
    if finish==1
        break;
    end
    %
end

% remove trees that is disconnected to the root tree
LL=sparse(paid_set);
[S, C] = graphconncomp(LL);
for im=1:N
    for in=1:N
        if paid_set(im,in)==1 && abs(C(im)-C(R))+abs(C(in)-C(R))>0
            paid_set(im,in)=0; paid_set(in,im)=0;
        end
    end
end
%

% Strong pruning
[prune_set]=Function_StrongPrune2(paid_set,NW,N,L,R);
%

% solution quality
[FinalPrize,TotalCost]=Function_solution(FinalPrize,TotalCost,1,1,prune_set,L,N,R,node_weight,Terminal,0);

% analysis
GW_set=prune_set;
Analysis_set=GW_set;
fprintf(['\n'])
[degree]=Function_OutputDegree(Analysis_set,N);
fprintf(['In GW_set\n'])
IncludedVertex1=zeros(N,1);
for i=1:N
    if degree(i)>0
        IncludedVertex1(i)=1;
    end
end
fprintf(['There are ',num2str(sum(IncludedVertex1)),' out of ',num2str(N),' vertices included.\n']);
removeV=N-sum(IncludedVertex1);
fprintf(['There are ',num2str(removeV),' out of ',num2str(N),' vertices removed.\n']);

save(['Result_',dataname,'_NW_',num2str(node_weight(1)),...
    '_GW_FitnessValue_',num2str(FinalPrize(1)),'_RemovedVertices_',num2str(removeV)]);



