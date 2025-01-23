clear;
clc;
close all;
e=exp(1);
for i=1:8
    X(i,1)=rand()*(5-1)+1;
end
es=0.0000000001;
delta_t=1;
L=2100;        %????????
N=round(L/delta_t);  %????????
ts=zeros(N,1);
sigma=0.5;
p=zeros(1,25);
p(1)=-0.5;
p(2)=-0.475;
p(3)=-0.45;
p(4)=-0.4;
p(5)=-0.375;
p(6)=-0.35;
p(7)=-0.325;
p(8)=-0.3;
p(9)=-0.275;
p(10)=-0.25;
p(11)=-0.225;
p(12)=-0.2;
p(13)=-0.175;
p(14)=-0.15;
p(15)=-0.125;
p(16)=-0.1;
p(17)=-0.08;
p(18)=-0.05;
p(19)=-0.01;
p(20)=-0.001;
p(21)=0.01;
p(22)=0.05;
p(23)=0.08;
p(24)=0.1;
p(25)=0.15;

D= [-2 1 0 0 0 0 0 0 ;...
    -1 -1 0 0 0 0 0 0 ;...
    1 0 -1 0 0 0 0 0 ;...
    1 0 0 -1 0 0 0 0 ;...
    1 0 0 0 -1 0 0 0 ;...
    0 0 0 0 1 -1 0 1 ;...
    0 0 0 0 0 0 -1 1 ;...
    0 0 0 0 0 0 -1 -1];

adjacent_network= [1 2 3 4 5 6;...
    2 1 3 4 5 6;...
    3 1 2 0 0 0;...
    4 1 2 0 0 0;...
    5 1 2 6 0 0;...
    6 1 2 5 7 8;...
    7 6 8 0 0 0;...
    8 6 7 0 0 0];

% 取十个正常样本
sample_num=10;
reference_data=zeros(8,sample_num);
q(1)=0.96^(1/abs(p(1)));
 E=[-2/5*q(1) 0 0 0 0 0 0 0 ;...
            0 -4/5 0 0 0 0 0 0 ;...
            0 0 -5/5 0 0 0 0 0 ;...
            0 0 0 -6/5 0 0 0 0 ;...
            0 0 0 0 -7/5 0 0 0 ;...
            0 0 0 0 0 -8/5 0 0 ;...
            0 0 0 0 0 0 -9/5 0 ;...
            0 0 0 0 0 0 0 -10/5];
J=D*E*inv(D);
for k=1:sample_num
    for i=1:N-1
        ts(i+1)=ts(i)+delta_t;
        eJ=e^(J*delta_t);
        for jj=1:8
            X(jj,i+1)=eJ(jj,:)*X(:,i)+sigma*normrnd(0,1)*delta_t;
        end
    end
   reference_data(:,k)=X(:,2000);
end


TT=200;  %模拟次数
single_sample_num=1;

for s=1:TT
    for t=1:25
        q(t)=0.96^(1/abs(p(t)));
        E=[-2/5*q(t) 0 0 0 0 0 0 0 ;...
            0 -4/5 0 0 0 0 0 0 ;...
            0 0 -5/5 0 0 0 0 0 ;...
            0 0 0 -6/5 0 0 0 0 ;...
            0 0 0 0 -7/5 0 0 0 ;...
            0 0 0 0 0 -8/5 0 0 ;...
            0 0 0 0 0 0 -9/5 0 ;...
            0 0 0 0 0 0 0 -10/5];
        J=D*E*inv(D);
        
        for k=1:single_sample_num
            for i=1:N-1
                ts(i+1)=ts(i)+delta_t;
                eJ=e^(J*delta_t);
                for jj=1:8
                    X(jj,i+1)=eJ(jj,:)*X(:,i)+sigma*normrnd(0,1)*delta_t;
                end
            end
            single_sample=X(:,2000);
        end
        
        for na=1:8
            edges_list=[];
            neighbour_list=adjacent_network(na,2:6);
            num=0;
            sPCC_num=[];
            node_num=[];
            for n1=2:6
                nei=adjacent_network(na,n1);
                if nei==0
                    break;
                end
                num=num+1;
                edges_list(num,:)=[na nei];
                center_mean=mean(reference_data(na,:));
                nei_mean=mean(reference_data(nei,:));
                center_var=var(reference_data(na,:));
                nei_var=var(reference_data(nei,:));
                delt_center_mean=single_sample(na)-center_mean;
                delt_nei_mean=single_sample(nei)-nei_mean;
                sPCC_num(num)=abs(delt_center_mean*delt_nei_mean)/sqrt(center_var*nei_var+es);
                node_num(num)=exp(-delt_nei_mean^2/(2*nei_var+es));
            end
            sPCC_num=sPCC_num/sum(sPCC_num);
            node_num=node_num/sum(node_num);
            sPCC_entropy=-(1/log(num))*sum(sPCC_num.*log(sPCC_num+es));
            node_entropy=-(1/log(num))*sum(node_num.*log(node_num+es));
            Entropy_matrix(na,t,s)=abs(node_entropy-sPCC_entropy);
            Entropy_list(na)=abs(node_entropy-sPCC_entropy);
        end
        Entropy_list=fillmissing(Entropy_list,'constant',0);
        Entropy_list=sort(Entropy_list,'descend');
        Entropy(s,t)=sum(Entropy_list(1:4));
        t,Entropy(s,t)
    end
    
end

E_CI=mean(Entropy);
plot(p(1:25),E_CI,'LineWidth',2,'color','red');






