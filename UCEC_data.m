-
clear;
clc;
close all;
[profile,pipi]= xlsread('UCEC_good.xlsx');

fid=fopen('UCEC_Gene_network.txt'); % 改成对应的文件名
adjacent_network={};
j=0;
while ~feof(fid)
    tline=fgetl(fid);
    j=j+1;
    adjacent_network{j}=regexp(tline, '\t', 'split');
end
fclose(fid);
total_node_num=j;
normal=profile(:,1:35);
case_mprofile=profile(:,36:463);

reference_sample_num=35;
case_mprofile=fillmissing(case_mprofile,'constant',0);
es=0.0000000001;
patients_num=[167,146,25,6,13,40,6,25];

tempcase(:,1,1:patients_num(1))=case_mprofile(:,1:167);     % Stage I 
tempcase(:,2,1:patients_num(2))=case_mprofile(:,168:313);   % Stage II
tempcase(:,3,1:patients_num(3))=case_mprofile(:,314:338);   % Stage III
tempcase(:,4,1:patients_num(4))=case_mprofile(:,339:344);   % Stage IV
tempcase(:,5,1:patients_num(5))=case_mprofile(:,345:357);
tempcase(:,6,1:patients_num(6))=case_mprofile(:,358:397);
tempcase(:,7,1:patients_num(7))=case_mprofile(:,398:403);
tempcase(:,8,1:patients_num(8))=case_mprofile(:,404:428);



psize=size(tempcase);
Entropy=zeros(psize(1),167,8);
stage=8;


for t=1:psize(2)
    for  s=1:patients_num(t)
         single_sample=tempcase(:,t,s);
         for na=1:total_node_num
             center=adjacent_network{na}{1};
             num=0;
             sPCC_num=[];
             node_num=[];
             for n=2:length(adjacent_network{na})
                 nei=adjacent_network{na}{n};
                 num=num+1;
                 center_mean=mean(normal(str2num(center),:));
                 nei_mean=mean(normal(str2num(nei),:));
                 center_var=var(normal(str2num(center),:));
                 nei_var=var(normal(str2num(nei),:));
                 delt_center_mean=single_sample(str2num(center))-center_mean;
                 delt_nei_mean=single_sample(str2num(nei))-nei_mean;
                 sPCC_num(num)=abs(delt_center_mean*delt_nei_mean)/sqrt(center_var*nei_var+es);
                 node_num(num)=exp(-delt_nei_mean^2/(2*nei_var+es));
             end
             if num<2
                Entropy(na,s,t)=0;
                continue;
             else
                 sPCC_num=sPCC_num/sum(sPCC_num);
                 node_num=node_num/sum(node_num);
                 sPCC_entropy=-(1/log(num))*sum(sPCC_num.*log(sPCC_num+es));
                 node_entropy=-(1/log(num))*sum(node_num.*log(node_num+es));
                 Entropy(na,s,t)=abs(node_entropy-sPCC_entropy);
             end
         end
    end
    t
end

save UCEC_entropy.mat;
Entropy=fillmissing(Entropy,'constant',0);
Entropy_size=size(Entropy);
case_result=zeros(167,8);
for t=1:Entropy_size(3)
    for case_num=1:patients_num(t)
        [sort_entropy,idx]=sort(Entropy(:,case_num,t),'descend');
        case_result(case_num,t)=mean(sort_entropy(1:0.05*Entropy_size(1)));
    end
    result(t)=mean(case_result(1:patients_num(t),t));
end

t=[1 2 3 4 5 6 7 8];
plot(t,result,'r','LineWidth',3);
set(gca,'XTick',1:8);
B={'IA'  'IB' 'IC' 'IIA' 'IIB' 'IIIA' 'IIIB' 'IV'};
set(gca,'XTickLabel',B);
% ylim([0.08,0.12])
% set(gca,'Ylim',[0.08 0.1]);
% set(gca,'YTick',[0.08:0.01:0.1]);
xlabel('Stages');
ylabel('Entropy');
%plot(t,aver_comidx,'r','LineWidth',3);
title('Average Entropy for UCEC');



%%

%%%%%%%%% selected marker genes %%%%%%%%%%%%%%%%%%%%%%%%%%%%
for na=1:total_node_num
    center=adjacent_network{na}{1};
    gene_name(na)=pipi(str2num(center));
end
Land_entropy_size=size(Entropy);
patient_num=[167,146,25,6,13,40,6,25];
selected_deltaH_genes=[];
fid=fopen('selected_UCEC_genes.txt','w');
selected_deltaH_genes=containers.Map;
for s=1:patient_num(5)
    [tmp_com_idx,idx]=sort(Entropy(:,s,5),'descend');
    %aver_com_idx(l)=mean(tmp_com_idx(1:floor(total_node_num*0.01)));
    tmp_genes=gene_name(idx(1:floor(Land_entropy_size(1)*0.05)));
    for i=1:length(tmp_genes)
        if isKey(selected_deltaH_genes,tmp_genes{i})==0
            selected_deltaH_genes(tmp_genes{i})=1;
        else
            selected_deltaH_genes(tmp_genes{i})=selected_deltaH_genes(tmp_genes{i})+1;
        end
    end
end
genes=selected_deltaH_genes.keys();
for i=1:length(genes)
    fprintf(fid,'%s\t%d\n',genes{i},selected_deltaH_genes(genes{i}));
end
fclose(fid);








