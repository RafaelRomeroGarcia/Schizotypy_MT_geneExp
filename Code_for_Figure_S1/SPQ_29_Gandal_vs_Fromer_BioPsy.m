clear

%Gandal differential schizophrenia expression data downloaded from
%https://science.sciencemag.org/content/suppl/2018/12/12/362.6420.eaat8127.DC1?_ga=2.229195356.1935284581.1565860151-2106640697.1565860151
%Table S1 converted from .xls to .csv
gandal_csv=importdata('aat8127_Table_S1.csv');
gandal_geneNames=gandal_csv.textdata(2:end,8);
gandal_difExp=gandal_csv.data(:,12);

%Fromer differential schizophrenia expression data downloaded from
%https://www.nature.com/articles/nn.4399#supplementary-information (Supplementary Data File 3) and converted from .xls to .csv 
fromer_csv=importdata('nn.4399-S5.csv');
fromer_geneNames=fromer_csv.textdata(3:end,2);
fromer_difExp=fromer_csv.data(:,1);

%Iterate across genes (each dot represent a gene)
x=[];
y=[];
num_valid=1;
for ig=1:numel(fromer_geneNames)
    %Look for the position of the gene on gandal list    
    gene_pos=find(strcmp(gandal_geneNames,fromer_geneNames{ig}));
    %gene_pos=find(not(cellfun('isempty',gene_pos_cell)));
    if not(isempty(gene_pos))
        x(num_valid)=gandal_difExp(gene_pos(1));
        y(num_valid)=fromer_difExp(ig);
        num_valid=num_valid+1;
    end
end

%Plot association
figure
hold on
scatter(x,y,20,[0 0 0],'filled');
scatter(x,y,10,[0.6 0.6 0.6],'filled');
f=polyfit(x,y,1);
plot(x,f(2)+f(1)*x,'linewidth',2,'color',[0 0 0]);
plot(x,f(2)+f(1)*x,'linewidth',1,'color',[0.2 0.2 0.2]);
set(gca,'fontsize',22);
set(gca,'TickDir','out');
set(gca,'LineWidth',3);
xlabel('Differential expression (Gandal)');
ylabel('Differential expression (Fromer)');

%Calculate correlation
[r p]=corrcoef(x,y);
display(['Correlation between Gandal and Fromer expression R2=' num2str(r(1,2)*r(1,2)) ' p= '  num2str(p(1,2))]);