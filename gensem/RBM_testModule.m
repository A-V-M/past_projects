
%
%
% Modified AVM/CSLB 2014


function [L DM_L L_hidpen L_hidbiases L_hidgenbiases secondOrder_DMs cc] = RBM_testModule(nLayers,num_units,epochs,LR,inp_data)

%clear all
%close all

%In the Science paper we use maxepoch=50, but it works just fine. 
%numhid=1000; numpen=500; numpen2=250; numopen=30;
%numhid=750; numpen=750; numpen2=750; numpen3=750; numpen4=750; numpen5=750; numopen=750;

fprintf(1,'Pretraining a deep autoencoder. \n');

% makebatches;

load McRae_feats_MoBo 
load animals_tools.mat
load cat_matrix.mat

%inp_indices = randperm(517); inp_indices = inp_indices(1:100);

%inp_data = McRae';
batchdata = inp_data;

[numcases numdims]=size(batchdata);

for nLayer = 1:nLayers

fprintf(1,'Pretraining Layer %d with RBM: %d-%d \n',nLayer,numdims,num_units(nLayer));
restart=1;
maxepoch = epochs(nLayer); numhid=num_units(nLayer);
epsilonw      = LR(nLayer);   % 0.1 Learning rate for weights 
epsilonvb     = LR(nLayer);   % 0.1 Learning rate for biases of visible units 
epsilonhb     = LR(nLayer);
rbmhidlinear; 
L{nLayer} = batchposhidprobs;
DM_L{nLayer} = squareform(pdist(L{nLayer},'cosine'));
L_hidpen{nLayer}=vishid; 
L_hidbiases{nLayer}=hidbiases; 
L_hidgenbiases{nLayer}=visbiases;

%save mnistvh vishid hidrecbiases visbiases;

end


% fprintf(1,'\nPretraining top layer with RBM: %d-%d \n',num_units(end),num_units(end));
% batchdata=batchposhidprobs;
% numhid=num_units(end); 
% restart=1;
% rbmhidlinear; 
% L{end+1} = batchposhidprobs;
% DM_L{end+1} = squareform(pdist(L{end},'cosine'));
% 
% hidtop=vishid; toprecbiases=hidbiases; topgenbiases=visbiases;
% save mnistpo hidtop toprecbiases topgenbiases;


%produce DMs
DMs = (pdist(inp_data,'cosine'))';
    
for nDM = 1:length(L)
    
DMs = [DMs squareform(DM_L{nDM})'];

end

secondOrder_DMs=(1-corr(DMs,'type','spearman'));
coords = mdscale(secondOrder_DMs,2);
figure;scatter(coords(:,1),coords(:,2),500,'.')

DM_labels{1} = 'InpLayer';

for n = 1:(nLayers+1); 
    
    DM_labels{end+1} = ['L' num2str(n)]; 
    
end

for n=1:length(coords);
    
    cc(n) = ccohesion_pw(squareform(DMs(:,n)),cl_all_categories);
    
    text(coords(n,1),coords(n,2),DM_labels{n},'FontWeight','bold');
end

%backprop; 

