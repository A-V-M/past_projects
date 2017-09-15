% test battery for 4-layer DBM

runs = 25;

obj_index = 1:517;


L = { L1{1} L2{1} L3{1} L4{1}};

layer_uncertainty = zeros(4,runs,517);

layer_uncertainty_cat = zeros(4,runs,517);

layer_uncertainty_dom = zeros(4,runs,517);

layer_condProb = zeros(4,runs,517,517); 

pCat = zeros(4,runs,517,24);

pDom = zeros(4,runs,517,2);

predBL = zeros(4,runs,517);

predCat = zeros(4,runs,517);

predDom = zeros(4,runs,517);

for q=1:runs
    tic
    for nLayer = 1:4;
        
        for obj = obj_index; 
            
            clear p 
            
            conc=obj; 
            inp_vec=(L{nLayer}(conc,:))>rand;
            dataMat = L{nLayer};
            
            for c = 1:517; 
                n = c; 
                cl = zeros(1,517); 
                cl(n) = 1; 
                p(c)= compute_cond_P(dataMat,(inp_vec),cl,1); 
                r(c)=corr(dataMat(c,:)',(inp_vec'));
            end; 
            
            condProb = p; 
            p = p(~isnan(p)); 
            p = p(p~=0);
            p=p./sum(p);
            
            uncertainty =-sum(p.*log(p)); 
            
            layer_condProb(nLayer,q,conc,:) = condProb; % conditional probability 
            
            [y i] = max(condProb);
            
            predBL(nLayer,q,conc) = i; % basic-level prediction
            
            layer_uncertainty(nLayer,q,conc) = uncertainty; % uncertainty
            
            pCat(nLayer,q,conc,:) = condProb*cat_matrix; % probability category
            
            pC = squeeze(pCat(nLayer,q,conc,:));
            
            pC = pC(~isnan(pC) & pC~=0); pC = pC./sum(pC);
                        
            uncertainty_cat = -sum(pC.*log(pC)); 
            
            layer_uncertainty_cat(nLayer,q,conc) = uncertainty_cat; % uncertainty category
                
            [y i] = max(pCat(nLayer,q,conc,:));
            
            predCat(nLayer,q,conc) = i; % predicted domain
            
            pDom(nLayer,q,conc,:) = condProb*domain_matrix;
            
            pD = squeeze(pDom(nLayer,q,conc,:));
            
            pD = pD(~isnan(pD) & pD~=0); pD = pD./sum(pD);
                        
            uncertainty_dom = -sum(pD.*log(pD)); 
            
            layer_uncertainty_dom(nLayer,q,conc) = uncertainty_dom; % uncertainty domain
                
            [y i] = max(pDom(nLayer,q,conc,:));
            
            predDom(nLayer,q,conc) = i;
            
                
        end;
        
        
        
        
    end
    
    disp(['Run: ' num2str(q)]);toc
    
    
end



for l=1:4; 
    for r=1:runs; 
        for o=1:517; 
            
            accuracy_layer(l,r,o) = (predBL(l,r,o) == o);
            accuracy_layer_cat(l,r,o) = (predCat(l,r,o) == find(cat_matrix(o,:)));
            accuracy_layer_dom(l,r,o) = (predDom(l,r,o) == find(domain_matrix(o,:)));
            
        end;
    end;
end


for n=1:4;

    average_layer_accuracy(n,:) = mean(squeeze(accuracy_layer(n,:,:)));
    average_layer_accuracy_cat(n,:) = mean(squeeze(accuracy_layer_cat(n,:,:)));
    average_layer_accuracy_dom(n,:) = mean(squeeze(accuracy_layer_dom(n,:,:)));
    
    average_layer_uncertainty(n,:) = mean(squeeze(layer_uncertainty(n,:,:)));
    average_layer_uncertainty_cat(n,:) = mean(squeeze(layer_uncertainty_cat(n,:,:)));
    average_layer_uncertainty_dom(n,:) = mean(squeeze(layer_uncertainty_dom(n,:,:)));
                
end






