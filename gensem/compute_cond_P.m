%compute conditional probability  p(C|f)

%for conc = 1:517

function [pCf pfC] = compute_cond_P(feat_matrix,feat_vector,cl,target_class)

%feat_matrix = L4;

%feat_vector = feat_matrix(conc,:) > 0.5;%abs(normrnd(0,1));
feat_vector = reshape(feat_vector,1,[]);

%cl = domain + 1;
%cl = kmeans(feat_matrix,2,'distance','cosine');

if length(find(cl == target_class)) > 1
    
    F_probs = mean(feat_matrix(find(cl == target_class),:));

else
    
    F_probs = (feat_matrix(find(cl == target_class),:)); 

end

if length(find(cl ~= target_class)) > 1
    
    not_F_probs = mean(feat_matrix(find(cl ~= target_class),:));

else
   
    not_F_probs = (feat_matrix(find(cl ~= target_class),:));

end





%L1 = 1- L1;prob = 1; 

P = zeros(length(feat_vector),1);
not_P = zeros(length(feat_vector),1);

for n=1:length(feat_vector); 
    
    if feat_vector(1,n) == 1; 
        
        P(n) = F_probs(n);
        not_P(n) = not_F_probs(n);
        
    else
        
        P(n) = 1 - F_probs(n); 
        not_P(n) = 1- not_F_probs(n);
        
    end
       
end

pfC = prod(P);
pfnC = prod(not_P);

%p(C) 

C = length(find(cl == target_class))/length(cl);
notC = length(find(cl ~= target_class))/length(cl);

pCf = (C * pfC)/[(C * pfC) + (notC * pfnC)];
%end

