function [stimulus_entropy equivocation I category_entropy equivocation_categories I_categories] = info_content2(layer_matrix,cat_matrix)

    mean_probs=mean(layer_matrix)'; 
   
    num_responses = size(layer_matrix,1);
    num_features = size(layer_matrix,2);
    num_categories = size(cat_matrix,2);
    
    pResponse = zeros(1,num_responses);
    inp = zeros(num_responses,num_features);
    
    pStimulus = ones(1,num_responses) .* (1/num_responses); % we assume here that all stimuli are equiprobable
    pCategory = ones(1,num_categories) .* (1/num_categories);
    
    stimulus_entropy = -sum(pStimulus.*log(pStimulus));
    category_entropy = -sum(pCategory.*log(pCategory));
            
% computing p(r)

for r=1:size(layer_matrix,1); 
    
    inp(r,:)=layer_matrix(r,:)>rand; 

    prob_vector = zeros(1,size(layer_matrix,2)); 
    
    prob_vector(inp(r,:) == 1) = mean_probs(inp(r,:) == 1); 
    
    prob_vector(inp(r,:) == 0) = 1 - mean_probs(inp(r,:) == 0); 
    
    pResponse(r) = prod(prob_vector); 
    
end

pResponse=(pResponse./sum(pResponse)); % p(r)

pResponse = repmat(pResponse',1,num_responses);

% computing p(r|S) for all stimuli and response combinations --> num
% responses X num stimuli

p_sR = zeros(num_responses,num_responses);
p_Rs = zeros(num_responses,num_responses);


for r=1:num_responses; 
    
    for c=1:num_responses; 
        n = c; 
        cl = zeros(1,517); 
        cl(n) = 1; 
        [p_sR(r,c) p_Rs(r,c)]= compute_cond_P(layer_matrix,(inp(r,:)),cl,1); 
    end;
    
end

%for individual stimuli
T1 = (pResponse .* p_sR);
T2 = log(p_sR);
T3 = T1 .* T2;
T3 = reshape(T3,1,[]);

equivocation = -sum(T3(~isnan(T3)));

I = stimulus_entropy - equivocation;



%for individual categories

pResponse = repmat(pResponse(:,1),1,num_categories);

cat_p_sR = p_sR * cat_matrix;

for n=1:size(cat_matrix,2); 
    cat_p_sR(:,n) = cat_p_sR(:,n) ./ sum(cat_p_sR(~isnan(cat_p_sR(:,n)),n)); 
end

T1 = (pResponse .* cat_p_sR);
T2 = log(cat_p_sR);
T3 = T1 .* T2;
T3 = reshape(T3,1,[]);

equivocation_categories = -sum(T3(~isnan(T3)));

I_categories = category_entropy - equivocation_categories;






