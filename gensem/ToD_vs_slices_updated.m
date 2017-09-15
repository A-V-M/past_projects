 
%set up vox indices -> only sig. semantic greater than visual within ventral
%stream

p_map = spm_vol('pMap_RSA_gen_model_testing_10mm_imask01_regAnalysis_C1_reggedOut.img ');p_map = spm_read_vols(p_map);
 
 beta_C1 = spm_vol('betaMap_RSA_gen_model_testing_10mm_imask01_regAnalysis_betaMap_[semPlusC1]C1.img'); beta_C1 = spm_read_vols(beta_C1);
 
 beta_categories = spm_vol('betaMap_RSA_gen_model_testing_10mm_imask01_regAnalysis_betaMap_[semPlusC1]categories.img'); beta_categories = spm_read_vols(beta_categories);
 
 beta_semnorms = spm_vol('betaMap_RSA_gen_model_testing_10mm_imask01_regAnalysis_betaMap_[semPlusC1]semantic_norms.img'); beta_semnorms = spm_read_vols(beta_semnorms);
 
 ventral_stream = spm_vol('rventral_stream_mask.img'); ventral_stream = spm_read_vols(ventral_stream);
 
 ccMap = spm_vol('ccMap_OBJ_MVPA_RSA_gen_model_testing_10mm_imask01.img'); ccMap = spm_read_vols(ccMap);
 
 p_vals = reshape(p_map,[],1);  
        
 beta_C1_vals = reshape(beta_C1,[],1); 
        
 beta_categories_vals = reshape(beta_categories,[],1); 
        
 beta_semnorms_vals = reshape(beta_semnorms,[],1); 
        
 ventral_stream_vals = reshape(ventral_stream,[],1);
        
 

savePath = '/work/imaging7/Functional/Object_PCA_MVPA/Object_MVPA_exp/RSA_searchlight_103/imask01_131concepts_searchlights/RDMs/compromiseRDMs/slice_analysis';

% Parameters %

measure_vals = average_layer_uncertainty([1 4],items_index)';

slice_range = [1:63];

lateral_slice_range = [1:53];

min_slice_range = 0

sem_pval = FDR(reshape(p_map,[],1),0.01)

min_vox = 20

nPerms = 0;

% set-up filter Map

indices = logical((p_vals < sem_pval) .* (abs(beta_C1_vals) < abs(beta_semnorms_vals))); indices = logical(indices .* (ventral_stream_vals > 0.9));
 
filterMap = reshape(indices,(size(p_map)));

% %
trendName = [ num2str(slice_range(1)) 'to' num2str(slice_range(end)) '_VS_UNCL1_L4_pFDR01_min' num2str(min_slice_range) '_minvox' num2str(min_vox)]


clear distSum_vals_mean distSum_vals p_vals r p norm_regress cc_vals slice_matrix* slice_windows* regress_distSum_vs_cat_nPerm
eval(['clear r_' trendName 'norm_regress_' trendName 'p_' trendName]);

% extract average mean within distance val for each concept for each slice
% [131 concepts x 63 slices]

distSum_slice_vals = cell(63,131);

num_of_voxs = zeros(1,63);
cc_vals = zeros(1,63);


for n = 1:63; 
        
        slice_vox_indices= logical(reshape(squeeze(filterMap(lateral_slice_range,n,:)),[],1));
        
        num_of_voxs(n) = length(find(slice_vox_indices));
        
        cc_slice = (reshape(squeeze(ccMap(lateral_slice_range,n,:)),[],1));
        
        cc_slice = cc_slice(~isnan(cc_slice));
        
        cc_vals(n) = mean(cc_slice);
        
       % slice_index = [slice_index ones(1,length(find(slice_vox_indices))).*n];
        
    for o=1:131
        
        distSum_vals = reshape(squeeze(distSum_cat_map(o,lateral_slice_range,n,:)),[],1);  
                
        vals = distSum_vals(slice_vox_indices);
        
        distSum_slice_vals{o,n} = [vals]; 
          
    end;
      
end

num_of_voxs_within_range = num_of_voxs(slice_range);
slice_range = slice_range(num_of_voxs_within_range > min_vox);

fprintf(1,'Slice range = %d - %d (%d slices)\n', slice_range(1),slice_range(end),length(slice_range));

%compute ToD (trend of dissociation)


c=combnk(slice_range,2); 
c=c(abs(c(:,1)-c(:,2))>=min_slice_range,:);

slice_corr_matrix = nan(size(measure_vals,2) ,63,63);
slice_p_matrix = nan(size(measure_vals,2), 63,63);

for n=1:length(c); 
    vals = cc_vals(c(n,1):c(n,2)); 
    ccVal_corr(n) = corr((1:length(vals))',vals');
end

tt=0;
for nSrange = 1:length(c)
tic    
        Srange = c(nSrange,1):c(nSrange,2);
        
        total_num_of_voxs(nSrange) = sum(num_of_voxs(Srange));
        
        slice_index=[]; tnum_voxs_srange = [];
        for n=Srange; slice_index = [slice_index ones(1,num_of_voxs(n)).*n]; tnum_voxs_srange = [tnum_voxs_srange ones(1,num_of_voxs(n)).*total_num_of_voxs(n)]; end

    for n=1:131; 
    
    dists=[];
    for iS = 1:length(Srange); dists = [dists distSum_slice_vals{n,Srange(iS)}'];end    
        
    Y = (dists'./norm(dists));
    
    %X = [(1:length(Srange))'./norm((1:length(Srange))) num_of_voxs(Srange)'./norm(num_of_voxs(Srange))];
    X = [slice_index'./norm(slice_index)];% tnum_voxs_srange'./norm(tnum_voxs_srange)];
    
    s=regstats(Y,X); 
    
    B = s.tstat.beta(2:end); 
    
    %p = s.tstat.pval(2:end); 
    
    %betas=[B p]; 
    
    regress_distSum_vs_cat(n) = B(1) * s.rsquare; 
    
    for nPerm=1:nPerms
        
        s=regstats(Y,X(randperm(length(X))));
        B = s.tstat.beta(2:end); 
        regress_distSum_vs_cat_nPerm(n,nPerm) = B(1) * s.rsquare; 
        
    end
    
    
    end

% correlate against measure
[r p]=corr(measure_vals,regress_distSum_vs_cat');

slice_corr_matrix(:,c(nSrange,1),c(nSrange,2)) = r; 
slice_corr_matrix(:,c(nSrange,2),c(nSrange,1)) = slice_corr_matrix(:,c(nSrange,1),c(nSrange,2));

slice_p_matrix(:,c(nSrange,1),c(nSrange,2)) = p;
slice_p_matrix(:,c(nSrange,2),c(nSrange,1)) = slice_p_matrix(:,c(nSrange,1),c(nSrange,2));

if nPerms > 0

[rPerm ] = corr(measure_vals,regress_distSum_vs_cat_nPerm);
p = sum(r < rPerm) / length(rPerm);

end

norm_regress = sum(regress_distSum_vs_cat);

t=toc;tt = tt+t;
fprintf(1,'%d [%d -> %d - r-val: %1.3f, %1.3f / p: %3.2f, %3.2f  -  %% done: %3.2f (%3.2f)]\n',nSrange,Srange(1),Srange(end),r(1),r(2),p(1),p(2),(nSrange/length(c)*100),tt );
%eval(['r_' trendName '(nSrange) = r; p_' trendName '(nSrange) = p;']);
%eval(['norm_regress_' trendName '(nSrange) = norm_regress;']);
%eval(['cc_lin_corr_' trendName '(nSrange) = cc_lin_corr;']);
        
end

eval(['slice_corr_matrix_' trendName ' = slice_corr_matrix; slice_p_matrix_' trendName ' = slice_p_matrix;']);


break

eval(['r = (r_' trendName '); p = (p_' trendName ');']);
eval(['norm_regress = (norm_regress_' trendName ');']);
eval(['cc_lin_corr = (cc_lin_corr_' trendName ');']);

[y i] = max(r) 
eval(['slice_range_' trendName ' = slice_range(slice_windows(i,:) == 1)'])

slI_r = find((r>0) .* (p<1));
sl_r = slice_windows .* repmat(r',1,length(slice_range));

slI_norm = find(isnan(norm_regress) == 0);
sl_norm = slice_windows .* repmat(norm_regress',1,length(slice_range));

slI_cc = find(isnan(cc_lin_corr) == 0);
sl_cc = slice_windows .* repmat(cc_lin_corr',1,length(slice_range));

slice_vals_r = zeros(1,63);
slice_vals_norm = zeros(1,63);
slice_vals_cc = zeros(1,63);

slice_vals_r(slice_range) = mean(sl_r(slI_r,:));
slice_vals_norm(slice_range) = mean(sl_norm(slI_norm,:));
slice_vals_cc(slice_range) = mean(sl_cc(slI_cc,:));

figure;bar([slice_vals_r; slice_vals_norm; slice_vals_cc]);

sliceMap_r = zeros(53,63,37);
sliceMap_norm = zeros(53,63,37);
sliceMap_cc = zeros(53,63,37); 

for n=1:63; 
    
    sliceMap_r(lateral_slice_range,n,:) = slice_vals_r(n); 
    sliceMap_norm(lateral_slice_range,n,:) = slice_vals_norm(n); 
    sliceMap_cc(lateral_slice_range,n,:) = slice_vals_cc(n);
    
    [y i] = max(r);
    
    
      
end

A = spm_vol('pMap_RSA_gen_model_testing_10mm_imask01_regAnalysis.img ');

cd(savePath);

sliceMap_vec =  reshape(sliceMap_r,[],1); p_map_vec = reshape(p_map,[],1); beta_C1_vec = reshape(beta_C1,[],1); beta_semnorms_vec = reshape(beta_semnorms,[],1);
ventral_stream_vec = reshape(ventral_stream,[],1); 

filteredMap = ((sliceMap_vec .* (p_map_vec < sem_pval)) .* (beta_C1_vec < beta_semnorms_vec)) .* (ventral_stream_vec > 0.9);
sliceMap = reshape(filteredMap,53,63,37);

A.fname = [ 'slices_cLengthXtrend_' trendName '.img '];
spm_write_vol(A, sliceMap);

sliceMap_vec =  reshape(sliceMap_norm,[],1); p_map_vec = reshape(p_map,[],1); beta_C1_vec = reshape(beta_C1,[],1); beta_semnorms_vec = reshape(beta_semnorms,[],1);
ventral_stream_vec = reshape(ventral_stream,[],1); 

filteredMap = ((sliceMap_vec .* (p_map_vec < sem_pval)) .* (beta_C1_vec < beta_semnorms_vec)) .* (ventral_stream_vec > 0.9);
sliceMap = reshape(filteredMap,53,63,37);

A.fname = [ 'slices_trend_' trendName '.img '];
spm_write_vol(A, sliceMap);


sliceMap_vec =  reshape(sliceMap_cc,[],1); p_map_vec = reshape(p_map,[],1); beta_C1_vec = reshape(beta_C1,[],1); beta_semnorms_vec = reshape(beta_semnorms,[],1);
ventral_stream_vec = reshape(ventral_stream,[],1); 

filteredMap = ((sliceMap_vec .* (p_map_vec < sem_pval)) .* (beta_C1_vec < beta_semnorms_vec)) .* (ventral_stream_vec > 0.9);
sliceMap = reshape(filteredMap,53,63,37);

A.fname = [ 'slices_cc_trend_' trendName '.img '];
spm_write_vol(A, sliceMap);

sliceMap = zeros(53,63,37);
[y i] = max(cc_lin_corr);
sliceMap(lateral_slice_range,slice_range(slice_windows(1,:)==1),:) = y;
sliceMap_vec =  reshape(sliceMap,[],1); p_map_vec = reshape(p_map,[],1); beta_C1_vec = reshape(beta_C1,[],1); beta_semnorms_vec = reshape(beta_semnorms,[],1);
ventral_stream_vec = reshape(ventral_stream,[],1); 

filteredMap = reshape(ccMap,[],1) .* (((sliceMap_vec .* (p_map_vec < sem_pval)) .* (beta_C1_vec < beta_semnorms_vec)) .* (ventral_stream_vec > 0.9));
sliceMap = reshape(filteredMap,53,63,37);


A.fname = [ 'slices_cc_MAX_trend_' trendName '.img '];
spm_write_vol(A, sliceMap);

cd('..');














