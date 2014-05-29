function LoadAlpha2BetaDoses
tic;
% prepare
fp = 'Z:\elw\MATLAB\bpx_analy\meta_data\';

fn = {'BPx_DiVj_DVHs_fx-1_a2bInf.mat'};
a2b_corr='BED';

CGobj = cell(length(fn),1);
screen_size=get(0,'ScreenSize');


% load data
for m = 1:length(fn)
    load(strcat(fp,fn{m}),'CGobj_current');
    CGobj{m} = CGobj_current;
end

for m = 1:length(fn)
    CG = CGobj{m};
    for n=1:length([CG.mGrp]);
        CG.mGrp(n).mBeta2AlphaCorrection = a2b_corr;
    end
    
    a2b_range = [.1:.1:30];
   
   % #pts x #a2b values
    a2b_doses = cell(length(CG.mGrp),length(a2b_range)+1);
    a2b_vols = cell(length(CG.mGrp),length(a2b_range)+1);
    
    %last entry is a2b=Inf
    a2b_doses(:,end) = {CG.mGrp.mDoseBins_LQ}';
    a2b_vols(:,end) = {CG.mGrp.mVolCum}';
    
    n_pt = CG.mNumInGrp;
    
    for k=1:length(a2b_range),
        
        tmpCG = CG;       
        tmpCG.mBeta2Alpha = 1/a2b_range(k);
        a2b_doses(:,k) = {tmpCG.mGrp.mDoseBins_LQ}';
        a2b_vols(:,k) ={tmpCG.mGrp.mVolCum}'; 
        clear tmpCG;
    end
   
     a2b_range = [a2b_range Inf];
  
     
      a2b_d05 = inf(n_pt,length(a2b_range));
      a2b_d35 = inf(n_pt,length(a2b_range));
      a2b_dmax = inf(n_pt,length(a2b_range));
      a2b_dmean = inf(n_pt,length(a2b_range));

    
        for l=1:length(a2b_range)%loop over all a/b dose bins
            cur_doses = [a2b_doses{:,l}];
            cur_vols = [a2b_vols{:,l}]; %501 x 125
        
            for m=1:n_pt
                vol = cur_vols(:,m);
                dose = cur_doses(:,m);
                a2b_dmean(m,l) = mean(dose);
                
                v05 = vol<5; %% < 5 cc
                d05_inds = find(v05);
                d05 = min(dose(d05_inds));
                a2b_d05(m,l) = d05;
            
                v35 = vol<3.5; %% < 5 cc
                d35_inds = find(v35);
                d35 = min(dose(d35_inds));
                a2b_d35(m,l) = d35;
            
                vol(~vol)=nan;
                nan_inds = find(isnan(vol));
                if ~isempty(nan_inds)
                    min_ind = nan_inds(1)-1;
                else
                    min_ind = length(vol);
                end
            %[~,min_ind] = min(vol);
               a2b_dmax(m,l) = dose(min_ind);
            end
    end

     save([fp 'BPx_a2b_dosebins.mat'],'a2b_dmax','a2b_d05','a2b_d35','a2b_range');
    toc;
end
end