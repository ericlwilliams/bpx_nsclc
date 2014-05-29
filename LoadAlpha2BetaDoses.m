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
    
    %last entry is a2b=Inf
    a2b_doses(:,end) = {CG.mGrp.mDoseBins_LQ}';
    

    
    for i=1:length(a2b_range),
    
        tmpCG = CG;       
        tmpCG.mBeta2Alpha = 1/a2b_range(i);
        a2b_doses(:,i) = {tmpCG.mGrp.mDoseBins_LQ}';
    end
   
     a2b_range = [a2b_range Inf];
  
     a2b_dmax = [cellfun(@(x) max(x),a2b_doses)];

     save([fp 'BPx_a2b_dosebins.mat'],'a2b_dmax','a2b_range');
    toc;
end
end