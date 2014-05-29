function BPxIsoEffect
tic;
% prepare
fp = 'Z:\elw\MATLAB\bpx_analy\meta_data\';

screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

   
do_print = true;
fig_loc = 'Z:\elw\MATLAB\bpx_analy\slides\figures\latest\';
 

%fn = {'BPx_DiVj_DVHs_fx-1_a2busc8p6.mat'};dose_calib='uscbed';a2b_str='8p6';

fn = {'BPx_DiVj_DVHs_fx-1_a2bInf.mat'};dose_calib='phys';a2b_str='\infty';
%fn = {'BPx_DiVj_DVHs_fx-1_a2b3.mat'};dose_calib='a2b3';a2b_str='3';


CGobj = cell(length(fn),1);
screen_size=get(0,'ScreenSize');

%scrsz = get(0,'ScreenSize');
%set(0,'DefaultFigurePosition',[scrsz(1)+scrsz(1)/4 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);

% load data
for m = 1:length(fn)
    load(strcat(fp,fn{m}),'CGobj_current');
    CGobj{m} = CGobj_current;
end

for m = 1:length(fn)
    
    CG = CGobj{m};
    CGgrp = CG.mGrp;
    %% Find TD50 for each fractionation
    
    fractions = unique([CGgrp.mFxNum]);
    num_unique_fractions = length(fractions);
    
    td50s = inf(num_unique_fractions,1);
    
    td_threshold = 0.1;
    %% for each fraction size, find TD50
    for i=1:num_unique_fractions
        
       cur_frac = fractions(i);
       
       cur_grp = [CGgrp([CGgrp.mFxNum]==cur_frac)];
       
       if length(cur_grp)<=1
           continue;
       end
       
       pttotal = ones(length(cur_grp),1);
       ptcomp = ones(length(cur_grp),1); 
       ptcomp([cur_grp.mFlgCensor])=0;
       doses = [cur_grp.mDoseTx]./100;
       
       [b,~,st]=glmfit(doses,[ptcomp pttotal],'binomial','link','logit');
       
       if b(2)<0 
           continue;
       end
       rp_doses = [0:max(doses)];
       rpb = zeros(length(rp_doses),1);
        
       dose_range = 1;
       while max(rpb)<td_threshold
        
        cur_doses = [0:(max(rp_doses)*dose_range)];
        [rpb,~,~] = glmval(b, cur_doses,'logit',st); % the responding function values at doses
        dose_range = dose_range+0.1;%increase range by 10%
       end
       
       if ~all(diff(rpb)) % not monotonically increasing
        continue;
       end
       
       cur_td50 = interp1(rpb,cur_doses,td_threshold);
       td50s(i) = cur_td50;
       
    end
   
    figure(1);
    plot(td50s./fractions',1./td50s,'*');
    
   
end