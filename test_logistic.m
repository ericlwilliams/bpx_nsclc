function test_logistic
tic;
% prepare

fp = 'Z:\elw\MATLAB\bpx_analy\meta_data\';

screen_size=get(0,'ScreenSize');

do_print = true;
fig_loc = 'Z:\elw\MATLAB\bpx_analy\slides\figures\latest\';

% use for complication info
fn = 'BPx_DiVj_DVHs_fx-1_a2bInf.mat';

CGobj = cell(length(fn),1);
screen_size=get(0,'ScreenSize');

% load data
load(strcat(fp,fn),'CGobj_current');
CG = CGobj_current;

LymanN = log10(CG.mLymanN);
CG.mLymanN = LymanN;

% survival/complication time
f2 = ~cellfun('isempty',{CG.mGrp.mDateComp}); % patients with no complication date
f3 = ~cellfun('isempty',{CG.mGrp.mDateLastFollowup}); % patients with no last follow up date
compdate = inf(CG.mNumInGrp,1);
lastfollowup = inf(CG.mNumInGrp,1);
compdate(f2) = ([CG.mGrp(f2).mDateComp] - [CG.mGrp(f2).mDateBaseline])' / 30;
lastfollowup(f3) = ([CG.mGrp(f3).mDateLastFollowup] - [CG.mGrp(f3).mDateBaseline])' / 30;
compdate = min( lastfollowup, compdate );
flgcensor = [CG.mGrp.mFlgCensor]';

pttotal = ones(CG.mNumInGrp,1);
ptcomp = ones(CG.mNumInGrp,1); 
ptcomp([CG.mGrp.mFlgCensor])=0;
    
log10a = -CG.mLymanN;
% gEUD log-likelihood    
euds = [CG.mGrp.mEUD]';
sa = CG.mKaplanMeierCompOverall;

loglikelihood = -inf(length(CG.mLymanN),1);
pvals = -inf(length(CG.mLymanN),1);
nr_loglikelihood = -inf(length(CG.mLymanN),1);
mm_loglikelihood = -inf(length(CG.mLymanN),1);

for k=1:size(CG.mLymanN,1) % loop over each logn value, calculate fit of gEUDs
    doses=euds(:,k);
    
    %% Newton-Raphson
    [b,~,s]=glmfit(doses,[ptcomp pttotal],'binomial','link','logit');
    pvals(k) = s.p(2);
    
    B0 = b(1);
    B1 = b(2);
    pr = exp(B0+B1*euds(:,k));
    
    pr = pr./(1+pr); % logistic probability
    pr(flgcensor) = 1-pr(flgcensor); % non-complication patients
    pr = log(pr); % log likelihood of each patients
    loglikelihood(k) = sum(pr); % loglikelihood of all
    
    %% MLE
    mle_data = [doses ptcomp];
    mle_test_b = b;
    [~,mle_llhd] = MLE_logistic(mle_data,CG,mle_test_b);
    mle_loglikelihood(k) = mle_llhd; % loglikelihood of all
    
     %% Mixed-Model
    mm_data = [doses ptcomp];
    mm_test_b = b;
    %mm_test_b = [0 0]';
    [~,mm_llhd] = MM_logistic(mm_data,CG,mm_test_b,(do_print && k==1));
    mm_loglikelihood(k) = mm_llhd; % loglikelihood of all
        
end


cur_fig=figure(1);clf reset;hold on;
set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);

h_ll=plot(log10a,loglikelihood,'LineWidth',2);
h_mle_ll=plot(log10a,mle_loglikelihood,'r','LineWidth',2);
h_mm_ll=plot(log10a,mm_loglikelihood,'g','LineWidth',2);
%ylim([-39 -29]);
xlabel('log_1_0(a)','FontSize',14);
ylabel('Log-likelihood','FontSize',14);
lgnd=legend([h_ll h_mle_ll h_mm_ll],'NR LogL','MLE LogL','MM LogL',...
    'Location','Best');
set(lgnd,'FontSize',14);
title('Logistic Regression Log-likelihoods','FontSize',14);
grid on;
 set(gca,'LineWidth',2);% grid lines visible on pdf export
 set(gca,'FontSize',12);

if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_loc,'test_logreg_llhds'],'-pdf','-painters');
end;

end
