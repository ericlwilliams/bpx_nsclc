function BPxAlpha2BetaAnaly
tic;
% prepare
fp = 'Z:\elw\MATLAB\bpx_analy\meta_data\';

screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

do_plot_km = false;
do_print = true;
fig_loc = 'Z:\elw\MATLAB\bpx_analy\slides\figures\latest\';

% use for complication info
fn = 'BPx_DiVj_DVHs_fx-1_a2bInf.mat';
%fn2 = 'BPx_a2b_dosebins.mat';
fn2 = 'BPx_a2b_dosebins.mat';

CGobj = cell(length(fn),1);
screen_size=get(0,'ScreenSize');


% load data
load(strcat(fp,fn),'CGobj_current');
CG = CGobj_current;

load(strcat(fp,fn2)); % a2b_max_doses, a2b_range


% survival/complication time
f2 = ~cellfun('isempty',{CG.mGrp.mDateComp}); % patients with no complication date
f3 = ~cellfun('isempty',{CG.mGrp.mDateLastFollowup}); % patients with no last follow up date
compdate = inf(CG.mNumInGrp,1);
lastfollowup = inf(CG.mNumInGrp,1);
compdate(f2) = ([CG.mGrp(f2).mDateComp] - [CG.mGrp(f2).mDateBaseline])' / 30;
lastfollowup(f3) = ([CG.mGrp(f3).mDateLastFollowup] - [CG.mGrp(f3).mDateBaseline])' / 30;
compdate = min( lastfollowup, compdate );
flgcensor = [CG.mGrp.mFlgCensor]';


a2b_var = a2b_dmax;


pttotal = ones(CG.mNumInGrp,1);
ptcomp = ones(CG.mNumInGrp,1); 
ptcomp([CG.mGrp.mFlgCensor])=0;
    
%a2b_max_doses = [cellfun(@(x) max(x),a2b_doses)];

sa=classKaplanMeierCurve(); % initialize a survivalanalysis obj

lr_pvals = Inf(size(a2b_var,2),1);
km_hrs = zeros(size(a2b_var,2),1);

pvals = Inf(size(a2b_var,2),1);

llhds = Inf(size(a2b_var,2),1);
llhds = Inf(size(a2b_var,2),1);
fig_ctr=1;
%loop over all a2b valus
for i=1:size(a2b_var,2)
    cur_a2b = a2b_range(i);
    %cur_a2b_max_doses = [a2b_max_doses(:,i)];
    %%new
    cur_a2b_var = [a2b_var(:,i)];
    
    [b,~,s]=glmfit(cur_a2b_var,[ptcomp pttotal],'binomial','link','logit');
    pvals(i) = s.p(2);
    
    B0 = b(1);
    B1 = b(2);
    pr = exp(B0+B1*cur_a2b_var);
    
    pr = pr./(1+pr); % logistic probability
    pr(flgcensor) = 1-pr(flgcensor); % non-complication patients
    pr = log(pr); % log likelihood of each patients
    llhds(i) = sum(pr); % loglikelihood of all
    
    
    %% KM split
    
    med_a2b_max_doses = median(cur_a2b_var);
    flg_below=cur_a2b_var<=med_a2b_max_doses;
    
    survivedate={compdate(flg_below); compdate(~flg_below)}; % survive time of each group
    fcensor={flgcensor(flg_below); flgcensor(~flg_below)}; % censor flag for each group
    sa.mSurvivalTime=survivedate;
    sa.mFlgCensor=fcensor;
    % compute survival curves and compare them
    sa=sa.fCalculateSurvivalCurve();
    sa=sa.fCombineSurvivalTime();
    sa=sa.fCompareSurvivalByLogrank();
    lr_pvals(i) = sa.mpValue;
    
%     if sum(~flgcensor(flg_below))
%         cox_beta=coxphfit(~flg_below,compdate,'baseline',0,'censoring',flgcensor);
%         km_hrs(i) = exp(cox_beta);
%     end
    cox_beta=coxphfit(~flg_below,compdate,'baseline',0,'censoring',flgcensor);
    km_hrs(i) = exp(cox_beta);

    %      % plot
if (i==1 || i==size(a2b_var,2) || mod(a2b_range(i),1)==0) && do_plot_km,
  
    cur_fig=figure(100+fig_ctr);clf reset;hold on;
    set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    
    h_km(1)=stairs(sa.mSurvivalTimeSorted{1},1-sa.mSurvivalCurve{1});
    plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1)),...
        1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'+');
    h_km(2)=stairs(sa.mSurvivalTimeSorted{2},1-sa.mSurvivalCurve{2},'r');
    plot(sa.mSurvivalTimeSorted{2}(sa.mCensorStatistics{2}(:,1)),...
        1-sa.mSurvivalCurve{2}(sa.mCensorStatistics{2}(:,1)),'r+');
    %         xticks = get(gca,'Xlim'); set(gca,'XTick',0:6:max(xticks));
    
    str_pval1 = ['Log-Rank p-value = ',num2str(sa.mpValue,2),10,...
        '\alpha/\beta = ',num2str(a2b_range(i))];
    text(.6,0.7,str_pval1,...
        'FontSize',16,'Units','normalized','BackgroundColor','w');
    lgnd=legend(h_km,...
        ['D$_{\rm{max}} \leq',num2str(med_a2b_max_doses,4),...
        '~\rm{Gy}_{',num2str(a2b_range(i)),'}$'],...
        ['D$_{\rm{max}} > ',num2str(med_a2b_max_doses,4),...
        '~\rm{Gy}_{',num2str(a2b_range(i)),'}$'],...
        'Location','SouthEast');
    set(lgnd,'FontSize',18);
    h=legend;
    set(h,'interpreter','latex');
    
    set(gca,'xminortick','on','yminortick','on');
    xlabel(['Months'],'fontsize',14);
    ylabel(['Probability of BPx'],'fontsize',14);
    title('Brachial Plexopathy Incidence','fontsize',14);
    
    
    if do_print,
     set(cur_fig,'Color','w');
     export_fig(cur_fig,[fig_loc,'bpx_a2b_km-',num2str(fig_ctr)],'-pdf');
      disp(['Saving ',fig_loc,'bpx_a2b_km-',num2str(fig_ctr),'.pdf...']);
    end;

    
    fig_ctr=fig_ctr+1;
    
end    
end

cur_fig=figure(1);clf reset;hold on;
set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
plot(a2b_range(1:end-1),pvals(1:end-1),'LineWidth',2);% excluding Inf
[min_pval,min_idx] = min(pvals(1:end-1));
text(.1,0.8,...
    ['Min p-value of ',num2str(min_pval),10,'at \alpha/\beta = ',num2str(a2b_range(min_idx))],...
    'FontSize',16,'Units','normalized','BackgroundColor','w');
set(gca,'FontSize',14);
ylabel('p-value (Logistic Regression)','FontSize',16);
xlabel('\alpha/\beta value [Gy]','FontSize',16);
title(['Logistic Regression',10,'BPx vs. Max BED for given \alpha/\beta'],'FontSize',16);
%grid on;

if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_loc,'bpx_a2b_pvals'],'-pdf');
    disp(['Saving ',fig_loc,'bpx_a2b_pvals.pdf...']);
end;

x_axis = a2b_range(1:end-1);
cur_fig=figure(2);clf reset;hold on;
set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
plot(x_axis,llhds(1:end-1),'LineWidth',2);% excluding Inf
[max_llhd,max_idx] = max(llhds(1:end-1));
lowCI68 = max_llhd - 0.5; % 68% confidence
plot([0 max(x_axis)],[lowCI68 lowCI68],'g--','LineWidth',2);


text(.1,0.4,...
    ['Max LogL of ',num2str(max_llhd),10,'at \alpha/\beta = ',num2str(a2b_range(max_idx))],...
    'FontSize',16,'Units','normalized','BackgroundColor','w');
set(gca,'FontSize',14);
ylabel('Log-Likelihood (Logistic Regression)','FontSize',16);
xlabel('\alpha/\beta value [Gy]','FontSize',16);
title(['Logistic regression',10,'BPx vs. Max BED for given \alpha/\beta'],'FontSize',16);
%grid on;

disp(['Logistic LLHD Physical dose: ',num2str(llhds(end))])

if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_loc,'bpx_a2b_llhds'],'-pdf');
    disp(['Saving ',fig_loc,'bpx_a2b_llhds.pdf...']);
end;

    
%% Logrank p-values
%  cur_fig=figure(3);clf reset;hold on;
% set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
% plot(a2b_range(1:end-1),lr_pvals(1:end-1),'LineWidth',2);% excluding Inf
% [min_pval,min_idx] = min(lr_pvals(1:end-1));
% text(.1,0.8,...
%     ['Min p-value of ',num2str(min_pval),10,'at \alpha/\beta = ',num2str(a2b_range(min_idx))],...
%     'FontSize',16,'Units','normalized','BackgroundColor','w');
% set(gca,'FontSize',14);
% ylabel('Logrank p-values','FontSize',16);
% xlabel('\alpha/\beta value','FontSize',16);
% title(['KM Logrank Analysis',10,'BPx vs. Max BED for given \alpha/\beta'],'FontSize',16);
% grid on;
% 
% if do_print,
%     set(cur_fig,'Color','w');
%     export_fig(cur_fig,[fig_loc,'bpx_a2b_lr_pvals'],'-png');
%     disp(['Saving ',fig_loc,'bpx_a2b_lr_pvals.png...']);
% end;

cur_fig=figure(3);clf reset;hold on;
set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);

x_axis = a2b_range(1:end-1);
[ax,h1,h2]=plotyy(x_axis,lr_pvals(1:end-1),x_axis,log(km_hrs(1:end-1)),@semilogy);hold on;
set(h1,'LineWidth',2);
set(h2,'LineWidth',2);
%set(ax(1),'YLim',[0.001 1]);
set(get(ax(1),'Ylabel'),'String','p-value','FontSize',16);
[min_pval,min_idx] = min(lr_pvals(1:end-1));
 text(.5,0.6,...
    ['Min p-value of ',num2str(min_pval,3),10,'at \alpha/\beta = ',num2str(a2b_range(min_idx))],...
     'FontSize',16,'Units','normalized','BackgroundColor','w');
%set(ax(2),'YLim',[ylim(1) max(log(km_hrs(1:end-1)))*1.2]);
set(get(ax(2),'Ylabel'),'String','ln(Hazard Ratio)','FontSize',16);
%semilogy(x_axis,repmat(0.05,length(x_axis)),'r--','LineWidth',1);
set(ax(1),'XLim',[min(x_axis) max(x_axis)]);
set(ax(2),'XLim',[min(x_axis) max(x_axis)]);
set(ax,'FontSize',14);

hold off; % grid on;
xlabel('\alpha/\beta value [Gy]','fontsize',15);


if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_loc,'bpx_a2b_lr_pvals_km_hr'],'-pdf');
    disp(['Saving ',fig_loc,'bpx_a2b_lr_pvals.pdf...']);
end;


end
