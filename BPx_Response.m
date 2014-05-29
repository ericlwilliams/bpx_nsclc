function BPx_Response
tic;
% prepare
if ~isunix,
    analy_loc = 'Z:/elw/MATLAB/bpx_analy/';
else
    analy_loc = '/Users/elw/Dropbox/eMacs/mskcc/analy/bpx_analy/';
end

fp = strcat(analy_loc,'meta_data/');
fig_loc = strcat(analy_loc,'slides/figures/latest/');

fig_ctr = 1;
screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

   
do_print = true;

%fn = 'BPx_DiVj_DVHs_fx-1_a2bInf.mat';dose_calib='phys';a2b_str='\infty';
fn = 'BPx_DiVj_DVHs_fx-1_a2b3.mat';dose_calib='a2b3';a2b_str='3';

if isequal(a2b_str,'\infty')
    dose_unit = 'Gy';
else
    dose_unit = ['Gy$_{',a2b_str,'}$'];
end

CGobj = cell(length(fn),1);
screen_size=get(0,'ScreenSize');

% load data
load(strcat(fp,fn),'CGobj_current');


% get preditors
doses = [CGobj_current.mGrp.mDoseBins_LQ];
vols = [CGobj_current.mGrp.mVolCum];
d05s = zeros(size(vols,2),1);

dmax = max(doses)';
%dmax2 = zeros(size(vols,2),1);

for j=1:size(vols,2) % for each patient
    vol = vols(:,j);
    
    vol = vol./max(vol);
    v05 = vol<0.05;
    d05_inds = find(v05);
    dose = doses(:,j);
    d05 = min(dose(d05_inds));
    d05s(j) = d05;
    
%     %tmp
%     vol(~vol)=nan;
%     nan_inds = find(isnan(vol));
%     if ~isempty(nan_inds)
%         min_ind = nan_inds(1)-1;
%     else
%         min_ind = length(vol);
%     end
%     %[~,min_ind] = min(vol);
%     dmax2(j) = dose(min_ind);
%
end

pttotal = ones(CGobj_current.mNumInGrp,1);
ptcomp = ones(CGobj_current.mNumInGrp,1); 
ptcomp([CGobj_current.mGrp.mFlgCensor])=0;


%% D_Max Response

cur_fig=figure(fig_ctr);
fig_ctr=fig_ctr+1;
set(gcf,'Position',ss_four2three);

% regression 
[b,dev,s]=glmfit(dmax,[ptcomp pttotal],'binomial','link','logit');

disp('Dmax model coefficients');
disp(['beta: ',num2str(b(2))]);
disp(['68% CI: ',num2str(b(2) + 0.9945*[-s.se(2) s.se(2)])]);

pr = exp(b(1)+b(2)*dmax);
pr = pr./(1+pr); % logistic probability
pr(~logical(ptcomp)) = 1-pr(~logical(ptcomp)); % non-complication patients
pr = log(pr); % log likelihood of each patients
loglikelihood = sum(pr); % loglikelihood of al
disp(['SE: ',num2str(s.se(2)),' p: ',num2str(s.p(2)),' Logl: ',num2str(loglikelihood)]);

pvalue = [s.p];
pval = pvalue(2); % the p-value corresponding to gEUD

%tmp
doses = (0:max(dmax))'; % doses (gEUD) whose RP probability will be computed
%doses = (0:90)';
[rpb,rplo,rphi] = glmval(b, doses,'logit',s); % the responding function values at doses

% plot
hold on;
plot(doses,rpb,'k','LineWidth',2); % responding function
plot(doses,rpb-rplo,'k','LineWidth',1); % low CI curve
plot(doses,rpb+rphi,'k','LineWidth',1); % high CI curve


flg=[CGobj_current.mGrp.mFlgCensor]; % censor flags of patients

[medianeud,~,~,binlow,binhigh,numcomp,numtotal,betainv84,betainv16] = EventObserved(flg,dmax,4);
prob = numcomp./numtotal;
% plot
errorbar(medianeud,prob,max(0,prob-betainv16),max(0,betainv84-prob),'k*','LineWidth',1);
errorbar_x(medianeud,prob,(medianeud-binlow),(binhigh-medianeud),'k*');

%tmp
ylim([0 0.5]);
xlim([0 max(binhigh)*1.05]);

%xlims = xlim;
%xmax = xlims(2);
%xlim([0 xmax]);
%ylim([0 0.5]);

loga_str = ['$p$-value~$= ',num2str(pval,'%3.2e'),'$'];
%str = sprintf(loga_str,num2str(loga),num2str(pval,2),structures{j});
% str = sprintf(loga_str,num2str(pval,'%3.2e'));
%text(0.65,0.85,str,'FontSize',24,'Units','normalized');
textbp(loga_str,'FontSize',24,'interpreter','latex');

set(gca,'xminortick','on','yminortick','on');
set(gca,'box','on');
set(gca,'FontSize',20);
xlabel(['D$_{\rm{max}}$~[',dose_unit,']'],'interpreter','latex','FontSize',24);
ylabel('Probability of Brachial Plexopathy','interpreter','latex','FontSize',24);

if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,...
        [fig_loc,'bpx_',dose_calib,'_lr_response_dmax'],'-pdf');
    disp(['Saving ',fig_loc,'bpx_',dose_calib,'_lr_response_dmax']);
end

%% D05 Response

cur_fig=figure(fig_ctr);
fig_ctr=fig_ctr+1;
set(gcf,'Position',ss_four2three);

% regression 
[b,~,s]=glmfit(d05s,[ptcomp pttotal],'binomial','link','logit');

disp('D_{05} model coefficients');
disp(['beta: ',num2str(b(2))]);
disp(['68% CI: ',num2str(b(2) + 0.9945*[-s.se(2) s.se(2)])]);

pr = exp(b(1)+b(2)*d05s);
pr = pr./(1+pr); % logistic probability
pr(~logical(ptcomp)) = 1-pr(~logical(ptcomp)); % non-complication patients
pr = log(pr); % log likelihood of each patients
loglikelihood = sum(pr); % loglikelihood of al
disp(['SE: ',num2str(s.se(2)),' p: ',num2str(s.p(2)),' Logl: ',num2str(loglikelihood)]);

pvalue = [s.p];
pval = pvalue(2); % the p-value corresponding to gEUD

%tmp
doses = (0:max(d05s))'; % doses (gEUD) whose RP probability will be computed
%doses = (0:90)';
[rpb,rplo,rphi] = glmval(b, doses,'logit',s); % the responding function values at doses

% plot
hold on;
plot(doses,rpb,'k','LineWidth',2); % responding function
plot(doses,rpb-rplo,'k','LineWidth',1); % low CI curve
plot(doses,rpb+rphi,'k','LineWidth',1); % high CI curve


flg=[CGobj_current.mGrp.mFlgCensor]; % censor flags of patients

[medianeud,~,~,binlow,binhigh,numcomp,numtotal,betainv84,betainv16] = EventObserved(flg,d05s,4);
prob = numcomp./numtotal;
% plot
errorbar(medianeud,prob,max(0,prob-betainv16),max(0,betainv84-prob),'k*','LineWidth',1);
errorbar_x(medianeud,prob,(medianeud-binlow),(binhigh-medianeud),'k*');

%tmp
ylim([0 0.5]);
xlim([0 max(binhigh)*1.05]);

%xlims = xlim;
%xmax = xlims(2);
%xlim([0 xmax]);
%ylim([0 0.5]);

loga_str = ['$p$-value~$= ',num2str(pval,'%3.2e'),'$'];
%str = sprintf(loga_str,num2str(loga),num2str(pval,2),structures{j});
% str = sprintf(loga_str,num2str(pval,'%3.2e'));
%text(0.65,0.85,str,'FontSize',24,'Units','normalized');
textbp(loga_str,'FontSize',24,'interpreter','latex');

set(gca,'xminortick','on','yminortick','on');
set(gca,'box','on');
set(gca,'FontSize',20);
xlabel(['D$_{5\%}$~[',dose_unit,']'],'interpreter','latex','FontSize',24);
ylabel('Probability of Brachial Plexopathy','interpreter','latex','FontSize',24);

if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,...
        [fig_loc,'bpx_',dose_calib,'_lr_response_d05s'],'-pdf');
    disp(['Saving ',fig_loc,'bpx_',dose_calib,'_lr_response_d05s']);
end

end