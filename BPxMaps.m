function BPxMaps
tic;
% prepare
fp = 'Z:\elw\MATLAB\bpx_analy\meta_data\';

screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];
ss_sixteen2nine = [0 0 screen_size(3) screen_size(4)];


   
do_print = true;
do_allfx = true;

do_inclfx = true;
do_1fx = false;
do_3fx = false;
do_4fx = false;
do_5fx = false;

fig_loc = 'Z:\elw\MATLAB\bpx_analy\slides\figures\latest\';
 
%fn = {'BPx_DiVj_DVHs_fx-1_a2bInf.mat'};dose_calib='phys';a2b_str='\infty';
fn = {'BPx_DiVj_DVHs_fx-1_a2b3.mat'};dose_calib='a2b3';a2b_str='3';


CGobj = cell(length(fn),1);

%scrsz = get(0,'ScreenSize');
%set(0,'DefaultFigurePosition',[scrsz(1)+scrsz(1)/4 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);

% load data
for m = 1:length(fn)
    load(strcat(fp,fn{m}),'CGobj_current');
    CGobj{m} = CGobj_current;
end



for m = 1:length(fn)
    
    
    %% KM curve displays
    
    % compose complication KM curves
     cm1 = colormap(jet(300)); cm1=cm1(1:256,:); %cm1(end,:) = 0.5;
    cm2 = colormap(jet(10));
    
    
    CG = CGobj{m};
    % survival/complication time
    f2 = ~cellfun('isempty',{CG.mGrp.mDateComp}); % patients with no complication date
    f3 = ~cellfun('isempty',{CG.mGrp.mDateLastFollowup}); % patients with no last follow up date
    compdate = inf(CG.mNumInGrp,1);
    lastfollowup = inf(CG.mNumInGrp,1);
    compdate(f2) = ([CG.mGrp(f2).mDateComp] - [CG.mGrp(f2).mDateBaseline])' / 30;
    lastfollowup(f3) = ([CG.mGrp(f3).mDateLastFollowup] - [CG.mGrp(f3).mDateBaseline])' / 30;
    compdate = min( lastfollowup, compdate );
    flgcensor = [CG.mGrp.mFlgCensor]';
    
     % log10(n) correction for figures
    LymanN = log10(CG.mLymanN);   
    CG.mLymanN = LymanN;
    CG.mBetaCumulativeThreshold = 0.2;
    

    
    
    
    
    
    if do_inclfx || do_allfx
    
    CGall = CG.fCrudeAtlas_DVH(-1);    
    %% DVH atlas

    
     %% Full DVH Atlas
    cur_fig=figure(2);  clf reset; hold on; % grid on;
    %set(gcf,'Position',ss_sixteen2nine);
    set(gcf,'Position',ss_four2three);
    %dose_step = 50;
    dose_step = 25;
    vol_step = 5;% 5% steps
    fontsz = 10;
    CGall.fAtlasFig_DVH(dose_step,vol_step,fontsz);
    
    if do_print,
       set(cur_fig,'Color','w');
       export_fig(cur_fig,[fig_loc,'bpx_dvh_atlas_allfx'],'-pdf');
    end;
    
    
    %% maps
    cur_fig=figure(8); clf reset;
    set(gcf,'Position',ss_sixteen2nine);
    
    CGall = CGall.fBetaCumulativeProbability_DVH();
    CGall.fProbabilityFig_DVH('');
    
    if do_print,
        set(cur_fig,'Color','w');
        export_fig(cur_fig,[fig_loc,'bpx_bcum_map_allfx'],'-pdf');
        disp(['Saving ',fig_loc,'bpx_bcum_map_allfx.pdf...']);
    end
    
    
     %% maps
    cur_fig=figure(13); clf reset;
    set(gcf,'Position',ss_sixteen2nine);
    
    CGall = CGall.fBetaInverseProbability_DVH();
    CGall.fLow68pctConfidenceFig_DVH();
    
    if do_print,
        set(cur_fig,'Color','w');
        export_fig(cur_fig,[fig_loc,'bpx_binv_map_allfx'],'-pdf');
        disp(['Saving ',fig_loc,'bpx_binv_map_allfx.pdf...']);
    end
    
    end
    
    
    
   if do_1fx || do_allfx
   
    CGone = CG.fCrudeAtlas_DVH(1);
   
     %% 1fx DVH Atlas
    cur_fig=figure(3);  clf reset; hold on; % grid on;
    set(gcf,'Position',ss_sixteen2nine);
    dose_step = 50;
    vol_step = 5;% 5% steps
    fontsz = 14;
    CGone.fAtlasFig_DVH(dose_step,vol_step,fontsz);
    
    if do_print,
       set(cur_fig,'Color','w');
       export_fig(cur_fig,[fig_loc,'bpx_dvh_atlas_1fx'],'-pdf');
    end;
    
        %% maps
    cur_fig=figure(9); clf reset;
    set(gcf,'Position',ss_sixteen2nine);
    
    CGone = CGone.fBetaCumulativeProbability_DVH();
    CGone.fProbabilityFig_DVH('');
    
    if do_print,
        set(cur_fig,'Color','w');
        export_fig(cur_fig,[fig_loc,'bpx_bcum_map_1fx'],'-pdf');
        disp(['Saving ',fig_loc,'bpx_bcum_map_1fx.pdf...']);
    end
    
    %% maps
    cur_fig=figure(14); clf reset;
    set(gcf,'Position',ss_sixteen2nine);
    
    CGone = CGone.fBetaInverseProbability_DVH();
    CGone.fLow68pctConfidenceFig_DVH();
    
    if do_print,
        set(cur_fig,'Color','w');
        export_fig(cur_fig,[fig_loc,'bpx_binv_map_1fx'],'-pdf');
        disp(['Saving ',fig_loc,'bpx_binv_map_1fx.pdf...']);
    end
   end
   
   if do_3fx || do_allfx
        
   CGthree = CG.fCrudeAtlas_DVH(3);
       
    %% DVH atlas nfx==3  
    cur_fig=figure(5);  clf reset; hold on; % grid on;
    set(gcf,'Position',ss_sixteen2nine);
    

    
    dose_step = 50;
    vol_step = 5;% 5% steps
    fontsz = 14;
    
    CGthree.fAtlasFig_DVH(dose_step,vol_step,fontsz);
    
    if do_print,
       set(cur_fig,'Color','w');
       export_fig(cur_fig,[fig_loc,'bpx_dvh_atlas_3fx'],'-pdf');
    end; 
         
    %% maps
    cur_fig=figure(10); clf reset;
    set(gcf,'Position',ss_sixteen2nine);
    
    CGthree = CGthree.fBetaCumulativeProbability_DVH();
    CGthree.fProbabilityFig_DVH('');
    
    if do_print,
        set(cur_fig,'Color','w');
        export_fig(cur_fig,[fig_loc,'bpx_bcum_map_3fx'],'-pdf');
        disp(['Saving ',fig_loc,'bpx_bcum_map_3fx.pdf...']);
    end
    
     %% maps
    cur_fig=figure(15); clf reset;
    set(gcf,'Position',ss_sixteen2nine);
    
    CGthree = CGthree.fBetaInverseProbability_DVH();
    CGthree.fLow68pctConfidenceFig_DVH();
    
    if do_print,
        set(cur_fig,'Color','w');
        export_fig(cur_fig,[fig_loc,'bpx_binv_map_3fx'],'-pdf');
        disp(['Saving ',fig_loc,'bpx_binv_map_3fx.pdf...']);
    end
    
   end;
   
   if do_4fx || do_allfx
   
       CGfour = CG.fCrudeAtlas_DVH(4);
       
       %% DVH atlas nfx==3 4
    cur_fig=figure(6);  clf reset; hold on; % grid on;
    set(gcf,'Position',ss_sixteen2nine);
    
   

    dose_step = 50;
    vol_step = 5;% 5% steps
    fontsz = 14;
    
    CGfour.fAtlasFig_DVH(dose_step,vol_step,fontsz);

    if do_print,
       set(cur_fig,'Color','w');
       export_fig(cur_fig,[fig_loc,'bpx_dvh_atlas_4fx'],'-pdf');
    end; 
    
        %% maps
    cur_fig=figure(11); clf reset;
    set(gcf,'Position',ss_sixteen2nine);
    
    CGfour = CGfour.fBetaCumulativeProbability_DVH();
    CGfour.fProbabilityFig_DVH('');
    
    if do_print,
        set(cur_fig,'Color','w');
        export_fig(cur_fig,[fig_loc,'bpx_bcum_map_4fx'],'-pdf');
        disp(['Saving ',fig_loc,'bpx_bcum_map_4fx.pdf...']);
    end
    
        cur_fig=figure(16); clf reset;
    set(gcf,'Position',ss_sixteen2nine);
    
    CGfour = CGfour.fBetaInverseProbability_DVH();
    CGfour.fLow68pctConfidenceFig_DVH();
    
    if do_print,
        set(cur_fig,'Color','w');
        export_fig(cur_fig,[fig_loc,'bpx_binv_map_4fx'],'-pdf');
        disp(['Saving ',fig_loc,'bpx_binv_map_4fx.pdf...']);
    end
    
   end;
   
   if do_5fx || do_allfx
       
       CGfive = CG.fCrudeAtlas_DVH(5);
       
       %% DVH atlas nfx==3 4
    cur_fig=figure(7);  clf reset; hold on; % grid on;
    set(gcf,'Position',ss_sixteen2nine);
    
    dose_step = 50;
    vol_step = 5;% 5% steps
    fontsz = 14;
    
    CGfive.fAtlasFig_DVH(dose_step,vol_step,fontsz);
    
    if do_print,
       set(cur_fig,'Color','w');
       export_fig(cur_fig,[fig_loc,'bpx_dvh_atlas_5fx'],'-pdf');
    end; 
        
    %% maps
    cur_fig=figure(12); clf reset;
    set(gcf,'Position',ss_sixteen2nine);
    
    CGfive = CGfive.fBetaCumulativeProbability_DVH();
    CGfive.fProbabilityFig_DVH('');
    
    if do_print,
        set(cur_fig,'Color','w');
        export_fig(cur_fig,[fig_loc,'bpx_bcum_map_5fx'],'-pdf');
        disp(['Saving ',fig_loc,'bpx_bcum_map_5fx.pdf...']);
    end
    
    
    cur_fig=figure(17); clf reset;
    set(gcf,'Position',ss_sixteen2nine);
    
    CGfive = CGfive.fBetaInverseProbability_DVH();
    CGfive.fLow68pctConfidenceFig_DVH();
    
    if do_print,
        set(cur_fig,'Color','w');
        export_fig(cur_fig,[fig_loc,'bpx_binv_map_5fx'],'-pdf');
        disp(['Saving ',fig_loc,'bpx_binv_map_5fx.pdf...']);
    end
    
   end;
    
end
end

    