function BPx_Anonymize

% prepare
fp = 'Z:\elw\MATLAB\bpx_analy\meta_data\';

%fn = {'BPx_DiVj_DVHs_fx-1_a2bInf.mat'};dose_calib='phys';a2b_str='\infty';
fn = {'BPx_DiVj_DVHs_fx-1_a2b3.mat'};


CGobj = cell(length(fn),1);
% load data
for m = 1:length(fn)
    load(strcat(fp,fn{m}),'CGobj_current');
   
    
     CGgrp = CGobj_current.mGrp;

    % encrypt mID 
   anon_mrns = cell(length(CGgrp),2);
   for j=1:length(CGgrp)
       
       
       tmp_id = CGgrp(j).mID;
        tmp_md5 = DataHash(tmp_id);
        anon_mrns{j,1} = tmp_id;
        anon_mrns{j,2} = tmp_md5;
       
        CGgrp(j).mID = tmp_md5;
    
   end
   
   [~,md5_inds] = sort({anon_mrns{:,2}});
   new_ids = mat2cell(md5_inds,1,ones(1,size(md5_inds,2)));
   anon_mrns(:,2) = new_ids';
   
   for k=1:length(CGgrp)
       CGgrp(k).mID = num2str(anon_mrns{k,2});
   end
   
   CGobj_current.mGrp = CGgrp;
   xlswrite(['Z:/elw/MATLAB/bpx_analy/meta_data/anon_data/BPx_mrn_encrypt.xlsx'],anon_mrns)
   anon_fn=['Z:/elw/MATLAB/bpx_analy/meta_data/anon_data/BPx_DiVj_DVHs_fx-1_a2b3_anon.mat'];
 
   
   save(anon_fn,'CGobj_current');

end
end
