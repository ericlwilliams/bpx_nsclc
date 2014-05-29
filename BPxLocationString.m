function locStr=BPxLocationString(pt_name,tar_loc)

locStr='';
%# Patients with non-standard DVH names
%anom_pts = {'Cohen','Lapid','Menaker','Pisano','Sutton','Underkoffler'};
anom_pts = {''};
matches = strcmp(anom_pts,pt_name);

if sum(matches)>0,
    if tar_loc(2)=='U',
        locStr='CW2SUP';
    else
        locStr='CW2INF';
    end
else %# not anomoly, just determine L/R
     if tar_loc(1)=='R',
        locStr='SD_BPX_R';
     elseif tar_loc(1)=='L',
        locStr='SD_BPX_L';
     elseif tar_loc(1)=='U',
         locStr='SD_BPX_U';
     else
         disp(['!! Unknown identifier: ',tar_loc,' for patient ',pt_name,10]);
     end
end
    
end
