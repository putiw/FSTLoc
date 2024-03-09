function roi = get_mt_range(subject,serverDir,hemi)


lcurv = read_curv(fullfile(serverDir,'/derivatives/freesurfer', subject,'surf', 'lh.curv'));

 mtrange = read_ROIlabel(fullfile(serverDir,'/derivatives/freesurfer', subject,'label/0localizer/',[hemi,'h.range.label']));
 
 if contains(hemi,'l')
     roi = mtrange;
 else
     roi = numel(lcurv)+mtrange;
 end
    

end