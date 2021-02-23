% Produce individual contrast images containing the experimental effects of
% interest on a 3D triangular mesh (canonical scalp surface). 
%
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

%
% $Id$

spm('Defaults', 'EEG')
[fname, sts] = spm_select(Inf, '^con*.*.\.mat$', 'Select contrast beta files');
nsubj = size(fname, 1);

gpos = []; r = 30;
for i = 1:nsubj
    load(deblank(fname(i,:)));
    
    if isfield(S.ch, 'xyzH')
        pos = S.ch.xyzH';
    else
        pos = S.ch.xyzC';
    end
    
    vinterp = spm_fnirs_roivert(pos, r);

    mdir = fileparts(which('spm_fnirs'));
    load(fullfile(mdir, 'gen_conimg', 'mask_iskull.mat'));
    
    vinterp = intersect(vinterp, vmask);
    
    if i == 1
        gpos = vinterp;
    else
        gpos = intersect(gpos, vinterp);
    end
end

for i = 1:nsubj
    fname_c = deblank(fname(i, :)); 
    load(fname_c); 
    [G, interpY] = spm_fnirs_3dscalp(S, fname_c, gpos);
end



