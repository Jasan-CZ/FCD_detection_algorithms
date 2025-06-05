function []=extension_ver20241030
% Detection Algorithm FCD, m-code
% EXTENSION images, MAC ver.2024-10-30
%
% INPUTS [2]:
% 1. T1w volume in "T1*.nii" or "T1*.nii.gz" format
% 2. wCORTEX_BinarMask.nii [made by CortThk2BinarMask.m]
% 
% A. Simple subtraction of data with templates
% B. Statistics with ADULT-template [z-score]
% "template_81zk_extension_mean" & "template_81zk_extension_std"
% C. Statistics with PEDIATRIC-template [z-score]
% "template_37hr_extension_mean" & "template_37hr_extension_std"
% D. Transformation maps back to structural space
%
% *************************************************************************
% 0. Definition of input constants and paths
% 1. Open T1w MR, automatic / selection by file
% 2. Segmentation
% 3. Normalisation
% 4. Smoothing 6mm + masking
% 5. Statistics [mean, std] (z-score maps) / Subtraction
%    - ADULT
%    - PEDIATRIC
% 6. Transformation to the structural space
% 7. Saving variables, Final clearing of data
% *************************************************************************

%% 0. Definition of input constants and paths
clearvars; clc;

% masks (BasalGanglias, Cerebellum, NaturallyBlurredStructures, Midline = All)
MaskAll='/Users/jasan/Documents/matlab/toolbox/fcd/fcd_masks/mask_mni_binar_All.nii';

% path to templates extension mean & std - ADULT / PEDIATRIC
path_to_adult_templates='/Users/jasan/Documents/matlab/toolbox/fcd/fcd_templates/adult/extension/';
path_to_pediatric_templates='/Users/jasan/Documents/matlab/toolbox/fcd/fcd_templates/pediatric/extension/';

% gauss kernel for smoothing
GaussKernel = 6;

% Subtraction & Z-score thresholds; cluster volume threshold
SubtractionIntensityThreshold = 0.45;
ZscoreIntensityThreshold = [2, 5]; %0=(0.1*MAX)
Threshold_for_small_clusters = 60; % prah malych clusteru v mm^3

%% 1. Open T1w MR volume
% automatical selection of T1-file
T1file = dir('T1*.nii');
T1 = [T1file.folder,'/',T1file.name];
[T1path, T1filename, T1ext] = fileparts(T1);
% selection by file
% T1 = spm_select(1,'.*','T1 MR v NII/Analyze formátu');
% if isempty(T1)
%     return
% end
% % pokud T1 ve formatu *.NII.GZ
% if contains(T1,'.nii.gz')
%     gunzip(T1);
%     T1=strrep(T1,'.nii.gz','.nii,1');
% end
% ziskani cesty, jmena souboru a pripony T1 ['cesta', 'jmeno', 'pripona,1']
% [T1path, T1filename, T1ext] = fileparts(T1);
% cd(T1path);

%% 2. Segmentation
clc;
disp('S T A R T   O F   A L G O R I T H M');
disp(T1file.name);
disp(' ... Segmentation');
if exist('matlabbatch','var'); clear matlabbatch; end
spm('defaults','fmri');
spm_jobman('initcfg');

matlabbatch{1}.spm.spatial.preproc.channel.vols = {T1};
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {'/Users/jasan/Documents/matlab/toolbox/spm12/tpm/TPM.nii,1'};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {'/Users/jasan/Documents/matlab/toolbox/spm12/tpm/TPM.nii,2'};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {'/Users/jasan/Documents/matlab/toolbox/spm12/tpm/TPM.nii,3'};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {'/Users/jasan/Documents/matlab/toolbox/spm12/tpm/TPM.nii,4'};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {'/Users/jasan/Documents/matlab/toolbox/spm12/tpm/TPM.nii,5'};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {'/Users/jasan/Documents/matlab/toolbox/spm12/tpm/TPM.nii,6'};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];
spm_jobman('run',matlabbatch);

%% 3. Normalisation
disp(' ... Normalisation');
if exist('matlabbatch','var'); clear matlabbatch; end
spm('defaults','fmri');
spm_jobman('initcfg');
matlabbatch{1}.spm.spatial.normalise.write.subj.def = {[T1path,'/y_',T1filename,'.nii']};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {[T1path,'/c1',T1filename,T1ext];[T1path,'/m',T1filename,T1ext]};
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                          78 76 125];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 1;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
spm_jobman('run',matlabbatch);

%% 4. Smoothing Gauss kernel 6mm + masking
disp(' ... Smoothing Gauss kernel 6mm');
if exist('matlabbatch','var'); clear matlabbatch; end
spm('defaults','fmri'); spm_jobman('initcfg');
matlabbatch{1}.spm.spatial.smooth.data = {[T1path,'/wc1',T1filename, T1ext,',1']};
matlabbatch{1}.spm.spatial.smooth.fwhm = [GaussKernel GaussKernel GaussKernel];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';
spm_jobman('run',matlabbatch);

%% 5. Statistics [mean, std] (z-score maps) / Subtraction
disp(' ... Subtraction, Statistics [z-score, ADULT and PEDIATRIC template], Masking');
if exist('matlabbatch','var'); clear matlabbatch; end
spm('defaults','fmri'); spm_jobman('initcfg');
% ADULT:
% ADULT subtraction
matlabbatch{1}.spm.util.imcalc.input = {
                                    [T1path,'/swc1',T1filename,T1ext,',1']
                                    [path_to_adult_templates,'template_81zk_extension_mean.nii,1']
                                    [path_to_adult_templates,'template_81zk_extension_std.nii,1']
                                    [MaskAll,',1']
                                    [T1path,'/wCORTEX_BinarMask.nii,1']
                                    };
matlabbatch{1}.spm.util.imcalc.outdir = {''};
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

matlabbatch{1}.spm.util.imcalc.output = 'wEXTENSION_subtraction_AdultTemplateCortThkMask';
matlabbatch{1}.spm.util.imcalc.expression = '(((i1-i2).*(((i1>0.5)+(i2>0.25))>0))-((i1-i2).*(((i1>0.5)+(i2>0.25))>0).*(i4>0))).*(i5>0)';
spm_jobman('run',matlabbatch);
% ADULT Z-score
matlabbatch{1}.spm.util.imcalc.output = 'wEXTENSION_z_score_AdultTemplateCortThkMask';
matlabbatch{1}.spm.util.imcalc.expression = '(((i1-i2)./i3).*(((i1>0.5)+(i2>0.25))>0)-(((i1-i2)./i3).*(((i1>0.5)+(i2>0.25))>0).*(i4>0))).*(i5>0)';
spm_jobman('run',matlabbatch);
% *************************************************************************
% PEDIATRIC:
% PEDIATRIC subtraction
matlabbatch{1}.spm.util.imcalc.input = {
                                    [T1path,'/swc1',T1filename,T1ext,',1']
                                    [path_to_pediatric_templates,'template_37hr_extension_mean.nii,1']
                                    [path_to_pediatric_templates,'template_37hr_extension_std.nii,1']
                                    [MaskAll,',1']
                                    [T1path,'/wCORTEX_BinarMask.nii,1']
                                    };
matlabbatch{1}.spm.util.imcalc.output = 'wEXTENSION_subtraction_PediatricTemplateCortThkMask';
matlabbatch{1}.spm.util.imcalc.expression = '(((i1-i2).*(((i1>0.5)+(i2>0.25))>0))-((i1-i2).*(((i1>0.5)+(i2>0.25))>0).*(i4>0))).*(i5>0)';
spm_jobman('run',matlabbatch);
% PEDIATRIC Z-score
matlabbatch{1}.spm.util.imcalc.output = 'wEXTENSION_z_score_PediatricTemplateCortThkMask';
matlabbatch{1}.spm.util.imcalc.expression = '(((i1-i2)./i3).*(((i1>0.5)+(i2>0.25))>0)-(((i1-i2)./i3).*(((i1>0.5)+(i2>0.25))>0).*(i4>0))).*(i5>0)';
spm_jobman('run',matlabbatch);

%% 6. Transformation to the structural space
disp(' ... Transformace do strukturalniho prostoru');
if exist('matlabbatch','var'); clear matlabbatch; end
spm('defaults','fmri');
spm_jobman('initcfg');
matlabbatch{1}.spm.spatial.normalise.write.subj.def = {[T1path,'/iy_',T1filename,'.nii']};
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -100 -70
                                                          78 110 125];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 1; % byla 4
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'str';

% ADULT subraction
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {[T1path,'/wEXTENSION_subtraction_AdultTemplateCortThkMask',T1ext];[T1path,'/wEXTENSION_subtraction_AdultTemplateCortThkMask',T1ext]};
spm_jobman('run',matlabbatch);
copyfile('strwEXTENSION_subtraction_AdultTemplateCortThkMask.nii','strEXTENSION_subtraction_AdultTemplateCortThkMask.nii');
delete('strwEXTENSION_subtraction_AdultTemplateCortThkMask.nii');
% ADULTsubtractionCLUSTERS = find_thresholded_clusters([T1path,'/strEXTENSION_subtraction_AdultTemplateCortThkMask',T1ext],SubtractionIntensityThreshold,Threshold_for_small_clusters);

% ADULT z-score
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {[T1path,'/wEXTENSION_z_score_AdultTemplateCortThkMask',T1ext];[T1path,'/wEXTENSION_z_score_AdultTemplateCortThkMask',T1ext]};
spm_jobman('run',matlabbatch);
copyfile(['strwEXTENSION_z_score_AdultTemplateCortThkMask',T1ext],'strEXTENSION_z_score_AdultTemplateCortThkMask.nii');
delete(['strwEXTENSION_z_score_AdultTemplateCortThkMask',T1ext]);
% ADULTz2scoreCLUSTERS = find_thresholded_clusters([T1path,'/strEXTENSION_z_score_AdultTemplateCortThkMask',T1ext],ZscoreIntensityThreshold(1),Threshold_for_small_clusters);
% ADULTz5scoreCLUSTERS = find_thresholded_clusters([T1path,'/strEXTENSION_z_score_AdultTemplateCortThkMask',T1ext],ZscoreIntensityThreshold(2),Threshold_for_small_clusters);

% PEDIATRIC subraction
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {[T1path,'/wEXTENSION_subtraction_PediatricTemplateCortThkMask',T1ext];[T1path,'/wEXTENSION_subtraction_PediatricTemplateCortThkMask',T1ext]};
spm_jobman('run',matlabbatch);
copyfile(['strwEXTENSION_subtraction_PediatricTemplateCortThkMask',T1ext],'strEXTENSION_subtraction_PediatricTemplateCortThkMask.nii');
delete(['strwEXTENSION_subtraction_PediatricTemplateCortThkMask',T1ext]);
% PEDIATRICsubtractionCLUSTERS = find_thresholded_clusters([T1path,'/strEXTENSION_subtraction_PediatricTemplateCortThkMask',T1ext],SubtractionIntensityThreshold,Threshold_for_small_clusters);

% PEDIATRIC z-score
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {[T1path,'/wEXTENSION_z_score_PediatricTemplateCortThkMask',T1ext];[T1path,'/wEXTENSION_z_score_PediatricTemplateCortThkMask',T1ext]};
spm_jobman('run',matlabbatch);
copyfile(['strwEXTENSION_z_score_PediatricTemplateCortThkMask',T1ext],'strEXTENSION_z_score_PediatricTemplateCortThkMask.nii');
delete(['strwEXTENSION_z_score_PediatricTemplateCortThkMask',T1ext]);
% PEDIATRICz2scoreCLUSTERS = find_thresholded_clusters([T1path,'/strEXTENSION_z_score_PediatricTemplateCortThkMask',T1ext],ZscoreIntensityThreshold(1),Threshold_for_small_clusters);
% PEDIATRICz5scoreCLUSTERS = find_thresholded_clusters([T1path,'/strEXTENSION_z_score_PediatricTemplateCortThkMask',T1ext],ZscoreIntensityThreshold(2),Threshold_for_small_clusters);

%% 7. Saving variables, Final clearing of data
clc;
save ('extension.mat');
% disp('Final clearing of data');
% delete('*.mat');
% delete('*T1*.nii');
% delete('wEXTENSION*.nii');
disp('E N D   O F   A L G O R I T H M');



%% functions

function [Clusters]=find_thresholded_clusters(PathToNIIFile, IntensityThreshold, VolumeThreshold) %napr. find_thresholded_clusters('path/to/file.nii')
CLUSTERvol=spm_vol([PathToNIIFile,',1']);
CLUSTERimg=spm_read_vols(CLUSTERvol);
CLUSTERimg_thresholded = CLUSTERimg > IntensityThreshold;
ConnectedComponents = bwconncomp(bwareaopen(CLUSTERimg_thresholded,VolumeThreshold));
LabelsMatrix = labelmatrix(ConnectedComponents);
Clusters.Count=num2str(max(LabelsMatrix,[],'all'));
Clusters.Volume=num2str(numel(find(LabelsMatrix>0)));
for i=1:max(LabelsMatrix,[],'all')
    disp([sprintf('%2.0d.',i),' ... ',sprintf('%4.0d voxel(s)',(numel(find(LabelsMatrix==i))))]);
    Clusters.Cluster(i)=(numel(find(LabelsMatrix==i)));
end
disp(['Find ',num2str(max(LabelsMatrix,[],'all')),' clusters with volume ',num2str(numel(find(LabelsMatrix>0))),' voxels/mm^3:']);