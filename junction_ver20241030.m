function []=junction_ver20241030
% Detection Algorithm FCD, m-code
% JUNCTION images (BLURRING), MAC ver.2024-10-30
% 
% INPUTS [2]:
% 1. T1w volume in "T1*.nii" or "T1*.nii.gz" format
% 2. wCORTEX_BinarMask.nii [made by CortThk2BinarMask.m]
%
% A. Simple subtraction of data with templates
% B. Statistics with ADULT-template [z-score]
% "template_81zk_junction_blurred5mm_mean" & "template_81zk_junction_blurred5mm_std"
% C. Statistics with PEDIATRIC-template [z-score]
% "template_37hr_junction_blurred5mm_mean" & "template_37hr_junction_blurred5mm_std"
% D. Transformation maps back to structural space
%
% *************************************************************************
%  0. Definition of input constants and paths
%  1. Open T1w MR, automatic / selection by file
%  2. Segmentation
%  3. Normalisation
%  4. Separation of GrayMatter
%  5. Separation of WhiteMatter
%  6. Thresholding and statistics
%  7. Junction (Blurring) calculating
%  8. Statistics [mean, std]
%  9. Transformation to the structural space
% 10. Saving variables, Final clearing of data
% *************************************************************************

%% 0. Definition of input constants and paths
clearvars;
lowerthreshold = 0.5;
upperthreshold = 0.5;

% masks (BasalGanglias, Cerebellum, NaturallyBlurredStructures, Midline = All)
MaskAll='/Users/jasan/Documents/matlab/toolbox/fcd/fcd_masks/mask_mni_binar_All.nii';

% path to templates extension mean & std - ADULT / PEDIATRIC
path_to_adult_templates='/Users/jasan/Documents/matlab/toolbox/fcd/fcd_templates/adult/junction/';
path_to_pediatric_templates='/Users/jasan/Documents/matlab/toolbox/fcd/fcd_templates/pediatric/junction/';

% gauss kernel for smoothing = 5mm
GaussKernel = 5;

%% 1. Open T1w MR, automatic / selection by file
% automatical selection of T1-file
T1file = dir('T1*.nii');
T1 = [T1file.folder,'/',T1file.name];
[T1path, T1filename, T1ext] = fileparts(T1);
% selection by file
%T1 = spm_select(1,'^T1fslroi.nii.gz','T1 MR v NII/Analyze formátu');
% T1 = spm_select(1,'.*','T1 MR v NII/Analyze formátu');
% if isempty(T1)
%     return
% end
% % pokud T1 ve formatu *.NII.GZ
% if contains(T1,'.nii.gz')
%     gunzip(T1);
%     T1=strrep(T1,'.nii.gz','.nii,1');
% end
% T1hdr = spm_vol(T1);
% % ziskani cesty, jmena souboru a pripony T1 ['cesta', 'jmeno', 'pripona,1']
% [T1path, T1filename, T1ext] = fileparts(T1);
% cd(T1path);
% T1img = spm_read_vols(T1hdr);

%% 2. Segmentation
clc;
disp('S T A R T   O F   A L G O R I T H M');
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
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
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
clc;
disp(' ... Normalisation');
if exist('matlabbatch','var'); clear matlabbatch; end
spm('defaults','fmri');
spm_jobman('initcfg');
matlabbatch{1}.spm.spatial.normalise.write.subj.def = {[T1path,'/y_',T1filename,'.nii']};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {[T1path,'/c1',T1filename,T1ext];[T1path,'/c2',T1filename,T1ext];[T1path,'/m',T1filename,T1ext]};
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                          78 76 125];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
spm_jobman('run',matlabbatch);

%% 4. Separation of GrayMatter
clc;
disp(' ... Separation of GrayMatter');
if exist('matlabbatch','var'); clear matlabbatch; end
spm('defaults','fmri'); spm_jobman('initcfg');
matlabbatch{1}.spm.util.imcalc.input = {[T1path,'/wc1',T1filename,T1ext];[T1path,'/wm',T1filename,T1ext]};
matlabbatch{1}.spm.util.imcalc.output = 'GrayMatter';
matlabbatch{1}.spm.util.imcalc.outdir = {''};
matlabbatch{1}.spm.util.imcalc.expression = '(i1>0.5).*i2';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
spm_jobman('run',matlabbatch);
% nacteni GrayMatter
GrayMatterhdr = spm_vol([T1path,'/GrayMatter',T1ext]);
GrayMatterimg = spm_read_vols(GrayMatterhdr);

%% 5. Separation of WhiteMatter
clc;
disp(' ... Separation of WhiteMatter');
if exist('matlabbatch','var'); clear matlabbatch; end
spm('defaults','fmri');
spm_jobman('initcfg');
matlabbatch{1}.spm.util.imcalc.input = {[T1path,'/wc2',T1filename,T1ext];[T1path,'/wm',T1filename,T1ext]};
matlabbatch{1}.spm.util.imcalc.output = 'WhiteMatter';
matlabbatch{1}.spm.util.imcalc.outdir = {''};
matlabbatch{1}.spm.util.imcalc.expression = '(i1>0.5).*i2';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
spm_jobman('run',matlabbatch);
% nacteni WhiteMatterr
WhiteMatterhdr = spm_vol([T1path,'\WhiteMatter',T1ext]);
WhiteMatterimg = spm_read_vols(WhiteMatterhdr);

%% 6. Thresholding and statistics
[GMmean, GMdev] = img_mean_dev(GrayMatterimg);
[WMmean, WMdev] = img_mean_dev(WhiteMatterimg);
disp(['Lower threshold: ',num2str(lowerthreshold)]);
disp(['Upper threshold: ',num2str(upperthreshold)]);
Tlowerthreshold = GMmean+(lowerthreshold*GMdev);
Tupperthreshold = WMmean-(upperthreshold*WMdev);
disp([num2str(Tupperthreshold),' > GM > ',num2str(Tlowerthreshold)]);

%% 7. Junction (Blurring) calculating
clc;
disp(' *** BLURRING *** ');
disp(' ... Thresholding Gray Matter');
if exist('matlabbatch','var'); clear matlabbatch; end
spm('defaults','fmri'); spm_jobman('initcfg');
matlabbatch{1}.spm.util.imcalc.input = {[T1path,'/GrayMatter',T1ext]};
matlabbatch{1}.spm.util.imcalc.output = 'GMthr';
matlabbatch{1}.spm.util.imcalc.outdir = {''};
matlabbatch{1}.spm.util.imcalc.expression = ['(i1>',num2str(sprintf('%f0',Tlowerthreshold)),').*(i1<',num2str(sprintf('%f0',Tupperthreshold)),')'];
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
spm_jobman('run',matlabbatch);
% maskovani Thresholded GM
if exist('matlabbatch','var'); clear matlabbatch; end
spm('defaults','fmri'); spm_jobman('initcfg');
matlabbatch{1}.spm.util.imcalc.input = {
                                        [T1path,'/GMthr',T1ext]
                                        [MaskAll,',1']
                                        };
matlabbatch{1}.spm.util.imcalc.expression = '(i1-(i2.*i1))>0';
matlabbatch{1}.spm.util.imcalc.output = 'GMthr_masked';
matlabbatch{1}.spm.util.imcalc.outdir = {''};
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
spm_jobman('run',matlabbatch);

% konvoluce s Gaussovym filtrem 5mm
clc;
disp(' ...Convolution GMthr_masked');
disp(['GaussKernel = ',num2str(GaussKernel)]);
if exist('matlabbatch','var'); clear matlabbatch; end
spm('defaults','fmri'); spm_jobman('initcfg');
matlabbatch{1}.spm.spatial.smooth.data = {[T1path,'/GMthr_masked.nii,1']};
matlabbatch{1}.spm.spatial.smooth.fwhm = [GaussKernel GaussKernel GaussKernel];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';
spm_jobman('run',matlabbatch);
% Ulozeni pod nazvem "wBLURRING_sXmm", kde X je filtr
disp('Saving JUNCTION');
copyfile('sGMthr_masked.nii',['wJUNCTION_s',num2str(GaussKernel),'mm.nii']);
delete('sGMthr_masked.nii');

%% 8. Statistics [mean, std]
% subtrakce [subject - template]
if exist('matlabbatch','var'); clear matlabbatch; end
spm('defaults','fmri'); spm_jobman('initcfg');
% ADULT:
% ADULT subtraction
matlabbatch{1}.spm.util.imcalc.input = {
                                    [T1path,'/wJUNCTION_s',num2str(GaussKernel),'mm',T1ext,',1']
                                    [path_to_adult_templates,'template_81zk_junction_blurred',num2str(GaussKernel),'mm_mean.nii,1']
                                    [path_to_adult_templates,'template_81zk_junction_blurred',num2str(GaussKernel),'mm_std.nii,1']
                                    [MaskAll,',1']
                                    [T1path,'/wCORTEX_BinarMask.nii,1']
                                    };
matlabbatch{1}.spm.util.imcalc.outdir = {''};
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

matlabbatch{1}.spm.util.imcalc.output = ['wJUNCTION_subtraction_s',num2str(GaussKernel),'mm_AdultTemplate'];
matlabbatch{1}.spm.util.imcalc.expression = '((i1-i2).*(i1>0.1))-((i1-i2).*(i1>0.1).*(i4>0))';
spm_jobman('run',matlabbatch);
% ADULT subtraction - CortThkMask
matlabbatch{1}.spm.util.imcalc.output = ['wJUNCTION_subtraction_s',num2str(GaussKernel),'mm_AdultTemplateCortThkMask'];
matlabbatch{1}.spm.util.imcalc.expression = '(((i1-i2).*(i1>0.1))-((i1-i2).*(i1>0.1).*(i4>0))).*(i5>0)';
spm_jobman('run',matlabbatch);
%delete('sGMthr_masked.nii');
% ADULT Z-score
matlabbatch{1}.spm.util.imcalc.output = 'wJUNCTION_z_score_AdultTemplate';
matlabbatch{1}.spm.util.imcalc.expression = '(((i1-i2)./i3).*(i1>0.1))-(((i1-i2)./i3).*(i1>0.1).*(i4>0))';
spm_jobman('run',matlabbatch);
% ADULT Z-score - CortThkMask
matlabbatch{1}.spm.util.imcalc.output = 'wJUNCTION_z_score_AdultTemplateCortThkMask';
matlabbatch{1}.spm.util.imcalc.expression = '((((i1-i2)./i3).*(i1>0.1))-(((i1-i2)./i3).*(i1>0.1).*(i4>0))).*(i5>0)';
spm_jobman('run',matlabbatch);
% *************************************************************************
% PEDIATRIC:
% PEDIATRIC subtraction
matlabbatch{1}.spm.util.imcalc.input = {
                                    [T1path,'/wJUNCTION_s',num2str(GaussKernel),'mm',T1ext,',1']
                                    [path_to_pediatric_templates,'template_37hr_junction_blurred',num2str(GaussKernel),'mm_mean.nii,1']
                                    [path_to_pediatric_templates,'template_37hr_junction_blurred',num2str(GaussKernel),'mm_std.nii,1']
                                    [MaskAll,',1']
                                    [T1path,'/wCORTEX_BinarMask.nii,1']
                                    };

matlabbatch{1}.spm.util.imcalc.output = ['wJUNCTION_subtraction_s',num2str(GaussKernel),'mm_PediatricTemplate'];
matlabbatch{1}.spm.util.imcalc.expression = '((i1-i2).*(i1>0.1))-((i1-i2).*(i1>0.1).*(i4>0))';
spm_jobman('run',matlabbatch);
% PEDIATRIC subtraction - CortThkMask
matlabbatch{1}.spm.util.imcalc.output = ['wJUNCTION_subtraction_s',num2str(GaussKernel),'mm_PediatricTemplateCortThkMask'];
matlabbatch{1}.spm.util.imcalc.expression = '(((i1-i2).*(i1>0.1))-((i1-i2).*(i1>0.1).*(i4>0))).*(i5>0)';
spm_jobman('run',matlabbatch);
% PEDIATRIC Z-score
matlabbatch{1}.spm.util.imcalc.output = 'wJUNCTION_z_score_PediatricTemplate';
matlabbatch{1}.spm.util.imcalc.expression = '(((i1-i2)./i3).*(i1>0.1))-(((i1-i2)./i3).*(i1>0.1).*(i4>0))';
spm_jobman('run',matlabbatch);
% PEDIATRIC Z-score - CortThkMask
matlabbatch{1}.spm.util.imcalc.output = 'wJUNCTION_z_score_PediatricTemplateCortThkMask';
matlabbatch{1}.spm.util.imcalc.expression = '((((i1-i2)./i3).*(i1>0.1))-(((i1-i2)./i3).*(i1>0.1).*(i4>0))).*(i5>0)';
spm_jobman('run',matlabbatch);

%% 9. Transformation to the structural space
disp(' ... Transformace do strukturalniho prostoru');
if exist('matlabbatch','var'); clear matlabbatch; end
spm('defaults','fmri');
spm_jobman('initcfg');
matlabbatch{1}.spm.spatial.normalise.write.subj.def = {[T1path,'/iy_',T1filename,'.nii']};
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -100 -70
                                                          78 110 125];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'str';

% ADULT subraction
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {[T1path,'/wJUNCTION_subtraction_s',num2str(GaussKernel),'mm_AdultTemplate',T1ext];[T1path,'/wJUNCTION_subtraction_s',num2str(GaussKernel),'mm_AdultTemplate',T1ext]};
spm_jobman('run',matlabbatch);
copyfile(['strwJUNCTION_subtraction_s',num2str(GaussKernel),'mm_AdultTemplate.nii'],['strJUNCTION_subtraction_s',num2str(GaussKernel),'mm_AdultTemplate.nii']);
delete(['strwJUNCTION_subtraction_s',num2str(GaussKernel),'mm_AdultTemplate.nii']);

% ADULT subraction - CortThkMask
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {[T1path,'/wJUNCTION_subtraction_s',num2str(GaussKernel),'mm_AdultTemplateCortThkMask',T1ext];[T1path,'/wJUNCTION_subtraction_s',num2str(GaussKernel),'mm_AdultTemplateCortThkMask',T1ext]};
spm_jobman('run',matlabbatch);
copyfile(['strwJUNCTION_subtraction_s',num2str(GaussKernel),'mm_AdultTemplateCortThkMask.nii'],['strJUNCTION_subtraction_s',num2str(GaussKernel),'mm_AdultTemplateCortThkMask.nii']);
delete(['strwJUNCTION_subtraction_s',num2str(GaussKernel),'mm_AdultTemplateCortThkMask.nii']);

% ADULT z-score
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {[T1path,'/wJUNCTION_z_score_AdultTemplate',T1ext];[T1path,'/wJUNCTION_z_score_AdultTemplate',T1ext]};
spm_jobman('run',matlabbatch);
copyfile('strwJUNCTION_z_score_AdultTemplate.nii','strJUNCTION_z_score_AdultTemplate.nii');
delete('strwJUNCTION_z_score_AdultTemplate.nii');

% ADULT z-score - CortThkMask
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {[T1path,'/wJUNCTION_z_score_AdultTemplateCortThkMask',T1ext];[T1path,'/wJUNCTION_z_score_AdultTemplateCortThkMask',T1ext]};
spm_jobman('run',matlabbatch);
copyfile(['strwJUNCTION_z_score_AdultTemplateCortThkMask.nii'],['strJUNCTION_z_score_AdultTemplateCortThkMask.nii']);
delete(['strwJUNCTION_z_score_AdultTemplateCortThkMask.nii']);

% PEDIATRIC subraction
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {[T1path,'/wJUNCTION_subtraction_s',num2str(GaussKernel),'mm_PediatricTemplate',T1ext];[T1path,'/wJUNCTION_subtraction_s',num2str(GaussKernel),'mm_PediatricTemplate',T1ext]};
spm_jobman('run',matlabbatch);
copyfile(['strwJUNCTION_subtraction_s',num2str(GaussKernel),'mm_PediatricTemplate.nii'],['strJUNCTION_subtraction_s',num2str(GaussKernel),'mm_PediatricTemplate.nii']);
delete(['strwJUNCTION_subtraction_s',num2str(GaussKernel),'mm_PediatricTemplate.nii']);

% PEDIATRIC subraction - CortThkMask
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {[T1path,'/wJUNCTION_subtraction_s',num2str(GaussKernel),'mm_PediatricTemplateCortThkMask',T1ext];[T1path,'/wJUNCTION_subtraction_s',num2str(GaussKernel),'mm_PediatricTemplateCortThkMask',T1ext]};
spm_jobman('run',matlabbatch);
copyfile(['strwJUNCTION_subtraction_s',num2str(GaussKernel),'mm_PediatricTemplateCortThkMask.nii'],['strJUNCTION_subtraction_s',num2str(GaussKernel),'mm_PediatricTemplateCortThkMask.nii']);
delete(['strwJUNCTION_subtraction_s',num2str(GaussKernel),'mm_PediatricTemplateCortThkMask.nii']);

% PEDIATRIC z-score
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {[T1path,'/wJUNCTION_z_score_PediatricTemplate',T1ext];[T1path,'/wJUNCTION_z_score_PediatricTemplate',T1ext]};
spm_jobman('run',matlabbatch);
copyfile('strwJUNCTION_z_score_PediatricTemplate.nii','strJUNCTION_z_score_PediatricTemplate.nii');
delete('strwJUNCTION_z_score_PediatricTemplate.nii');

% PEDIATRIC z-score - CortThkMask
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {[T1path,'/wJUNCTION_z_score_PediatricTemplateCortThkMask',T1ext];[T1path,'/wJUNCTION_z_score_PediatricTemplateCortThkMask',T1ext]};
spm_jobman('run',matlabbatch);
copyfile(['strwJUNCTION_z_score_PediatricTemplateCortThkMask.nii'],['strJUNCTION_z_score_PediatricTemplateCortThkMask.nii']);
delete(['strwJUNCTION_z_score_PediatricTemplateCortThkMask.nii']);

%% 10. Saving variables, Final clearing of data
% delete('*.mat');
% delete('GM*.nii');
% delete('WhiteMatter.nii');
% delete('GrayMatter.nii');
% delete('w*.nii');
% delete('*T1*.nii');
clc;
disp('E N D   O F   A L G O R I T H M');