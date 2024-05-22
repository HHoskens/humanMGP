%% COMPILE PHENOTYPE FILES
close all; clear all; clc
addpath(genpath('/mnt/BHServer4/FaceBase_3/Pipeline/Aid_Functions/'));

loadPath = '/mnt/BHServer4/FaceBase_3/Data/';
studyPath = '/mnt/BHServer4/FaceBase_3/Analysis/HumanMGP/';
cd(studyPath)

parpool('LocalSingle',8); maxNumCompThreads(8)

%% LOAD DATA
%  ATLAS
[v,f] = read_ply_hh([loadPath 'Images/Atlas/Dense_5k/dense_5k_atlas.ply']);
atlas = shape3D; atlas.Vertices = v; atlas.Faces = f; atlas.SingleColor = [.9 .9 .9]; atlas.Material = 'Dull';

% Center Vertices
tmp_avg = mean(atlas.Vertices,1);
atlas.Vertices = atlas.Vertices - tmp_avg;

% Quick check
scan1 = clone(atlas);
v = viewer(scan1); v.SceneLightVisible = true;



% META DATA
meta = readtable([studyPath 'Data/Harry Data Compiled.xlsx']);
head(meta)

% Add gen IIDs
opts = detectImportOptions('/mnt/BHServer/Hallgrimsson/Collaborations/FaceBase/Data/Covariates/Raw_Imports/Pitt/TDFN-Export_landmarks_measurements.csv');
opts.VariableTypes(1:2) = {'char'};
metakey = readtable('/mnt/BHServer/Hallgrimsson/Collaborations/FaceBase/Data/Covariates/Raw_Imports/Pitt/TDFN-Export_landmarks_measurements.csv',opts);

meta.IID = meta.SubjectID;
[m,indm] = ismember(meta.IID,metakey.StudyID); tabulate(m)
meta.IID(m) = metakey.FBID(indm(m));



% COVARIATE DATA
covs = meta(:,[1:14 end]);
indiv_id = covs.SubjectID;
image_id = erase(covs.fileName,'.mat');
gen_iid = covs.IID;

% Rename
[~,idx] = ismember({'Sex','Age_yrs_','height','weight'},covs.Properties.VariableNames);
covs.Properties.VariableNames(idx) = {'Sex','Age','Height','Weight'};

% Add imageID, age2, age3 and centroid size
covs.ImageID = image_id;
covs.IID = gen_iid;
covs.Age2 = covs.Age.^2;
covs.Age3 = covs.Age.^3;
%covs.CSize = csize';



% GENOMIC DATA         
famTANZ = readtable([loadPath 'Genetics/Tanz/QC/Prune/Spritz_imputed_qc_prune.fam'],'FileType','text');
famTANZ.Properties.VariableNames(1:2) = {'FID','IID'};
fam3DFN = readtable([loadPath 'Genetics/Pitt/QC/Prune/Marazita_imputed_qc_prune.fam'],'FileType','text');
fam3DFN.Properties.VariableNames(1:2) = {'FID','IID'}; fam3DFN.FID = cellstr(num2str(fam3DFN.FID)); fam3DFN.IID = cellstr(num2str(fam3DFN.IID)); 


%%
%
% === 3D FACIAL NORMS ===
%
%
%% Load Shapes
% SHAPE DATA
files = dir('/mnt/BHServer4/FaceBase_3/Data/Images/Pitt/2_Registered/Dense_5k/Scans/Original/*.ply');

lms = nan(atlas.nVertices,3,length(files));
parfor i=1:length(files)
    [v,f] = read_ply_hh([files(i).folder '/' files(i).name]);
    [vr,~] = read_ply_hh([replace(files(i).folder,'Original','Reflected') '/' files(i).name]);

    % align orig and reflected copy
    data = cat(3,v,vr);
    iter = 3; scale = false; reflect = false; display = false;
    [vsym,~,~,~] = GeneralizedProcrustesAnalysis(data,atlas,iter,scale,reflect,display);
    % symmetrize
    lms(:,:,i) = mean(vsym,3);
end

nLM = size(lms,1);% 5629
nID = size(lms,3);% 2267

% Align full sample
iter = 3; scale = true; reflect = false; display = false;
[lms,avg,csize,~] = GeneralizedProcrustesAnalysis(lms,atlas,iter,scale,reflect,display);

%% Match Pheno and Geno data
[m1,idx1] = ismember(erase({files.name},'.ply'),covs.SubjectID); tabulate(m1)

% Reduce
sym3DFN = lms(:,:,m1);
cov3DFN = covs(idx1(m1),:);
n3DFN = sum(m1);% 2343

% Add centroid size
cov3DFN.CSize = csize(m1)';

%% Quick check
data = permute(sym3DFN,[2 1 3]);
data = reshape(data,atlas.nVertices*3,size(sym3DFN,3));

space = shapePCA;
space.RefScan = clone(atlas);
getAverage(space,data);
getModel(space,data);
space.stripPercVar(98)

% pc to lm
pc3DFN = space.AvgVec + space.EigVec*space.Tcoeff';
pc3DFN = reshape(pc3DFN,3,atlas.nVertices,n3DFN);

% get z-scores 
outM = sqrt(sum(((space.Tcoeff-repmat(space.AvgCoeff',space.n,1))./repmat(space.EigStd',space.n,1)).^2,2))';
zM = (outM-mean(outM))/std(outM);

OutlierIndex = find(abs(zM)>3);
cprintf([0 0 1],[num2str(length(OutlierIndex)) ' outliers were detected\n'])

% Manual check of outliers
scan = clone(atlas); v = viewer(scan); v.SceneLightVisible = true;
for i=1:length(OutlierIndex), scan.Vertices = pc3DFN(:,:,OutlierIndex(i)); scan.VertexSize = 20; scan.Visible = true; pause; end

%% Export
% writetable(cov3DFN,[studyPath 'Data/COV_3DFN_dense.csv'],'Delimiter',',');
% 
% pc = space.Tcoeff;
% eigvec = space.EigVec;
% eigval = space.EigVal;
% avg = space.AvgVec';
% 
% % add IDs
% pc = [cov3DFN.IID num2cell(pc)];
% writecell(pc,[studyPath 'Data/PC_3DFN_dense.csv'],'Delimiter',',');
% writematrix(eigvec,[studyPath 'Data/EIGVEC_3DFN_dense.csv'],'Delimiter',',');
% writematrix(eigval,[studyPath 'Data/EIGVAL_3DFN_dense.csv'],'Delimiter',',');
% writematrix(avg,[studyPath 'Data/AVG_3DFN_dense.csv'],'Delimiter',',');
% 
% save([studyPath 'Data/PCASPACE_3DFN_dense.csv'],'space','cov3DFN','-v7.3')

%% 
% === TANZANIA ===

%% Match Pheno and Geno data
idx1 = ismember(covs.Dataset,'Tanzania'); tabulate(idx1)
idx2 = ismember(covs.IID,famTANZ.IID); tabulate(idx2)
% Remove duplicates
idx3 = ismember(covs.Method,'Harry'); tabulate(idx3)

idx = idx1 & idx2 & idx3;

% Reduce
lmTANZ = lms(:,:,idx);
covTANZ = covs(idx,:);
nTANZ = sum(idx);% 3478

%% Align and Symmetrize data
left =[40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65];
right = [14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39];

idx_orig = 1:sAtlas.nVertices;
idx_refl = 1:sAtlas.nVertices; idx_refl([left right]) = idx_refl([right left]); 

orig = lmTANZ;
refl = orig; refl(:,1,:) = -1 * refl(:,1,:); refl = refl(idx_refl,:,:);

iter = 3; scale = true; reflect = false; display = false;  
totalLM = cat(3,orig,refl);
[alignedTANZ,avgTANZ,csizeTANZ,~] = GeneralizedProcrustesAnalysis(totalLM,sAtlas,iter,scale,reflect,display);

orig = alignedTANZ(:,:,1:nTANZ);
refl = alignedTANZ(:,:,nTANZ+1:end);
symTANZ = (orig+refl)/2;

% Add centroid size
covTANZ.CSize = csizeTANZ(1:nTANZ)';

%% Quick check
data = permute(symTANZ,[2 1 3]);
data = reshape(data,sAtlas.nVertices*3,size(symTANZ,3));

space = shapePCA;
space.RefScan = clone(sAtlas);
getAverage(space,data);
getModel(space,data);
space.stripPercVar(98)

% get z-scores 
outM = sqrt(sum(((space.Tcoeff-repmat(space.AvgCoeff',space.n,1))./repmat(space.EigStd',space.n,1)).^2,2))';
zM = (outM-mean(outM))/std(outM);

OutlierIndex = find(abs(zM)>3);
cprintf([0 0 1],[num2str(length(OutlierIndex)) ' outliers were detected\n'])

% Manual check of outliers
scan = shape3D; v = viewer(scan); v.SceneLightVisible = true;
for i=1:length(OutlierIndex), scan.Vertices = symTANZ(:,:,OutlierIndex(i)); scan.VertexSize = 20; scan.Visible = true; pause; end


%% Export
% writetable(covTANZ,[studyPath 'Data/COV_TANZ.csv'],'Delimiter',',');
% 
% sym = permute(symTANZ,[2 1 3]);
% sym = reshape(sym,3*sAtlas.nVertices,nTANZ);
% sym = num2cell(sym');
% % add IDs
% sym = [covTANZ.IID sym];
% 
% writecell(sym,[studyPath 'Data/LM_TANZ.csv'],'Delimiter',',');

%% END
