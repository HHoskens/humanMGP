%% COMPILE PHENOTYPE FILES
close all; clear all; clc
addpath(genpath('/mnt/BHServer4/FaceBase_3/Pipeline/Aid_Functions/'));

loadPath = '/mnt/BHServer4/FaceBase_3/Data/';
studyPath = '/mnt/BHServer4/FaceBase_3/Analysis/HumanMGP/';
cd(studyPath)

parpool('LocalSingle',8); maxNumCompThreads(8)

%% LOAD DATA
%  ATLAS
% Sparse
v = readcell([loadPath 'Images/Atlas/Sparse_65_ears/landmarks_65_list.xlsx']);
sAtlas = shape3D; sAtlas.Vertices = cell2mat(v(2:end,2:end)); sAtlas.VertexSize = 20;

% Dense
[v,f] = read_ply_hh([loadPath 'Images/Atlas/Dense_2k_ears/dense_2k_ears_atlas.ply']);
dAtlas = shape3D; dAtlas.Vertices = v; dAtlas.Faces = f; dAtlas.SingleColor = [.9 .9 .9]; dAtlas.Material = 'Dull';

% Center Vertices
tmp_avg = mean(dAtlas.Vertices,1);
dAtlas.Vertices = dAtlas.Vertices - tmp_avg;
sAtlas.Vertices = sAtlas.Vertices - tmp_avg;

% Quick check
scan1 = clone(dAtlas);
scan2 = clone(sAtlas);
v = viewer(scan1); v.SceneLightVisible = true;
scan2.RenderAxes = v.RenderAxes; scan2.Visible = true;



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



% SHAPE DATA
lms = table2array(meta(:,15:end-1));
nLM = size(lms,2)/3;% 65
nID = size(lms,1);% 31724

% Reshape
lms = reshape(lms',3,nLM,size(lms,1));
lms = permute(lms,[2 1 3]);

% Make sure they are aligned
iter = 3; scale = false; reflect = false; display = false;
[lms,avg,csize,~] = GeneralizedProcrustesAnalysis(lms,sAtlas,iter,scale,reflect,display);



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
covs.CSize = csize';



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
%% Match Pheno and Geno data
idx1 = ismember(covs.Dataset,'Pittsburgh'); tabulate(idx1)
idx2 = ismember(covs.IID,fam3DFN.IID); tabulate(idx2)

idx = idx1 & idx2;

% Reduce
lm3DFN = lms(:,:,idx);
cov3DFN = covs(idx,:);
n3DFN = sum(idx);% 4559

%% Align and Symmetrize data
left =[40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65];
right = [14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39];

idx_orig = 1:sAtlas.nVertices;
idx_refl = 1:sAtlas.nVertices; idx_refl([left right]) = idx_refl([right left]); 

orig = lm3DFN;
refl = orig; refl(:,1,:) = -1 * refl(:,1,:); refl = refl(idx_refl,:,:);

iter = 3; scale = false; reflect = false; display = false;  
totalLM = cat(3,orig,refl);
[aligned3DFN,avg3DFN,~,~] = GeneralizedProcrustesAnalysis(totalLM,sAtlas,iter,scale,reflect,display);

orig = aligned3DFN(:,:,1:n3DFN);
refl = aligned3DFN(:,:,n3DFN+1:end);
sym3DFN = (orig+refl)/2;

% Quick check
scan = shape3D; v = viewer(scan); v.SceneLightVisible = true;
for i=1:nID, scan.Vertices = sym3DFN(:,:,i); scan.VertexSize = 20; scan.Visible = true; pause; end

%% Export
% writetable(cov3DFN,[studyPath 'Data/COV_3DFN.csv'],'Delimiter',',');
% 
% sym = permute(sym3DFN,[2 1 3]);
% sym = reshape(sym,3*sAtlas.nVertices,n3DFN);
% sym = num2cell(sym');
% % add IDs
% sym = [cov3DFN.IID sym];
% 
% writecell(sym,[studyPath 'Data/LM_3DFN.csv'],'Delimiter',',');

%% 
% === TANZANIA ===

%% Match Pheno and Geno data
idx1 = ismember(covs.Dataset,'Tanzania'); tabulate(idx1)
idx2 = ismember(covs.IID,famTANZ.IID); tabulate(idx2)

idx = idx1 & idx2;

% Reduce
lmTANZ = lms(:,:,idx);
covTANZ = covs(idx,:);
nTANZ = sum(idx);% 6967

%% Align and Symmetrize data
left =[40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65];
right = [14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39];

idx_orig = 1:sAtlas.nVertices;
idx_refl = 1:sAtlas.nVertices; idx_refl([left right]) = idx_refl([right left]); 

orig = lmTANZ;
refl = orig; refl(:,1,:) = -1 * refl(:,1,:); refl = refl(idx_refl,:,:);

iter = 3; scale = false; reflect = false; display = false;  
totalLM = cat(3,orig,refl);
[alignedTANZ,avgTANZ,~,~] = GeneralizedProcrustesAnalysis(totalLM,sAtlas,iter,scale,reflect,display);

orig = alignedTANZ(:,:,1:nTANZ);
refl = alignedTANZ(:,:,nTANZ+1:end);
symTANZ = (orig+refl)/2;

% Quick check
scan = shape3D; v = viewer(scan); v.SceneLightVisible = true;
for i=1:nTANZ, scan.Vertices = symTANZ(:,:,i); scan.VertexSize = 20; scan.Visible = true; pause; end

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
