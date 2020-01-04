%%

%Outputs the Raw (filtered) and inferred traces for CalmAn data generated with drgCalmAn_script

function [raw,inferred]=drgGetCAtraces(Y,Ag,Cg,bg,fg,Cng,options)
 
no_images = size(Cg,2);
bg = double(bg);
Cg = double(Cg);
fg = double(fg);
nA = full(sqrt(sum(Ag.^2))');
[K,~] = size(Cg);

Ag = Ag/spdiags(nA,0,K,K);    % normalize spatial components to unit energy


Cg = bsxfun(@times,Cg,nA(:)); %spdiags(nA,0,K,K)*C;

nROIs = size(Ag,2);     % number of ROIs
nb = size(fg,1);     % number of background components
%nA = full(sum(A.^2))';  % energy of each row
%Y_r = spdiags(nA,0,nr,nr)\(A'*Y- (A'*A)*C - (A'*full(b))*f) + C;
 
AY = mm_fun(Ag,Y);
Y_rg = (AY- (Ag'*Ag)*Cg - full(Ag'*double(bg))*fg) + Cg;


[~,Dfg] = extract_DF_F(Y,Ag,Cg,[],options,AY);

 
Dfmat=repmat(Dfg,1,no_images);
raw=Y_rg./Dfmat;
inferred=Cg./Dfmat;
 

