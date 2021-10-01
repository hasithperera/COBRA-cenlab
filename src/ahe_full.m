% full 1D processing script
% Stop at r1 <.08
% need to run at the analysis folder

%dependancies
%   1. add anal_00l to the path
%   2. ahe_process_1D is needed in the path
%   3. ahe_ref folder should be available in the pwd for this to work

% How to use
% change the sample_id and the file variable to select correct files to
% process



%data source

addpath('C:\Users\HasithPerera\Documents\MATLAB\COBRA\Data')
addpath('C:\Users\HasithPerera\Documents\MATLAB\COBRA\Data\se_capped_org')

sample_id = 104;
scale = 12;

name= split(pwd(),"_");
root_tag = 'ahe';
if strcmp(name{length(name)},num2str(sample_id))== 0
    system(strcat("mkdir ",root_tag,'_',num2str(sample_id)));
    system(strcat("copy ahe_ref ",root_tag,'_ ',num2str(sample_id)))
    % %se capped files
    % file = strcat('cleaned_Se_capped_1ML_FeSe_ISTO_light_',num2str(sample_id),'_00L_CTR_R_Data_out.dat');
    
    
    folder_loc = strcat(root_tag,"_",num2str(sample_id));
    cd(folder_loc)
    system('pwd')
    system('mkdir old')
%     system('move bfcr*.dat old/')
    system('move Peakedge*.dat old/')
    
    if exist('old/zstr.dat')== 0
        %save the old zstr file with peak edges
        system('move zstr.dat old/zstr.dat')
    end
    clc
    fprintf('Files cleaned\n')
    
end

%Se Si capped files
%file = strcat('cleaned_Se_capped_1ML_FeSe_STO_light_',num2str(sample_id),'_00L_CTR_R_Data_out.dat');

%se capped files
file  = strcat("cleaned_Se_capped_1ML_FeSe_ISTO_light_",num2str(sample_id),"_00L_CTR_R_Data_out.dat")

%%
qz = 0.0001:0.025:4.5001;
d = load(file);

lamda = 12.398/20;
C_0 = 3.905;  %for STO
L = qz;

sin_th = lamda * L ./ (2 * C_0);

data = interp1(d(:,1),abs(d(:,2)),qz,'pchip')'

%before correction
semilogy(data)
hold on


%correction for sin theta
data = data.*sin_th';


slct = find(qz < min(d(:,1)));
data(slct) = 0;

% semilogy(d(:,1),abs(d(:,2)),'.-')
% hold on
% semilogy(qz,data)

semilogy(data)

% hold on
% load 'org_data/data1.dat'
hold on
semilogy(data*scale)

data = data*scale;

% sum(d(:,2)<0)
close()

save data1.dat data -ascii
save qz.dat qz -ascii

% ahe_sim
ahe_process_1D
cd ..
load('ahe_7/zstr.dat')
hold on
plot(zstr)