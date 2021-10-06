
%cmd to clean the S#, in thespec file
% sed -i 's/#S,//g' rod_se_capped.csv
% single cmd to get form spec file
%  cat <file> | grep "rodscan"| sed 's/  /,/g' | sed 's/ /,/g' | sed 's/#S,//g'

rod_file = './Data/rod_info/rod_Se_capped.csv';
rod_info = load(rod_file);

%start with rod10 scan
start_indx = find((rod_info(:,2)==1 & rod_info(:,3)==0))
    
qz = 0.0001:0.025:4.5001;
lamda = 12.398/20;
C_0 = 3.905;  %for STO
L = qz;

sin_th = lamda * L ./ (2 * C_0);

sample_name = 'Se_capped'
root_file = 'cleaned_Se_capped_1ML_FeSe_ISTO_light_';
file_end = 'L_CTR_R_Data_out.dat';

clc

for i = 1:length(start_indx)
  
    fprintf('%d: %d\n',i,rod_info(start_indx(i),1))

    %store notebook info
    
    %data(i).sample = 'Se_capped'
    %data(i).qz = qz;
    data(i).temp =input("Sample temp:")
    data(i).exp_id = i;
    clc
    fprintf('3D Sequence\n')
    rod_info(start_indx(i)+(0:8),1)
    if(input("3D scan sequence correct")==1)
    data(i).COBRA3 =rod_info(start_indx(i)+(0:8),1)
    else
        data(i).COBRA3 = []
    end
    data(i).rod00 = input("Scan ID of rod 00:");
    
    data(i).is_actual_rod00 = input("rod 00 actual? (1/0):");
    data(i).light = input('Light 0- NO,1-UV,2-:wht,3:wht-cts');
end
% 
% tmp_i = input(sprintf('enter selection for scan start(1-%d):',length(start_indx)));


% %save a mat file with the data

sample.rod_info = rod_info;
sample.fileformat = 'cleaned_Se_capped_1ML_FeSe_ISTO_light_%d_%d%dL_CTR_R_Data_out.dat';
sample.data_loc = 'C:\Users\HasithPerera\Documents\MATLAB\COBRA-cenlab\Data\processed\';
sample.lamda = 12.398/20;
sample.C_0 = 3.905;
sample.tag = sample_name;
sample.qz = qz;
% sample.experiments = data;

save(sprintf('./Data/%s_sample.mat',sample_name),'sample')

