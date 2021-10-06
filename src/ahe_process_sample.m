
%se capped
sample_file = '../Data/se_capped_sample.mat';
save_loc = '../Data/3D/'
load(sample_file)

rod_info = sample.rod_info;
known_scale =[10 12 12 12 12]; %based on 1D experimental values from the 1D reconstruction

for exp_id = [1:length(sample.experiments)]
    exp_struct = sample.experiments;
    
    data = zeros([length(sample.qz),10]);
    k = 1
    for scanid = [exp_struct(exp_id).rod00,exp_struct(exp_id).COBRA3']
        scan_file = get_filename(scanid,sample);
        dat_1D = get_data(scan_file,sample,known_scale(exp_id));
        data(:,k) = dat_1D;
        semilogy(data(:,k),'.')
        hold on
        
        k = k + 1;
        
    end 
    save(strcat(save_loc,sample.tag,'_',num2str(exp_struct(exp_id).COBRA3(1)),"_3D_data1.dat"),'data','-ascii')
%     break
    %     scan_id = exp_struct(exp_id).rod00;
    
%     break
    %     exp_struct(exp_id).temp;
end




%local function definitions
%There should not be any commands below the function definitions

% sample file structure

%       rod_info: [49×9 double]
%      fileformat: string
%        data_loc: string
%           lamda: float
%             C_0: float
%             tag: string
%              qz: array
%     experiments: struct

function file_name = get_filename(scan_id,sample)
% convers a scan id to a file name

indx = find(sample.rod_info(:,1)==scan_id);
file_name = strcat(sample.data_loc,'\',sprintf(sample.fileformat,sample.rod_info(indx,1:3)));
end


function scaled_data = get_data(file,sample,scale)
% read raw data from file
% use a scaling factor : only way to decide if this is correct is to look
% at the construction and make the substrate look natural as possible.


% Need to handel the scaling factor in the case of a 3D scan
% use the one used in the 1D reconstruction (current option)


qz = sample.qz;
lamda = sample.lamda;
C_0 = sample.C_0;  %for STO
L = qz;

d = load(file);

sin_th = lamda * L ./ (2 * C_0);
data_1D = interp1(d(:,1),abs(d(:,2)),qz,'pchip')'

%correction for sin theta
data_1D = data_1D.*sin_th';

slct = find(qz < min(d(:,1)));
data_1D(slct) = 0;

scaled_data = data_1D*scale;

% semilogy(data_1D)
% hold on
% semilogy(data_1D)
% hold on
% semilogy(data_1D*scale)

end
