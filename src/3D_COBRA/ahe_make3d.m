% make data files for 3D COBRA

% use the exp and sample structs as input files to
% make 3D COBRA files

exp_con = load('se_capped_exp.mat').data
exp_dat = load('se_capped_sample.mat').sample


rod_info = exp_dat.rod_info;

qz = exp_dat.qz;
data_final = zeros(length(qz),10);
sample_tag = 'Se_capped_new'
file_loc = "C:\Users\HasithPerera\Documents\MATLAB\COBRA-cenlab\Data\3D\"

for exp_id=1:5
    
    %output file named using 10 rod
    out_name = sprintf('%s%s_%d_3D_data1.dat',file_loc,sample_tag,exp_con(exp_id).COBRA3(1))
    
    
    %process rods 00 to 33
    for j=1:10
        if j==1
            scan_id = exp_con(exp_id).rod00;
        else
            scan_id = exp_con(exp_id).COBRA3(j-1)
        end
        rod_index =find( rod_info(:,1)==scan_id)
        if length(rod_index)>1
            warning("More than one scan info found")
        end
        scan_info = rod_info(rod_index,:)
        file_name = strcat(exp_dat.data_loc,sprintf(exp_dat.fileformat,scan_info(1:3)))
        
        %%load data
        raw_CTR = load(file_name);

        
        CTR_sig = interp1(raw_CTR(:,1),abs(raw_CTR(:,2)),qz,'pchip');
                
        lamda = exp_dat.lamda
        L = qz;
        C_0 = 3.905;
        
        sin_th = lamda * L ./ (2 * C_0);
        
        %correction for sin theta
        data = CTR_sig.*sin_th;
        
        
        %clip to measurements
        slct = find(qz < scan_info(4))
        data(slct) = 0;
        
        slct = find(qz > scan_info(5))
        data(slct) = 0;
        
        data_final(:,j)=data;
        
%         semilogy(qz,data,'o')
%         hold on 
%         semilogy(raw_CTR(:,1),raw_CTR(:,2)) 
%         
        
        
    end
%     out_name
        save(out_name,'data_final','-ascii') 
%         save data1.dat data -ascii
    
end
