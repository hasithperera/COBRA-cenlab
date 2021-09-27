% Preprocessing - BG removal

% enter sc_num = with the scan ID
% f_spec = name of the file
% shift: change the center of the signal window
% sig: signal window shape 
addpath(genpath('D:\Research\COBRA\COBRA CTR Data Files'))
%[scan_id, tag] = ahe_importfile_ms('D:\Research\COBRA\COBRA CTR Data Files\SPEC File Folder\ahe_scan_data_03302021.csv');
%[scan_id, tag] = ahe_importfile_ms('E:\Research\COBRA\COBRA CTR Data Files\SPEC File Folder\ahe_scan_data_04052021.csv');

[scan_id, tag] = ahe_importfile_ms('D:\Research\COBRA\COBRA CTR Data Files\SPEC File Folder\ahe_scan_data_04052021.csv');

sc_num = [11]; %enter scan ID
select_scans = mod(find(scan_id == sc_num),114)
if(exist('k'))
fprintf("%d-%d\n",select_scans(k),sc_num(k))
end
retry = str2num(input("retry (1/0):","s"));
if (exist('k'))
    if (retry)
        k = k
    else 
        k = k + 1
    end
else
    k = 1;
end
fprintf("Scan id: %d\n",sc_num(k))
f_spec = 'sb2Te3_20QL_capped_1ML_FeSe_ISTO_04052021_1';
shift = [0,0];
sig = [8,20];
window = [24,40];

state = 0;
i = select_scans(k);


file_name = ahe_fun(scan_id(i),strcat(tag(i,:)),f_spec,shift,sig,window,1);
files = ls(strcat(file_name,"*"));

%remove previous files with the same names
if size(files,1)>0
    fprintf("removing old files.\n");
    for j=1:size(files,1)
    delete(strcat('../Data/',files(j,:)))
    end
end

%preocess bg removal
ahe_fun(scan_id(i),strcat(tag(i,:)),f_spec,shift,sig,window,0);
ahe_CTR_correction(file_name)
title(strcat("#",num2str(scan_id(i)),"\_",strcat(tag(i,:))))

% plot_tag = "sb2Te3_";
plot_tag = "Se_Si_capped_";
saveas(gcf,strcat("../Data/plots/",plot_tag,num2str(scan_id(i)),"_",strcat(tag(i,:),".png")))

%% old process - ahe week 1
% % ahe_fun(36,'22','Se_capped_1ML_FeSe_ISTO_light_04032021_1',shift)
% % while i<size(scan_id,1)
% %     fprintf("%d:%s %d,of%d\n",scan_id(i),tag(i,:),i,size(scan_id,1))
% %     if state == 0
% %         shift = [0,0];
% %     end
% %     ahe_fun(scan_id(i),strcat(tag(i,:)),'Se_capped_1ML_FeSe_ISTO_light_04032021_1',shift);
% %  
% %     clc
% %     usr_in = input(strcat("Is scan id:",num2str(i), " ok (n)>"),"s")
% %     if usr_in == 'n'
% %         log_ahe = fopen("../Data/ahe_log.csv",'a+');
% %         fprintf(log_ahe,"Shift needed in scan_id %d \n",scan_id(i));
% %         
% %         fclose(log_ahe);
% %         while (1)
% %             try
% %                 subplot(231)
% %                 hold on
% %                 fprintf(2,"Select starting edge\n")
% %                 
% %                 shift = str2num(input("x,y shift>","s"))
% %                 if length(shift)==2
% %                     
% %                     ahe_fun(scan_id(i),strcat(tag(i,:)),'Se_capped_1ML_FeSe_ISTO_light_04032021_1',shift);
% %                     shift
% %                     usr_in = input(strcat("Is scan id:",num2str(i), " ok (n)>"),"s")
% %                     if usr_in == 'n'
% %                         state = 1;
% %                     else
% %                         i = i+1;
% %                         break
% %                     end
% %                 end
% %             catch
% %                 
% %             end
% %         end
% %         fprintf(2,"change location")
% %         
% %         
% %     else
% %         i = i+1;
% %     end
% %     
% %     
% % end
