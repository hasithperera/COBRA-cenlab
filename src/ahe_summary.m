
%Se capped
a ='C:\Users\HasithPerera\Documents\MATLAB\COBRA\Cheng_Cen_EPSCoR_Project_Light_Induced_FeSe_STO\Analysis for Se_capped FeSe_SrTiO3\00L Analysis\';
id_a = [7,46,101,104];
sample_a = 1
%id_a = [7]


%Se_Si capped
b = 'C:\Users\HasithPerera\Documents\MATLAB\COBRA\Cheng_Cen_EPSCoR_Project_Light_Induced_FeSe_STO\Analysis for Se_Si_capped\';
% = [182,184,186,274,352,358];
id_b = [9,149,176,182,184,186,274,352,358];
sample_b = 2;




%change folder and id based on the setting you need to process
folder_loc = a;
id = id_a;
sample = sample_a;

physical_scale = 0.14462963;

%settings
export_VI_data = 0; %set to 1 to export ahe_ED.dat for each folder
manual_peak_select = 0;
output_filename = 'Peakedge_5_raw.dat'


%Manual peak selection
% In this mode select 6 point to detect the 3 peaks
% this is done in 2 ranges for the film and the substrate

peak_extract = 0;
peakedge_file = 'Peakedge_5_raw.dat';


visualize_selection = 1;
% subplot_x = 2;
% subplot_y = 2;

% sample 2
% subplot_x = 3;
% subplot_y = 3;

summary_location = 'C:\Users\HasithPerera\Documents\MATLAB\COBRA-cenlab\Data\'

% peak extraction
% ahe_peak(zstr,peakedge,step_size)
% this function contains functions from Dr. Zhou
% Return : CoM,CoArea
% once the data is there export them to an array

for i = 1:length(id)
    fprintf("Processing scan:%d\n",id(i))
    
    %clear old data
    clear('dat')
    clear('peakedge')
    
    
    dat = load(strcat(folder_loc,'ahe_',num2str(id(i)),'\zstr.dat'));
    
    % Manual peak selection
    if manual_peak_select == 1
        
        
        %load other scan line Se_capped/ahe_7 as ref
        
        if exist('dat_old')==0
            dat_old = load('C:\Users\HasithPerera\Documents\MATLAB\COBRA\Cheng_Cen_EPSCoR_Project_Light_Induced_FeSe_STO\Analysis for Se_capped FeSe_SrTiO3\00L Analysis\ahe_7\zstr.dat');
        end
        
        plot(dat_old)
        hold on
        
        plot(dat)
        %ylim([-.0005 .02])
        %         xlim([50 150])
        %         xlim([220 inf])
        %         break
        pk_edge = ginput(16);
        
        %         %selecting the top film
        %         xlim([250 450])
        %         pk_2 = ginput(6);
        %         pk_edge = [pk_edge;pk_2];
        
        %save peak selection data to Peakedge.dat
        cd(strcat(folder_loc,'ahe_',num2str(id(i))))
        save(output_filename,'pk_edge','-ascii');
        clf()
        
    end
    
    % save scale corrected data
    if export_VI_data == 1
        fprintf('Exporting ahe_ED.dat for labVIEW\n')
        cd(strcat(folder_loc,'ahe_',num2str(id(i))))
        step_size=0.14462963; %based on s3tpram file line 104
        z_scale = [1:length(dat)]*step_size;
        
        data = [z_scale' ,dat];
        save ahe_ED.dat data -ascii -tabs
    end
    %
    
    %peak extraction
    if peak_extract ==1
        peakedge = load(strcat(folder_loc,'ahe_',num2str(id(i)),'\',peakedge_file));
        if contains(peakedge_file,'raw')
            peakedge = format_rawpeak(peakedge)
        end
        [CoM,CoArea] = ahe_peakdetect(dat,peakedge,physical_scale,0)
        
        %save to file
        
        data_CoM = [sample,id(i),CoM];
        data_CoA = [sample,id(i),CoArea];
        
        cd(summary_location)
        
        save ahe_CoM.dat data_CoM -ascii -tabs -append
        save ahe_CoA.dat data_CoA -ascii -tabs -append
    end
    
    % code to show selection plots
    if visualize_selection == 1
        
%         subplot(subplot_y,subplot_x,i)

        peakedge = load(strcat(folder_loc,'ahe_',num2str(id(i)),'\',peakedge_file));
        
        if contains(peakedge_file,'raw')
            %correction for raw peakdata
            peakedge = format_rawpeak(peakedge);
        end
         [CoM,CoA] = ahe_peakdetect(dat,peakedge,physical_scale,1);
        for pk = 1:2:length(CoA)
            plot(CoA(pk)*[1 1],[0 .08],':black')
        end
       
        %         xlim([30 60])
        ylim([0 0.08])
        title(strcat('scan id:',num2str(id(i))))
    end
end

if visualize_selection ==1
%     ahe_plot()
end


