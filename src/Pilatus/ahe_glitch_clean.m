data_path = "../Data/";
log_file = "../Data/ml_log";

files = dir(data_path+"sb2Te3_*.dat");

if(exist('j'))
    kk = 1;
end

plot_clean_only = 1
plot_tag="Se_Si";
manual_select = 1

for i = kk:length(files)
    file = data_path +files(i).name
    out_name = "cleaned_"+files(i).name
    
    if (manual_select)
        
        data = load(file)
        semilogy(data(:,1),abs(data(:,2)),".-")
        title(files(i).name)
        
        selected = 0
        
        while(selected <=0)
            fprintf("Press enter to continue after glitch selection\n")
            [x,y] = ginput();
            
            glitch_ids= zeros([1,length(x)]);
            for j = 1:length(x)
                %find closest point to the data
                [~,indx] = min(abs(data(:,1)-x(j)));
                glitch_ids(j)=indx;
            end
            disp(glitch_ids)
            hold on
            semilogy(data(glitch_ids,1),abs(data(glitch_ids,2)),"ro",'Linewidth',2)
            selected = str2num(input("confirm selection(1/0):","s"));
            
        end
        %remove from data
        fprintf("old length:%d\n",length(data))
        data(glitch_ids,:)=[];
        semilogy(data(:,1),abs(data(:,2)),"r.",'Linewidth',2)
        hold off
        fprintf("New length:%d\n",length(data))
        dlmwrite(data_path + out_name,abs(data),'delimiter',' ','precision','%.7e');
        log_glitch(file,glitch_ids)
    else
        %use preprocessed files to save plots
        if plot_clean_only
            data = load(data_path+out_name)
            data2 = load(data_path+file)
        end
        scan_data = split(out_name,"_")
        plot_name = strcat(scan_data(end-5)," ",scan_data(end-4))
        subplot(212)
        semilogy(data(:,1),abs(data(:,2)),"r.-")
        title("cleaned")
        hold on
        subplot(211)
        semilogy(data2(:,1),abs(data2(:,2)),".-b")
        title(plot_name)
        break
        saveas(gcf,strcat("../Data/plots/",plot_tag,plot_name,".png"))
        
    end
    
    
   % break; %testing single iterration
end


function log_glitch(file_id,indx)
log_file = "../Data/ml_log.csv";
fp = fopen(log_file,"a+");
fprintf(fp,"%s,",file_id)
for j=1:length(indx)
    
    fprintf(fp,"%d",indx(j))
    if j~= length(indx)
        fprintf(fp,",")
    end
end
fprintf(fp,"\n")
fclose(fp)
end

function clean_plot()

end
