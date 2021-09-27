function write_root_filename = ahe_fun(ahe_scan_num,ahe_exp_tag,spec_file,shift,sig,window,debug)
%AHE_FUN Summary of this function goes here

%   Detailed explanation goes here
% runs Pilatus data file integration by H.Z. 03/01/2011
% Updated correction files from 33-BM February 2011

% ***setup choices*****
    PLOT = 1;  %0=no plots, 1 = all plots, 2 = 1D plots only
    auto_fit = 2; %0=manual, 1&6=automatic integration (no user control), using
                                %  1(both linefits),  2(only line along x),  3(only line along y)')
                                %  4(both curvefits), 5(only curve along x), 6(only curve along y)')
        
    poly_fit_order = [6 6]; % <=6, orders of polynomial fits along x and y, respectively.
    poly_fit_error_message_showing_option = 0; % 1: show message
    
    c_lat = 3.905; %c0 lattice constant of the substrate
    
    log_out = 1;  % 1 = log scale, 0 = linear
    ERROR_MINIMUM = 0; %1= enforces minimum error, 0 = purely statistical error
        error_min = 0.03;  %set fractional minimum error
    BCK_SUB= 1; %0=no corrections; 1=do corrections (dark/efficiency)
    SB_COR = 0; %0=no scanbar correction; 1=Correct time using scanbar
        SB_PERSEC = 514590;%111191;%
    E_ph = 20;  % nominal photon energy in KeV 
    %   (the program will use this ONLY if the photon Energy is NOT in the specfile)
    
    zinger_opt = 1;  % 1= zinger correction/0 = no zinger correction

    Data_type=1; %1:CTR, 2: RAXR, 3:time series

    RAXR_Escan_direction = 0; % 0: high to low E , 1: low to high E
    
    otheroption={'33ID'};   %for filter correction 33ID,  6ID14keV  6ID17keV  6ID13p05keV 6ID16keV 
    
    overlap_plot={''};
    
%**Data format******
    File_type = 2;   % 1=SPE files for PI-Roper CCD, 2 = Tiff files for Pilatus
    
    in_situ_plot_option = 1;
    Qmax_ctr = 7;  % for plotting A^-1
    
% ***** file information********
    %root_filename = 'sb2Te3_20QL_capped_1ML_FeSe_ISTO_04052021_1';%enter filename root here % SPEC file file name
    root_filename = spec_file;
    scan_number = [ahe_scan_num];
    TUNE = 0; % 0=integrate data, 1=show tune parameters for integration areas
    ahe_root = split(root_filename,"_");
    write_root_filename = strcat('../Data/',strjoin({ahe_root{1:end-2},num2str(scan_number),ahe_exp_tag},"_"),'L');%'hem001_rosso1_diw_L760a_';%enter output filename root her
 
    % ****Integration window controls********   
    xbar_ccd = 237+shift(1); %   center of signal region
    ybar_ccd = 86+shift(2);%86; %110; % 168;%

    delx_signal = sig(1);%8; %325;%50;%half width of signal region
    dely_signal = sig(2); %windows y 20; 30;
    
    % Caution: always keep the background window larger than the signal window!
    delx_back = window(1);%345; %half width of background region
    dely_back = window(2);% 

%default num
%     delx_back = 24;%345; %half width of background region
%     dely_back = 30;% 
   
    ymin=1;  % Slit window ROI range top Y
    ymax=195; % Slit window ROI range bottom Y
    xmin=1; % Slit window ROI range left X
    xmax=487; % S0lit window ROI range right X
   
    slit_ROI = [ymin ymax xmin xmax];
    
    outputfile1 = [write_root_filename,'_CTR_R.ipg']; %enter output filename1 suffix here    
    % file_number(jjjj),LL,epoch,energy,monitor,Signal_0d,mse_sSignal_0d,Signal,sigSignal,Signaly,sigSignaly,Background,ccd_num_back_ccd_px,time,time_correct,xbar,ybar,HH,KK,LL,theta,ttheta,chi,phi,filter,alfa,beta,azimuth
    outputfile2 = [write_root_filename, '_CTR_P.ipgx']; %enter output filename2 suffix here
    %   file_number(jjjj),LL,epoch,energy,monitor,Signal_0d,sSignal_0d,Background,Backgroundy,fwhmx,fwhmy,chi2_back,chi2_backy,regions;

    if (debug==1)
        disp(outputfile1)
       return; 
    end
    
    if (TUNE == 0 )
        log_ahe = fopen("../Data/ahe_log.csv",'a+');
        fprintf(log_ahe,"%s,%s,%d,%d,%d,%d,%d,%d,%d\n",write_root_filename,ahe_exp_tag,scan_number,xbar_ccd,ybar_ccd,delx_signal,dely_signal, delx_back,dely_back);
        fclose(log_ahe);
    end
    
    if min(poly_fit_order) <=0 || max(poly_fit_order)>6 || sum(poly_fit_order~=round(poly_fit_order))
        error(' Adjust the ''poly_fit_order'' parameter!! ')
    end

% ****Image correction, including dark counts, Pilatus efficiency correction
    global pilatus_dark pilatus_efficiency

if BCK_SUB == 1
    
   % dark_image_number = [1:5];
   % Dark_count_time = 100; %dark counts counting time
    
   % for i=1:length(dark_image_number)
   % dark_image=['S6_pilatus_test_0708.000.2061.00194_0000' num2str(dark_image_number(i)) '.tif']; % background dark counts image
   % dark_image_array(:,:,i) = imread(dark_image);
   % end 
   
   % pilatus_dark = sum(sum(sum(double(dark_image_array))))/(size(dark_image_array,1)*size(dark_image_array,2)*size(dark_image_array,3))/Dark_count_time; %pilatus dark count (cts/pixel/s)
     pilatus_dark = 0;
    % flatfield_image_number = [67:86]; % first image to last image for flat field images
     efficiency_count_time = 20; %flat field image counting time
    
%     for i=1:length(flatfield_image_number)
%     flatfield_image=['S6_pilatus_test_0708.000.2061.15140_000' num2str(flatfield_image_number(i)) '.tif'];
%     pilatus_efficiency_array(:,:,i) = imread(flatfield_image);
%     end 
%     pilatus_efficiency = double(pilatus_efficiency_array(ymin:ymax,xmin:xmax,:))-efficiency_count_time*pilatus_dark;
%     pilatus_efficiency = mean(pilatus_efficiency,3)/mean(mean(mean(pilatus_efficiency,3)));
%     pilatus_efficiency = ones(ymax-ymin+1,xmax-xmin+1);

     flatfield_image=['flat_field_correction.tif'];
     pilatus_efficiency_array(:,:,1) = imread(flatfield_image);
     pilatus_efficiency = double(pilatus_efficiency_array(ymin:ymax,xmin:xmax,:))-efficiency_count_time*pilatus_dark+0.1;
     pilatus_efficiency = mean(pilatus_efficiency,3)/mean(mean(mean(pilatus_efficiency,3)));
    
% run the integration function   

% try
    Pilatus_fit_function10_ZM_12IDD(root_filename,outputfile1,outputfile2...
    ,xbar_ccd,ybar_ccd,delx_signal,dely_signal,delx_back,dely_back...
    ,PLOT,auto_fit,TUNE,log_out,ERROR_MINIMUM,error_min,File_type,E_ph,BCK_SUB,SB_COR,SB_PERSEC,zinger_opt...
    ,scan_number,Data_type, in_situ_plot_option, Qmax_ctr, poly_fit_order, poly_fit_error_message_showing_option...
    ,otheroption,overlap_plot,RAXR_Escan_direction,c_lat,slit_ROI);
% catch
%     fprintf(2,"User terminated due to error\n")
%     return
% end

end

