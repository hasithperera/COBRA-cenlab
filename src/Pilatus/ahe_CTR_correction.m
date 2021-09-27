function ahe_CTR_correction(file_name)
%AHE_CTR_CORRECTION Summary of this function goes here
%   Detailed explanation goes here

%% Plotdata_Pilatus
%  Pilatus Image Data Plot Program
%  04/12/2011 H.Z.
%
% general inputs
alat =  3.905;
dplot = 1;	% 0 = linear; 1=semilogy; 2=normI (semilogy)
x_opt = 2; % 1 = data vs. Q(A-1), 2 = data vs. L
filter_corr = 1;% 0 = no filter correction; 1 = filter
% data specific inputs
indir='';%'C:\Hua Zhou Project at ANL and BNL\Paul Fenter CTR programs\1D COBRA programs\Eom's Data Feb. 2011 and Analysis\Eom_02162011\Feb11\Eom_pilatus\';
file1_root = [strcat(file_name,'_CTR_R')] ;
file2_root = ['sb2Te3_20QL_capped_1ML_FeSe_ISTO_36_10_aL_CTR_R'];
file3_root = ['QCl_Thin_CTR_10112018_5uA_2min_CTR_R'];
file4_root = ['QCl_Thin_CTR_10112018_10uA_9min_CTR_R'];


file1 = [file1_root,'.ipg'] ;
file2 = [file2_root,'.ipg'];
file3 = [file3_root,'.ipg'];
file4 = [file4_root,'.ipg'];

% file2 = 'STO_SIO_STO_symmetry_m1m2L_CTR_R.ipg';
% file3 = 'STO_SIO_STO_symmetry_21L_CTR_R.ipg';
% file4 = 'LCO214_032015B_8UC_00L_115CTR_R.ipg';

file1=[indir,file1];file2=[indir,file2];file3=[indir,file3];file4=[indir,file4];
% plot title
plot_title = [];
% plot option
nplots = [1 0 0 0 ];     % which data to plot
fileopt = [1 1 1 1 ];    % 1 = plot I/Mon with EBs; 2= plot data w/o norm & w/o EBs; 3=plot I/Mon w/o EBs
%
%   1             2           3            4            5           6            7                8          9                 10          11              12          13             14           15                   16    17                  18            19                  20          21        22        23          24           25        26          27
%   file_number   L           epoch        Energy       Monitor  Signal_Best  sigSignal_Best   Signal_0d  mse_sSignal_0d    Signal      sigSignal       Signaly     sigSignaly     Background   ccd_num_back_ccd_px  time  time_corr_factor    xbar          ybar                HH          KK        LL         fwhm_x      fwhm_y        atten   accept_option  specscan#\n
%xcol = 2;
xcol = 2;
ycol = 6;
sycol= 7;
mcol = 5;
filterline=25;
trans=30;
energycol = 4;      % if filter_corr ==1
%  **********
dataopt = [1 1 1 1];    % Integrated area from: 1 = peakfit, 2 = numerical, 3 = numerical&Deadtime
xscale = [1.00 1.00 1.0000 1.00];     % xscale for each data set
xoffset = 1*[0.0 0.0 0.0 0.0];
yscale = [1 1 1 1];   % yscale for each data set
yoffset = [0 0 0 0];
sym1 = 'b.-';           % symbol for each data set
sym2 = 'rs-';
sym3 = 'k<-';
sym4 = 'gd-';
mica_fudge = 6; % =1 for most crystal, = 2 for mica-like ctrs (with only even order reflections)

%  *******************************
if dplot ==2
    x_opt =1;
end
if x_opt == 1
    xscale =xscale*2*pi/alat;
end
% *******************************
%  axes labels
%plot_title = 'tio2: 0.001MRb ';
y_txt = 'Intensity';%
set_axis_limits = 0;
if set_axis_limits ==1
    axis_length = [0 6 1e-5 3];
    %    axis_length = [11 12 1e-9 .01];
end
%
if or(x_opt == 1,dplot==2)
    x_txt = 'Q (A-1)';
else
    x_txt = 'L (rlu)'; %%'Azimuthal Angle';    % %% ' xsamx (mm)';%'Epoch'; 'scan #';%
end
if dplot ==2
    y_txt = 'normIntensity';
end

% **********************************************************
%   load data files
%  *********************************************************
if nplots(1) == 1
    t1 = load(file1);	%red triangle (right facing)
    t1 = sortrows(t1,xcol);%
end
if nplots(2) == 1
    t2 = load(file2);	%red triangle (right facing)
    t2  = sortrows(t2,xcol);%
end
if nplots(3) == 1
    t3 = load(file3);	%blue triangles (left facing)
    t3 = sortrows(t3,xcol);%cal2;
end
if nplots(4) == 1
    t4 = load(file4);
    t4 = sortrows(t4,xcol);
end

ATTEN=1;
% *************************************************
if filter_corr == 1
    FILTER = t1(:,filterline);
    ATTEN = ones(size(FILTER));
    
    for jj=1:length(FILTER)
        
        
        %             if FILTER(jj) == 0
        %                 ATTEN(jj) = 1;
        %             elseif FILTER(jj) == 1
        %                 ATTEN(jj) = 1/3.62;
        %             elseif FILTER(jj) == 2
        %                 ATTEN(jj) = (1/3.62)^2;
        %             elseif FILTER(jj) == 3
        %                 ATTEN(jj) = (1/3.62)^3;
        %             elseif FILTER(jj) == 4
        %                 ATTEN(jj) = (1/3.62)^4;
        %             elseif FILTER(jj) == 5
        %                 ATTEN(jj) = (1/3.62)^5;
        %             elseif FILTER(jj) == 6
        %                 ATTEN(jj) = (1/3.62)^6;
        %             elseif FILTER(jj) == 7
        %                 ATTEN(jj) = (1/3.62)^7;
        %             elseif FILTER(jj) == 8
        %                 ATTEN(jj) = (1/3.62)^8;
        %             elseif FILTER(jj) == 9
        %                 ATTEN(jj) = (1/3.62)^8*(1/3.62);
        %             elseif FILTER(jj) == 10
        %                 ATTEN(jj) = (1/3.62)^8*(1/3.62)^2;
        %             elseif FILTER(jj) == 11
        %                 ATTEN(jj) = (1/3.62)^8*(1/3.62)^3;
        %             elseif FILTER(jj) == 12
        %                 ATTEN(jj) = (1/3.62)^8*(1/3.62)^4;
        %             elseif FILTER(jj) == 13
        %                 ATTEN(jj) = (1/3.62)^8*(1/3.62)^5;
        %             end
        %
        %      absorptivity = 10.13;
        absorptivity = 4;
        % absorptivity = 19.15;
        % absorptivity = 2.5;
        if FILTER(jj) == 0
            ATTEN(jj) = 1;
        elseif FILTER(jj) == 1
            ATTEN(jj) = 1/2.25;
            ATTEN(jj) = 1/absorptivity;
        elseif FILTER(jj) == 2
            ATTEN(jj) = (1/2.3)^2;
            ATTEN(jj) = (1/absorptivity)^2;
        elseif FILTER(jj) == 3
            ATTEN(jj) = (1/absorptivity)^2*(1/2.33);
            ATTEN(jj) = (1/absorptivity)^3;
        elseif FILTER(jj) == 4
            ATTEN(jj) = (1/2.65)^4;
            ATTEN(jj) = (1/absorptivity)^4;
        elseif FILTER(jj) == 5
            ATTEN(jj) = (1/2.59)^5;
            ATTEN(jj) = (1/absorptivity)^5;
        elseif FILTER(jj) == 6
            ATTEN(jj) = (1/2.62)^6;
            ATTEN(jj) = (1/absorptivity)^6;
        elseif FILTER(jj) == 7
            ATTEN(jj) = (1/2.589)^7;
            ATTEN(jj) = (1/absorptivity)^7;
        elseif FILTER(jj) == 8
            ATTEN(jj) = (1/2.552)^8;
            ATTEN(jj) = (1/absorptivity)^8;
        elseif FILTER(jj) == 9
            ATTEN(jj) = (1/absorptivity)^8*(1/2.8);
            ATTEN(jj) = (1/absorptivity)^8*(1/absorptivity);
        elseif FILTER(jj) == 10
            ATTEN(jj) = (1/absorptivity)^8*(1/2.8)^2;
            ATTEN(jj) = (1/absorptivity)^8*(1/absorptivity)^2;
        elseif FILTER(jj) == 11
            ATTEN(jj) = (1/absorptivity)^8*(1/2.6)^3;
            ATTEN(jj) = (1/absorptivity)^8*(1/absorptivity)^3;
        elseif FILTER(jj) == 12
            ATTEN(jj) = (1/absorptivity)^8*(1/2.79)^4;
            ATTEN(jj) = (1/absorptivity)^8*(1/absorptivity)^4;
        elseif FILTER(jj) == 13
            ATTEN(jj) = (1/absorptivity)^8*(1/absorptivity)^5;
        elseif FILTER(jj) == 15
            ATTEN(jj) = (1/absorptivity)^8*(1/absorptivity)^7;
        end
        
    end
    
end
%  *******************************

if nplots(1) == 1
    if fileopt(1) ==1    % use XXXR files
        q1  =  xscale(1)*(t1(:,xcol)-xoffset(1));
        R1 =   yscale(1)*t1(:,ycol)./t1(:,mcol) + yoffset(1);
        s1  =  yscale(1)*t1(:,sycol)./t1(:,mcol);%0.02*R1;%
        nR1  = R1.*(q1.*sin(q1*alat/2/mica_fudge)).^2;  % normalized reflectivity
        ns1  = s1.*(q1.*sin(q1*alat/2/mica_fudge)).^2;
    elseif fileopt(1) ==2  % use XXXP files
        q1  =  xscale(1)*(t1(:,xcol)-xoffset(1));
        R1 =   yscale(1)*t1(:,ycol) + yoffset(1);
        nR1  = R1.*(q1.*sin(q1*alat/2/mica_fudge)).^2;  % normalized reflectivity
    elseif fileopt(1) ==3    % use XXXR files
        q1  =  xscale(1)*(t1(:,xcol)-xoffset(1));
        R1 =   yscale(1)*t1(:,ycol)./t1(:,mcol) + yoffset(1);
        nR1  = R1.*(q1.*sin(q1*alat/2/mica_fudge)).^2;  % normalized reflectivity
    end
    
    %  if filter_corr == 1
    %      FILTER = t1(:,25);
    %      ATTEN=FILTER;
    %      R1=R1./ATTEN;
    %      nR1=nR1./ATTEN;
    %  end
    if filter_corr == 1
        ATTEN = t1(:,trans);
    end
    % FILTER = t1(:,27);
    % ATTEN=FILTER;
    R1=R1./ATTEN;
    s1=s1./ATTEN;
    nR1=nR1./ATTEN;
    ns1=ns1./ATTEN;
end

% filter correction
% for jjj = 1:length(t1(:,1))
%     if t1(jjj,filterline) > 0
%         if filter_corr == 1
%             filternum = t1(jjj,filterline);
%             filtertemp = dec2bin(filternum);
%             filterndx = [0 0 0 0];
%             for kkk=1:length(filtertemp)
%                 filterndx(kkk) = str2num(filtertemp(length(filtertemp)-kkk+1));
%             end
%             t_filter1 = sum(filter_thick(:).*filterndx(:));
%             f_cor1 = exp(t_filter1/LAM(1));
%         end
%         R1(jjj) = R1(jjj)*f_cor1;
%         nR1(jjj) = nR1(jjj)*f_cor1;
%         if fileopt(1) ==1
%             s1(jjj) = s1(jjj)*f_cor1;
%             ns1(jjj) = ns1(jjj)*f_cor1;
%         end
%     end
% end
%
if nplots(2) == 1
    if fileopt(2) ==1
        q2  =  xscale(2)*(t2(:,xcol)-xoffset(2));
        R2 =   yscale(2)*t2(:,ycol)./t2(:,mcol) + yoffset(2);
        s2  =  yscale(2)*t2(:,sycol)./t2(:,mcol);
        nR2  = R2.*(q2.*sin(q2*alat/2/mica_fudge)).^2;%/10;  % normalized reflectivity
        ns2  = s2.*(q2.*sin(q2*alat/2/mica_fudge)).^2;%/10;
    elseif fileopt(2) ==2
        q2  =  xscale(2)*(t2(:,xcol)-xoffset(2));
        R2 =   yscale(2)*t2(:,ycol) + yoffset(2);
        nR2  = R2.*(q2.*sin(q2*alat/2/mica_fudge)).^2;
    elseif fileopt(2) ==3
        q2  =  xscale(2)*(t2(:,xcol)-xoffset(2));
        R2 =   yscale(2)*t2(:,ycol)./t2(:,mcol) + yoffset(2);
        nR2  = R2.*(q2.*sin(q2*alat/2/mica_fudge)).^2;
    end
    if filter_corr == 1
        ATTEN = t2(:,trans);
    end
    R2=R2./ATTEN;
    s2=s2./ATTEN;
    nR2=nR2./ATTEN;
    ns2=ns2./ATTEN;
end


if nplots(3) == 1
    if fileopt(3) ==1
        q3  =  xscale(3)*(t3(:,xcol)-xoffset(3));
        R3 =   yscale(3)*t3(:,ycol)./t3(:,mcol) + yoffset(3);
        s3  =  yscale(3)*t3(:,sycol)./t3(:,mcol);
        nR3  = R3.*(q3.*sin(q3*alat/2/mica_fudge)).^2;
        ns3  = s3.*(q3.*sin(q3*alat/2/mica_fudge)).^2;
    elseif fileopt(3) ==2
        q3  =  xscale(3)*(t3(:,xcol)-xoffset(3));
        R3 =   yscale(3)*t3(:,ycol) + yoffset(3);
        nR3  = R3.*(q3.*sin(q3*alat/2/mica_fudge)).^2;
    elseif fileopt(3) ==3
        q3  =  xscale(3)*(t3(:,xcol)-xoffset(3));
        R3 =   yscale(3)*t3(:,ycol)./t3(:,mcol) + yoffset(3);
        nR3  = R3.*(q3.*sin(q3*alat/2/mica_fudge)).^2;
    end
    
    if filter_corr == 1
        ATTEN = t3(:,trans);
    end
    R3=R3./ATTEN;
    s3=s3./ATTEN;
    nR3=nR3./ATTEN;
    ns3=ns3./ATTEN;
end

if nplots(4) == 1
    if fileopt(4) == 1
        q4  =  xscale(4)*(t4(:,xcol)-xoffset(4));
        R4 =   yscale(4)*t4(:,ycol)./t4(:,mcol) + yoffset(4);
        s4  =  yscale(4)*t4(:,sycol)./t4(:,mcol);
        nR4  = R4.*(q4.*sin(q4*alat/2/mica_fudge)).^2;
        ns4  = s4.*(q4.*sin(q4*alat/2/mica_fudge)).^2;
        for jjj = 1:length(t4(:,1))
            if t4(jjj,filterline) == 1
                R4(jjj) = R4(jjj);
                s4(jjj) = s4(jjj);
                nR4(jjj) = nR4(jjj);
                ns4(jjj) = ns4(jjj);
            end
            if t4(jjj,filterline) == 2
                R4(jjj) = R4(jjj);
                s4(jjj) = s4(jjj);
                nR4(jjj) = nR4(jjj);
                ns4(jjj) = ns4(jjj);
            end
        end
    elseif fileopt(4) ==2
        q4  =  xscale(4)*(t4(:,xcol)-xoffset(4));
        R4 =   yscale(4)*t4(:,ycol) + yoffset(4);
        nR4  = R4.*(q4.*sin(q4*alat/2/mica_fudge)).^2;
    elseif fileopt(4) == 3
        q4  =  xscale(4)*(t4(:,xcol)-xoffset(4));
        R4 =   yscale(4)*t4(:,ycol)./t4(:,mcol) + yoffset(4);
        nR4  = R4.*(q4.*sin(q4*alat/2/mica_fudge)).^2;
    end
    
    if filter_corr == 1
        ATTEN = t4(:,trans);
    end
    R4=R4./ATTEN;
    s4=s4./ATTEN;
    nR4=nR4./ATTEN;
    ns4=ns4./ATTEN;
end

%
figure(3);
if   dplot == 0  % linear plots
    if nplots(1) == 1
        hold off
        plot(q1,R1,sym1)
        if fileopt(1) ==1
            hold on
            errorbar(q1,R1,s1,sym1)
            hold off
        end
    end
    if nplots(2) == 1
        if nplots(1) ==0
            hold off
        else
            hold on
        end
        plot(q2,R2,sym2)
        if fileopt(2) ==1
            hold on
            errorbar(q2,R2,s2,sym2)
            hold off
        end
    end
    if nplots(3) == 1
        hold on
        plot(q3,R3,sym3)
        if fileopt(3) ==1
            hold on
            errorbar(q3,R3,s3,sym3)
            hold off
        end
    end
    if nplots(4) == 1
        hold on
        plot(q4,R4,sym4)
        if fileopt(4) ==1
            hold on
            errorbar(q4,R4,s4,sym4)
            hold off
        end
    end
    % *****
elseif   dplot == 1  % semilogy plots
    if nplots(1) == 1
        hold off
        semilogy(q1,R1,sym1)
        if fileopt(1) ==1
            hold on
            plot(q1,R1,sym1) % Guichuan has included this line to remove the error bar
            errorbar(q1,R1,s1,sym1) %comented to remove the error bar
            hold off
        end
    end
    if nplots(2) == 1
        if nplots(1) ==0
            hold off
        else
            hold on
        end
        semilogy(q2,R2,sym2)
        if fileopt(2) == 1
            hold on
            plot(q2,R2,sym2)
            errorbar(q2,R2,s2,sym2)
            hold off
        end
    end
    if nplots(3) == 1
        hold on
        semilogy(q3,R3,sym3)
        if fileopt(3) ==1
            hold on
            plot(q3,R3,sym3)
            errorbar(q3,R3,s3,sym3)
            hold off
        end
    end
    if nplots(4) == 1
        hold on
        semilogy(q4,R4,sym4)
        if fileopt(4) ==1
            hold on
            plot(q4,R4,sym4)
            errorbar(q4,R4,s4,sym4)
            hold off
        end
    end
    
elseif dplot ==2  % semilogy of normalized data
    if nplots(1) == 1
        hold  off
        semilogy(q1,nR1,sym1)
        if fileopt(1) ==1
            hold on
            errorbar(q1,nR1,sym1)
            hold off
        end
    end
    if nplots(2) == 1
        if nplots(1) ==0
            hold off
        else
            hold on
        end
        semilogy(q2,nR2,sym2)
        if fileopt(2) ==1
            hold on
            errorbar(q2,nR2,ns2,sym2)
            hold off
        end
    end
    if nplots(3) == 1
        hold on
        semilogy(q3,nR3,sym3)
        if fileopt(3) ==1
            hold on
            errorbar(q3,nR3,ns3,sym3)
            hold off
        end
    end
    if nplots(4) == 1
        hold on
        semilogy(q4,nR4,sym4)
        if fileopt(4) ==1
            hold on
            errorbar(q4,nR4,ns4,sym4)
            hold off
        end
    end
end
% *************************
%  plot titles
title(plot_title);
if set_axis_limits ==1
    AXIS(axis_length)
end

xlabel(x_txt);
ylabel(y_txt);
% ************************
%axis tight;

%% correction
Output_datafile = file1(1:end-4);

if nplots(1) == 1
    Data1(:,1)=q1;
    Data1(:,2)=R1;
    Data1(:,3)=s1;
    Data1(:,4)=t1(:,26);
    Data1(:,5)=t1(:,27);
    Data1(:,6)=t1(:,4);
end
if nplots(2) == 1
    Data2(:,1)=q2;
    Data2(:,2)=R2;
    Data2(:,3)=s2;
    Data2(:,4)=t2(:,26);
    Data2(:,5)=t2(:,27);
    Data2(:,6)=t2(:,4);
end
if nplots(3) == 1
    Data3(:,1)=q3;
    Data3(:,2)=R3;
    Data3(:,3)=s3;
    Data3(:,4)=t3(:,26);
    Data3(:,5)=t3(:,27);
    Data3(:,6)=t3(:,4);
end
if nplots(4) == 1
    Data4(:,1)=q4;
    Data4(:,2)=R4;
    Data4(:,3)=s4;
    Data4(:,4)=t4(:,26);
    Data4(:,5)=t4(:,27);
    Data4(:,6)=t4(:,4);
end
Clat = 3.905;
Max_fp = 0;

%good for all our smaples
sample_size = 3; % mm in unit based on sample size
Vertical_beam_size = 100; % micron meter in unit
Alaph = 5; % angle of incidence for non-specualr CTRs

CTR_corr = 1; % 0 = no CTR correction; 1 = correction for 00L rod;  2=correction for non-specualr rods
% the correction includes activation area for 00L rod, and intersection
% corrections for all rod.

if CTR_corr == 0
    if nplots(1) == 1
        Data_out=Data1(:,1:3);
        save([file1_root,'_Data_out_nocorrection.dat'], 'Data_out','-ASCII');
    end
    if nplots(2) == 1
        Data_out=Data2(:,1:3);
        save([file2_root,'_Data_out_nocorrection.dat'], 'Data_out','-ASCII');
    end
    if nplots(3) == 1
        Data_out=Data3(:,1:3);
        save([file3_root,'_Data_out_nocorrection.dat'], 'Data_out','-ASCII');
    end
    if nplots(4) == 1
        Data_out=Data4(:,1:3);
        save([file4_root,'_Data_out_nocorrection.dat'], 'Data_out','-ASCII');
    end
    %    save Data_out.dat Data_out -ASCII
    %
elseif CTR_corr == 1
    
    Data1(:,2)=Data1(:,2).*((12.398./Data1(:,6)).*Data1(:,1)/2/Clat).^2;  % activation area and CTR intersection corrections
    Data1(:,3)=Data1(:,3).*((12.398./Data1(:,6)).*Data1(:,1)/2/Clat).^2;  % activation area and CTR intersection corrections
    
    for i=1:size(Data1,1)
        if ((12.398/Data1(i,6))*Data1(i,1)/2/Clat) <= (Vertical_beam_size/1000/sample_size)
            Max_fp=i;
        end
    end
    if Max_fp == 0
        Data_out=Data1(:,1:3);
    else
        Data1(1:Max_fp,2)=Data1(1:Max_fp,2)./((12.398./Data1(1:Max_fp,6)).*Data1(1:Max_fp,1)/2/Clat)*(Vertical_beam_size/1000/sample_size);  % X-ray footprint corrections
        Data1(1:Max_fp,3)=Data1(1:Max_fp,3)./((12.398./Data1(1:Max_fp,6)).*Data1(1:Max_fp,1)/2/Clat)*(Vertical_beam_size/1000/sample_size);  % X-ray footprint corrections
        Data_out=Data1(:,1:3);
    end
    save([ Output_datafile,'_Data_out.dat'], 'Data_out','-ASCII');
    %    save Data_out.dat Data_out -ASCII
    
elseif CTR_corr == 2
    
    Data1(:,2)=Data1(:,2)./(1-(cos(Data1(:,5)/180*pi)).^2.*(sin(Data1(:,4)/180*pi)).^2).*cos(Data1(:,4)/180*pi).*sin((Data1(:,5)-Alaph)/180*pi);
    Data1(:,3)=Data1(:,3)./(1-(cos(Data1(:,5)/180*pi)).^2.*(sin(Data1(:,4)/180*pi)).^2).*cos(Data1(:,4)/180*pi).*sin((Data1(:,5)-Alaph)/180*pi);
    Data_out=Data1(:,1:3);
    %    save Data_out.dat Data_out -ASCII
    save([ Output_datafile,'_Data_out.dat'], 'Data_out','-ASCII');
    
end
%
% semilogy(Data_out(:,1),Data_out(:,2),'>-');
% hold on;
% semilogy(Data_out(:,1),Data_out(:,3),'>-');
% title(['Now Scan: ',Output_datafile]);
%
close all
semilogy(Data_out(:,1),Data_out(:,2),'.-b')
hold on
semilogy(Data_out(:,1),abs(Data_out(:,2)),'o-r')
fprintf("No of negatives: %d\n",sum(Data_out(:,2)<0))

end

