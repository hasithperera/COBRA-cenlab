
%This script combines all the data from
%different peak detection runs and construct all peak locations
%Using a Unitcell definition for plotting the location
%and the error bar is based on the last few substrate layers

%experimental conditions are added manually from the notebook


%film
data_src = '../Data/peak_4/'

load(strcat(data_src,'ahe_CoA.dat'))
substrate_CA = zeros([size(ahe_CoA,1),3]);
film_peaks = ahe_CoA(:,3:2:end);

%STO substrate
data_src = '../Data/peak_5/'

load(strcat(data_src,'ahe_CoA.dat'))
substrate_CA = zeros([size(ahe_CoA,1),3]);
STO_peaks = ahe_CoA(:,3:2:end);

loc_error = 0.0145;


full_sample = [STO_peaks(:,1:end-1) film_peaks];


%% experimental conditions

%Sampeid
%scanid
%Temp
%light : 0,1 - UV ,2 - White 3- white cts
%se capped
exp_conditions = [  2,9,300,0; %1
    2,149,10,0; %2
    2,176,300,2;%3
    2,182,10,0; %4
    2,184,10,2; %5
    2,186,10,3  %6
    2,274,19,1  %7
    2,352,22.5,1%8
    2,358,22.5,1%9
    1,7,300,0;  %10
    1,46,10,0;  %11
    1,101,30,1; %12
    1,104,80,1; %13
    ]
exp_id = 10
lbl = {} 
h = []
k = 1


for exp_id = [4:8]
    unit_cell = [-8:0,1:2];
    peak_distance = diff(full_sample(exp_id,[1:2:17,18:20]));
    h(k) = errorbar(unit_cell,peak_distance,loc_error*ones(1,length(unit_cell)),'.:')
    
    hold on
    lbl{k} = sprintf('Temp:%dK, light: %d sample:%d',exp_conditions(exp_id,3),exp_conditions(exp_id,4),exp_conditions(exp_id,1))
    k = k +1;
%     if k == 1
%         
%     hold on
%     end
end

xlabel('Z (peak)')


ahe_plot()
legend(h,lbl)

% 
