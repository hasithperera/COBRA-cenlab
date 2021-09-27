

load ahe_CoA.dat
load ahe_CoM.dat

substrate_CA = zeros([size(ahe_CoA,1),3]);
substrate_CM = zeros([size(ahe_CoA,1),3]);

%find sustrate
peaks = ahe_CoA;
substrate = [];
for ii = 1:3
    tmp = peaks(:,2*ii+1+4) - peaks(:,2*ii+1);
    substrate(:,ii) = tmp;
end
substrate_CA = substrate;

peaks = ahe_CoM;
substrate = [];
for ii = 1:3
    tmp = peaks(:,2*ii+1+4) - peaks(:,2*ii+1);
    substrate(:,ii) = tmp;
end
substrate_CM = substrate;


%Sampeid 
%scanid 
%Temp 
%light : 0,1 - UV ,2 - White 3- white cts
%se capped
exp_conditions = [  
                    2,9,300,0;
                    2,149,10,0;
                    2,176,300,2;
                    2,182,10,0;
                    2,184,10,2;
                    2,186,10,3
                    2,274,19,1
                    2,352,22.5,1
                    2,358,22.5,0
                    1,7,300,0;
                    1,46,10,0;
                    1,101,30,1;
                    1,104,80,1; 
                    ]
% substrate_2 = ahe_CoM(:,5) - ahe_CoM(:,3)
% substrate_2(indx)
%% STO layer distance

% STO_latice = 3.905;
% plot([1 length(ahe_CoM)],STO_latice*[1 1])
% hold on
% for k =1    
%     sample_id = 1 
%     indx = find(ahe_CoM(:,1)==sample_id);
%     plot(indx,substrate(indx,k),'or')
%     
%     sample_id = 2 
%     indx = find(ahe_CoM(:,1)==sample_id);
%     plot(indx,substrate(indx,k),'oblack')
% 
% end
% ylim([3.89,3.93])
% ylabel('STO Lattice distance (A)')
% legend('known', 'Se capped','Se Si capped')

%Se capped

% %temp 
% 
% smpl = find(exp_conditions(:,1)==1)
% plot(exp_conditions(smpl,3),substrate_CA(smpl,1),'o')
% hold on
% plot(exp_conditions(smpl,3),substrate_CM(smpl,1),'o')
% legend('COA','COM')
% xlabel('Temperature (K)')
% ylabel('Distance (A)')

% Se si capped
%temp 

smpl = find(exp_conditions(:,1)==2)
plot(exp_conditions(smpl,3),substrate_CA(smpl,1),'o')
hold on
for i = 1:length(exp_conditions(smpl,3))
    tt = strcat(num2str(exp_conditions(smpl(i),2)));
    if exp_conditions(smpl(i),4) == 0
    plot(exp_conditions(smpl(i),3),substrate_CA(smpl(i),1),'*black')
    text(exp_conditions(smpl(i),3)-5,substrate_CA(smpl(i),1),tt)
    end
end
% hold on
% plot(exp_conditions(smpl,3),substrate_CM(smpl,1),'or')
% legend('COA','COM')
xlabel('Temperature (K)')
xlim([0 310])

ylabel('Distance (A)')
ylim([3.902 3.93])
%% film distances
% 
% film = ahe_CoM(:,15:2:end);
% 
% d1 = film(:,2)-film(:,1);
% d2 = film(:,3)-film(:,2);
% 
% 
% plot(d2)
% 
% uv = find(exp_conditions(:,4));
% smpl = find(exp_conditions(:,1)==2);
% 
% 
% plot(d1,'o')
% hold on
% plot(d1(smpl),'o')
