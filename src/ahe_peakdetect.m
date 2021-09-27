function [CoM,CoArea] = ahe_peakdetect(zstr,Peakedge,step_size,get_plot)
%AHE_PEAKDETECT Summary of this function goes here
%   This extracts the peaks using the implementation from Dr Zhou,
%   Refer for detailed help

%% CoArea

%step_size=0.14462963; %0.180787; % just for example. You need to adjust that based on the true step size from your CTR and cobra analysis.
step_size_density=10;
CoArea=[];
Z_original=1*step_size:step_size:length(zstr)*step_size;
Z_original=Z_original';
Z_dense=1*step_size/step_size_density:step_size/step_size_density:length(zstr)*step_size;
Z_dense=Z_dense';
zstr_dense=interp1(Z_original,zstr,Z_dense,'PCHIP');

for i=1:length(Peakedge)
    for j=1:length(Z_dense)
        if abs(Peakedge(i,1)*step_size-Z_dense(j))<step_size/step_size_density/2
            Peakedge(i,2)=j;
        end
    end
end


for jj=1:2:length(Peakedge)-1
    ED=zstr_dense(Peakedge(jj,2):Peakedge(jj+1,2));
    
    y=ED';
    x=Z_dense(Peakedge(jj,2):Peakedge(jj+1,2));
    Total_area(jj)=sum(y);
    
    for kk=Peakedge(jj,2)+1:Peakedge(jj+1,2)
        Area_integral(kk-Peakedge(jj,2))=sum(zstr_dense(Peakedge(jj,2):kk));
        Difference(kk-Peakedge(jj,2))=abs(Area_integral(kk-Peakedge(jj,2))-Total_area(jj)/2);
    end
    
    [mindiff ind]=min(Difference);
    CoArea(jj)=Z_dense(Peakedge(jj,2)+ind);
end
% CoArea=CoArea';


%% CoM

CoM=[];
Z_original=1*step_size:step_size:length(zstr)*step_size;
Z_triple=1*step_size/10:step_size/10:length(zstr)*step_size;

zstr_triple=interp1(Z_original,zstr,Z_triple,'PCHIP');

for i=1:length(Peakedge)
    for j=1:length(Z_triple)
        if abs(Peakedge(i,1)*step_size-Z_triple(j))<step_size/10/2
            Peakedge(i,2)=j;
        end
    end
end

%ahe - plotting
if get_plot==1
    plot(Z_triple,zstr_triple)
    hold on
end


for jj=1:2:length(Peakedge)-1
    ED=zstr_triple(Peakedge(jj,2):Peakedge(jj+1,2));
    pk_range = Peakedge(jj,2):Peakedge(jj+1,2);
    
    %ahe - plotting
    if get_plot==1
        plot(Z_triple(pk_range),zstr_triple(pk_range))
    end
    
    y=ED';
    x=Z_triple(Peakedge(jj,2):Peakedge(jj+1,2));
    Top=0;
    Bottom=0;
    for n=1:length(x)
        Top=Top+(x(n)*y(n));
        Bottom=Bottom+y(n);
    end
    CoM(jj)=Top/Bottom;
    % CoM=CoM;
end

end

