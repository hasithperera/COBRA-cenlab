
load("./2 and 3 data/ahe_CoA.dat")
load("./2 and 3 data/ahe_CoM.dat")


substrate_CA = zeros([size(ahe_CoA,1),3]);
substrate_CM = zeros([size(ahe_CoA,1),3]);

%find sustrate
peaks = ahe_CoA;
substrate = [];
for ii = 1:3
    tmp = peaks(:,2*ii+1+4) - peaks(:,2*ii+1);
    substrate(:,ii) = tmp;
  
end

% for s = 1:length(substrate)
%    plot([1 2 3],substrate(s,:),'o') 
%    hold on
% end