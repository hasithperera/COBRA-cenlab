
clc
%fprintf('s3smo')
s3smo
clc
%fprintf('s3tths')
s3tths
clc
%fprintf('s3tbl')
s3tbl
clc
%fprintf('s3ndata1')
s3ndata1
clc
%fprintf('s3ttbf')
s3ttbf
clc
%fprintf('s3pht')
s3pht
clc
%fprintf('s3bft')
s3bft
s3tdbft

plot(bfcr)

ginput(1)
clc
%fprintf('s3tdbft')

s3tdbft

clear('bfcr1')

load('bfcr1.dat')


%final step before itr
save edorg.dat bfcr1 -ascii
s3itr_imp
