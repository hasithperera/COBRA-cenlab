%s1ths Simulates experimental data and 
%optimizes the parameters to fit the experimental results.


%AHE modification
%210621 modified the script to do the simulation without user inputs for
%education purposes and to see how each paramtr in the s3tparam file
%changes the CTR


clear all
lnt= '0';%input('lognum default-0   ','s');

if lnt=='0'
    load lognum.txt
    lg=lognum+1;
    
else
    lgr=str2num(lnt);
    slg=num2str(lgr(1));
    eval(['!copy /A /Y s3tatbl' slg '.m s3tatbl.m '])
    eval(['!copy /A /Y s3tparam' slg '.m s3tparam.m '])
    eval(['!copy /A /Y s3param' slg '.m s3param.m '])
    eval(['load vpr' slg '.dat'])
    eval(['save vpr.dat vpr' slg ' -ascii'])
    try
        eval(['load ndata' slg '_1.dat'])
        eval(['save ndata.dat ndata' slg '_1 -ascii'])
        eval(['load otbl' slg '_1.dat'])
        eval(['save otbl.dat otbl' slg '_1 -ascii'])
    catch
    end
end
s3tatbl
save par
s3param
load par

global mask0 fn0 dat2

load qz.dat
nqz=round(skzr*snk+1);

load data2.dat
[nl nc]=size(data2);
data2(1:nl,nc+1:srod)=0; %
tsp=round((nl-1)/(skzr*snk));
kk=0;
for ii=1:tsp:nl
    kk=kk+1;
    dat2(kk,:)=data2(ii,:);
    qzp(kk)=qz(ii);
end
qzp=qzp.';
[nqzp dum]=size(qzp);
fn0=0;
im=sqrt(-1);


flg0= '0'; %input('use stored fit parameters, yes-1, no - 0')
if flg0==1
    load vpr.dat
else
    vpr=svpr;
end

flms= '1';%input('simulation - 1; real data - 0    ');
mask0=ones(nqzp,srod);
if flms==0
    mask0=zeros(nqzp,srod);
    out=tbl1;
    for ii=1:nc
        sti=ceil(out(ii,4)/tsp);
        eni=floor(out(ii,5)/tsp);
        mask0(sti:eni,ii)=1;
    end
end

nbrz=floor(skzr);
ners=round(ner/tsp);
mask0(1:ners,:)=0;
for ii=1:nbrz
   mask0(ii*snk-ners+1:ii*snk+ners+1,:)=0;
end

options=optimset('MaxFunEvals',100);

[mfv nfv]=size(sfv);


flg=0;
clear xx;
for ii=1:nfv
   if sfv(ii)==1
      if flg==0
         xx0=vpr(ii);
         flg=1;
      else
         xx0=[xx0 vpr(ii)];
      end
   end
end

fn0=0;

flg1= '0' %input('perform fit, yes - 1, no - 0')
if flg1==1   
   xx=fminsearch('s3tssq',xx0,options,vpr,srul);
else
   xx=xx0;
end



sv=s3tscfactor(xx,vpr,srul);
svt(:,:)=sv(:,:,1)+sv(:,:,2);
siv(:,:)=sv(:,:,2);

kk=0;
for uu=0:sn0
    for vv=0:uu
        kk=kk+1;
        
        %removed the figure from plotting - ahe
        figure(kk)
        plot(qzp,log10(dat2(:,kk)),'k',qzp,log10(abs(svt(:,kk)).^2),'r');
    end
end

jj=1;
for ii=1:nfv
   if sfv(ii)==0
      xt(ii)=vpr(ii);
   else
      xt(ii)=xx(jj);
      jj=jj+1;
   end
end
xt
ss=s3tssq(xx,vpr,srul)
sm=sum(sum(mask0).');
r1=sqrt(ss/sm)/sum(sum(sqrt(dat2).*mask0/sm).')

svtr=real(svt(:,1:nc));
svti=imag(svt(:,1:nc));
sivr=real(siv(:,1:nc));
sivi=imag(siv(:,1:nc));
mask00=mask0(:,1:nc);
save mask0.dat mask00 -ascii
save vpr.dat xt -ascii
save svtro.dat svtr -ascii
save svtio.dat svti -ascii
save sivro.dat sivr -ascii
save sivio.dat sivi -ascii
save svtrorg.dat svtr -ascii %save for safe keeping.
save svtiorg.dat svti -ascii

flg0='0';% input('record log yes-1 no-0');

if flg0==1
    itr=0;
    load lognum.txt
    lg=lognum+1;

    slg=num2str(lg);
    save logitr.txt itr -ascii
    save lognum.txt lg -ascii
    eval(['save vpr' slg '.dat xt -ascii'])    
    eval(['!copy /A /Y s3tatbl.m s3tatbl' slg '.m '])
    eval(['!copy /A /Y s3tparam.m s3tparam' slg '.m '])
    eval(['!copy /A /Y s3param.m s3param' slg '.m '])
    
    fid=fopen('log.txt', 'a');
    fprintf(fid,'%s\n',' ');
    fprintf(fid,'%s',[slg '   ']);
    clk=clock;
    clks=datestr(clk,31);
    fprintf(fid,'%s',clks);
    fprintf(fid,'%s\n','  s3tths');
    
    fprintf(fid,'%s',['sfv=']);
    fprintf(fid,'%2d',sfv);
    fprintf(fid,'%s\n',' ');
    
    fprintf(fid,'%s',['vpr=']);
    fprintf(fid,'%8.3e\t',xt);
    fprintf(fid,'%s\n',' ');
    
    fprintf(fid,'%s',['r1=']);
    fprintf(fid,'%8.3e\n',r1);
    
    fclose(fid);
end
    ahe_1 = qzp;
    ahe_2 = log10(abs(svt(:,kk)).^2);
return 



