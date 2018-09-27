% GENVDATA Generation of simulated vowel data
%

Ntot = 20000;

Pa = 0.25;
mu_a = [730 1090];
std_a = [35  20
         20  230];
Na = Pa*Ntot;
a = randn(Na,2) * std_a + repmat(mu_a,Na,1);
	
Pe = 0.3;
mu_e = [530 1840];
std_e = [120  25
         25  190];
Ne = Pe*Ntot;
e = randn(Ne,2) * std_e + repmat(mu_e,Ne,1);

Pi = 0.25;
mu_i = [270 2290];
std_i = [50  5
         5  190];
Ni = Pi*Ntot;
i = randn(Ni,2) * std_i + repmat(mu_i,Ni,1);

Po = 0.15;
mu_o = [570 840];
std_o = [40  20
         20 140];
No = Po*Ntot;
o = randn(No,2) * std_o + repmat(mu_o,No,1);

Py = 0.05;
mu_y = [440 1020];
std_y = [80 40
         40 130];
Ny = Py*Ntot;
y = randn(Ny,2) * std_y + repmat(mu_y,Ny,1);

allvow = [a; e; i; o; y];
allvow = allvow(randperm(Ntot),:);

save vowels.mat a e i o y allvow

clf;
%set(gca,'xlim',[0 4000],'ylim',[0 4000],'dataaspectratio',[1 1 1e-6]);
set(gca,'xlim',[0 3000],'ylim',[0 3000],'dataaspectratio',[1 1 1e-6]);
%set(gca,'xlim',[0 1500],'ylim',[0 4000]);
hold on; grid on;
plot(a(:,1),a(:,2),'y+');
plot(e(:,1),e(:,2),'r+');
plot(i(:,1),i(:,2),'g+');
plot(o(:,1),o(:,2),'b+');
plot(y(:,1),y(:,2),'m+');
