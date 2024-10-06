%% The accuracy of Li(2017) is verified by utilizing the Ostrander model 
clc;
clear;
close all;
angle1=0:1:35;

% Ostrander model       
Vp1=3048;  Vs1=1244;  den1=2.400;  psb1=0.4;
Vp2=2428;  Vs2=1625;  den2=2.140;  psb2=0.1;
y=(Vs2+Vs1)^2/(Vp1+Vp2)^2;
angle2=asind(Vp1/Vp2*sind(angle1));
angle=0.5*(angle2+angle1);%Average of angle of incidence and angle of transmission
 
%Li equation: vp,psb,den
A_theta=1/2*(secd(angle).*secd(angle))-4*y*(sind(angle).*sind(angle));
B_theta=2*(sind(angle).*sind(angle))*(1+2*y*y-3*y);
C_theta=1/2-2*y*sind(angle).*sind(angle);
R_Li=2*A_theta*(Vp2-Vp1)/(Vp1+Vp2)+2*B_theta*(psb2-psb1)/(psb2+psb1)+2*C_theta*(den2-den1)/(den2+den1);
   
% zoeppritz equation
% calling the zoeppritz function
coef=zoeppritz(den1,Vp1,Vs1,den2,Vp2,Vs2,1,1,1,angle);
coef2=zoeppritz(den1,Vp1,Vs1,den2,Vp2,Vs2,1,1,2,angle);
A_coef=real(coef).*cos(imag(coef));
wucha1=sum(abs(R_Li-A_coef))./sum(abs(A_coef));%calculation error
wucha2=sum(abs(R_Li-A_coef))./sum(abs(A_coef));

% Figure 1
figure ('Name','Figure 1！！Comparison of reflection coefficient')
plot(angle,R_Li,'g-','LineWidth',1)
hold on
plot(angle,A_coef,'k','LineWidth',1)
set(gca,'yLim',[-0.51 0.01]);%Y-axis data display range
set(gca,'yTick',[-0.5:0.1:0],'Fontname', 'Times New Roman','Fontsize',10);%Set the coordinate scale to be displayed
set(gca,'xLim',[0 35]);
set(gca,'xTick',[0:10:35],'Fontname', 'Times New Roman','Fontsize',10);
set(gca,'XGrid','off');
set(gca,'YGrid','off');
ylabel('Reflection coefficient')
xlabel('Angle,^{\circ}');
set(gca,'xaxislocation','top')
legend('Li','Zoeppritz','Fontname','Times New Roman','Fontsize',10)


%% forward record and inversion result
%This program is used to extract the Well Control AVO operator and verify it with the forward record
%Further inversion using the Well Control AVO operator
clc;
clear;

%% Model data
load('vp.mat')
load('ps.mat')
load('den.mat')
load('vs.mat')
[iSample_num,iTrace_num]=size(ps);

% Angle range
angle=1:2:35;
angle=angle';
len_angle=length(angle); 

%% Calculate the reflection coefficient！！ref_vp;ref_ps;ref_den��ref_zop
ref_vp=zeros(iSample_num,iTrace_num);
ref_ps=zeros(iSample_num,iTrace_num);
ref_den=zeros(iSample_num,iTrace_num);
ref_zoep=zeros(iSample_num,iTrace_num,len_angle);
for j=1:iTrace_num
    for i=1:iSample_num-1
        ref_vp(i+1,j)=(vp(i+1,j)-vp(i,j))/(vp(i,j)+vp(i+1,j));
        ref_ps(i+1,j)=(ps(i+1,j)-ps(i,j))/(ps(i,j)+ps(i+1,j));
        ref_den(i+1,j)=(den(i+1,j)-den(i,j))/(den(i,j)+den(i+1,j));
        coef=zoeppritz(den(i,j),vp(i,j),vs(i,j),den(i+1,j),vp(i+1,j),vs(i+1,j),1,1,1,angle);
        ref_zoep(i+1,j,1:len_angle)=real(coef).*cos(imag(coef));
    end
end
ref_vppsden=[ref_vp;ref_ps;ref_den];

%% Li equation！！Vp,Ps,Den
k=1/4;
A_theta=0.5.*secd(angle).*secd(angle)-4.*k.*sind(angle).*sind(angle);
B_theta=2.*(1+2*k^2-3*k).*sind(angle).*sind(angle);
C_theta=0.5-2.*k.*sind(angle).*sind(angle);
A_matrix_Li=2*[A_theta B_theta C_theta];

%% Wavelet
wavelet_length=50;
dt=0.002;
fmain=30;
[wavelet,tw]=ricker(dt,fmain,wavelet_length*dt);
wavelet_length=length(wavelet);
G_wavelet=operaterG(wavelet,iSample_num);

%% Forward record 
% Reflection coefficient of Li equation
for j=1:iTrace_num
    for ang=1:1:len_angle
        aa=A_matrix_Li(ang,:)*[ref_vp(:,j)';ref_ps(:,j)';ref_den(:,j)'];
        ref_Li(:,j,ang)=aa';
    end
end

% Forward record of Zoeppritz equation (traditional forward framework)
sei_vp=zeros(iSample_num,iTrace_num);
sei_ps=zeros(iSample_num,iTrace_num);
sei_den=zeros(iSample_num,iTrace_num);
for j=1:iTrace_num
    for ang=1:1:len_angle
        sei_zoep(:,j,ang)=G_wavelet*ref_zoep(:,j,ang);
        sei_zoep_20(:,j,ang)=rnoise(sei_zoep(:,j,ang),5);
        sei_zoep_50(:,j,ang)=rnoise(sei_zoep(:,j,ang),2);
        sei_Li_1(:,j,ang)=G_wavelet*ref_Li(:,j,ang);
    end
    sei_vp(:,j)=G_wavelet*ref_vp(:,j);
    sei_ps(:,j)=G_wavelet*ref_ps(:,j);
    sei_den(:,j)=G_wavelet*ref_den(:,j);
end 
sei_zoep_no=sei_zoep;

% Simplified forward framework 
for j=1:iTrace_num
    for ang=1:1:len_angle
        aa=A_matrix_Li(ang,:)*[sei_vp(:,j)';sei_ps(:,j)';sei_den(:,j)'];
        sei_Li_2(:,j,ang)=aa;
    end
end

% Forward record of Li equation (traditional forward framework)
Gp_Li=kron(A_matrix_Li,G_wavelet);
sei_xkp_3spare=Gp_Li*ref_vppsden;
A_matrix_xkp_spare=angleterm_spare_matrix_wu(A_matrix_Li,iSample_num);
G_wavelet_spare=Gwavelet_spare_matrix_wu(G_wavelet,len_angle);
Gp_Li=G_wavelet_spare*A_matrix_xkp_spare;
sei_xkp_4spare=Gp_Li*ref_vppsden;

for j=1:iTrace_num
    a_wc=reshape(sei_xkp_3spare(:,j),iSample_num,len_angle);
    b_wc=reshape(sei_xkp_4spare(:,j),iSample_num,len_angle);
    for ang=1:1:len_angle
        sei_Li_3(:,j,ang)=a_wc(:,ang);
        sei_Li_4(:,j,ang)=b_wc(:,ang);
    end
end 

% If it is a multi-channel (2D) model, take one of the multi-angle seismic data
jj=1;
sei_zoep_single=zeros(iSample_num,iTrace_num);
sei_Li_1_single=zeros(iSample_num,iTrace_num);
sei_Li_2_single=zeros(iSample_num,iTrace_num);
sei_Li_3_single=zeros(iSample_num,iTrace_num);
sei_Li_4_single=zeros(iSample_num,iTrace_num);

for ang=1:1:len_angle
    sei_zoep_single(:,ang)=sei_zoep(:,jj,ang);
    sei_zoep_single_no(:,ang)=sei_zoep_no(:,jj,ang);
    sei_Li_1_single(:,ang)=sei_Li_1(:,jj,ang);
    sei_Li_2_single(:,ang)=sei_Li_2(:,jj,ang);
    sei_Li_3_single(:,ang)=sei_Li_3(:,jj,ang);    
    sei_Li_4_single(:,ang)=sei_Li_4(:,jj,ang);
end
sei_vp_single=sei_vp(:,jj);
sei_ps_single=sei_ps(:,jj);
sei_den_single=sei_den(:,jj);

%% Well Control AVO 
sei_atr=[sei_vp_single';sei_ps_single';,sei_den_single'];
[s,v,d]=svd(sei_atr'); 
v(v <= 0.001) = 0; 
v_1=1./v;
v_1(v_1==-inf)= 0;  v_1(v_1==inf)= 0;  
ginv_sei_atr=d*(v_1')*s';
yz_1=sei_atr'*ginv_sei_atr;
yz_2=ginv_sei_atr*sei_atr';

% matrix division  
A_matrix_wc_1=sei_zoep_single'/sei_atr;
% Generalized Inverse-Row Full Rank
A_matrix_wc_2=sei_zoep_single'*(sei_atr'*inv(sei_atr*sei_atr'));
% Singular Value Decomposition
A_matrix_wc_3=sei_zoep_single'*ginv_sei_atr';
% iterative algorithm
max_iter_irls=5; 
max_iter_cgls=10; 
Wr = ones(size(sei_zoep_single));%Wr=1;
Wx = ones(size(A_matrix_Li'));%Wx=1;
muu=0;
[A_matrix_wc_4,J_6] = mirls(A_matrix_Li', sei_zoep_single, sei_atr', muu, max_iter_cgls, max_iter_irls, 1e-6,1e-6);
A_matrix_wc_4=A_matrix_wc_4';

% Forward record of simplified forward framework 
Gp_wc=kron(A_matrix_wc_3,G_wavelet);
sei_wja_spare=Gp_wc*ref_vppsden;
for j=1:iTrace_num
    a_wc=reshape(sei_wja_spare(:,j),iSample_num,len_angle);
    for ang=1:1:len_angle
        sei_wc(:,j,ang)=a_wc(:,ang);
    end
end 

sei_wc_single=zeros(iSample_num,iTrace_num);
for ang=1:1:len_angle
    sei_wc_single(:,ang)=sei_wc(:,jj,ang);
end
sei_wc_single_20=[];
sei_wc_single_50=[];
for ang=1:1:len_angle    
    sei_wc_single_20(:,ang)=rnoise(sei_wc_single(:,ang),5);
    sei_wc_single_50(:,ang)=rnoise(sei_wc_single(:,ang),2);
    for j=1:iTrace_num
        sei_wja_20(:,j,ang)=rnoise(sei_wc(:,j,ang),5);
        sei_wja_50(:,j,ang)=rnoise(sei_wc(:,j,ang),2);
    end
end

for j=1:iTrace_num
    for ang=1:1:len_angle
        aa=A_matrix_wc_3(ang,:)*[ref_vp(:,j)';ref_ps(:,j)';ref_den(:,j)'];
        ref_wja(:,j,ang)=aa';
    end
end



%% Figure in the article
num_start=10;
num_end=330;
sei_zoep_ed=sei_zoep_single(num_start:num_end,:);
sei_zoep_single_no_ed=sei_zoep_single_no(num_start:num_end,:);
sei_Li_spa_ed=sei_Li_3_single(num_start:num_end,:);
sei_Li_syn_ed=sei_Li_2_single(num_start:num_end,:);
sei_wc_ed=sei_wc_single(num_start:num_end,:);
sei_cancha_ed=sei_Li_spa_ed-sei_zoep_ed;
[ns ls]=size(sei_zoep_ed);
siz=10;
y_tic=[55 110 165 220 275];
y_lab=[0.1 0.2 0.3 0.4 0.5];
x_tic=min(angle):5:max(angle);

% Figure3:Accuracy verification of the simplified forward framework
figure('Name','Figure3:Accuracy verification of the simplified forward framework')
subplot(1,2,1)
plotseis(sei_Li_spa_ed,1:ns,angle',1,1,1,1,[0 0 0]);
ylabel('Time,s')
xlabel('Angle,^{\circ}');
set(gca,'xTick',x_tic,'xlim',([min(angle)-2,max(angle)+2]),'yTick',y_tic,'yTicklabel',y_lab,'Fontname', 'Times New Roman','Fontsize', siz);
subplot(1,2,2)
plotseis(sei_Li_syn_ed,1:ns,angle',1,1,1,1,[0 0 0]);
ylabel('Time,s')
xlabel('Angle,^{\circ}');
set(gca,'xTick',x_tic,'xlim',([min(angle)-2,max(angle)+2]),'yTick',y_tic,'yTicklabel',y_lab,'Fontname', 'Times New Roman','Fontsize', siz);

% Figure4:Feasibility validation of Well Control AV
figure('Name','Figure4:Feasibility validation of Well Control AVO')
subplot(2,2,1)
plotseis(sei_zoep_ed,1:ns,angle',1,1,1,1,[0 0 0]);
ylabel('Time,s')
xlabel('Angle,^{\circ}');
set(gca,'xTick',x_tic,'xlim',([min(angle)-2,max(angle)+2]),'yTick',y_tic,'yTicklabel',y_lab,'Fontname', 'Times New Roman','Fontsize', siz);
subplot(2,2,2)
plotseis(sei_Li_spa_ed,1:ns,angle',1,1,1,1,[0 0 0]);
ylabel('Time,s')
xlabel('Angle,^{\circ}');
set(gca,'xTick',x_tic,'xlim',([min(angle)-2,max(angle)+2]),'yTick',y_tic,'yTicklabel',y_lab,'Fontname', 'Times New Roman','Fontsize', siz);
subplot(2,2,3)
plotseis(sei_wc_ed,1:ns,angle',1,1,1,1,[0 0 0]);
ylabel('Time,s')
xlabel('Angle,^{\circ}');
set(gca,'xTick',x_tic,'xlim',([min(angle)-2,max(angle)+2]),'yTick',y_tic,'yTicklabel',y_lab,'Fontname', 'Times New Roman','Fontsize', siz);

% Figure5a！！Error between forward record of Well Control AVO and forward record of Zoeppritz equation
figure('name','Figure5a！！Error between forward record of Well Control AVO and forward record of Zoeppritz equation');
plotseis(sei_wc_ed-sei_zoep_single_no_ed,1:ns,angle',1,1,1,1,[0 0 0]);
set(gca,'xTick',x_tic,'xlim',([min(angle)-2,max(angle)+2]),'yTick',y_tic,'yTicklabel',y_lab,'Fontname', 'Times New Roman','Fontsize', siz);

%Figure5b！！Error between forward record of theoretical approximation formula and forward record of Zoeppritz equation
figure('name','Figure5b！！Error between forward record of theoretical approximation formula and forward record of Zoeppritz equation')
plotseis(sei_cancha_ed,1:ns,angle',1,1,1,1,[0 0 0]);
set(gca,'xTick',x_tic,'xlim',([min(angle)-2,max(angle)+2]),'yTick',y_tic,'yTicklabel',y_lab,'Fontname', 'Times New Roman','Fontsize', siz);

% Figure5c！！Reflection coefficient at 0.4 seconds
siz_line=12;
time_num=220;
time_num=time_num+num_start;
r_1=reshape(ref_zoep(time_num,:,:),[1,len_angle]);
r_2=reshape(ref_wja(time_num,:,:),[1,len_angle]);
r_3=reshape(ref_Li(time_num,:,:),[1,len_angle]);
figure('name','Figure5c！！Reflection coefficient at 0.4 seconds')
plot(angle,r_1,'k--',angle,r_2,'r',angle,r_3,'b')
ylabel('Reflectivity')
xlabel('Angle,^{\circ}');
set(gca,'xTick',[0,10,20,30],'xlim',([min(angle)-1,max(angle)+1]),'XAxisLocation', 'top','ylim',([0.03,0.055]),'yTick',[0.03,0.04,0.05],'Fontname', 'Times New Roman','Fontsize', siz_line);
legend({'Zoeppritz','New','Li'},'Location','southwest','FontSize',8)

% Figure2！！1D model data
y_tic=[];
y_lab=[];
figure('name','Figure2！！1D model data')
subplot(1,3,1)
imagesc(vp(num_start:num_end,:))    
set(gca,'xlim',([0.5,1]),'xTick',[],'yTick',y_tic,'yTicklabel',y_lab,'Fontname', 'Times New Roman','Fontsize', siz);
colormap(flipud(jet))
colorbar('Ticks',[3600,4600,5600],'Fontname', 'Times New Roman','Fontsize', siz)
subplot(1,3,2)
imagesc(ps(num_start:num_end,:))  
set(gca,'xlim',([0.5,1]),'xTick',[],'yTick',y_tic,'yTicklabel',y_lab,'Fontname', 'Times New Roman','Fontsize', siz);
colormap(flipud(jet))
colorbar('Ticks',[0.19,0.21,0.23],'Fontname', 'Times New Roman','Fontsize', siz)
subplot(1,3,3)
imagesc(den(num_start:num_end,:))    
set(gca,'xlim',([0.5,1]),'xTick',[],'yTick',y_tic,'yTicklabel',y_lab,'Fontname', 'Times New Roman','Fontsize', siz);
colormap(flipud(jet))
colorbar('Ticks',[2370,2540,2710],'Fontname', 'Times New Roman','Fontsize', siz)

%% Inversion
for j=1:iTrace_num
 sei_for_inv(:,j)=reshape(sei_zoep(:,j,:),iSample_num*len_angle,1);
end
magnitude_seis=max(abs(sei_for_inv(:)));
sei_for_inv=sei_for_inv/magnitude_seis;
Gp_Li=Gp_Li/magnitude_seis;
Gp_wc=Gp_wc/magnitude_seis;

% Constructing low-frequency models
fvp_min=3;      fvp_max=6;
fps_min=3;      fps_max=6;
fden_min=3;     fden_max=6;
E_matrix=Reverse_extension_matrix(iSample_num);
F=dftmtx(3*iSample_num);
h_win2_vp = f_win_wu(0,0,fvp_min,fvp_max,3*iSample_num,dt);
Win_matrix=filter_vec2mat(h_win2_vp);
L1_vp=real(E_matrix'*F^-1*Win_matrix*F*E_matrix)/3;
h_win2_ps = f_win_wu(0,0,fps_min,fps_max,3*iSample_num,dt);
Win_matrix=filter_vec2mat(h_win2_ps);
L1_ps=real(E_matrix'*F^-1*Win_matrix*F*E_matrix)/3;
h_win2_den = f_win_wu(0,0,fden_min,fden_max,3*iSample_num,dt);
Win_matrix=filter_vec2mat(h_win2_den);
L1_den=real(E_matrix'*F^-1*Win_matrix*F*E_matrix)/3;
LP=Blocking_matrix_f_win_3(L1_vp,L1_ps,L1_den);
ivp2=[];ips2=[];iden2=[];
for j=1:1:size(vp,2)
    ivp2(:,j)=L1_vp*vp(:,j);
    ips2(:,j)=L1_ps*ps(:,j);
    iden2(:,j)=L1_den*den(:,j);
end
ivp=ivp2; ips=ips2; iden=iden2;

tic

for j=1:iTrace_num
    fprintf('trace %dth\n',j);
    
    % No noise  
    gamma_wc=0.001;
    mu_vp_wc =1.0;
    mu_ps_wc = 1.5;
    mu_den_wc = 1.3;
    gamma_Li=0.005;
    mu_vp_Li = 1.3;
    mu_ps_Li = 1.5;
    mu_den_Li =1.3;

    % gamma_wc=0.001;
    % mu_vp_wc =1.5;
    % mu_ps_wc = 1.5;
    % mu_den_wc = 1.3;
    % gamma_Li=0.005;
    % mu_vp_Li = mu_vp_wc ;
    % mu_ps_Li = mu_ps_wc ;
    % mu_den_Li =mu_den_wc ;

    C=ones(iSample_num,iSample_num);
    C=tril(C);
    LL=Make_L_mat(iSample_num);
    L_matrix=kron(eye(3,3),LL);
    
    % Low frequency model
    Xivp=Xi_matrix(ivp(:,j),iSample_num);
    Xips=Xi_matrix(ips(:,j),iSample_num);
    Xiden=Xi_matrix(iden(:,j),iSample_num);
    
    Cxi_vppsden=[Xivp;Xips;Xiden];
    CXi_wc=[mu_vp_wc*Xivp;mu_ps_wc*Xips;mu_den_wc*Xiden];
    Cp_wc=kron([mu_vp_wc 0 0;0 mu_ps_wc 0;0 0 mu_den_wc],C);
    CXi_Li=[mu_vp_Li*Xivp;mu_ps_Li*Xips;mu_den_Li*Xiden];
    Cp_Li=kron([mu_vp_Li 0 0;0 mu_ps_Li 0;0 0 mu_den_Li],C);
    
    % Parameter de-correlation
    V=decorrelation([ref_vp(:,j) ref_ps(:,j) ref_den(:,j)],iSample_num);

    A_wc=[Gp_wc*V;LP*Cp_wc*V];
    b_wc=[sei_for_inv(:,j);CXi_wc];
    A_li=[Gp_Li*V;Cp_Li*V];
    b_li=[sei_for_inv(:,j);CXi_Li];
    R_test_wc=zeros(iSample_num*3,1);
    R_test_Li=zeros(iSample_num*3,1);
    
    % Algorithmic solution
    max_iter_cgls=200;
    max_iter_irls=7;
    [R_test_wc,J] = mirls(0.5*L_matrix*Cxi_vppsden(:), b_wc, A_wc, gamma_wc, max_iter_cgls, max_iter_irls, 1e-6,1e-6);
    [R_test_Li,J] = mirls(0.5*L_matrix*Cxi_vppsden(:), b_li, A_li, gamma_Li, max_iter_cgls, max_iter_irls, 1e-6,1e-6);
    r_test_wja=reshape(V*R_test_wc,iSample_num,3);
    r_test_xkp=reshape(V*R_test_Li,iSample_num,3);

    % Integration of reflection coefficients
    for m=1:iSample_num
        vp_inv_wc(m,j)=ivp(1,1)*exp(2.0*sum_matrix(r_test_wja(:,1),m));
        ps_inv_wc(m,j)=ips(1,1)*exp(2.0*sum_matrix(r_test_wja(:,2),m));
        den_inv_wc(m,j)=iden(1,1)*exp(2.0*sum_matrix(r_test_wja(:,3),m));
        vp_inv_Li(m,j)=ivp(1,1)*exp(2.0*sum_matrix(r_test_xkp(:,1),m));
        ps_inv_Li(m,j)=ips(1,1)*exp(2.0*sum_matrix(r_test_xkp(:,2),m));
        den_inv_Li(m,j)=iden(1,1)*exp(2.0*sum_matrix(r_test_xkp(:,3),m));    
    end

%% Figure of Inversion result 
y_tic=[55 110 165 220 275];
y_lab=[0.1 0.2 0.3 0.4 0.5]; 
siz=12;
num_start=17;
num_end=330;
[ns ls]=size(vp(num_start:num_end,:));
yt_sei=1:1:ns;
x_tic_vp=[3500 5500];
x_tic_ps=[0.18 0.24];
x_tic_den=[2300 2700];

% Figure9a！！Inversion results of model data
figure('Name','Figure9a！！Inversion results of model data')
subplot(1,3,1)
plot(vp(num_start:num_end,:),yt_sei,'k')
hold on
plot(ivp(num_start:num_end,:),yt_sei,'g--')
hold on 
plot(vp_inv_wc(num_start:num_end,:),yt_sei,'r')
hold on 
plot(vp_inv_Li(num_start:num_end,:),yt_sei,'b')
set(gca,'xTick',x_tic_vp,'XAxisLocation','top','YDir', 'reverse','yTick',y_tic,'yTicklabel',y_lab,'Fontname', 'Times New Roman','Fontsize', siz);
axis([2500,6500,0,ns])
subplot(1,3,2)
plot(ps(num_start:num_end,:),yt_sei,'k')
hold on
plot(ips(num_start:num_end,:),yt_sei,'g--')
hold on 
plot(ps_inv_wc(num_start:num_end,:),yt_sei,'r')
hold on 
plot(ps_inv_Li(num_start:num_end,:),yt_sei,'b')
set(gca,'xTick',x_tic_ps,'XAxisLocation','top','YDir', 'reverse','yTick',y_tic,'yTicklabel',[],'Fontname', 'Times New Roman','Fontsize', siz);
axis([0.16,0.27,0,ns])
subplot(1,3,3)
plot(den(num_start:num_end,:),yt_sei,'k')
hold on
plot(iden(num_start:num_end,:),yt_sei,'g--')
hold on 
plot(den_inv_wc(num_start:num_end,:),yt_sei,'r')
hold on 
plot(den_inv_Li(num_start:num_end,:),yt_sei,'b')
set(gca,'xTick',x_tic_den,'XAxisLocation','top','YDir', 'reverse','yTick',y_tic,'yTicklabel',[],'Fontname', 'Times New Roman','Fontsize', siz);
axis([2150,2900,0,ns])

%% Error analysis of inversion results
figure('name','Figure10a！！Error analysis of inversion results')
size=12;
vp_eor_xkp=vp_inv_Li(num_start:num_end,:)-vp(num_start:num_end,:);
vp_eor_wja=vp_inv_wc(num_start:num_end,:)-vp(num_start:num_end,:);
ps_eor_xkp=ps_inv_Li(num_start:num_end,:)-ps(num_start:num_end,:);
ps_eor_wja=ps_inv_wc(num_start:num_end,:)-ps(num_start:num_end,:);
den_eor_xkp=den_inv_Li(num_start:num_end,:)-den(num_start:num_end,:);
den_eor_wja=den_inv_wc(num_start:num_end,:)-den(num_start:num_end,:);
subplot(131)
h2 = histogram(vp_eor_xkp); h2.FaceColor = [0 0 1];  h2.EdgeColor = 'b';  h2.BinEdges = [-600:100:600];
hold on
h1 = histogram(vp_eor_wja); h1.FaceColor = [1 0 0];  h1.EdgeColor = 'r';  h1.BinEdges = [-600:100:600];
set(gca,'xTick',-500:500:500,'ylim',[0 180],'yTick',0:80:200,'Fontname', 'Times New Roman','Fontsize', siz); 
subplot(132)
h2 = histogram(ps_eor_xkp); h2.FaceColor = [0 0 1];  h2.EdgeColor = 'b';  h2.BinEdges = [-0.018:0.003:0.015];
hold on
h1 = histogram(ps_eor_wja); h1.FaceColor = [1 0 0];  h1.EdgeColor = 'r';  h1.BinEdges = [-0.018:0.003:0.015];
set(gca,'xTick',-0.01:0.01:0.01,'ylim',[0 150],'yTick',0:60:160,'Fontname', 'Times New Roman','Fontsize', siz); 
subplot(133)  
h2 = histogram(den_eor_xkp); h2.FaceColor = [0 0 1];  h2.EdgeColor = 'b';  h2.BinEdges = [-105:15:90];
hold on
h1 = histogram(den_eor_wja); h1.FaceColor = [1 0 0];  h1.EdgeColor = 'r';  h1.BinEdges = [-105:15:90];
set(gca,'xTick',-80:80:80,'ylim',[0 140],'yTick',0:55:160,'Fontname', 'Times New Roman','Fontsize', siz); 

end

toc
