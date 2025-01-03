clear
clc
%% ion data & flux data & end member
data = importdata('Data_1.xlsx'); %load the data
ion_data = data.data;
K = ion_data(:,1);
CaMg = ion_data(:,2);
Na = ion_data(:,3);
n = 100000;
Nsam = length(K);
sigma = ion_data(:,4);
Cl = ion_data(:,5);
SO4 = ion_data(:,6);
Re = ion_data(:,7);  %pmol/L
discharge = ion_data(:,8);
Area = ion_data(:,9);
Re_C = ion_data(:,10);
Runoff = ion_data(:,11);  %Runoff_all(mm/yr)
Re_error = ion_data(:,12)./2;   %pmol/L
Re_nor_l = ion_data(:,13); 
Re_nor_u = ion_data(:,14);
Re_C_err = ion_data(:,15)./2;

member = importdata('end_member_1.xlsx');
end_member = member.data;
Nm = length(end_member);

F_sil = NaN(n,Nsam);
F_carb = NaN(n,Nsam);
F_eva = NaN(n,Nsam);
MemK_sil = NaN(n,Nsam);
MemCaMg_sil = NaN(n,Nsam);
MemNa_sil = NaN(n,Nsam);
MemCaMg_eva = NaN(n,Nsam);
MemNa_eva = NaN(n,Nsam);

%% random pd
for i = 1:Nm
    pd = makedist('Uniform','lower',end_member(i,1),'upper',end_member(i,2));
    eval(['pd_',member.textdata{i+1,1}, '= pd;']);
end

%% PART ONE inverse model for end-member mixing analysis
for i = 1:n
    for j = 1:Nm
        eval(['pdi = pd_',member.textdata{j+1,1},';']); 
        eval([member.textdata{j+1,1}, '= random(pdi);']); 
    end
    Mem = [K_sil,     0,         0    
           CaMg_sil,CaMg_carb, CaMg_eva;
           Na_sil,    0,        Na_eva];      
    Mar = [K';CaMg';Na'];
    F = Mem\Mar;
    %verification
    neg = find(F < 0);
    F(neg) = 1000;
    F_sum = sum(F,1); %sum absolute values of each row
	F_indi = find(F_sum < 2); %find indices where sum(x)=sum(abs(x))
    
    F_sil(i,F_indi) = F(1,F_indi); %mixing fractions     
    F_carb(i,F_indi) = F(2,F_indi);
    F_eva(i,F_indi) = F(3,F_indi);
    Class = ~isnan(F_sil(i,:));
    for pp = 1:length(Class)
        if Class(pp) == 1
           MemK_sil(i,pp) = Mem(1,1);
           MemCaMg_sil(i,pp) = Mem(2,1);
           MemNa_sil(i,pp) = Mem(3,1);
           MemCaMg_eva(i,pp) = Mem(2,3);
           MemNa_eva(i,pp) = Mem(3,3);
        end
    end
end
F_sil_mean = mean(F_sil,"omitnan");
F_carb_mean = mean(F_carb,"omitnan");
F_eva_mean = mean(F_eva,"omitnan");   %proportion of catio derived from evaporite (eq/eq)
%{
K_sil_cc = mean(F_sil.*MemK_sil,"omitnan");   %eq/eq
CaMg_sil_cc = mean(F_sil.*MemCaMg_sil,"omitnan");
Na_sil_cc = mean(F_sil.*MemNa_sil,"omitnan");
CaMg_eva_cc = mean(F_eva.*MemCaMg_eva,"omitnan");
Na_eva_cc = mean(F_eva.*Na_eva,"omitnan");
%}
Na_sil_cs = F_sil.*MemNa_sil;  %eq/eq
cation_sil = F_sil.*MemK_sil + F_sil.*MemCaMg_sil + F_sil.*MemNa_sil;   

%% silicate calculation
cation_sil_cc = cation_sil'.*sigma.*discharge.*(1E-3*365*3600*24/10000/2)./Area;  % (10^5 mol/km2/yr)
%{
cation_sil_m = mean(cation_sil_cc,2,"omitnan");   % mean value of data
cation_sil_s = std(cation_sil_cc,0,2,"omitnan");  % std value of data
cation_sil_median = median(cation_sil_cc,2,"omitnan");  % median value of data
cation_sil_p = prctile(cation_sil_cc,[25,50,75],2);   % 25th,50th,75th value of data (J_CO2_silicate)
cation_sil_ml = zeros(Nsam,1);
cation_sil_medianl = zeros(Nsam,1);
cation_sil_sl = zeros(Nsam,1);
%}
pro = 1;
for row = 1:Nsam
    if sum(isnan(cation_sil_cc(row,:))) ~= n
       figure (1)
       subplot(3,4,pro),histogram(cation_sil_cc(row,:),30);
       xlabel('J-CO2-silicate (mol/km2/yr)');
       ylabel('Frequency');
       title(['Sample QH' num2str(row)]);
       hold on;
       pro = pro+1;
       ication = cation_sil_cc(row,:)';
       pdi = fitdist(ication,'Lognormal');
       %cation_sil_ml(row,1) = exp(pdi.mu+pdi.sigma^2/2);  % mean value of lognormal distribution
       %cation_sil_medianl(row,1) =exp(pdi.mu);  % median value of lognormal distribution
       %cation_sil_sl(row,1) = ((exp(pdi.sigma^2)-1)*exp(2*pdi.mu+pdi.sigma^2))^0.5; % std value of lognormal distribution
    end
end
%% Re source portioning
Re_OC = NaN(n,15);
Re_OC_pt = NaN(n,15);
J_OC = NaN(n,15);
J_OC_r_m = NaN(1,15);
J_OC_r_s = NaN(1,15);
for row = 1:Nsam
    Re_e = Re(row,1) + randn(n,1).*Re_error(row,1);   %normal distrrbution
    Re_sil_ccs = Na_sil_cs(:,row)*sigma(row,1).*(2*10^(-4)+(5*10^(-3)-2*10^(-4)).*rand(n,1)); %uniform distribution 
    Re_so4 = SO4(row,1).*(2*10^(-4)+(4*10^(-3)-2*10^(-4)).* rand(n,1));  %uniform distribution 
    Re_OC(:,row) = Re_e - Re_sil_ccs - Re_so4;  %Re contributed by OC
    Re_OC_pt(:,row) = (Re_e - Re_sil_ccs - Re_so4)./Re_e;  %ratio of Re contributed by OC (fc)
end
Re_OC_mean = mean(Re_OC,"omitnan");
Re_OC_std = std(Re_OC,"omitnan");
Re_OC_pt_mean = mean(Re_OC_pt,"omitnan");
Re_OC_pt_std = std(Re_OC_pt,"omitnan");

%% CO2 yield
for row = 1:Nsam
    Re_nor = Re_nor_l(row,1) + (Re_nor_u(row,1)-Re_nor_l(row,1))*rand(n,1);
    Re_C_r = Re_C(row,1) + randn(n,1).*Re_C_err(row,1); %normal distrrbution
    F_gra = 0.7 + (1-0.7)*rand(n,1);
    J_OC_r = Re_nor.*Re_OC_pt(:,row).*(185/100000)./Re_C_r.*Runoff(row,1).*F_gra;  %tC/km2/yr
    J_OC(:,row) = J_OC_r;  %tC/km2/yr
    if sum(isnan(J_OC_r)) ~= n
       figure (1)
       subplot(3,4,pro),histogram(J_OC_r,30,'FaceColor','r');  %normal distribuion test
       xlabel('J-CO2-OC (tC/km2/yr)');
       ylabel('Frequency');
       title(['Sample QH' num2str(row)]);
       hold on;
       pro = pro+1;
       pdir = fitdist(J_OC_r,'Normal');
       J_OC_r_m(1,row) = pdir.mu;
       J_OC_r_s(1,row) = pdir.sigma;
    end
end 
J_OC_mean =  mean(J_OC,"omitnan");   %mean value
J_OC_std = std(J_OC,"omitnan");  %std value

%% Net carbon budget for sample QH07
J_pyrite = 0.66 +(1.188-0.66).*rand(n,1);
OC_bio = 0.27 +(0.33-0.27).*rand(n,1);
Silicate_sink = (cation_sil_cc(7,:))'*12/1000000;
Net = J_OC(:,7) + J_pyrite - OC_bio - Silicate_sink;
figure (2)
histogram(Net, 30,'FaceColor','c');
xlabel('Net carbon budget (tC/km2/yr)');
ylabel('Frequency');
title('Sample QH07');
Net_prc = prctile(Net,[25,50,75]);   % 25th,50th,75th value of data



