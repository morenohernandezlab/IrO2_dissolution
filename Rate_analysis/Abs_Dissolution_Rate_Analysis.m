%% This script will calculate the absolute dissolution rates of IrO2 nanocrystals from LP-TEM dissolution trajectories
% Matteo Fratarcangeli
% Ivan A. Moreno-Hernandez
% August 2024
%
% for information on preparing the video_info file, see https://github.com/morenohernandezlab/IrO2_Dissolution

clear all
clc
load video_info_IrO2_800C_FeCl_final_2024_11_15_rates_M_drive.mat

%% Look at dimensions
clearvars -except video_info
v = 1; %select video to analyze
report_fit = [video_info(v).fit(:,2),video_info(v).fit(:,3:9)];


report_dim = [report_fit(:,1),report_fit(:,2)+report_fit(:,3),report_fit(:,4)+report_fit(:,7),report_fit(:,5)+report_fit(:,6),report_fit(:,8)];
% report_dim structure: [time, (110) dim, (111) side 1, (111) side 2, (001)]

offset = 8;

report_dim = report_dim(offset:end,:);
report_dim(:,1) = report_dim(:,1)-report_dim(1,1);
video_info_good(v).offset = offset;
%
hold on
tiledlayout(2,2)
for i = 1:4
    nexttile
plot(report_dim(:,1),report_dim(:,1+i))
end
hold off

%% Fit data
k_rates = linspace(0,1,31);
t_array = [.1,1,2,4];
s_array = [2,2.1,2.2,2.3,2.4,2.5,2.75,3,4,8,16];
[ka,ta,sa] = ndgrid(k_rates,t_array,s_array);
trial_values = [ka(:),ka(:),ka(:),ka(:),ta(:),sa(:)];
time_end = report_dim(end,1); %Seconds
M_IrO2 = 52.0025; %Molar density (mol/L) of IrO2 material
M_s = 20; %Molar density (mol/L) of saturated solution, estimate
tol = 4;
opts = odeset('AbsTol',10^-tol,'RelTol',10^-tol);
dim_initial = report_dim(1,2:5);
V_initial = (2*dim_initial(1))*(2*dim_initial(1))*dim_initial(4);

for i = 1:size(trial_values,1) %make par
i;
k_rates = trial_values(i,1:4); %Reaction rate for etchin in nm/s
t_rad = trial_values(i,5); %Radiolysis time constant in s, set max to 10 s
s = trial_values(i,6); % Radius of liquid pocket in nm, usually 1.1 to 10
ode_param = [M_IrO2,M_s,dim_initial,k_rates,t_rad,s,V_initial];
[time,dimensions] = ode45(@(t,y) myode3(t,y,ode_param),[0,time_end], [dim_initial],opts);
dimensions = max(dimensions,0);
interp_dimensions = interp1(time,dimensions,report_dim(:,1));
clf
hold on
plot(report_dim(:,1),interp_dimensions(:,1))
plot(report_dim(:,1),report_dim(:,2))

hold off 
drawnow
error_fit(i,:) = [i,sum(((interp_dimensions(:,[1])-report_dim(:,2))).^2,'all')];


end
%
error_no_fail = [];
for i = 1:size(error_fit,1)
if error_fit(i,1) ~= 0
error_no_fail = [error_no_fail;error_fit(i,:)];
end
end
[~,min_indx] = min(error_no_fail(:,2));
best_fit = trial_values(error_no_fail(min_indx,1),:); %[t_a(:),k_L_a(:),k_W_a(:),s_a(:)];
ode_param = [M_IrO2,M_s,dim_initial,best_fit(1:4),best_fit(6),best_fit(6),V_initial];
[time,dimensions] = ode45(@(t,y) myode3(t,y,ode_param),[0,time_end], dim_initial,opts);
interp_dimensions = interp1(time,dimensions,report_dim(:,1));
clf
subplot(2,2,1)
for k = 1:4
    nexttile
hold on
plot(report_dim(:,1),interp_dimensions(:,k))
plot(report_dim(:,1),report_dim(:,k+1))

hold off 
end
drawnow


% Refine fits
for u = 1:10
for j = 1:4
p = .4;
counts = 0;
while p > .0005
k_rates = linspace(1-p,1+p,9);

trial_values = k_rates(:);
error_fit = [];
for i = 1:size(trial_values,1) %make par
k_rates = best_fit(1:4);
k_rates(j) = k_rates(j)*trial_values(i);
t_rad = best_fit(5);
s = best_fit(6);
ode_param = [M_IrO2,M_s,dim_initial,k_rates,t_rad,s,V_initial];
[time,dimensions] = ode45(@(t,y) myode3(t,y,ode_param),[0,time_end], [dim_initial],opts);
interp_dimensions = interp1(time,dimensions,report_dim(:,1));
error_fit(i,:) = [i,sum(((interp_dimensions(:,j)-report_dim(:,j+1))).^2,'all')];
end
error_no_fail = [];
for i = 1:size(error_fit,1)
if error_fit(i,1) ~= 0
error_no_fail = [error_no_fail;error_fit(i,:)];
end
end

[~,min_indx] = min(error_no_fail(:,2));
best_fit(j) = best_fit(j)*trial_values(min_indx);



if min_indx == 5
    p = p*.9; %/2
end
counts = counts +1;
if counts > 5
p = p*.9; %/2
end

end

end
% See best fit
ode_param = [M_IrO2,M_s,dim_initial,best_fit(1:4),best_fit(5),best_fit(6),V_initial];

[time,dimensions] = ode45(@(t,y) myode3(t,y,ode_param),[0,time_end], dim_initial,opts);
%dimensions = max(dimensions,0);
interp_dimensions = [];
interp_dimensions = interp1(time,dimensions,report_dim(:,1));
clf
for k = 1:4
    nexttile
hold on
plot(report_dim(:,1),interp_dimensions(:,k))
plot(report_dim(:,1),report_dim(:,k+1))

hold off 
end
drawnow
end

% save best fit

%[(110) , (111) side 1, (111) side 2, (001), time constant for radiolysis, confinement]

video_info(v).kinetic_best_fit = best_fit;
video_info(v).kinetic_best_fit



%%
function dydt = myode3(t,y,ode_param)


%ode_param = [M_RuO2,M_s,dim_initial,k_rates,t_rad,s];
%ode_param = [M_RuO2,M_s,dim_initial,k_rates,t_rad,s,V_initial];
M_IrO2 = ode_param(1); %Molar density (mol/L) of IrO2 material
M_s = ode_param(2); %Molar density (mol/L) of saturated solution
dim_initial = ode_param(3:6);
k1 = ode_param(7);
k2 = ode_param(8);
k3 = ode_param(9);
k4 = ode_param(10);
t_rad = ode_param(11); %Radiolysis time constant in s
s = ode_param(12); % Radius of liquid pocket in nm, usually 1.1 to 10
V_initial = ode_param(13);
volume_shape = (2*y(1))*(2*y(1))*y(4);
v_ratio =  volume_shape./V_initial;
v_ratio = min(v_ratio,1);
%v_change = max(v_change,0)
V_t = V_initial*v_ratio; %volume(shp3d);

f = (V_initial-V_t)./(s*V_initial-V_t);
%(1-((M_IrO2/M_s)*f))

dydt(1) = -k1*(1-((M_IrO2/M_s)*f))*(1-exp(-t./t_rad)); 
dydt(2) = -k2*(1-((M_IrO2/M_s)*f))*(1-exp(-t./t_rad)); 
dydt(3) = -k3*(1-((M_IrO2/M_s)*f))*(1-exp(-t./t_rad)); 
dydt(4) = -k4*(1-((M_IrO2/M_s)*f))*(1-exp(-t./t_rad)); 


dydt = dydt';
end






