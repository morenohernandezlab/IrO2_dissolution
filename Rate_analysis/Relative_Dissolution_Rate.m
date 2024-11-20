%% This script will calculate the relative dissolution rates of IrO2 nanocrystals from LP-TEM dissolution trajectories
% Matteo Fratarcangeli
% Ivan A. Moreno-Hernandez
% August 2024
%
% for information on preparing the video_info file, see https://github.com/morenohernandezlab/IrO2_Dissolution

clear all
clc

load video_info_IrO2_800C_FeCl_final_2024_05_16_rates_M_drive.mat

%% Analyze and adjust fits
clf
v = 1; %select video to analyze

start_index = 1;

fit_trajectory = video_info(v).fit;
fit_time = video_info(v).fit(:,2);
fit_L110 = video_info(v).fit(:,3)*2;
fit_L111_side1 = video_info(v).fit(:,5) + video_info(v).fit(:,8);
fit_L111_side2 = video_info(v).fit(:,6) + video_info(v).fit(:,7);
fit_L001 = video_info(v).fit(:,9);

dim_etched_L110 = fit_L110(start_index)-fit_L110(start_index:end);
dim_etched_L111_side1 = fit_L111_side1(start_index)-fit_L111_side1(start_index:end);
dim_etched_L111_side2 = fit_L111_side2(start_index)-fit_L111_side2(start_index:end);
dim_etched_L001 = fit_L001(start_index)-fit_L001(start_index:end);

dim_etched_L110_line = linspace(min(dim_etched_L110),max(dim_etched_L110),1001);

p_111_side1 = polyfit(dim_etched_L110,dim_etched_L111_side1,1)
p_111_side2 = polyfit(dim_etched_L110,dim_etched_L111_side2,1)
p_001 = polyfit(dim_etched_L110,dim_etched_L001,1)

clf
subplot(3,3,1)

plot(fit_time,fit_L110)

subplot(3,3,2)
plot(fit_time,fit_L111_side1)

subplot(3,3,3)
plot(fit_time,fit_L111_side2)

subplot(3,3,4)
plot(fit_time,fit_L001)


subplot(3,3,5)
hold on
plot(dim_etched_L110,dim_etched_L111_side1,'o')
plot(dim_etched_L110_line,dim_etched_L110_line*p_111_side1(1)+p_111_side1(2))

hold off
box on
subplot(3,3,6)
hold on

plot(fit_L110(start_index)-fit_L110(start_index:end),fit_L111_side2(start_index)-fit_L111_side2(start_index:end),'o')

plot(dim_etched_L110_line,dim_etched_L110_line*p_111_side2(1)+p_111_side2(2))
hold off
box on
subplot(3,3,7)
hold on
plot(fit_L110(start_index)-fit_L110(start_index:end),fit_L001(start_index)-fit_L001(start_index:end),'o')
plot(dim_etched_L110_line,dim_etched_L110_line*p_001(1)+p_001(2))

hold off
box on

video_info(v).L111_side1_relative_rate_fit = [start_index,p_111_side1];
video_info(v).L111_side2_relative_rate_fit = [start_index,p_111_side2];
video_info(v).L001_relative_rate_fit = [start_index,p_001];

%% Calculate difference in activation barrier
clearvars -except video_info
kb = 8.617333262*10^-5; %eV/K
T = 298.13; %K
c = 1000;%eV to meV

for i = 1:size(video_info,2)

L111_rel_rate_avg = (video_info(i).L111_side1_relative_rate_fit(:,2)+video_info(i).L111_side2_relative_rate_fit(:,2))/2;
video_info(i).delta_act_barrier = -kb*T*log(L111_rel_rate_avg)*c;
delta_activation_barrier(i) = video_info(i).delta_act_barrier;

end

delta_activation_barrier