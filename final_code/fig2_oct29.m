%% Plot Fig 2
% Katharine Z Coyte April 2021

%%
figure

plot_selection = 0;
plot_log=1;
community_measure =1;
low_g = 1;
high_g = 72;

pre_selection = [89,116,164]/255;
post_selection = [203,137,99]/255;

Pm_ix = 1;

m_ix=1;

%% all species


delta_R_all = squeeze(abs(all_robustness_all(1, :,Pm_ix,:,1)) - abs(all_robustness_all(1, :,Pm_ix,:,2)))';
delta_E_all = squeeze(abs(all_robustness_all(1, :,Pm_ix,:,2)) - abs(all_robustness_all(2, :,Pm_ix,:,2)))';

pf_before_pla = squeeze(all_pf_subset(1, :,Pm_ix,:,2))';
pf_after_pla = squeeze(all_pf_subset(2, :,Pm_ix,:,2))';


%% Subset

delta_R_background = squeeze(abs(all_robustness_subset(1, :,Pm_ix,:,1)) - abs(all_robustness_subset(1, :,Pm_ix,:,2)))';
delta_E_background = squeeze(abs(all_robustness_subset(1, :,Pm_ix,:,2)) - abs(all_robustness_subset(2, :,Pm_ix,:,2)))';

%% Focal

delta_R_focal = squeeze(abs(all_robustness_resistant(1, :,Pm_ix,:,1)) - abs(all_robustness_resistant(1, :,Pm_ix,:,2)))';
delta_E_focal = squeeze(abs(all_robustness_resistant(1, :,Pm_ix,:,2)) - abs(all_robustness_resistant(2, :,Pm_ix,:,2)))';


%% Robustness

subplot(1,3,1)

my_gray = [0.25,0.25,0.25];

% DeltaR all
shadedErrorBar(gamma_space,delta_R_all,{@mean,@std},'lineprops', {'color', my_gray, 'LineWidth', 2},'patchSaturation',0.13)
hold on

% DeltaR background
shadedErrorBar(gamma_space,delta_R_background,{@mean,@std},'lineprops', {'color', my_gray, 'LineWidth', 2, 'LineStyle', '--'},'patchSaturation',0.13)
hold on

% DeltaR focal
shadedErrorBar(gamma_space,delta_R_focal,{@mean,@std},'lineprops', {'color', my_gray, 'LineWidth', 2, 'LineStyle', ':'},'patchSaturation',0.13)
hold on

xlabel("Mobility")
ylabel("\Delta R")


%% Selection all

subplot(1,3,2)

% DeltaE all
shadedErrorBar(gamma_space,delta_E_all,{@mean,@std},'lineprops', {'color', my_gray, 'LineWidth', 2},'patchSaturation',0.13)
hold on

% DeltaE background
shadedErrorBar(gamma_space,delta_E_background,{@mean,@std},'lineprops', {'color', my_gray, 'LineWidth', 2, 'LineStyle', '--'},'patchSaturation',0.13)

% DeltaE focal
shadedErrorBar(gamma_space,delta_E_focal,{@mean,@std},'lineprops', {'color', my_gray, 'LineWidth', 2, 'LineStyle', ':'},'patchSaturation',0.13)


% set(gca,'Yscale','log')
xlabel("Mobility")
ylabel("\Delta E")

%% Plasmid frequency

subplot(1,3,3)

shadedErrorBar(gamma_space,pf_before_pla,{@mean,@std},'lineprops', {'color', pre_selection, 'LineWidth', 2},'patchSaturation',0.13)
hold on
ylim([0,1])

shadedErrorBar(gamma_space,pf_after_pla,{@mean,@std},'lineprops', {'color', post_selection, 'LineWidth', 2},'patchSaturation',0.13)
hold on
ylim([0,1])


xlabel("Mobility")
ylabel("Plasmid frequency")

%%

hold on

if plot_log
    for i = 1:3
        subplot(1,3,i)
        set(gca,'Xscale','log')
        
    end
end

for i = 1:3
    subplot(1,3,i)
    axis square
end