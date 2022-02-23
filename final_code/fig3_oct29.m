%% Plot Fig 3
% Katharine Z Coyte April 2021


%%
%

my_fig = figure




plot_log=1;
community_measure =1;

Pm_ix=9;



pre_selection = [152,152,152]/255;
post_selection = [152,152,152]/255;

CT=cbrewer('div','RdBu',100);



%% all species

m_ix=1;

% Calculate difference as proportion of susceptible
all_robustness_all(all_robustness_all>0)=0;
all_robustness_subset(all_robustness_subset>0)=0;
all_robustness_resistant(all_robustness_resistant>0)=0;

robust_sus_all = squeeze(all_robustness_all(m_ix, :, :, :, 1));
robust_plasmid_all = squeeze(all_robustness_all(m_ix, :, :, :, 2));
diff_robust_all = (abs(robust_sus_all) - abs(robust_plasmid_all));%./abs(robust_sus_all);
diff_robust_all = mean(diff_robust_all, 3);

robust_sus_subset = squeeze(all_robustness_subset(m_ix, :, :, :, 1));
robust_plasmid_subset = squeeze(all_robustness_subset(m_ix, :, :, :, 2));
diff_robust_subset = (abs(robust_sus_subset) - abs(robust_plasmid_subset));%./abs(robust_sus_subset);
diff_robust_subset = mean(diff_robust_subset, 3);

robust_sus_resistant = squeeze(all_robustness_resistant(m_ix, :, :, :, 1));
robust_plasmid_resistant = squeeze(all_robustness_resistant(m_ix, :, :, :, 2));
diff_robust_resistant = (abs(robust_sus_resistant) - abs(robust_plasmid_resistant));%./abs(robust_sus_resistant);
diff_robust_resistant = mean(diff_robust_resistant, 3);


my_clims = [-0.5,0.5];

subplot(2,3,1)
imagesc(diff_robust_all');
hold on 
%contour(real(diff_robust_all'), [0,0], 'LineWidth',1, 'LineColor','w')

caxis(my_clims)
colormap(CT)
ylabel("% Positive interactions")
xlabel("Conjugation rate")

xticks(linspace(1,length(gamma_space),length(gamma_space)))
xticklabels(gamma_space)
yticks(linspace(1,length(PmSpace),length(PmSpace)))
yticklabels(PmSpace)
set(gca,'YDir','normal')

h = colorbar;
ylabel(h, '\Delta stability')
title("Whole community")
axis square

subplot(2,3,2)
imagesc(diff_robust_subset');

hold on 
%contour(real(diff_robust_subset'),[0,0], 'LineWidth',1, 'LineColor','k')

my_clims = [-0.1,0.1];

caxis(my_clims)
colormap(CT)
ylabel("% Positive interactions")
xlabel("Conjugation rate")

xticks(linspace(1,length(gamma_space),length(gamma_space)))
xticklabels(gamma_space)
yticks(linspace(1,length(PmSpace),length(PmSpace)))
yticklabels(PmSpace)
set(gca,'YDir','normal')

h = colorbar;
ylabel(h, '\Delta stability')
title("Background community")
axis square


subplot(2,3,3)
imagesc(diff_robust_resistant');

hold on 
%contour(real(diff_robust_resistant'),[0,0], 'LineWidth',1, 'LineColor','k')

my_clims = [-1,1];
%my_clims = [-0.5,0.5];
caxis(my_clims)
colormap(CT)
ylabel("% Positive interactions")
xlabel("Conjugation rate")

xticks(linspace(1,length(gamma_space),length(gamma_space)))
xticklabels(gamma_space)
yticks(linspace(1,length(PmSpace),length(PmSpace)))
yticklabels(PmSpace)
set(gca,'YDir','normal')

h = colorbar;
ylabel(h, '\Delta stability')
title("Focal species")
axis square

%%



% all species

m_ix=1;

% Calculate difference as proportion of susceptible
all_robustness_all(all_robustness_all>0)=0;
all_robustness_subset(all_robustness_subset>0)=0;
all_robustness_resistant(all_robustness_resistant>0)=0;


robust_plasmid_pre_all = squeeze(all_robustness_all(1, :, :, :, 2));
robust_plasmid_post_all = squeeze(all_robustness_all(2, :, :, :, 2));
diff_robust_all = (abs(robust_plasmid_pre_all) - abs(robust_plasmid_post_all));%./abs(robust_plasmid_pre_all);
diff_robust_all = mean(diff_robust_all, 3);

robust_plasmid_pre_subset = squeeze(all_robustness_subset(1, :, :, :, 2));
robust_plasmid_post_subset = squeeze(all_robustness_subset(2, :, :, :, 2));
diff_robust_subset = (abs(robust_plasmid_pre_subset) - abs(robust_plasmid_post_subset));%./abs(robust_plasmid_pre_subset);
diff_robust_subset = mean(diff_robust_subset, 3);

robust_plasmid_pre_resistant = squeeze(all_robustness_resistant(1, :, :, :, 2));
robust_plasmid_post_resistant = squeeze(all_robustness_resistant(2, :, :, :, 2));
diff_robust_resistant = (abs(robust_plasmid_pre_resistant) - abs(robust_plasmid_post_resistant));%./abs(robust_plasmid_pre_resistant);
diff_robust_resistant = mean(diff_robust_resistant, 3);

my_clims = [-0.3,0.3];

subplot(2,3,4)
imagesc(diff_robust_all');
caxis(my_clims)
colormap(CT)
ylabel("% Positive interactions")
xlabel("Conjugation rate")

xticks(linspace(1,length(gamma_space),length(gamma_space)))
xticklabels(gamma_space)
yticks(linspace(1,length(PmSpace),length(PmSpace)))
yticklabels(PmSpace)
set(gca,'YDir','normal')

h = colorbar;
ylabel(h, 'Selection effect')
axis square

subplot(2,3,5)
imagesc(diff_robust_subset');

my_clims = [-0.3,0.3];
caxis(my_clims)
colormap(CT)
ylabel("% Positive interactions")
xlabel("Conjugation rate")

xticks(linspace(1,length(gamma_space),length(gamma_space)))
xticklabels(gamma_space)
yticks(linspace(1,length(PmSpace),length(PmSpace)))
yticklabels(PmSpace)
set(gca,'YDir','normal')

h = colorbar;
ylabel(h, 'Selection effect')
axis square

subplot(2,3,6)
imagesc(diff_robust_resistant');

hold on 
%contour(real(diff_robust_resistant'), [0,0], 'LineWidth',1, 'LineColor','k')



my_clims = [-0.05,0.05];
caxis(my_clims)
colormap(CT)
ylabel("% Positive interactions")
xlabel("Conjugation rate")

xticks(linspace(1,length(gamma_space),length(gamma_space)))
xticklabels(gamma_space)
yticks(linspace(1,length(PmSpace),length(PmSpace)))
yticklabels(PmSpace)
set(gca,'YDir','normal')

h = colorbar;
ylabel(h, 'Selection effect')
axis square


