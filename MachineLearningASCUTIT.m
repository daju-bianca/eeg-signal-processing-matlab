%% CLUSTERING NESUPERVIZAT - K-MEANS EEG ASCUTIT
% Sursa 1: MathWorks. Cluster Analysis Example
%Sursa 2:AI
cale_folder_ascutit = 'C:\Users\Bianca\Documents\CercetareIRA\Stimul_ascutit (1)\';
fisier_date = fullfile(cale_folder_ascutit, 'REZULTATE_Entropie_ASCUTIT.csv');
tabel_date = readtable(fisier_date);

canale_eeg = {'Fp1','Fp2','F3','F4','C3','C4','P3','P4','O1','O2','Fz','Pz'};
X = tabel_date{:, canale_eeg};
fprintf('Date incarcate: %d subiecti x %d canale\n', size(X,1), size(X,2));
% 2. NORMALIZARE
X_norm = normalize(X);
% 3. PCA -
[coeff, score, ~, ~, explained] = pca(X_norm);
fprintf('PC1: %.1f%% | PC2: %.1f%% | Total: %.1f%%\n', ...
        explained(1), explained(2), explained(1)+explained(2));

% 4. ELBOW METHOD 
figure('Name','Elbow Method - ASCUTIT','Position',[100 100 500 400]);
sumd_vec = zeros(1,6);
for k = 1:6
    [~,~,sumd] = kmeans(X_norm, k, 'Replicates', 10, 'Display', 'off');
    sumd_vec(k) = sum(sumd);
end
plot(1:6, sumd_vec, '-o', 'LineWidth', 2, 'MarkerFaceColor','b', 'MarkerSize', 8);
xlabel('Numar clustere (k)');
ylabel('Suma distante intra-cluster');
title('Elbow Method - Numar optim clustere (ASCUTIT)');
xline(2,'--r','k=2','LineWidth',1.5);
grid on;

% 5. K-MEANS k=2
rng('default');  
rng(42);       
[idx_clustere, centroizi, sumdist] = kmeans(X_norm, 2, ...
    'Distance',   'sqeuclidean', ...
    'Replicates', 10, ...
    'Display',    'off');

fprintf('\n=== REZULTATE CLUSTERING ===\n');
fprintf('Cluster 1: %d subiecti\n', sum(idx_clustere==1));
fprintf('Cluster 2: %d subiecti\n', sum(idx_clustere==2));

% 6. SILHOUETTE SCORE
s = silhouette(X_norm, idx_clustere);
fprintf('Silhouette Score mediu: %.4f\n', mean(s));

% 7. PURITATE CLUSTERE
media_c1 = mean(X(idx_clustere==1, :));
media_c2 = mean(X(idx_clustere==2, :));
fprintf('\nMedia SampEn per cluster:\n');
fprintf('%-6s | %-12s | %-12s | %s\n', 'Canal','Cluster 1','Cluster 2','Diferenta');
fprintf('%s\n', repmat('-',1,45));
for j = 1:length(canale_eeg)
    fprintf('%-6s | %-12.4f | %-12.4f | %.4f\n', ...
            canale_eeg{j}, media_c1(j), media_c2(j), media_c1(j)-media_c2(j));
end
figure('Name','Clustering Nesupervizat - ASCUTIT', ...
       'Position',[50 50 1400 500]);

% Subplot 1 - PCA
subplot(1,3,1);
gscatter(score(:,1), score(:,2), idx_clustere, ...
         [0.2 0.5 0.8; 0.8 0.3 0.2], 'os', 8);
hold on;
centroizi_pca = centroizi * coeff(:,1:2);
plot(centroizi_pca(:,1), centroizi_pca(:,2), 'kx', ...
     'MarkerSize',15, 'LineWidth',3);
hold off;
xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
title('Clustere K-Means pe PCA (ASCUTIT)');
legend({'Cluster 1','Cluster 2','Centroizi'}, 'Location','best','FontSize',8);
grid on;

% Subplot 2 - Media SampEn 
subplot(1,3,2);
bar([media_c1; media_c2]', 'grouped');
set(gca,'XTickLabel', canale_eeg, 'XTickLabelRotation', 45);
legend({'Cluster 1','Cluster 2'}, 'Location','northeast');
ylabel('SampEn mediu');
title('SampEn mediu per canal (ASCUTIT)');
grid on;

% Subplot 3 - Silhouette
subplot(1,3,3);
silhouette(X_norm, idx_clustere);
title(sprintf('Silhouette Score = %.3f', mean(s)));
sgtitle('Clustering Nesupervizat K-Means - SampEn EEG ASCUTIT', ...
        'FontSize',13, 'FontWeight','bold');
saveas(gcf, fullfile(cale_folder_ascutit, 'CLUSTERING_NESUPERVIZAT_ASCUTIT1.png'));
SampEn_Global = mean(X, 2); 

TabelClustering = tabel_date(:, {'Fisier'});
TabelClustering.SampEn_Mediu = SampEn_Global; 
TabelClustering.Cluster = idx_clustere;
TabelClustering.Cluster_Label = repmat({''}, height(TabelClustering), 1);
TabelClustering.Cluster_Label(idx_clustere==1) = {'Posibil_Sanatos'};
TabelClustering.Cluster_Label(idx_clustere==2) = {'Posibil_Bolnav'};
% Salvăm tabelul în format CSV
writetable(TabelClustering, fullfile(cale_folder_ascutit, 'REZULTATE_Entropie_ASCUTIT.csv'));
disp(TabelClustering);