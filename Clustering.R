#-------------------------------------------------------------------------------
## Reproducible & generalisable clustering analysis script
# Covers:
#   - Data simulation and preprocessing
#   - Optimal cluster number: elbow, silhouette, gap statistic, NbClust
#   - K-means clustering
#   - Hierarchical clustering (agglomerative + divisive)
#   - Soft/fuzzy clustering
#   - Model-based clustering
#   - Semi-supervised model-based clustering
#   - FlexMix (mixture of regression models)
#   - DBSCAN (density-based)
#   - Spectral clustering
#   - UMAP (dimensionality reduction + cluster overlay)
#   - t-SNE (for comparison)
#   - Cluster validation and characterization
#   - Reusable pipeline function
#-------------------------------------------------------------------------------

#--------------------------------------------
## Step 1: Setup
rm(list = ls())
set.seed(123)

required_pkgs<-c("MASS", "dplyr", "tidyr", "ggplot2", "factoextra","cluster", "NbClust", "mclust", "flexmix",
                   "dbscan", "kernlab", "umap", "Rtsne","fclust", "e1071", "pheatmap", "RColorBrewer","ggrepel", "patchwork", "tibble", "corrplot")

is_installed<-required_pkgs %in% rownames(installed.packages(all.available=TRUE))
if(any(is_installed == FALSE)){
  install.packages(required_pkgs[!is_installed],repos = "http://cran.us.r-project.org")
}
invisible(lapply(required_pkgs, library, character.only = TRUE))

#--------------------------------------------
## Step 2: Simulating multi-cluster dataset

# 4 well-separated clusters in 8-dimensional space (mimics an omics / environmental mixture dataset)
n_per<-120
K_true<-4

centres<-list(c(0,0,1,1,-1,-1,0,0),
  c(3,3,0,0,1,1,2,2),
  c(-3,2,2,-2,0,1,-1,1),
  c(1,-3,-1,2,2,-1,1,-2))

sim_list<-lapply(seq_len(K_true), function(k){
  X <- MASS::mvrnorm(n_per, mu = centres[[k]],Sigma = diag(0.8, 8) + matrix(0.1, 8, 8))
  as.data.frame(X)
})

sim_df<-do.call(rbind, sim_list)
colnames(sim_df)<-paste0("V", seq_len(8))
true_labels<-rep(seq_len(K_true), each = n_per)
sim_df$true_k<-factor(true_labels)

# Scale features (essential for most clustering methods)
X_scaled<-scale(sim_df[, paste0("V", 1:8)])

#--------------------------------------------
## Step 3: EDA

# Correlation matrix of features
corrplot::corrplot(cor(X_scaled), method = "color",title = "", mar = c(0,0,1,0))

# PCA overview
pca_res<-prcomp(X_scaled, scale. = FALSE)
pca_df<-as.data.frame(pca_res$x[, 1:2])
pca_df$true_k<-sim_df$true_k

ggplot(pca_df, aes(x = PC1, y = PC2, color = true_k)) +
  geom_point(alpha = 0.6, size = 2) +
  labs(title = "PCA: true cluster structure",
       color = "True cluster") +
  theme_bw(base_size = 14)

#--------------------------------------------
## Determining the optimal number of clusters

# Approach 1 - elbow method (within-cluster sum of squares)
factoextra::fviz_nbclust(X_scaled, FUNcluster = kmeans,method = "wss", k.max = 10) +
  labs(title = "Elbow method: WSS by K") +
  theme_bw()

# Approach 2 - average silhouette width
factoextra::fviz_nbclust(X_scaled, FUNcluster = kmeans,method = "silhouette", k.max = 10) +
  labs(title = "Silhouette method") +
  theme_bw()

# Approach 3 - gap statistic
set.seed(123)
gap_stat<-cluster::clusGap(X_scaled, FUN = kmeans,nstart = 25, K.max = 10, B = 50)
factoextra::fviz_gap_stat(gap_stat) +
  labs(title = "Gap statistic") +
  theme_bw()

cat("Gap statistic optimal K:", cluster::maxSE(gap_stat$Tab[,"gap"],gap_stat$Tab[,"SE.sim"],method="Tibs2001SEmax"), "\n")

# Approach 4 -  NbClust: uses majority vote across many indices

nb_res<-NbClust::NbClust(data= X_scaled,
  distance="euclidean",
  min.nc=2, max.nc = 8,
  method="kmeans",
  index="alllong")

cat("\nNbClust best K (majority vote):",nb_res$Best.nc["Number_clusters", ] %>% table() %>% which.max() %>% names(),"\n")

# Approach 5 -  BIC via mclust (for model-based)
mclust_bic<-mclust::mclustBIC(X_scaled, G = 1:10)
plot(mclust_bic, main = "mclust BIC by K and model")
cat("mclust BIC-optimal G:", substr(names(summary(mclust_bic))[1],5,5), "\n")

# Summary table
sil_vals<-sapply(2:10, function(k){
  km<-kmeans(X_scaled, centers = k, nstart = 25)
  s<-cluster::silhouette(km$cluster, dist(X_scaled))
  mean(s[, 3])
})
wss_vals<-sapply(2:10, function(k)
  kmeans(X_scaled, centers = k, nstart = 25)$tot.withinss)

opt_df<-data.frame(K = 2:10, WSS = wss_vals, Silhouette = sil_vals)
cat("\nElbow/Silhouette table\n")
round(opt_df, 4)

#--------------------------------------------
## k-means clustering

K_sel<-cluster::maxSE(gap_stat$Tab[,"gap"],gap_stat$Tab[,"SE.sim"],method="Tibs2001SEmax")   # set based on optimal K selection above

set.seed(123)
km_fit <- kmeans(X_scaled, centers = K_sel, nstart = 50, iter.max = 300)

cat("K-means cluster sizes:\n")
table(km_fit$cluster)
cat("Total WSS:", round(km_fit$tot.withinss, 2), "\n")
cat("Between/Total SS ratio:", round(km_fit$betweenss/km_fit$totss, 4), "\n")

# Silhouette
km_sil <- cluster::silhouette(km_fit$cluster, dist(X_scaled))
cat("Average silhouette width:", round(mean(km_sil[,3]), 4), "\n")
fviz_silhouette(km_sil) + theme_bw() + labs(title = "K-means: silhouette plot")

# Cluster plot (PCA space)
fviz_cluster(km_fit, data = X_scaled,
             palette = "jco", ellipse.type = "convex",
             repel = TRUE, ggtheme = theme_bw()) +
  labs(title = "K-means cluster plot (PCA space)")

# Agreement with true labels
km_ari <- mclust::adjustedRandIndex(km_fit$cluster, true_labels)
cat("K-means Adjusted Rand Index (ARI) vs. true labels:", round(km_ari, 4), "\n")

#--------------------------------------------
## Hierachical clustering

dist_mat<-dist(X_scaled, method = "euclidean")

# 1- Agglomerative (Ward.D2)
hclust_ward<-hclust(dist_mat, method = "ward.D2")

# Dendrogram
fviz_dend(hclust_ward, k = K_sel,
          cex = 0.4, palette = "jco",
          rect = TRUE, rect_fill = TRUE,
          main = "Hierarchical (Ward.D2): dendrogram") +
  theme_bw()

hc_labels<-cutree(hclust_ward, k = K_sel)
hc_sil<-cluster::silhouette(hc_labels, dist_mat)
cat("Hierarchical (Ward.D2) average silhouette:",round(mean(hc_sil[,3]), 4), "\n")
cat("Adjusted Rand Index (ARI) vs. true labels:",round(mclust::adjustedRandIndex(hc_labels, true_labels), 4), "\n")

# 2 - Comparing linkage methods
linkages<-c("ward.D2","complete","average","single")
link_sil<-sapply(linkages, function(l) {
  hc<-hclust(dist_mat, method = l)
  lab<-cutree(hc, k = K_sel)
  mean(cluster::silhouette(lab, dist_mat)[, 3])
})

cat("\nSilhouette by linkage\n")
round(sort(link_sil, decreasing = TRUE), 4)

# 3 - Divisive clustering (DIANA)
diana_fit<-cluster::diana(X_scaled, metric = "euclidean")
diana_lab<-cutree(as.hclust(diana_fit), k = K_sel)
cat("DIANA average silhouette:",round(mean(cluster::silhouette(diana_lab, dist_mat)[,3]), 4), "\n")

#--------------------------------------------
## Soft / fuzzy clustering

# fuzzy C-means
fcm_fit <- e1071::cmeans(X_scaled, centers = K_sel,
                         iter.max = 200, m = 2, # m=2 is standard fuzziness
                         method = "cmeans")

cat("Fuzzy C-means cluster sizes (hard assignment):\n")
table(fcm_fit$cluster)

# Membership matrix (first 6 rows)
cat("\nMembership probabilities (first 6 obs)\n")
round(head(fcm_fit$membership), 3)

# Distribution of maximum membership probability
max_memb <- apply(fcm_fit$membership, 1, max)
hist(max_memb, main = "Fuzzy C-means: max membership probability per observation",
     xlab = "Max membership", col = "#A6DDCE", breaks = 20)
abline(v = 0.5, lty = 2, col = "red")
cat("Observations with max membership < 0.6 (ambiguous):",sum(max_memb < 0.6), "\n")

# Fuzziness index (partition coefficient; 1 = crisp, 1/K = fully fuzzy)
pc <- sum(fcm_fit$membership^2) / nrow(X_scaled)
cat("Partition coefficient (PC):", round(pc, 4),"(closer to 1 = crisper clusters)\n")

# Hard assignment ARI
fcm_ari<-mclust::adjustedRandIndex(fcm_fit$cluster, true_labels)
cat("Fuzzy C-means (hard assign) Adjusted Rand Index (ARI):", round(fcm_ari, 4), "\n")

# Visualizing soft assignments in PCA space
pca_fcm<-as.data.frame(pca_res$x[, 1:2])
pca_fcm$hard_cluster<-factor(fcm_fit$cluster)
pca_fcm$max_memb<-max_memb

ggplot(pca_fcm, aes(x = PC1, y = PC2,color = hard_cluster, size = max_memb)) +
  geom_point(alpha = 0.6) +
  scale_size_continuous(range = c(0.5, 4),name = "Max membership") +
  labs(title = "Fuzzy C-means: PCA plot (size = certainty)", color = "Cluster") +
  theme_bw(base_size = 14)

#--------------------------------------------
## Model-based clustering

mclust_fit<-mclust::Mclust(X_scaled, G = K_sel)

cat("\nmclust: selected model:", mclust_fit$modelName, "\n")
cat("BIC:", round(mclust_fit$bic, 2), "\n")
cat("Cluster sizes:\n"); print(table(mclust_fit$classification))
cat("Adjusted Rand Index (ARI) vs. true labels:",round(mclust::adjustedRandIndex(mclust_fit$classification, true_labels), 4),"\n")

# Plotting mclust diagnostics
plot(mclust_fit, what = "BIC",  main = "mclust: BIC by model")
plot(mclust_fit, what = "classification",main = "mclust: classification (PC space)")
plot(mclust_fit, what = "uncertainty",main = "mclust: uncertainty")

# Uncertainty (complement of max posterior probability)
mclust_uncertainty <- 1 - apply(mclust_fit$z, 1, max)
cat("Mean uncertainty:", round(mean(mclust_uncertainty), 4), "\n")
hist(mclust_uncertainty, main = "mclust: classification uncertainty",
     xlab = "Uncertainty", col = "#A6DDCE", breaks = 20)

# Posterior probabilities (first 6 rows)
cat("\nPosterior probabilities (first 6 obs)\n")
round(head(mclust_fit$z), 3)

#--------------------------------------------
## Semi-supervized model-based clustering

# Scenario: 20% of observations have known labels; rest are unlabelled

known_idx<-sample(seq_len(nrow(X_scaled)),size = round(0.2 * nrow(X_scaled)))
class_labels<-rep(NA, nrow(X_scaled))
class_labels[known_idx]<-true_labels[known_idx]

ss_fit<-mclust::MclustSSC(X_scaled, class = class_labels, G = K_sel)

cat("\nSemi-supervised mclust results\n")
cat("Model:", ss_fit$modelName, "\n")
cat("Adjusted Rand Index (ARI) vs. true labels:",round(mclust::adjustedRandIndex(ss_fit$classification, true_labels), 4),"\n")
cat("Compare: unsupervised ARI =",round(mclust::adjustedRandIndex(mclust_fit$classification, true_labels), 4),"\n")

plot(ss_fit, what = "classification",main = "Semi-supervised mclust: classification")

#--------------------------------------------
## Flexmix (mixture of regression / latent class models): fits mixture models where each component 
## has its own regression. It is useful when clusters differ in predictor-outcome relationships

# Creating a response variable for illustration
flex_df<-as.data.frame(X_scaled)
flex_df$response <- 2 * flex_df$V1 - 1.5 * flex_df$V2 + rnorm(nrow(flex_df), 0, 1) + rep(c(0, 2, -2, 1), each = n_per)

# Fitting FlexMix: mixture of linear regressions
set.seed(123)
flex_fit<-flexmix::flexmix(response ~ V1 + V2 + V3,data=flex_df,k=K_sel)

cat("\nFlexMix summary\n")
summary(flex_fit)
cat("Log-likelihood:", round(logLik(flex_fit), 2), "\n")
cat("BIC:", round(BIC(flex_fit), 2), "\n")
cat("Cluster sizes:\n"); print(table(flexmix::clusters(flex_fit)))
cat("Adjusted Rand Index (ARI):", round(mclust::adjustedRandIndex(flexmix::clusters(flex_fit), true_labels), 4), "\n")

# Parameters per component
cat("\nComponent-specific regression parameters\n")
flexmix::parameters(flex_fit)

# Posterior probabilities
flex_post<-flexmix::posterior(flex_fit)
cat("Mean max posterior:", round(mean(apply(flex_post, 1, max)), 4), "\n")

# BIC-based K selection for FlexMix
flex_bic<-sapply(2:6, function(k) {
  f<-tryCatch(
    flexmix::flexmix(response ~ V1 + V2 + V3, data = flex_df, k = k),
    error = function(e) NULL)
  if (is.null(f)) return(NA)
  BIC(f)
})
cat("\n-FlexMix BIC by K\n")
data.frame(K = 2:6, BIC = round(flex_bic, 2))

#--------------------------------------------
## DBSCAN (density-based; no K required)

# Advantages: finds arbitrary shapes, handles noise/outliers
# Key parameters: eps (neighbourhood radius), minPts (min neighbours)

# Estimating eps via k-nearest neighbour distance plot
dbscan::kNNdistplot(X_scaled, k = 5)
abline(h = 1.5, lty = 2, col = "red") # adjust based on plot

db_fit<-dbscan::dbscan(X_scaled, eps = 1.5, minPts = 8)

cat("DBSCAN cluster sizes (0 = noise):\n")
table(db_fit$cluster)
cat("Noise points:", sum(db_fit$cluster == 0), "\n")

# Silhouette (excluding noise)
non_noise <- db_fit$cluster != 0
if(length(unique(db_fit$cluster[non_noise])) > 1){
  db_sil<-cluster::silhouette(db_fit$cluster[non_noise],dist(X_scaled[non_noise, ]))
  cat("DBSCAN average silhouette (non-noise):",round(mean(db_sil[,3]), 4), "\n")
}

# Visualizing
pca_db<-as.data.frame(pca_res$x[, 1:2])
pca_db$cluster<-factor(db_fit$cluster)
ggplot(pca_db, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("0"="grey70", setNames(RColorBrewer::brewer.pal(8,"Dark2"),as.character(1:8))[seq_len(max(db_fit$cluster))])) +
  labs(title = "DBSCAN: cluster assignments (grey = noise)",color = "Cluster") +
  theme_bw(base_size = 14)

#--------------------------------------------
## Special clustering: for non-convex clusters; uses graph Laplacian eigen decomposition

spec_fit<-kernlab::specc(X_scaled, centers = K_sel)
spec_lab<-as.integer(spec_fit)

cat("Spectral cluster sizes:\n"); print(table(spec_lab))
spec_sil  <- cluster::silhouette(spec_lab, dist_mat)
cat("Spectral average silhouette:", round(mean(spec_sil[,3]), 4), "\n")
cat("Adjusted Rand Index (ARI) vs. true labels:",round(mclust::adjustedRandIndex(spec_lab, true_labels), 4), "\n")

#--------------------------------------------
## UMAP: non-linear dimensionality reduction preserving local structure
# Not a clustering method itself, but used to visualize cluster structure

set.seed(123)
umap_config<-umap::umap.defaults
umap_config$n_neighbors<-15
umap_config$min_dist<-0.1
umap_config$n_components<-2

umap_fit<-umap::umap(X_scaled, config = umap_config)
umap_df<-as.data.frame(umap_fit$layout)
colnames(umap_df)<-c("UMAP1","UMAP2")
umap_df$true_k<-sim_df$true_k
umap_df$kmeans_k<-factor(km_fit$cluster)
umap_df$mclust_k<-factor(mclust_fit$classification)

# True labels on UMAP
p_umap_true<-ggplot(umap_df, aes(x=UMAP1, y=UMAP2, color=true_k)) +
  geom_point(alpha=0.6, size=2) +
  labs(title="UMAP: true clusters", color="True K") +
  theme_bw(base_size=13)

# K-means labels on UMAP
p_umap_km<-ggplot(umap_df, aes(x=UMAP1, y=UMAP2, color=kmeans_k)) +
  geom_point(alpha=0.6, size=2) +
  labs(title="UMAP: K-means labels", color="K-means") +
  theme_bw(base_size=13)

# mclust labels on UMAP
p_umap_mc<-ggplot(umap_df, aes(x=UMAP1, y=UMAP2, color=mclust_k)) +
  geom_point(alpha=0.6, size=2) +
  labs(title="UMAP: mclust labels", color="mclust") +
  theme_bw(base_size=13)

p_umap_true + p_umap_km + p_umap_mc

# UMAP with fuzzy membership overlaid
umap_df$max_memb<-max_memb
ggplot(umap_df, aes(x=UMAP1, y=UMAP2,color=factor(fcm_fit$cluster), size=max_memb)) +
  geom_point(alpha=0.6) +
  scale_size_continuous(range=c(0.5,4)) +
  labs(title="UMAP: fuzzy C-means membership certainty",color="Cluster", size="Max membership") +
  theme_bw(base_size=14)

#--------------------------------------------
## t-SNE (comparison with UMAP)

set.seed(123)
tsne_fit<-Rtsne::Rtsne(X_scaled, dims=2, perplexity=30,check_duplicates=FALSE, verbose=FALSE)
tsne_df<-as.data.frame(tsne_fit$Y)
colnames(tsne_df)<-c("tSNE1","tSNE2")
tsne_df$true_k<-sim_df$true_k

ggplot(tsne_df, aes(x=tSNE1, y=tSNE2, color=true_k)) +
  geom_point(alpha=0.6, size=2) +
  labs(title="t-SNE: true cluster structure", color="True K") +
  theme_bw(base_size=14)

#--------------------------------------------
## Cluster validation and comparison

# Collecting all hard assignments
all_labels<-data.frame(true=true_labels,
  kmeans=km_fit$cluster,
  hclust=hc_labels,
  diana=diana_lab,
  fcm_hard=fcm_fit$cluster,
  mclust=mclust_fit$classification,
  mclust_ss=ss_fit$classification,
  flexmix=flexmix::clusters(flex_fit),
  spectral=spec_lab)

# ARI against true labels
ari_df <- data.frame(method = names(all_labels)[-1],
  ARI= sapply(names(all_labels)[-1], function(m)
    round(mclust::adjustedRandIndex(all_labels[[m]], true_labels), 4)))
cat("\nAdjusted Rand Index (ARI) vs. true labels\n")
ari_df[order(-ari_df$ARI), ]

# Average silhouette widths
sil_df <- data.frame(method=c("kmeans","hclust","diana","fcm","mclust","spectral"),
  avg_sil=c(mean(cluster::silhouette(km_fit$cluster,dist_mat)[,3]),
    mean(cluster::silhouette(hc_labels,dist_mat)[,3]),
    mean(cluster::silhouette(diana_lab,dist_mat)[,3]),
    mean(cluster::silhouette(fcm_fit$cluster,dist_mat)[,3]),
    mean(cluster::silhouette(mclust_fit$classification, dist_mat)[,3]),
    mean(cluster::silhouette(spec_lab,dist_mat)[,3])))
cat("\nAverage silhouette widths\n")
round(sil_df[order(-sil_df$avg_sil),"avg_sil" ], 4)

ggplot(sil_df, aes(x = reorder(method, avg_sil), y = avg_sil)) +
  geom_col(fill = "#A6DDCE", color = "black") +
  coord_flip() +
  labs(title = "Cluster validation: average silhouette by method",x = "", y = "Average silhouette width") +
  theme_bw(base_size = 14)

#--------------------------------------------
## Cluster characterization

best_labels<-km_fit$cluster # replace with your best method
char_df<-as.data.frame(X_scaled)
char_df$cluster<-factor(best_labels)

# Mean feature values per cluster
cluster_means<-char_df %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean)) %>%
  tidyr::pivot_longer(-cluster, names_to="feature", values_to="mean")

ggplot(cluster_means, aes(x=feature, y=mean, fill=cluster)) +
  geom_col(position="dodge", color="black") +
  labs(title="Mean feature values per cluster (scaled)",x="Feature", y="Mean (scaled)") +
  theme_bw(base_size=13) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

# Heatmap of cluster means
means_mat<-char_df %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean)) %>%
  tibble::column_to_rownames("cluster") %>%
  as.matrix()

pheatmap::pheatmap(t(means_mat),
                   cluster_cols = TRUE, cluster_rows = TRUE,
                   main = "Cluster mean profiles",
                   color = colorRampPalette(
                     rev(RColorBrewer::brewer.pal(9,"RdBu")))(50))

# ANOVA: is each feature different across clusters?
anova_res <- sapply(paste0("V",1:8), function(v) {
  summary(aov(as.formula(paste(v,"~ cluster")), data=char_df))[[1]][1,"Pr(>F)"]
})
cat("\nNOVA p-values: feature differences across clusters\n")
sort(anova_res)

#--------------------------------------------
## Reusable pipeline

run_clustering_pipeline<-function(data, K = NULL, K_max = 10, scale_data = TRUE,
                                    methods = c("kmeans","hclust","mclust","fuzzy"),
                                    seed = 123){
  set.seed(seed)
  X<-if(scale_data) scale(data) else as.matrix(data)
  
  # Determine K if not provided (silhouette)
  if(is.null(K)){
    sil_v <- sapply(2:K_max, function(k)
      mean(cluster::silhouette(kmeans(X, centers=k, nstart=25)$cluster, dist(X))[,3]))
    K<-(2:K_max)[which.max(sil_v)]
    message("Auto-selected K = ", K, " (max silhouette)")
  }
  
  res<-list(K = K)
  d_mat<-dist(X)
  
  if("kmeans"%in% methods){
    fit<- kmeans(X, centers=K, nstart=50, iter.max=300)
    res$kmeans<-fit$cluster
    res$km_sil<-round(mean(cluster::silhouette(fit$cluster, d_mat)[,3]), 4)
    message("K-means average silhouette: ", res$km_sil)
  }
  if("hclust"%in% methods){
    hc<-hclust(d_mat, method="ward.D2")
    res$hclust<-cutree(hc, k=K)
    res$hc_sil<-round(mean(cluster::silhouette(res$hclust, d_mat)[,3]), 4)
    message("Hierarchical average silhouette: ", res$hc_sil)
  }
  if("mclust"%in% methods){
    mc<-mclust::Mclust(X, G=K)
    res$mclust<-mc$classification
    res$mc_bic<-round(mc$bic, 2)
    message("mclust BIC: ", res$mc_bic)
  }
  if("fuzzy" %in% methods){
    fc<-e1071::cmeans(X, centers=K, m=2, iter.max=200)
    res$fuzzy<-fc$cluster
    res$fc_pc<-round(sum(fc$membership^2)/nrow(X), 4)
    message("Fuzzy PC: ", res$fc_pc)
  }
  
  # UMAP visualization
  umap_r<-umap::umap(X)
  umap_df<-as.data.frame(umap_r$layout)
  colnames(umap_df)<-c("UMAP1","UMAP2")
  
  if("kmeans" %in% methods) {
    umap_df$cluster<-factor(res$kmeans)
    p<-ggplot(umap_df, aes(UMAP1, UMAP2, color=cluster)) +
      geom_point(alpha=0.6) +
      labs(title=paste("UMAP — K-means (K=",K,")",sep="")) +
      theme_bw()
    print(p)
  }
  
  return(res)
}

pipe_res<-run_clustering_pipeline(data=sim_df[, paste0("V",1:8)],K=4,methods=c("kmeans","hclust","mclust","fuzzy"))
