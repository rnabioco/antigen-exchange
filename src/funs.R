
# Utility ----

#' Save Seurat object and meta.data
#' 
#' @param ob_in Seurat object to save.
#' @param prfx Prefix to use for saved files. If set to NULL, the name of the
#' object is used.
#' @param ob_dir Directory to save files. The default is the current project
#' directory.
#' @export
save_objs <- function(ob_in, prfx = NULL, ob_dir = "") {
  
  if (is.null(prfx)) {
    prfx <- deparse(substitute(ob_in))
  }
  
  ob_in %>%
    qs::qsave(file.path(ob_dir, str_c(prfx, ".qs")))
  
  if ("Seurat" %in% class(ob_in)) {
    ob_in@meta.data %>%
      tibble::as_tibble(rownames = "cell_id") %>%
      vroom::vroom_write(file.path(ob_dir, str_c(prfx, ".tsv.gz")))
  }
}

#' Plot color palette
#' 
#' @param cols_in Character vector containing colors to plot
#' @param ... Additional arguments to pass to ggplot
#' @return ggplot object
#' @export
plot_color_palette <- function(cols_in, ...) {
  
  if (is.null(names(cols_in))) {
    cols_in <- set_names(cols_in)
  }
  
  col_df <- tibble::tibble(color = factor(names(cols_in), names(cols_in)))
  
  res <- ggplot2::ggplot(col_df, aes(x = "color", fill = color)) +
    ggplot2::geom_bar(...) +
    ggplot2::scale_fill_manual(values = cols_in) +
    ggplot2::theme_void()
  
  res
}

#' Create function to generate color gradient
#' 
#' @param cols_in Character vector containing colors
#' @return Color gradient function
#' @export
create_col_fun <- function(cols_in) {
  
  create_gradient <- function(cols_in, n = NULL) {
    
    if (is.null(n)) {
      n <- length(cols_in)
    }
    
    colorRampPalette(cols_in)(n)
  }
  
  function(n = NULL) {
    create_gradient(cols_in, n)
  }
}

#' Capitalize first character without modifying other characters
#' 
#' @param string Charactor string to modify
#' @param exclude Word to exclude from output
#' @return Character string
#' @export
str_to_title_v2 <- function(string, exclude = "cell") {
  
  cap_first_chr <- function(string, exclude) {
    chrs <- string %>%
      str_split(pattern = "") %>%
      unlist()
    
    if (any(chrs %in% LETTERS) || string == exclude) {
      return(string)
    }
    
    chrs[1] <- str_to_upper(chrs[1])
    
    res <- chrs %>%
      reduce(str_c)
    
    res
  }
  
  res <- string %>%
    map_chr(~ {
      .x %>%
        str_split(pattern = " ") %>%
        unlist() %>%
        map_chr(cap_first_chr, exclude = exclude) %>%
        reduce(str_c, sep = " ")
    })
  
  res
}

#' Create list of cluster column names
#' 
#' @param rslns Clustering resolutions.
#' @param prefix Prefix to add to column names.
#' @return list of column names.
#' @export
set_clst_names <- function(rslns, prefix = "res_") {
  
  res <- rslns %>%
    set_names() %>%
    map(~ list(
      clst = str_c("RNA_snn_res.", .x),
      type = str_c(prefix, .x)
    ))
  
  res
}

#' Calculate a pseudo count for a given vector
#' 
#' @param x Vector of values to use for calculation
#' @param frac The pseudo count is calculated by multiplying the smallest
#' non-zero value by frac.
#' @return Values with pseudo count added
#' @export
add_pseudo <- function(x, frac = 0.5) {
  
  pseudo <- min(x[x > 0]) * frac
  
  x + pseudo
}

#' Identify cell types represented in all samples
#' 
#' @param so_in Seurat object
#' @param type_clmn meta.data column containing cell types
#' @param sample_clmn meta.data column containing sample names
#' @param n_cells Minimum number of cells
#' @return Vector with cell types that are represented by > n_cells for all
#' samples in sample_clmn
#' @export
get_cell_types <- function(so_in, type_clmn, sample_clmn, n_cells = 3) {
  meta <- so_in@meta.data
  
  sams <- meta %>%
    pull(sample_clmn) %>%
    unique()
  
  res <- meta %>%
    group_by(!!sym(type_clmn), !!sym(sample_clmn)) %>%
    summarize(n = n(), .groups = "drop") %>%
    filter(n > n_cells) %>%
    group_by(!!sym(type_clmn)) %>%
    filter(all(sams %in% !!sym(sample_clmn))) %>%
    pull(type_clmn) %>%
    unique()
  
  res
}

#' Format p values for labels
#' 
#' modified from djvdj
.format_pvalue <- function(p, digits = 1, cutoffs = NULL, show_decimal = 0.1) {
  
  # Set p label based on vector of cutoffs
  if (!is.finite(p)) return(as.character(NA))
  
  if (!is.null(cutoffs)) {
    if (any(duplicated(cutoffs))) {
      cli::cli_abort("Cutoff values for p_label must be unique.")
    }
    
    # Set default labels when not provided by user
    if (is.null(names(cutoffs))) {
      cutoffs <- sort(cutoffs, decreasing = TRUE)
      
      names(cutoffs) <- purrr::imap_chr(
        cutoffs, ~ paste0(rep("*", .y), collapse = "")
      )
    }
    
    cutoffs <- sort(cutoffs)
    p_label <- as.character(NA)
    
    for (val in names(cutoffs)) {
      if (p < cutoffs[val]) {
        p_label <- val
        
        break()
      }
    }
    
    # Treat "value" as a keyword that will allow user to display actual
    # p-value for a certain cutoff
    # All custom labels need to be wrapped in quotes for parsing
    if (!identical(p_label, "value")) {
      if (!is.na(p_label)) p_label <- paste0("\'", p_label, "\'")
      
      return(p_label)
    }
  }
  
  # Format p-value label
  # label_scientific will round 0.095 to 0.1 when digits = 1
  if (round(p, digits + 1) >= show_decimal) return(as.character(round(p, 1)))
  
  p <- scales::label_scientific(digits = digits)(p)
  
  ex <- str_extract_all(p, "[+\\-][0-9]+$")
  
  p <- sub(paste0("\\", ex, "$"), "", p)
  
  ex <- as.numeric(ex)
  ex <- as.character(ex)
  
  p <- sub("e", "*x*10^", p)
  p <- paste0(p, ex)
  
  p
}

# Processing helpers ----

#' Classify cell types based on cluster mean expression
#' 
#' @param so_in Seurat object.
#' @param feats List of features to use for classifying clusters
#' @param filt Expression to use for filtering clusters, e.g. Cd3e < 0.1
#' @param type Cell type label to use for cells identified by filtering
#' expression
#' @param clst_col meta.data column containing cell clusters to use for
#' calculating mean expression.
#' @param type_col meta.data column to add cell type label.
#' @param summary_fn Function to use for summarizing marker gene expression.
#' @return Seurat object containing new cell type classifications.
#' @export
classify_markers <- function(so_in, feats, filt, type_label, clst_col, type_col,
                             summary_fn = mean) {
  clsts <- so_in %>%
    FetchData(c(feats, clst_col, type_col)) %>%
    group_by(!!sym(clst_col)) %>%
    summarize(across(all_of(feats), summary_fn), .groups = "drop") %>%
    filter({{filt}}) %>%
    pull(clst_col) %>%
    as.character()
  
  res <- so_in %>%
    mutate_meta(
      mutate,
      !!sym(type_col) := ifelse(
        !!sym(clst_col) %in% clsts,
        type_label,
        !!sym(type_col)
      )
    )
  
  res
}

#' Generate shortened project names from matrix path
#' 
#' @param str_in Input string
#' @param extra_path Portion of path to remove from str_in
#' @return Shortened name
#' @export
shorten_names <- function(str_in, extra_path = "/outs/(filtered|raw)_feature_bc_matrix$") {
  res <- str_in %>%
    str_remove(extra_path) %>%
    basename() %>%
    str_extract("^[a-zA-Z0-9_]+") %>%
    str_remove("^GEX_")
  
  res
}

#' Wrapper to create Seurat object
#' 
#' @param mat_dir Directory containing matrix generated by Cell Ranger.
#' @param proj_name Project name to include in meta.data table.
#' @param hash_ids Name of cell hashing antibodies included in matrix.
#' @param adt_count_min If CITE-seq was performed, this option will remove
#' antibodies where the sum total counts is less than adt_count_min.
#' @param gene_min Minimum number of detected genes for cell.
#' @param gene_max Maximum number of detected genes for cell.
#' @param mito_max Maximum percentage of mitochondrial reads for cell.
#' @param mt_str String to use for identifying mitochondrial genes.
#' @param rna_assay Name of RNA assay if multiple assays are being added to the
#' object (e.g. if CITE-seq data is included).
#' @param adt_assay Name of ADT assay for Seurat object.
#' @return Seurat object
#' @export
create_sobj <- function(mat_dir, proj_name = "SeuratProject", hash_ids = NULL, adt_count_min = 0,
                        gene_min = 250, gene_max = 5000, mito_max = 20, mt_str = "^mt-",
                        rna_assay = "Gene Expression", adt_assay = "Antibody Capture") {
  
  # Load matrices
  mat_list <- Seurat::Read10X(mat_dir)
  rna_mat  <- mat_list
  
  # Create Seurat object using gene expression data
  if (is_list(mat_list)) {
    rna_mat <- mat_list[[rna_assay]]
  }
  
  res <- rna_mat %>%
    Seurat::CreateSeuratObject(
      project   = proj_name,
      min.cells = 5
    )
  
  # Add antibody capture data to Seurat object
  if (is_list(mat_list)) {
    adt_mat <- mat_list[[adt_assay]]
    
    # Double check that cells match for both assays
    if (!identical(colnames(res), colnames(adt_mat))) {
      adt_mat <- adt_mat[, colnames(res)]
      
      warning("Not all cells are shared between RNA and ADT assays.")
    }
    
    # Remove ADT features that have low total counts and likely failed or
    # were omitted
    n_feats    <- nrow(adt_mat)
    count_sums <- rowSums(as.matrix(adt_mat))
    
    adt_mat <- adt_mat[count_sums >= adt_count_min, ]
    
    if (n_feats != nrow(adt_mat)) {
      warning("Some ADT features were removed due to low counts (<", adt_count_min, ").")
    }
    
    res[["ADT"]] <- Seurat::CreateAssayObject(adt_mat)
  }
  
  # Calculate percentage of mitochondrial reads
  res <- res %>%
    Seurat::PercentageFeatureSet(
      pattern  = mt_str, 
      col.name = "pct_mito"
    )
  
  # Add QC classifications to meta.data
  res <- res %>%
    mutate_meta(
      mutate,
      qc_class = case_when(
        pct_mito     > mito_max ~ "high_mito_reads",
        nFeature_RNA > gene_max ~ "high_gene_count",
        nFeature_RNA < gene_min ~ "low_gene_count",
        TRUE ~ "pass"
      )
    )
  
  res
}

#' Create Seurat object with separate assay for virus counts
#' 
#' @param mat_dir Directory containing matrix generated by Cell Ranger.
#' @param proj_name Project name to include in meta.data table.
#' @param gene_min Minimum number of detected genes for cell to be included in
#' the object.
#' @param gene_max Maximum number of detected genes for cell to be included in
#' the object.
#' @param mito_max Maximum percentage of mitochondrial reads for a cell to be
#' included in the object.
#' @param cells_min Minimum number of cells a feature must be detected in to be
#' included in the object.
#' @param virus_min A cell with more viral reads than virus_min will be
#' included in the object regardless of other QC cutoffs with the exception of
#' gene_max.
#' @param virus_str String to use for identifying viral genes.
#' @param virus_assay Name of assay to use for viral counts.
#' @return Seurat object
#' @param mt_str String to use for identifying mitochondrial genes.
#' @export
create_virus_obj <- function(mat_dir, proj_name = "SeuratProject", gene_min = 250, gene_max = 5000,
                             mito_max = 15, cells_min = 5, virus_min = Inf, virus_str = "^CHIKV",
                             virus_assay = "CHIKV", mt_str = "^mt-") {
  
  rna_assay <- "RNA"

  # Load matrix
  mat <- Seurat::Read10X(mat_dir)
  
  v_feats <- grep(virus_str, rownames(mat), value = TRUE)
  
  # Create Seurat object containing expression data
  rna_mat <- mat[!rownames(mat) %in% v_feats, ]
  
  res <- rna_mat %>%
    Seurat::CreateSeuratObject(
      project   = proj_name,
      assay     = rna_assay,
      min.cells = cells_min
    )
  
  # Add virus assay
  # set drop = FALSE to keep as matrix when there is only one feature
  v_mat <- mat[v_feats, , drop = FALSE]
  
  rownames(v_mat) <- rownames(v_mat) %>%
    str_replace_all("_", "-")
  
  res[[virus_assay]] <- Seurat::CreateAssayObject(v_mat)
  
  # Add mito percentage to meta.data
  res <- res %>%
    Seurat::PercentageFeatureSet(
      assay    = rna_assay, 
      pattern  = mt_str, 
      col.name = "pct_mito"
    )
  
  # Add viral counts to meta.data
  names(v_feats) <- v_feats %>%
    str_replace_all("_", "-") %>%
    str_c(str_to_lower(virus_assay), "_", .)
  
  res <- res %>%
    AddMetaData(
      metadata = FetchData(., names(v_feats)),
      col.name = v_feats
    )
  
  # Add QC classifications to meta.data
  # Some viruses (e.g. CHIKV) inhibit host transcription by degrading the large
  # subunit of RNA polymerase II, so cells where the virus is replicating may
  # have low host mRNA counts and/or a high percentage of mitochondrial reads
  v_counts   <- sym(str_c("nCount_", virus_assay))
  v_pct      <- sym(str_c("pct_", virus_assay))
  rna_counts <- sym(str_c("nCount_", rna_assay))
  
  res <- res %>%
    djvdj::mutate_meta(
      mutate,
      qc_class = case_when(
        nFeature_RNA  > gene_max  ~ "high_gene_count",
        !!v_counts    > virus_min ~ "pass",
        pct_mito      > mito_max  ~ "high_mito_reads",
        nFeature_RNA  < gene_min  ~ "low_gene_count",
        TRUE ~ "pass"
      ),
      
      !!v_pct := !!v_counts / (!!v_counts + !!rna_counts)
    )
  
  res
}

#' Wrapper to normalize and scale Seurat object
#' 
#' @param sobj_in Seurat object.
#' @param rna_assay Name of RNA assay in object.
#' @param adt_assay Name of ADT assay in object.
#' @param cc_scoring Score cell cycle genes using cc.genes included in Seurat.
#' @param regress_vars Variables to regress out when scaling data.
#' @param rna_method Method to use with NormalizeData for RNA assay.
#' @param adt_method Method to use with NormalizeData for ADT assay.
#' @param scale_data Scale data after normalization.
#' @return Seurat object
#' @export
norm_sobj <- function(sobj_in, rna_assay = "RNA", adt_assay = "ADT", cc_scoring = FALSE,
                      regress_vars = NULL, rna_method = "LogNormalize", adt_method = "CLR",
                      scale_data = TRUE) {
  
  # Normalize counts
  res <- sobj_in %>%
    Seurat::NormalizeData(
      assay                = rna_assay,
      normalization.method = rna_method
    )
  
  # Score cell cycle genes
  if (cc_scoring) {
    s.genes <- cc.genes$s.genes %>%
      str_to_title()
    
    g2m.genes <- cc.genes$g2m.genes %>%
      str_to_title()
    
    res <- res %>%
      Seurat::CellCycleScoring(
        s.features   = s.genes,
        g2m.features = g2m.genes
      )
  }
  
  # Scale data
  # By default variable features will be used
  if (scale_data) {
    res <- res %>%
      Seurat::FindVariableFeatures(
        selection.method = "vst",
        nfeatures        = 2000
      ) %>%
      Seurat::ScaleData(vars.to.regress = regress_vars)
  }
  
  # Normalize ADT data
  if (!is.null(adt_assay) && adt_assay %in% names(res)) {
    res <- res %>%
      Seurat::NormalizeData(
        assay                = adt_assay,
        normalization.method = adt_method
      )
    
    if (scale_data) {
      res <- res %>%
        Seurat::ScaleData(assay = adt_assay)
    }
  }
  
  res
}

#' Wrapper to run doubletFinder
#' 
#' @param sobj_in Seurat object with low quality cells removed
#' @param assay Name of assay in object.
#' @param dbl_rate Expected doublet rate for experiment. If set to NULL,
#' expected doublet rate will be estimated based on number of captured cells.
#' @param gene_min Minimum number of detected genes to use for filtering low
#' quality cells
#' @param mito_max Maximum percentage of mitochondrial reads to use for
#' filtering low quality cells
#' @param prep Should pre-processing be performed on Seurat object
#' (normalization, clustering, etc.)
#' @param PCs Number of principal components to use
#' @param pN Number generated artificial doublets, expressed as a proportion of the merged real-artificial data
#' @param reuse.pANN meta.data column containing previously generated pANN
#' results.
#' @param rsln Resolution to use for clustering when prep = TRUE.
#' @param clust_column If prep = FALSE, meta.data column containing cell
#' clusters to use for estimating homotypic doublets
#' @param mito_clmn meta.data column containing the percentage of mitochondrial
#' counts to use for filtering low quality cells
#' @param ... Additional arguments to pass to doubletFinder_v3.
#' @return Seurat object with doublet classifications add to meta.data
#' @export
run_doubletFinder <- function(sobj_in, assay = "RNA", dbl_rate = NULL, mito_max = 15, gene_min = 250,
                              prep = TRUE, PCs = 1:40, pN = 0.25, reuse.pANN = FALSE, rsln = 1,
                              clust_column = "seurat_clusters", mito_clmn = "pct_mito", ...) {
  
  # Remove low quality cells
  res <- sobj_in %>%
    subset(!!sym(mito_clmn) <= mito_max & nFeature_RNA >= gene_min)
  
  # Preprocess object
  if (prep || (is.logical(reuse.pANN) && !reuse.pANN)) {
    clust_column <- "seurat_clusters"
    
    res <- sobj_in %>%
      norm_sobj(
        rna_assay  = assay,
        rna_method = "LogNormalize",
        scale_data = TRUE
      ) %>%
      cluster_RNA(
        assay      = assay,
        resolution = rsln,
        dims       = PCs
      )
  }
  
  # pK identification
  sweep_res   <- paramSweep_v3(res, PCs = PCs)
  sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
  bcmvn       <- find.pK(sweep_stats)
  
  pK <- bcmvn %>%
    filter(BCmetric == max(BCmetric)) %>%
    pull(pK) %>%
    as.character() %>%
    as.double()
  
  # Estimate expected doublet rate based on number of captured cells
  # assumes a capture rate of 57% and multiplet rate of 4.6e-6
  # assumptions are described here: https://satijalab.org/costpercell/
  if (is.null(dbl_rate)) {
    dbl_rate <- 0.0000046
    cap_rate <- 0.57
    
    dbl_rate <- (ncol(sobj_in) / cap_rate) * dbl_rate
  }
  
  # Homotypic doublet proportion estimate
  # homotypic doublets involve cells with the same transcriptional profile and
  # will not be detected by doubletFinder
  # homotypic doublet rate is estimated based on the distribution of cells
  # within cell clusters
  clsts       <- pull(res@meta.data, clust_column)
  htypic_prop <- modelHomotypic(clsts)
  nExp        <- round(dbl_rate * length(Cells(res)))
  nExp_adj    <- round(nExp * (1 - htypic_prop))
  
  # Run doubletFinder
  dbl_clmn <- str_c("DF.classifications", pN, pK, nExp_adj, sep = "_")
  
  res <- res %>%
    doubletFinder_v3(
      PCs        = PCs,
      pK         = pK,
      pN         = pN,
      nExp       = nExp_adj,
      reuse.pANN = reuse.pANN,
      ...
    ) %>%
    FetchData(dbl_clmn)
  
  # Add doublet classifications to original object
  res <- sobj_in %>%
    AddMetaData(res, "dbl_class")
  
  res
}

#' Run PCA and UMAP for gene expression data
#' 
#' @param sobj_in Seurat object.
#' @param assay Name of assay in object.
#' @param dims Dimensions to use for UMAP.
#' @param prefix Prefix to add to reduction keys and meta.data columns.
#' @param pca_meta Should PC-1 and PC-2 coordinates be added to meta.data table.
#' @param umap_meta Should UMAP coordinates be added to meta.data table.
#' @param ... Additional arguments to pass to RunPCA.
#' @return Seurat object
#' @export
run_UMAP_RNA <- function(sobj_in, assay = "RNA", dims = 1:40, prefix = "",
                         pca_meta = TRUE, umap_meta = TRUE, ...) {
  
  # Reduction keys
  pca_name  = str_c(prefix, "pca")
  pca_key   = str_c(prefix, "PC_")
  umap_name = str_c(prefix, "umap")
  umap_key  = str_c(prefix, "UMAP_")
  
  # Run PCA and UMAP
  # By default only variable features are used for PCA
  res <- sobj_in %>%  
    RunPCA(
      assay          = assay,
      reduction.name = pca_name,
      reduction.key  = pca_key,
      ...
    ) %>%
    RunUMAP(
      assay          = assay,
      dims           = dims,
      reduction.name = umap_name,
      reduction.key  = umap_key
    )
  
  # Add PCA to meta.data
  if (pca_meta) {
    pca_columns = str_c(pca_key, c(1, 2))
    
    res <- res %>%
      AddMetaData(
        metadata = FetchData(., pca_columns),
        col.name = pca_columns
      )
  }
  
  # Add UMAP to meta.data
  if (umap_meta) {
    umap_columns = str_c(umap_key, c(1, 2))
    
    res <- res %>%
      AddMetaData(
        metadata = Embeddings(., reduction = umap_name),
        col.name = umap_columns
      )
  }
  
  res
}

#' Run PCA, cluster, run UMAP, cluster cells
#' 
#' @param sobj_in Seurat object.
#' @param assay Name of assay in object.
#' @param resolution Resolution to use for clustering.
#' @param dims Dimensions to use for clustering.
#' @param prefix Prefix to add to reduction keys and meta.data columns.
#' @param pca_meta Should PC-1 and PC-2 coordinates be added to meta.data table.
#' @param umap_meta Should UMAP coordinates be added to meta.data table.
#' @param rerun_pca Re-run PCA and UMAP even if there is already a reduction
#' with the same name.
#' @param ... Additional arguments to pass to RunPCA.
#' @return Seurat object
#' @export
cluster_RNA <- function(sobj_in, assay = "RNA", resolution = 0.6, dims = 1:40, 
                        prefix = "", pca_meta = TRUE, umap_meta = TRUE, rerun_pca = TRUE,
                        ...) {
  
  # Use FindNeighbors to construct a K-nearest neighbors graph based on the euclidean distance in 
  # PCA space, and refine the edge weights between any two cells based on the
  # shared overlap in their local neighborhoods (Jaccard similarity).
  # Use FindClusters to apply modularity optimization techniques such as the Louvain algorithm 
  # (default) or SLM, to iteratively group cells together
  
  res <- sobj_in
  
  # Run PCA and UMAP
  # Data must be scaled
  umap_name <- str_c(prefix, "umap")
  
  if (!umap_name %in% names(sobj_in@reductions) || rerun_pca) {
    res <- res %>%
      run_UMAP_RNA(
        assay     = assay,
        prefix    = prefix,
        dims      = dims,
        pca_meta  = pca_meta,
        umap_meta = umap_meta,
        ...
      )
  }
  
  # Create nearest neighbors graph and find clusters
  pca_name <- str_c(prefix, "pca")
  
  res <- res %>%
    FindNeighbors(
      assay     = assay,
      reduction = pca_name,
      dims      = dims
    ) %>%
    FindClusters(
      resolution = resolution,
      verbose    = FALSE
    )
  
  res
}

#' Subset Seurat objects based on cell type
#' 
#' @param sobj_in Seurat object.
#' @param ... Arguments to pass to subset.
#' @param cluster_data Scale and cluster data after subsetting.
#' @param assay Assay to use for cell clustering.
#' @param var_p p-value cutoff for variable features.
#' @param rsln Resolution to use for cell clustering.
#' @param dims PCA components to use for cell clustering and UMAP.
#' @param regress_vars Variables to regress out when scaling data.
#' @return Seurat object
#' @export
subset_sobj <- function(sobj_in, ..., assay = "RNA", var_p = 0.05, rsln = 1,
                        dims = 1:40, regress_vars = NULL) {
  
  # Subset object
  res <- sobj_in %>%
    subset(...)
  
  # Find variable features with M3Drop
  # counts cannot be log-normalized
  counts <- res %>%
    GetAssayData(
      slot  = "counts",
      assay = assay
    ) %>%
    M3DropConvertData(
      is.log    = FALSE,
      is.counts = TRUE
    )
  
  var_genes <- counts %>%
    M3DropFeatureSelection(
      mt_threshold  = var_p,
      suppress.plot = TRUE
    )
  
  VariableFeatures(res) <- rownames(var_genes)
  
  # Cluster cells
  res <- res %>%
    ScaleData(
      assay = assay,
      vars.to.regress = regress_vars
    ) %>%
    cluster_RNA(
      assay      = assay,
      resolution = rsln,
      dims       = dims
    )
  
  res
}

#' Create labeller function to add cell n labels
#' 
#' @param sobj_in Seurat object or data.frame containing plot data.
#' @param lab_col meta.data column containing cell groups.
#' @param nm Should cell group be included in label.
#' @param sep Separator to use for labels.
#' @return Labeller function
#' @export
get_nlab_fun <- function(sobj_in, lab_col = NULL, nm = TRUE, l = "\n(", r = ")") {
  
  dat <- sobj_in
  
  if ("Seurat" %in% class(sobj_in)) {
    dat <- sobj_in@meta.data
  }
  
  if (!is.null(lab_col)) {
    dat <- dat %>%
      group_by(!!sym(lab_col))
  }
  
  labs <- dat %>%
    summarize(n = n(), .groups = "drop") %>%
    distinct() %>%
    arrange(desc(n)) %>%
    mutate(
      n = format(n, big.mark = ",", scientific = FALSE),
      n = str_trim(n),
      n = str_c("n = ", n)
    )
  
  if (nm && !is.null(lab_col)) {
    labs <- labs %>%
      mutate(n = str_c(!!sym(lab_col), l, n, r))
  }
  
  if (is.null(lab_col)) {
    return(pull(labs, "n"))
  }
  
  labs <- set_names(
    x  = pull(labs, "n"),
    nm = pull(labs, lab_col)
  )
  
  res <- function(x) labs[x]
  
  res
}

#' Perform k-means clustering on meta.data variable
#' 
#' @param dat data.frame with single column containing data to use for clustering.
#' @param k Number of clusters.
#' @param dat_clmn Column containing data to cluster
#' @param out_clmn Name of output column containing cell classifications.
#' @param clst_nms Labels to use for cell clusters.
#' @return data.frame containing cell clusters
#' @export
.run_km <- function(dat, k = 2, dat_clmn = "data", out_clmn = "km_cluster",
                    clst_nms = NULL) {
  
  # Data column name
  if (!is.null(colnames(dat))) {
    dat_clmn <- colnames(dat)
  }
  
  # K-means clustering
  res <- dat %>%
    stats::kmeans(centers = k)
  
  # Format results data.frame
  res <- res$cluster %>%
    data.frame()
  
  colnames(res) <- out_clmn
  
  if (!identical(rownames(dat), rownames(res))) {
    stop("Input and results rownames do not match.")
  }
  
  res <- dplyr::bind_cols(res, dat)
  
  # Add cluster names
  if (!is.null(clst_nms)) {
    if (length(clst_nms) != k) {
      stop("Must provide same number of cluster names as k.")
    }
    
    nms <- res %>%
      dplyr::group_by(!!sym(out_clmn)) %>%
      dplyr::summarize(mn = mean(!!sym(dat_clmn)), .groups = "drop") %>%
      dplyr::arrange(mn) %>%
      dplyr::pull(out_clmn)
    
    clst_nms <- purrr::set_names(clst_nms, nms)
    
    res <- res %>%
      tibble::rownames_to_column() %>%
      dplyr::mutate(
        !!sym(out_clmn) := clst_nms[as.character(!!sym(out_clmn))]
      ) %>%
      tibble::column_to_rownames()
  }
  
  res
}

#' Cluster meta.data variable using gaussian mixture model
#' 
#' @param dat data.frame with single column containing data to use for clustering.
#' @param k Number of clusters.
#' @param out_clmn Name of output column containing cell classifications.
#' @param clst_nms Labels to use for cell clusters.
#' @param prob Probability cutoff to use for classifying cells.
#' @param quiet Suppress output messages.
#' @return data.frame containing cell clusters
#' @export
.run_gmm <- function(dat, k = 2, out_clmn = "gmm_cluster", clst_nms = c("low", "high"),
                     prob = 0.5, quiet = TRUE) {
  
  if (length(clst_nms) != k) {
    stop("Must provide same number of cluster names as k.")
  }
  
  # Data column name
  dat_clmn <- "data"
  
  if (!is.null(colnames(dat))) {
    dat_clmn <- colnames(dat)
  }
  
  # Fit GMM for ova signal
  quiet_EM <- quietly(~ mixtools::normalmixEM(., k = k))
  
  if (!quiet) {
    quiet_EM <- mixtools::normalmixEM
  }
  
  set.seed(42)
  
  mdl <- dat %>%
    dplyr::pull(dat_clmn) %>%
    quiet_EM()
  
  if (quiet) {
    mdl <- mdl$result
  }
  
  # New column names
  comp_nms <- colnames(mdl$posterior)
  
  if (mdl$mu[1] > mdl$mu[2]) {
    clst_nms <- rev(clst_nms)
  }
  
  post              <- as.data.frame(mdl$posterior)
  colnames(post)    <- clst_nms
  names(comp_nms)   <- clst_nms
  names(mdl$mu)     <- clst_nms
  names(mdl$sigma)  <- clst_nms
  names(mdl$lambda) <- clst_nms
  
  # Format results data.frame
  clmns <- c("mu", "sigma", "lambda")
  
  clmns <- purrr::set_names(
    stringr::str_c(out_clmn, "_", clmns),
    clmns
  )
  
  res <- dplyr::bind_cols(dat, post)
  
  res <- res %>%
    tibble::rownames_to_column() %>%
    dplyr::mutate(
      !!sym(out_clmn) := if_else(
        !!sym(clst_nms[2]) >= prob,
        clst_nms[2],
        clst_nms[1]
      ),
      !!sym(clmns[["mu"]])     := mdl$mu[!!sym(out_clmn)],
      !!sym(clmns[["sigma"]])  := mdl$sigma[!!sym(out_clmn)],
      !!sym(clmns[["lambda"]]) := mdl$lambda[!!sym(out_clmn)],
      .before                   = !!sym(dat_clmn)
    ) %>%
    dplyr::select(-all_of(clst_nms)) %>%
    tibble::column_to_rownames()
  
  # Check that results match input data
  if (!identical(rownames(dat), rownames(res))) {
    stop("Input and results rownames do not match.")
  }
  
  res
}

#' Cluster meta.data variable
#' 
#' @param sobj_in Seurat object.
#' @param data_column meta.data column containing data to use for clustering.
#' @param k Number of clusters.
#' @param grp_column meta.data column containing cell labels to use for
#' dividing data. Clusters will be identified independently for each group.
#' @param filt Cell group present in grp_column to use for filtering cells before clustering.
#' All other cells will be labeled "other".
#' @param data_slot Slot to pull data from.
#' @param clust_column Name of meta.data column to output cell classifications.
#' @param clust_names Labels to use for cell clusters.
#' @param return_sobj Return a Seurat object. If FALSE a data.frame is
#' returned.
#' @param method Method to use for clustering, can be either "km" or "gmm".
#' @return Seurat object with cell classifications added to meta.data.
#' @export
cluster_signal_old <- function(sobj_in, data_column, k = 2, grp_column = NULL,
                           filt = NULL, data_slot = "counts",
                           clust_column = "clust",
                           clust_names = c("low", "high"),
                           return_sobj = TRUE, method = "gmm") {

  # Select method
  .funs <- list(
    km  = .run_km,
    gmm = .run_gmm
  )
  
  if (!method %in% names(.funs)) {
    stop("Must select one of the following methods: ", str_c(names(.funs), collapse = ", "))
  }
  
  .fun <- .funs[[method]]
  
  # Filter Seurat object
  so_flt <- sobj_in
  
  if (!is.null(filt)) {
    so_flt <- so_flt %>%
      subset(!!sym(grp_column) == filt)
  }
  
  # Split meta.data by grp_column
  so_flt <- list(so_flt)
  
  if (!is.null(grp_column)) {
    so_flt <- so_flt[[1]] %>%
      Seurat::SplitObject(grp_column)
  }
  
  # Cluster signal
  res <- so_flt %>%
    imap_dfr(~ {
      .x <- .x %>%
        Seurat::FetchData(data_column, slot = data_slot) %>%
        
        .fun(
          k        = k,
          out_clmn = clust_column,
          clst_nms = clust_names
        ) %>%
        
        tibble::rownames_to_column()
      
      if (!is.null(grp_column)) {
        .x <- .x %>%
          dplyr::mutate(!!sym(grp_column) := .y, .before = !!sym(clust_column))
      }  
      
      .x %>%  
        tibble::column_to_rownames()
    })
  
  # Return data.frame
  if (!return_sobj) {
    return(res)
  }
  
  # Add clusters to meta.data for input object
  res <- res %>%
    dplyr::select(-all_of(c(data_column, grp_column)))
  
  res <- sobj_in %>%
    Seurat::AddMetaData(res)
  
  # Add "other" label for cells not included in the comparison
  res <- res %>%
    djvdj::mutate_meta(mutate, !!sym(clust_column) := replace_na(!!sym(clust_column), "other"))
  
  res
}

cluster_signal <- function(sobj_in, data_column, k = 2, grp_column = NULL,
                           filt = NULL, data_slot = "counts",
                           clust_column = "clust",
                           clust_names = c("low", "high"),
                           return_sobj = TRUE, method = "gmm") {

  # Select method
  .funs <- list(
    km  = .run_km,
    gmm = .run_gmm
  )
  
  if (!method %in% names(.funs)) {
    stop("Must select one of the following methods: ", str_c(names(.funs), collapse = ", "))
  }
  
  .fun <- .funs[[method]]
  
  # Filter Seurat object
  so_flt <- sobj_in
  
  if (!is.null(filt)) {
    so_flt <- so_flt %>%
      subset(!!sym(grp_column) == filt)
  }
  
  # Split meta.data by grp_column
  so_flt <- list(so_flt)
  
  if (!is.null(grp_column)) {
    so_flt <- so_flt[[1]] %>%
      Seurat::SplitObject(grp_column)
  }
  
  # Cluster signal
  res <- so_flt %>%
    imap_dfr(~ {
      .x <- .x %>%
        Seurat::FetchData(data_column, slot = data_slot) %>%
        
        .fun(
          k        = k,
          out_clmn = clust_column,
          clst_nms = clust_names
        ) %>%
        
        tibble::rownames_to_column()
      
      if (!is.null(grp_column)) {
        .x <- .x %>%
          dplyr::mutate(!!sym(grp_column) := .y, .before = !!sym(clust_column))
      }  
      
      .x %>%  
        tibble::column_to_rownames()
    })
  
  # Return data.frame
  if (!return_sobj) {
    return(res)
  }
  
  # Add clusters to meta.data for input object
  res <- res %>%
    dplyr::select(-all_of(c(data_column, grp_column)))
  
  res <- sobj_in %>%
    Seurat::AddMetaData(res)
  
  # Add "other" label for cells not included in the comparison
  res <- res %>%
    djvdj::mutate_meta(mutate, !!sym(clust_column) := replace_na(!!sym(clust_column), "other"))
  
  res
}

#' Classify cell type based on module score
#' 
#' @param so_in Seurat object.
#' @param feats List of features to use for calculating module score. List name
#' will be used as cell type label.
#' @param prefix Prefix to add to meta.data column containing average module scores.
#' @param cutoff Minimum average score to use for identifying positive
#' clusters.
#' @param clst_col meta.data column containing cell clusters to use for
#' averaging module scores.
#' @param type_col meta.data column to add cell type label.
#' @return Seurat object containing new cell type classifications.
#' @export
classify_mod_score <- function(so_in, feats, prefix, cutoff = 1, clst_col, type_col = "cell_type") {
  
  nm   <- str_c(prefix, "_score")
  clmn <- str_c(nm, "1")
  
  res <- so_in %>%
    AddModuleScore(
      features = feats,
      name     = nm
    ) %>%
    mutate_meta(~ {
      .x %>%
        rename(!!sym(nm) := !!sym(clmn)) %>%
        group_by(!!sym(clst_col)) %>%
        mutate(
          !!sym(nm)       := mean(!!sym(nm)),
          !!sym(type_col) := if_else(!!sym(nm) > cutoff, names(feats), !!sym(type_col))
        )
    })
  
  res
}

#' Calculate feature fold change using median signal from control group
#' 
#' @param sobj_in Seurat object.
#' @param feat Name of feature to use for calculating fold changes.
#' @param data_slot Slot to pull data from.
#' @param count_pseudo Pseudo count to add to the median counts used for
#' calculating fold changes. This is useful if the median for control cells is
#' often 0.
#' @param fc_pseudo Should a pseudo count be added to calculated fold changes,
#' this is useful if the results will be log transformed.
#' @param fc_column Name of new column to add to meta.data table.
#' @param grp_column Name of meta.data column containing cell groups.
#' @param control_grps Name of control groups to use for calculating fold changes.
#' @return Seurat object
#' @export
calc_feat_fc <- function(sobj_in, feat = "adt_ovalbumin", data_slot = "counts", count_pseudo = 0,
                         fc_pseudo = FALSE, fc_column = "ova_fc", grp_column = "cell_type",
                         control_grps = c("B cell", "T cell")) {
  
  is_sobj <- "Seurat" %in% class(sobj_in)
  
  # Add data for calculations to meta.data
  meta <- sobj_in
  
  if (is_sobj) {
    res <- sobj_in %>%
      AddMetaData(FetchData(
        object = ., 
        vars   = feat,
        slot   = data_slot
      ))
    
    meta <- res@meta.data %>%
      rownames_to_column(".cell_id")
  }
    
  # Calculate fold changes
  meta <- meta %>%
    mutate(con_grp = if_else(                                               # Identify control groups
      !!sym(grp_column) %in% control_grps,
      TRUE,
      FALSE
    )) %>%
    
    group_by(con_grp, .add = TRUE) %>%                                                   # Calculate median counts for control
    mutate(con_med = ifelse(
      con_grp,
      median(!!sym(feat)),
      NA
    )) %>%
    ungroup() %>%
    
    mutate(
      con_med = replace_na(con_med, max(con_med, na.rm = TRUE)),  # Fill in NAs with median control counts
      con_med = con_med + count_pseudo,
      
      !!sym(fc_column) := !!sym(feat) / con_med
    ) %>%
    dplyr::select(-con_grp, -con_med)
  
  # Add pseudo count so fold changes can be log transformed
  # pseudo count is smallest non-zero value divided by 2
  if (0 %in% pull(meta, fc_column) && fc_pseudo) {
    meta <- meta %>%
      mutate(
        pseudo = ifelse(!!sym(fc_column) > 0, !!sym(fc_column), NA), 
        pseudo = min(pseudo, na.rm = TRUE) * 0.5,
        
        !!sym(fc_column) := !!sym(fc_column) + pseudo
      )
  }
  
  if (is_sobj) {
    meta <- meta %>%
      column_to_rownames(".cell_id")
    
    res@meta.data <- meta
    
  } else {
    res <- meta
  }
  
  res
}

#' Export counts and meta.data tables
#' 
#' @param sobj_in Seurat object
#' @param assays Assays to include in counts matrix
#' @param gene_prefix Prefix to add to gene names for counts matrix
#' @param columns meta.data columns to export
#' @param out_dir Output directory
#' @param file_prefix Prefix to add to output files
#' @return Counts and meta.data tables
#' @export
export_matrices <- function(sobj_in, assays = "RNA", gene_prefix = "", columns,
                            out_dir, file_prefix = "") {
  
  if (length(gene_prefix) == 1) {
    gene_prefix <- rep(gene_prefix, length(assays))
  }
  
  # Format count matrices
  counts <- map2_dfr(assays, gene_prefix, ~ {
    sobj_in %>%
      Seurat::GetAssayData(assay = .x, "counts") %>%
      tibble::as_tibble(rownames = "gene_symbol") %>%
      dplyr::mutate(gene_symbol = stringr::str_c(.y, gene_symbol))
  })
  
  # Write count matrix
  counts_out <- file.path(out_dir, str_c(file_prefix, "count_matrix.tsv.gz"))
  
  counts %>%
    readr::write_tsv(counts_out)
  
  # Write meta.data table
  meta_out <- file.path(out_dir, str_c(file_prefix, "metadata.tsv.gz"))
  
  sobj_in@meta.data %>%
    tibble::as_tibble(rownames = "cell_id") %>%
    dplyr::select(all_of(columns)) %>%
    readr::write_tsv(meta_out)
}


# Plotting ----

#' Create legend guide
#' @param nrow Number of rows for legend.
#' @param ncol Number of columns for legend.
#' @param size Legend glyph size.
#' @param shape Legend glyph shape.
#' @param stroke Legend glyph stroke size.
#' @param color Legend glyph color.
#' @param rev Reverse legend order.
#' @return guide_legend
#' @export
legd_gd <- function(nrow = NULL, ncol = NULL, size = 4, shape = 21,
                    stroke = 0.4, color = "black", rev = FALSE, ...) {
  opts <- list(
    ...,
    size   = size,
    shape  = shape,
    stroke = stroke,
    alpha  = 1
  )
  
  if (!is.null(color)) {
    opts <- append(opts, list(color = color))
  }
    
  guide_legend(
    nrow         = nrow,
    ncol         = ncol,
    reverse      = rev,
    override.aes = opts
  )
}

#' Create colorbar guide
#' 
#' @param wdth Colorbar width.
#' @param color Colorbar frame color.
#' @param ln_wdth Colorbar frame line width.
#' @param ttl Colorbar title.
#' @param direction Direction of colorbar.
#' @param ttl_pos Position of colorbar title.
#' @param ... Additional arguments to pass to guide_colorbar().
#' @return guide_colorbar
#' @export
bar_gd <- function(wdth = unit(0.15, "cm"), color = "black", ln_wdth = 0.3, ttl = waiver(),
                   direction = "vertical", ttl_pos = "top", ...) {
  ht <- NULL
  
  if (direction == "horizontal") {
    ht   <- wdth
    wdth <- NULL
  }
  
  guide_colorbar(
    barwidth        = wdth,
    barheight       = ht,
    frame.colour    = color,
    frame.linewidth = ln_wdth,
    title           = ttl,
    direction       = direction,
    title.position  = ttl_pos,
    ...
  )
}

#' Create boxplots showing feature signal
#' 
#' @param sobj_in Seurat object.
#' @param feats Features to plot.
#' @param grp_column meta.data column containing groups to plot.
#' @param box_colors Colors to use for boxplots.
#' @param median_pt Size of point used to mark median.
#' @param panels_n_row Number of rows for final figure.
#' @param panels_n_col Number of columns for final figure.
#' @param ... Additional arguments to adjust the boxplot theme.
#' @return ggplot object
#' @export
create_feat_boxes <- function(sobj_in, feats, grp_column = "subtype", box_colors = NULL,
                              median_pt = 1, panels_n_row = 4, panels_n_col = 10, plot_vln = FALSE, ...) {
  
  # Check for empty inputs
  if (is_empty(feats)) {
    res <- ggplot() +
      geom_blank()
    
    return(res)
  }
  
  # Boxplot data
  # What to arrange boxplots by median and 3rd quartile
  bx_dat <- sobj_in %>%
    FetchData(c(feats, grp_column)) %>%
    as_tibble(rownames = "cell_id")
  
  bx_labs <- get_nlab_fun(bx_dat, lab_col = grp_column, l = " (")
  
  bx_dat <- bx_dat %>%  
    pivot_longer(cols = c(-cell_id, -!!sym(grp_column))) %>%
    mutate(type_name = str_c(!!sym(grp_column), "_", name)) %>%
    group_by(type_name) %>%
    mutate(
      up_qt = boxplot.stats(value)$stats[4],
      med   = median(value),
      max   = max(value)
    ) %>%
    ungroup() %>%
    arrange(desc(med), desc(up_qt), desc(max)) %>%
    mutate(
      name      = fct_relevel(name, feats),
      type_name = fct_inorder(type_name)
    )
  
  # Set plot theme
  res <- bx_dat %>%
    ggplot(aes(type_name, value, fill = !!sym(grp_column))) +
    
    facet_wrap(
      ~ name,
      nrow   = panels_n_row,
      ncol   = panels_n_col,
      scales = "free_x"
    ) +
    
    labs(y = "Counts") +
    
    theme_minimal_hgrid() +
    base_theme +
    theme(
      legend.position    = "bottom",
      legend.title       = element_blank(),
      legend.key.height  = unit(0.45, "cm"),
      axis.title.x       = element_blank(),
      axis.text.x        = element_blank(),
      axis.line.x        = element_blank(),
      axis.ticks.x       = element_blank(),
      panel.grid.major.y = element_blank()
    ) +
    theme(...)
  
  if (!is.null(box_colors)) {
    res <- res +
      scale_fill_manual(values = box_colors, labels = bx_labs) +
      scale_color_manual(values = box_colors)
  }
  
  # Return violin plots
  if (plot_vln) {
    res <- res +
      geom_violin(size = 0.2, color = "black", draw_quantiles = c(0.25, 0.75)) +
      stat_summary(geom = "point", fun = median, size = 0.5)
    
    return(res)
  }
  
  # Return boxplots
  res <- res +
    stat_summary(geom = "point", shape = 22, fun = median, size = 0) +
    stat_summary(geom = "point", shape = 22, fun = median, size = median_pt + 1, color = "black") +
    geom_boxplot(
      size           = 0.3,
      width          = 0.6,
      fatten         = 0,
      outlier.colour = "grey85",
      outlier.alpha  = 1,
      outlier.size   = 0.1,
      show.legend    = FALSE
    ) +
    stat_summary(
      aes(color = !!sym(grp_column)),
      geom        = "point",
      shape       = 22,
      fun         = median,
      size        = median_pt,
      stroke      = 0.5,
      fill        = "#ffffff",
      show.legend = FALSE
    ) +
    guides(fill = legd_gd(nrow = 3, shape = 22)) +
    labs(y = "log normalized counts") +
    theme(
      panel.background = element_rect(fill = fade_0),
      strip.text       = element_text(size = txt_pt2)
    )
  
  # Fix plot columns
  axis_buf <- 0.5           # This is to account for the width of y-axis
  n_feats <- length(feats)
  
  if (n_feats < panels_n_col) {
    frac <- (n_feats + axis_buf) / panels_n_col
    
    res <- plot_grid(
      res, NULL,
      nrow = 1,
      rel_widths = c(frac, 1 - frac)
    )
  }
  
  # Fix plot rows
  legd_buf <- 0.45          # This is to account for the height of legend    
  n_row <- ceiling(n_feats / panels_n_col)
  
  if (n_row < panels_n_row) {
    frac <- (n_row + legd_buf) / panels_n_row
    
    res <- plot_grid(
      res, NULL,
      ncol = 1,
      rel_heights = c(frac, 1 - frac)
    )
  }
  
  res
}

#' Create GO plot
#' 
#' @param GO_df data.frame of GO results.
#' @param plot_colors Plot colors.
#' @param n_terms Number of terms to label.
#' @return ggplot object
#' @export
create_bubbles <- function(go_df, plot_colors = NULL, n_terms = 15) {
  
  # Check for empty inputs
  if (is_empty(go_df) || nrow(go_df) == 0) {
    res <- ggplot() +
      geom_blank()
    
    return(res)
  }
  
  # Shorten GO terms and database names
  go_nms <- c(
    "GO:BP" = "Biological\nProcess",
    "GO:CC" = "Cellular\nComponent",
    "GO:MF" = "Molecular\nFunction",
    "KEGG"  = "KEGG"
  )
  
  go_dat <- go_df %>%
    mutate(
      term_id = str_remove(term_id, "(GO|KEGG):"),
      term_id = str_c(term_id, " ", term_name),
      term_id = str_to_lower(term_id),
      term_id = str_trunc(term_id, 40, "right"),
      source  = go_nms[source]
    )
  
  # Reorder database names
  plot_lvls <- unname(go_nms)
  plot_lvls <- plot_lvls[plot_lvls %in% go_dat$source]
  plot_lvls <- unique(c(plot_lvls, go_dat$source))
  
  go_dat <- go_dat %>%
    mutate(source = fct_relevel(source, plot_lvls))
  
  # Extract top terms for each database
  top_go <- go_dat %>%
    group_by(source) %>%
    arrange(p_value) %>%
    dplyr::slice(1:n_terms) %>%
    ungroup()
  
  # Create bubble plots
  go_cols <- c("#D55E00", "#0072B2", "#009E73", "#6A51A3")
  
  if (!is.null(plot_colors)) {
    go_cols <- plot_colors
  }
  
  res <- go_dat %>%
    ggplot(aes(1.25, -log10(p_value), size = intersection_size)) +
    geom_point(aes(color = source), alpha = 0.5, show.legend = TRUE) +
    
    geom_text_repel(
      aes(2, -log10(p_value), label = term_id),
      data         = top_go,
      size         = 2.3,
      direction    = "y",
      hjust        = 0,
      segment.size = NA
    ) +
    
    guides(color = FALSE) +
    xlim(1, 8) +
    labs(y = "-log10(p-value)") +
    
    scale_color_manual(values = go_cols) +
    base_theme +
    theme(
      legend.position = "bottom",
      axis.title.x    = element_blank(),
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank(),
      axis.line.x     = element_blank()
    ) +
    facet_wrap(~ source, scales = "free", nrow = 1)
  
  res
}

#' Run gprofiler
#' 
#' @param gene_list List of input genes.
#' @param genome Genome to use for identifying GO terms.
#' @param gmt_id GMT ID to use instead of genome name.
#' @param p_max Maximum adjusted p-value for GO term to be included in results.
#' @param GO_size Minimum term size for GO term.
#' @param intrsct_size Minimum number of genes that intersect GO term.
#' @param order_query Run an ordered query where more significant genes are at
#' the top of the list.
#' @param dbases GO databases to query.
#' @param file_path File path to write results.
#' @param overwrite Overwrite existing file.
#' @param ... Additional arguments to pass to gprofiler2.
#' @return tibble containing GO terms.
#' @export
run_gprofiler <- function(gene_list, genome = NULL, gmt_id = NULL, p_max = 0.05, GO_size = 10, intrsct_size = 10, 
                          order_query = FALSE, dbases = c("GO:BP", "GO:MF", "KEGG"), file_path = NULL,
                          overwrite = FALSE, ...) {
  
  # Check for empty gene list
  if (is_empty(gene_list)) {
    return(as_tibble(NULL))
  }
  
  # Check arguments
  if (is.null(genome) && is.null(gmt_id)) {
    stop("Must specifiy genome or gmt_id.")
  }
  
  if (!is.null(file_path) && file.exists(file_path) && !overwrite) {
    warning("Loading existing file: ", file_path, ". Set overwrite = TRUE to rerun.")
    
    return(read_tsv(file_path))
  }
  
  # Use gmt id
  if (!is.null(gmt_id)) {
    genome <- gmt_id
    dbases <- NULL
  }
  
  # Run gProfileR
  res <- gene_list %>%
    gost(
      organism       = genome,
      sources        = dbases,
      domain_scope   = "annotated",
      significant    = TRUE,
      user_threshold = p_max,
      ordered_query  = order_query,
      ...
    )
  
  # Format and sort output data.frame
  res <- as_tibble(res$result)
  
  if (!is_empty(res)) {
    
    if (!is.null(GO_size)) {
      res <- res %>%
        filter(term_size > GO_size)
    }
      
    if (!is.null(intrsct_size)) {
      res <- res %>%
        filter(intersection_size > intrsct_size)
    }
    
    res <- res %>%
      arrange(source, p_value)
  }
  
  # Write table
  if (!is.null(file_path) && nrow(res) > 0) {
    res %>%
      select(-parents) %>%
      write_tsv(file_path)    
  }
  
  res
}

#' Find markers with presto
#' 
#' @param sobj_in Seurat object.
#' @param grp_column meta.data column containing labels to use for grouping cells.
#' @param exclude_grp Name of cell group to exclude from comparison.
#' @param groups_use Vector of cell groups to use for comparison.
#' @param filt_group Filter results to only include groups specified by filt_group.
#' @param fc_range Vector of length 2 containing minimum and maximum log fold change for markers.
#' @param auc_min Minimum AUC for markers.
#' @param p_max Maximum adjusted p-value for filtering markers.
#' @param pct_in_min Minimum percentage of cells within the group that express each marker.
#' @param pct_out_max Maximum percentage of cells outside the group that express each marker.
#' @param ... Additional arguments to pass to wilcoxauc.
#' @return tibble
#' @export
# find_markers <- function(sobj_in, grp_column = NULL, exclude_grp = NULL, groups_use = NULL, 
#                          filt_group = NULL, fc_range = c(0.25, Inf), auc_min = 0.5, p_max = 0.05, pct_in_min = 50, 
#                          pct_out_max = 100, ...) {
#   
#   if (!is.null(exclude_grp) && is.null(groups_use)) {
#     sobj_in <- sobj_in %>%
#       subset(subset = !!sym(grp_column) != exclude_grp)
#   }
#   
#   res <- sobj_in %>%
#     wilcoxauc(
#       group_by   = grp_column,
#       groups_use = groups_use,
#       ...
#     ) %>%
#     as_tibble() %>%
#     filter(
#       padj    < p_max,
#       logFC   > fc_range[1],
#       logFC   < fc_range[2],
#       auc     > auc_min,
#       pct_in  > pct_in_min,
#       pct_out < pct_out_max,
#       group %in% filt_group
#     ) %>%
#     arrange(desc(logFC)) %>%
#     rename(gene = feature)
#   
#   res
# }

#' Run FindConservedMarkers for all cell identities
#' 
#' @param sobj_in Seurat object.
#' @param ident_1 Identity class to define markers for.
#' @param idnet_2 Second identity class for comparison. If NULL use all other
#' cells for comparison.
#' @param grp_var Grouping variable for finding markers.
#' @param p_max Maximum p-value for filtering marker genes.
#' @param min_fc Minimum log2 fold change for filtering marker genes.
#' @param fc_range Vector of length 2 containing minimum and maximum log2 fold
#' change for markers. Use this argument to select positive or negative
#' markers.
#' @param file_path File path to write results.
#' @param filt_regex Remove genes from results that match the provided regular
#' expression
#' @param overwrite Overwrite existing file.
#' @return tibble containing marker genes.
#' @export
find_conserved_markers <- function(sobj_in, ident_1 = NULL, ident_2 = NULL, grp_var, p_max = 0.05,
                                   fc_range = c(log2(1.25), Inf), file_path = NULL, filt_regex = NULL,
                                   overwrite = FALSE) {
  
  if (is.null(ident_1) && !is.null(ident_2)) {
    stop("Must specify ident_1 when ident_2 is given.")
  }
  
  if (!is.null(file_path) && file.exists(file_path) && !overwrite) {
    warning("Loading existing file: ", file_path, ". Set overwrite = TRUE to rerun.")
    
    return(read_tsv(file_path))
  }
  
  if (length(fc_range) != 2) {
    stop("fc_range must be a vector of length 2")
  }
  
  if (!fc_range[1] < fc_range[2]) {
    stop("The first value of fc_range must be smaller than the second value")
  }
  
  if (is.null(ident_1)) {
    ident_1 <- unique(Idents(sobj_in))
  }
  
  # Set FC filtering parameters
  only_pos <- FALSE
  fc_lim   <- 0
  
  if (all(fc_range > 0)) {
    only_pos <- TRUE
    fc_lim   <- fc_range[1]
    
  } else if (all(fc_range < 0)) {
    fc_lim <- abs(fc_range[2])
  }
  
  # Find markers
  res <- ident_1 %>%
    map_dfr(~ {
      mks <- sobj_in %>%
        FindConservedMarkers(
          ident.1         = .x,
          ident.2         = ident_2,
          grouping.var    = grp_var,
          logfc.threshold = fc_lim,
          only.pos        = only_pos
        ) %>%
        rownames_to_column("gene")
      
      # Filter and format output table
      fc_clmns <- grep("_avg_log2FC$", colnames(mks), value = TRUE)
      
      if (nrow(mks) > 0) {
        mks <- mks %>%
          rowwise() %>%
          mutate(
            avg_log2FC = mean(!!!syms(fc_clmns)),
            ident_1    = .x,
            ident_2    = ident_2
          ) %>%
          filter(
            across(
              all_of(fc_clmns),
              ~ .x > fc_range[1] & .x < fc_range[2]
            ),
            max_pval < p_max
          ) %>%
          ungroup()
      }
      
      mks
    })
  
  # Remove genes that match filt_regex
  if (!is.null(filt_regex)) {
    res <- res %>%
      filter(!grepl(filt_regex, gene))
  }
  
  # Write table
  if (!is.null(file_path) && nrow(res) > 0) {
    write_tsv(res, file_path)
  }
  
  res
}

#' Find marker genes and create summary figures
#' 
#' @param sobj_in Seurat object.
#' @param marker_fn Function to use for finding markers.
#' @param n_genes Number of marker genes to plot.
#' @param box_x Group to use for x-axis of boxplots.
#' @param box_ncol Number of columns for arranging boxplots.
#' @param box_cols Colors for boxplots.
#' @param go_genome Genome to use for finding GO terms.
#' @param go_cols Colors for GO plots.
#' @param rel_h Relative height of boxplot and GO term panels.
#' @param plot_vln Should violin plots be generated instead of boxplots?
#' @param show_go Should panels showing GO terms be included?
#' @param file_dir Directory to save output files.
#' @param file_prefix Prefix to add to output files.
#' @param overwrite Should gene lists be overwritten if they already exist?
#' @param ... Arguments to pass to marker_fn.
#' @return ggplot object
#' @export
create_marker_fig <- function(sobj_in, marker_fn, n_genes = 40, box_x, box_nrow = 4, box_ncol = 10, box_cols = NULL,
                              go_genome = "mmusculus", n_terms = 10, go_cols = NULL, rel_h = c(1, 0.5), plot_vln = FALSE,
                              show_go = TRUE, file_dir = NULL, file_prefix = NULL, overwrite = FALSE, ...) {
  
  # Markers
  m_file <- file_dir
  
  if (!is.null(m_file)) {
    m_file <- str_c(file_prefix, "markers.tsv", sep = "_")
    m_file <- here(file_dir, m_file)
  }
  
  marks <- sobj_in %>%
    marker_fn(
      ...,
      file_path = m_file,
      overwrite = overwrite
    )
  
  N_MARKERS <<- nrow(marks)
  
  top_marks <- marks %>%
    head(n_genes) %>%
    pull(gene)
  
  # Boxplots
  boxes <- sobj_in %>%
    create_feat_boxes(
      feats        = top_marks,
      grp_column   = box_x,
      box_colors   = box_cols,
      panels_n_row = box_nrow,
      panels_n_col = box_ncol,
      plot_vln     = plot_vln
    ) +
    theme(plot.margin = unit(c(0.2, 0.2, 1, 0.2), "cm"))
  
  # GO terms
  go_file <- file_dir
  
  if (!is.null(go_file)) {
    go_file <- str_c(file_prefix, "GO.tsv", sep = "_")
    go_file <- here(file_dir, go_file)
  }
  
  go_df <- marks %>%
    pull(gene) %>%
    run_gprofiler(
      genome      = go_genome,
      order_query = TRUE,
      file_path   = go_file,
      overwrite   = overwrite
    )
  
  N_GO_TERMS <<- nrow(go_df)
  
  # Return boxplots if show_go = FALSE
  if (!show_go) {
    return(boxes)
  }
  
  # Create dot plots
  go_dot <- go_df %>%
    create_bubbles(
      n_terms     = n_terms,
      plot_colors = go_cols
    ) +
    theme(legend.title = element_text(size = txt_pt2))
  
  # Create final figure
  res <- plot_grid(
    boxes, go_dot,
    rel_heights = rel_h,
    ncol        = 1
  )
  
  res
}

#' Add p-value labels to plot
#' 
#' @param gg_in ggplot object.
#' @param x x start coordinate for line.
#' @param xend x end coordinate for line.
#' @param y y coordinate for line and label.
#' @param p_val p-value to add to label.
#' @param prefix Prefix to add to p-value label.
#' @param format_p Format p-value for label.
#' @param digits Number of digits to use for formatting p-value.
#' @param add_line Include line to show comparison.
#' @param scale_line Value to use for scaling length of line.
#' @param line_size Line size.
#' @param line_col Line color.
#' @param line_type Line type.
#' @param ... Additional arguments to pass to geom_text.
#' @return ggplot object.
#' @export
add_pvals <- function(gg_in, x, xend, y, p_val, prefix = "", format_p = TRUE,
                      digits = 2, add_line = TRUE, scale_line = 0.1, line_size = 0.1,
                      line_col = "black", line_type = 1, ...) {
  
  p   <- p_val
  prs <- TRUE
  
  if (format_p) {
    p <- p %>%
      scales::scientific(digits = digits) %>%
      str_split("e") %>%
      unlist()
    
    expt <- as.numeric(p[2])
    p[2] <- as.character(expt)  # Converting to numeric and back removes extra zeros and "+"
    
    if (expt < -digits) {
      p <- p %>%
        str_c(collapse = "x10\"^\"") %>%
        str_c("\"", prefix, ., "\"")
    } else {
      prs <- FALSE
      p   <- round(p_val, digits = digits)
      p   <- str_c(prefix, p)
    }
  }
  
  # Set label position
  dat <- gg_in$data
  
  x_nm <- as_name(gg_in$mapping$x)
  
  x_lvls <- dat %>%
    arrange(!!sym(x_nm)) %>%
    pull(x_nm) %>%
    unique()
  
  if (is.character(x)) {
    x <- grep(x, x_lvls)
  }
  
  if (is.character(xend)) {
    xend <- grep(xend, x_lvls)
  }
  
  x_mid   <- median(c(x, xend))
  x_len   <- xend - x
  x_scale <- x_len * scale_line
  x       <- x + x_scale
  xend    <- xend - x_scale
  
  # Add label
  res <- gg_in +
    geom_text(
      aes(x_mid, y + (y * 0.07)),
      label         = p,
      parse         = prs,
      check_overlap = TRUE,
      ...
    )
  
  # Add line
  if (add_line) {
    res <- res +
      geom_segment(
        aes(x = x, xend = xend, y = y, yend = y),
        size     = line_size,
        color    = line_col,
        linetype = line_type
      )
  }
  
  res
}


# Figures ----

#' UMAPs for different clustering resolutions
plot_clust_rslns <- function(df_in, c_clmns = NULL, t_clmns = NULL, r_clmns = NULL, c_cols = get_cols(130),
                             t_cols = get_cols(40), r_cols = c("white", "#E69F00"), ...) {
  
  # UMAP helper
  plot_umap <- function(..., pt_size = 0.01) {
    plot_features(
      ...,
      feature = "value",
      size    = pt_size,
      n_label = "none"
    ) +
      facet_wrap(~ name, nrow = 1) +
      guides(color = guide_legend(override.aes = list(size = 3))) +
      umap_theme +
      theme(
        aspect.ratio     = 0.9,
        panel.background = element_rect(color = fade_1, fill = NA, size = 0.5),
        legend.title     = element_blank()
      )
  }
  
  umaps <- list()
  
  # Cluster UMAPs
  if (!is.null(c_clmns)) {
    u <- df_in %>%
      pivot_longer(cols = all_of(c_clmns)) %>%
      group_by(name) %>%
      mutate(
        n    = n_distinct(value),
        name = str_remove(name, "^RNA_snn_res\\."),
        name = str_c(name, " (", n, ")")
      ) %>%
      plot_umap(plot_colors = c_cols, ...) +
      theme(legend.position = "none")
    
    umaps <- append(umaps, list(u))
  }
  
  # Cell type UMAPs
  if (!is.null(t_clmns)) {
    u <- df_in %>%
      pivot_longer(cols = all_of(t_clmns)) %>%
      plot_umap(plot_colors = t_cols, ...)
    
    umaps <- append(umaps, list(u))
  }
  
  # Correlation UMAPS
  if (!is.null(r_clmns)) {
    u <- df_in %>%
      pivot_longer(cols = all_of(r_clmns)) %>%
      
      plot_umap(plot_colors = r_cols, ...) +
      
      guides(color = bar_gd(ttl = "r")) +
      theme(legend.title = element_text())
    
    umaps <- append(umaps, list(u))
  }
  
  # Create final figure
  res <- plot_grid(
    plotlist = umaps,
    ncol     = 1,
    align    = "v",
    axis     = "rl"
  )
  
  res
}

#' Create panel with UMAP and bars
create_umap_bars <- function(df_in, fill, grps = "orig.ident", filt = TRUE, plot_clrs, plot_lvls = NULL,
                             order_bars_by_n = TRUE, size = 0.0001, list_out = TRUE, legd_rows = 10,
                             ttl = NULL, n_label = FALSE, flip_bars = TRUE, umap_aspect_ratio = NULL,
                             bar_aspect_ratio = NULL, legd_labs = NULL, plot_theme = base_theme, ...) {
  
  bar_lvls <- plot_lvls
  
  if (order_bars_by_n) {
    bar_lvls <- df_in %>%
      group_by(!!sym(fill)) %>%
      summarize(n = n()) %>%
      arrange(desc(n)) %>%
      pull(fill)
  }
  
  # Set legend labels
  tot_lab   <- get_nlab_fun(df_in)
  type_labs <- legd_labs
  
  if (is.null(legd_labs)) {
    type_labs <- get_nlab_fun(df_in, fill, l = " (")
  }
  
  # Create UMAP
  type_umap <- df_in %>%
    plot_features(
      feature     = fill,
      size        = size,
      plot_colors = plot_clrs,
      plot_lvls   = rev(plot_lvls),
      ...
    ) +
    
    ggtitle(ttl) +
    umap_theme +
    theme(legend.position = "none")
  
  if (n_label) {
    type_umap <- type_umap + 
      geom_text(
        aes(Inf, Inf),
        label         = tot_lab,
        check_overlap = TRUE,
        color         = "black",
        hjust         = 2,
        vjust         = 2,
        size          = txt_pt2 / .pt
      )
  }
  
  if (!is.null(umap_aspect_ratio)) {
    type_umap <- type_umap +
      theme(aspect.ratio = umap_aspect_ratio)
  }
  
  # Create bar graphs
  bar_df <- df_in
  
  bar_df <- bar_df %>%
    filter({{filt}})
  
  if (!is.null(bar_lvls)) {
    bar_df <- bar_df %>%
      mutate(!!sym(fill) := fct_relevel(!!sym(fill), bar_lvls))
  }
  
  type_bars <- bar_df %>%
    ggplot(aes(!!sym(grps), fill = !!sym(fill))) +
    geom_bar(
      position  = "fill",
      color     = "black",
      size      = 0.2,
      alpha     = 0.9,
      key_glyph = draw_key_point
    ) +
    
    scale_fill_manual(values = plot_clrs, labels = type_labs) +
    
    guides(fill = legd_gd(shape = 22, nrow = legd_rows)) +
    
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    
    labs(y = "fraction of cells") +
    plot_theme +
    theme(
      legend.position   = "top",
      legend.title      = element_blank(),
      legend.key.height = unit(18, "pt"),
      legend.text       = element_text(hjust = 0),
      axis.title.x      = element_blank(),
      axis.text.x       = element_text(angle = 45, hjust = 1)
    )
  
  if (!is.null(bar_aspect_ratio)) {
    type_bars <- type_bars +
      theme(aspect.ratio = bar_aspect_ratio)
  }
  
  # Flip bar graphs
  if (flip_bars) {
    type_bars <- type_bars +
      coord_flip() +
      labs(x = "fraction of cells") +
      theme(
        legend.position = "right",
        legend.title    = element_blank(),
        axis.title.y    = element_blank()
      )
    
    if (list_out) {
      legd <- get_legend(type_bars)
      
      type_bars <- type_bars +
        theme(
          legend.position = "none",
          axis.text.x     = element_text(angle = 0, hjust = 0.5),
          axis.title.x    = element_text(size = txt_pt2)
        )
      
      res <- list(type_umap, type_bars, legd)
      
      return(res)
    }
  }
  
  # Create final figure
  if (list_out) {
    legd <- get_legend(type_bars)
    
    type_bars <- type_bars +
      theme(legend.position = "none")
    
    res <- list(type_umap, legd, type_bars)
    
    return(res)
  }
  
  if (flip_bars) {
    res <- plot_grid(
      type_umap, type_bars,
      ncol        = 1,
      rel_heights = c(1, 0.5)
    )
    
    return(res)
  }
  
  res <- plot_grid(
    type_umap, type_bars,
    nrow       = 1,
    rel_widths = c(1, 0.6)
  )
  
  res
}

#' Create boxplots comparing expression of provided gene
#' 
#' @param df_in data.frame containing data for plotting.
#' @param type Named vector specifying cell type to filter for before
#' generating plots. The name should indicate the column to use for filtering.
#' @param gene Gene to plot.
#' @param plt_clrs Named vector specifying plot colors.
#' @param pval_out File path to use for saving calculated p-values.
#' @param p_lab_pos Vector containing y coordinates to use for positioning
#' p-values in the plot area.
#' @return ggplot object
#' @export
create_gene_boxes <- function(df_in, type, gene, plt_ttl = type, plt_clrs, pval_out = NULL, p_lab_pos) {
  
  type_clmn <- names(type)
  type      <- unname(type)
  lvls      <- names(plt_clrs)
  
  dat <- df_in %>%
    filter(!!sym(type_clmn) == type) %>%
    mutate(
      treat_chikv = if_else(treatment == "CHIKV", chikv_grp, treatment),
      treat_chikv = fct_relevel(treat_chikv, lvls)
    )
  
  v_sub <- get_nlab_fun(dat)
  
  # Mxra8 CHIKV-high boxes
  res <- dat %>%
    create_boxes(
      x         = "treat_chikv",
      y         = gene,
      plot_clrs = plt_clrs
    ) +
    box_theme +
    labs(title = plt_ttl, subtitle = v_sub, y = str_c(gene, " expression"))
  
  # Calculate p-values
  p_vals <- dat %>%
    calc_p_vals(
      data_column = gene,
      type_column = "treat_chikv"
    )
  
  if (!is.null(pval_out)) {
    p_vals %>%
      write_csv(pval_out)
  }
  
  comps <- lvls %>%
    combn(2, simplify = FALSE)
  
  comps <- list(
    x    = map_chr(comps, pluck, 1),
    xend = map_chr(comps, pluck, 2),
    y    = p_lab_pos
  )
  
  comps %>%
    pwalk(~ {
      v_args <- list(...)[c("x", "xend")]
      
      p <- p_vals %>%
        filter(`Cell type 1` %in% v_args & `Cell type 2` %in% v_args) %>%
        pull(p_adj)
      
      res <<- res %>%
        add_pvals(
          ...,
          p_val    = p,
          size     = 10 / .pt,
          line_col = "grey75"
        )
    })
  
  res
}

#' Create boxplots
#'
#' @param df_in data.frame containing data to plot.
#' @param x x-axis variable.
#' @param y y-axis variable.
#' @param plot_clrs Fill color for boxplots.
#' @param type Type of plot to generate, must be either "boxplot" or "violin".
#' @param ... Additional arguments to pass to ggplot2.
#' @return ggplot object
#' @export
create_boxes <- function(df_in, x, y, fill = x, plot_clrs = NULL, type = "boxplot",
                         add_n_labs = TRUE, ...) {
  
  if (!type %in% c("boxplot", "violin")) {
    stop("type must be either 'boxplot' or 'violin'.")
  }
  
  clrs <- plot_clrs
  
  if (add_n_labs) {
    x_labs <- df_in %>%
      get_nlab_fun(x)
    
    clrs        <- plot_clrs[names(x_labs())]
    names(clrs) <- x_labs(names(clrs))
    
    df_in <- df_in %>%
      mutate(!!sym(x) := x_labs(as.character(!!sym(x))))
  }
  
  df_in <- df_in %>%
    mutate(!!sym(x) := fct_relevel(!!sym(x), names(clrs)))
  
  res <- df_in %>%
    ggplot(aes(!!sym(x), !!sym(y), fill = !!sym(fill))) +
    scale_fill_manual(values = clrs)
  
  if (type == "boxplot") {
    res <- res +
      geom_boxplot(
        ...,
        color        = "black",
        size         = 0.5,
        outlier.size = 0.2,
        alpha        = 0.8
      )
  }
  
  if (type == "violin") {
    res <- res +
      geom_violin(
        ...,
        color          = "black",
        size           = 0.5,
        alpha          = 0.8,
        draw_quantiles = c(0.25, 0.75)
      ) +
      stat_summary(
        geom = "point",
        fun  = median,
        size = 2
      )
  }
  
  res
}

#' Create boxplots showing signal for provided features
#' 
#' @param sobj_in Seurat object.
#' @param feats Features to plot.
#' @param grp_column meta.data column containing groups to plot.
#' @param box_colors Colors to use for boxplots.
#' @param median_pt Size of point used to mark median.
#' @param panels_n_row Number of rows for final figure.
#' @param panels_n_col Number of columns for final figure.
#' @param show_x_labels Should x-axis labels be shown?
#' @param plot_theme Base plot theme to use for boxplots.
#' @param ... Additional arguments to adjust boxplot theme.
#' @return ggplot object
#' @export
create_feat_boxes <- function(sobj_in, feats, grp_column = "subtype",
                              box_colors = NULL, median_pt = 1, panels_n_row = 4,
                              panels_n_col = 10, plot_vln = FALSE,
                              show_x_labels = FALSE, plot_theme = base_theme,
                              ...) {
  
  # Check for empty inputs
  if (is_empty(feats)) {
    res <- ggplot() +
      geom_blank() +
      theme(panel.background = element_blank())
    
    return(res)
  }
  
  # Boxplot data
  # What to arrange boxplots by median and 3rd quartile
  bx_dat <- sobj_in %>%
    FetchData(c(feats, grp_column)) %>%
    as_tibble(rownames = "cell_id")
  
  legd_labs <- get_nlab_fun(bx_dat, lab_col = grp_column, l = " (")
  
  bx_dat <- bx_dat %>%  
    pivot_longer(cols = c(-cell_id, -!!sym(grp_column))) %>%
    mutate(type_name = str_c(!!sym(grp_column), "_", name)) %>%
    group_by(type_name) %>%
    mutate(
      up_qt = boxplot.stats(value)$stats[4],
      med   = median(value),
      max   = max(value)
    ) %>%
    ungroup() %>%
    arrange(desc(med), desc(up_qt), desc(max)) %>%
    mutate(
      name      = fct_relevel(name, feats),
      type_name = fct_inorder(type_name)
    )
  
  # Get x-axis labels
  x_labs <- bx_dat %>%
    distinct(!!sym(grp_column), type_name)
  
  x_labs <- set_names(
    x_labs[[grp_column]],
    x_labs$type_name
  )
  
  # Set plot theme
  res <- bx_dat %>%
    ggplot(aes(type_name, value, fill = !!sym(grp_column))) +
    
    facet_wrap(
      ~ name,
      nrow   = panels_n_row,
      ncol   = panels_n_col,
      scales = "free_x"
    ) +
    
    labs(y = "Counts") +
    scale_x_discrete(labels = x_labs) +
    
    theme_minimal_hgrid() +
    plot_theme +
    theme(
      legend.position    = "bottom",
      legend.title       = element_blank(),
      legend.key.height  = unit(0.45, "cm"),
      axis.title.x       = element_blank(),
      panel.grid.major.y = element_blank()
    )
  
  if (!show_x_labels) {
    res <- res +
      theme(
        axis.text.x  = element_blank(),
        axis.line.x  = element_blank(),
        axis.ticks.x = element_blank()
      )
  }
  
  res <- res +
    theme(...)
  
  if (!is.null(box_colors)) {
    res <- res +
      scale_fill_manual(values = box_colors, labels = legd_labs) +
      scale_color_manual(values = box_colors)
  }
  
  # Return violin plots
  if (plot_vln) {
    res <- res +
      geom_violin(size = 0.2, color = "black", draw_quantiles = c(0.25, 0.75)) +
      stat_summary(geom = "point", fun = median, size = 0.5)
    
    return(res)
  }
  
  # Return boxplots
  res <- res +
    stat_summary(geom = "point", shape = 22, fun = median, size = 0) +
    stat_summary(geom = "point", shape = 22, fun = median, size = median_pt + 1, color = "black") +
    geom_boxplot(
      size           = 0.3,
      width          = 0.6,
      fatten         = 0,
      outlier.colour = "grey85",
      outlier.alpha  = 1,
      outlier.size   = 0.1,
      show.legend    = FALSE
    ) +
    stat_summary(
      aes(color = !!sym(grp_column)),
      geom        = "point",
      shape       = 22,
      fun         = median,
      size        = median_pt,
      stroke      = 0.5,
      fill        = "#ffffff",
      show.legend = FALSE
    ) +
    guides(fill = legd_gd(nrow = 3, shape = 22)) +
    labs(y = "log normalized counts") +
    theme(
      panel.background = element_rect(fill = fade_0),
      strip.text       = element_text(size = txt_pt2)
    )
  
  # Fix plot rows
  n_feats  <- length(feats)
  legd_buf <- 0.45          # This is to account for the height of legend
  n_row <- ceiling(n_feats / panels_n_col)
  
  if (n_row < panels_n_row) {
    frac <- (n_row + legd_buf) / panels_n_row
    
    res <- plot_grid(
      res, NULL,
      ncol = 1,
      rel_heights = c(frac, 1 - frac)
    )
  }
  
  res
}

#' Create GO plot
#' 
#' @param GO_df data.frame of GO results.
#' @param plot_colors Plot colors.
#' @param n_terms Number of terms to label.
#' @param plot_theme Base plot theme to use for GO plots.
#' @return ggplot object
#' @export
create_bubbles <- function(go_df, plot_colors = NULL, n_terms = 15,
                           plot_theme = base_theme) {
  
  # Check for empty inputs
  if (is_empty(go_df) || nrow(go_df) == 0) {
    res <- ggplot() +
      geom_blank() +
      theme(panel.background = element_blank())
    
    return(res)
  }
  
  # Shorten GO terms and database names
  go_nms <- c(
    "GO:BP" = "Biological\nProcess",
    "GO:CC" = "Cellular\nComponent",
    "GO:MF" = "Molecular\nFunction",
    "KEGG"  = "KEGG"
  )
  
  go_dat <- go_df %>%
    mutate(
      term_id = str_remove(term_id, "(GO|KEGG):"),
      term_id = str_c(term_id, " ", term_name),
      term_id = str_to_lower(term_id),
      term_id = str_trunc(term_id, 40, "right"),
      source  = go_nms[source]
    )
  
  # Reorder database names
  plot_lvls <- unname(go_nms)
  plot_lvls <- plot_lvls[plot_lvls %in% go_dat$source]
  plot_lvls <- unique(c(plot_lvls, go_dat$source))
  
  go_dat <- go_dat %>%
    mutate(source = fct_relevel(source, plot_lvls))
  
  # Extract top terms for each database
  top_go <- go_dat %>%
    group_by(source) %>%
    arrange(p_value) %>%
    dplyr::slice(1:n_terms) %>%
    ungroup()
  
  # Create bubble plots
  go_cols <- c("#D55E00", "#0072B2", "#009E73", "#6A51A3")
  
  if (!is.null(plot_colors)) {
    go_cols <- plot_colors
  }
  
  res <- go_dat %>%
    ggplot(aes(1.25, -log10(p_value), size = intersection_size)) +
    geom_point(aes(color = source), alpha = 0.5, show.legend = TRUE) +
    
    geom_text_repel(
      aes(2, -log10(p_value), label = term_id),
      data         = top_go,
      size         = 2.3,
      direction    = "y",
      hjust        = 0,
      segment.size = NA
    ) +
    
    guides(color = FALSE) +
    xlim(1, 8) +
    labs(y = "-log10(p-value)") +
    
    scale_color_manual(values = go_cols) +
    plot_theme +
    theme(
      aspect.ratio    = 0.9,
      legend.position = "bottom",
      axis.title.x    = element_blank(),
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank(),
      axis.line.x     = element_blank()
    ) +
    facet_wrap(~ source, scales = "free", nrow = 1)
  
  res
}

#' Format labels for plotting
#' 
#' @param so_in Seurat object or data.frame.
#' @param lab_clmn meta.data column to use for generating new labels.
#' @param char_vec Named character vector specifying new labels. Names should
#' correspond to old labels present in the Seurat object.
#' @param expr_vec Named character vector with strings that should be converted
#' to expressions to use for the new labels. Names should correspond to old
#' labels present in the Seurat object.
#' @param add_n Should the number of cells be added to each new label.
#' @param sep Separator to use when add_n is TRUE. 
#' @return Named vector containing new labels.
#' @export
format_labels <- function(so_in, lab_clmn, char_vec = NULL, expr_vec = NULL,
                          add_n = TRUE, l = " (", r = ")") {
  
  if ("Seurat" %in% class(so_in)) {
    so_in <- so_in@meta.data
  }
  
  char_labs <- so_in[[lab_clmn]] %>%
    as.character() %>%
    unique() %>%
    sort() %>%
    set_names()
  
  # Get n labels
  if (add_n) {
    char_labs <- so_in %>%
      get_nlab_fun(
        lab_col = lab_clmn,
        l       = l,
        r       = r
      )
    
    char_labs <- char_labs()
  }
  
  # Adjust labels based on char_vec
  if (!is.null(char_vec)) {
    char_vec %>%
      iwalk(~ {
        lbs <- char_labs[str_detect(names(char_labs), .y)]
        nms <- names(lbs)
        
        char_labs[nms] <<- str_replace(lbs, .y, .x)
      })
  }
  
  # Adjust labels based on expr_vec
  if (!is.null(expr_vec)) {
    expr_vec <- expr_vec %>%
      map_chr(str_remove, "\"$")
    
    expr_labs <- NULL
    
    expr_vec %>%
      iwalk(~ {
        nms <- names(char_labs)
        nms <- char_labs[str_detect(nms, .y)]
        
        expr_labs[names(nms)] <<- unname(nms) %>%
          str_replace(.y, .x) %>%
          str_c("\"") %>%
          parse(text = .)
      })
    
    expr_labs <- expr_labs[names(expr_labs) != ""]
    char_labs <- char_labs[!names(char_labs) %in% names(expr_labs)]
    char_labs <- c(char_labs, expr_labs)
  }
  
  char_labs
}

#' Calculate pairwise p-values for cell subtypes
#' 
#' @param sobj_in Seurat object or data.frame.
#' @param sample_name Sample name to add to results table.
#' @param data_column meta.data column containing data to use for test.
#' @param type_column meta.data column containing cell subtypes to test.
#' @return tibble
#' @export
calc_p_vals <- function(sobj_in, sample_name = NULL, data_column, type_column) {
  
  # Helper to calc p-values for give cell types
  calc_p <- function(x, y, df_in, data_column, type_column) {
    
    x_dat <- df_in %>%
      filter(!!sym(type_column) == x) %>%
      pull(data_column)
    
    y_dat <- df_in %>%
      filter(!!sym(type_column) == y) %>%
      pull(data_column)
    
    res <- wilcox.test(
      x        = x_dat, 
      y        = y_dat,
      conf.int = TRUE,
      exact    = FALSE
    ) %>%
      tidy() %>%
      mutate(
        group1 = x,
        group2 = y
      )
    
    res
  }
  
  # Get object meta.data
  if ("Seruat" %in% class(sobj_in)) {
    sobj_in <- sobj_in@meta.data %>%
      as_tibble(rownames = "cell_id")
  }
  
  p_data <- sobj_in %>%
    mutate("cell_type" = as.character(!!sym(type_column)))
  
  # Calculate additional stats for table
  p_stats <- p_data %>%
    group_by(cell_type) %>%
    summarize(
      n_cells = n(),
      med     = median(!!sym(data_column)),
      .groups = "drop"
    ) %>%
    mutate(frac_cells = n_cells / sum(n_cells))
  
  # Create data.frame with all combinations of cell types for comparison
  # Flip group columns for B and T cells so group1 includes all cell types
  c_types <- p_data %>%
    pull(cell_type) %>%
    unique()
  
  c_comps <- gtools::combinations(
    n = length(c_types),
    r = 2,
    v = c_types
  ) %>%
    as_tibble(.name_repair = "unique") %>%
    rename(x = ...1, y = ...2)
  
  # Calculate pairwise wilcox test for all combinations of cell types
  res <- c_comps %>%
    pmap_dfr(
      calc_p,
      df_in       = p_data,
      data_column = data_column,
      type_column = type_column
    )
  
  # Apply multiple testing correction
  res <- res %>%
    mutate(p_adj = p.adjust(p.value, method = "bonf"))
  
  # Add medians to data.frame
  res <- p_stats %>%
    rename(med_1 = med) %>%
    right_join(res, by = c(cell_type = "group1"))
  
  res <- p_stats %>%
    select(cell_type, med_2 = med) %>%
    right_join(res, by = c("cell_type" = "group2"))
  
  # Format final table
  final_cols <- c(
    "Cell type 1"  = "cell_type.y",
    "n cells 1"    = "n_cells",
    "frac cells 1" = "frac_cells",
    "Median 1"     = "med_1",
    "Cell type 2"  = "cell_type",
    "Median 2"     = "med_2",
    
    "method",    "alternative",
    "statistic", "estimate",
    "conf.low",  "conf.high",
    "p.value",   "p_adj"
  )
  
  if (!is.null(sample_name)) {
    res <- res %>%
      mutate(Sample = sample_name)
    
    final_cols <- c("Sample", final_cols)
  }
  
  res <- res %>%
    select(all_of(final_cols)) %>%
    arrange(`Cell type 1`, p_adj)
  
  res
}

#' Create gene figure with UMAP and boxplot panels
create_gene_fig <- function(df_in, grps, feat, feat_clrs = c("white", "#D7301F"), box_clrs,
                            size = 0.0001, stroke = 0.75) {
  
  # Create UMAP
  u <- df_in %>%
    arrange(!!sym(feat)) %>%
    ggplot(aes(UMAP_1, UMAP_2, fill = !!sym(feat))) +
    geom_point_trace(size = size, stroke = stroke) +
    scale_fill_gradientn(colours = feat_clrs) +
    guides(fill = guide_colorbar(
      title.position = "top",
      barheight      = unit(6, "pt"),
      barwidth       = unit(120, "pt"),
      ticks          = FALSE
    )) +
    theme_void() +
    theme(
      plot.margin     = margin(t = 15),
      aspect.ratio    = 0.7,
      legend.position = "top",
      legend.title    = element_text(size = 18)
    )
  
  # Set boxplot order
  lvls <- df_in %>%
    group_by(!!sym(grps)) %>%
    summarize(
      stats = list(boxplot.stats(!!sym(feat))),
      stats = map(stats, pluck, "stats"),
      med   = map_dbl(stats, pluck, 3),
      q3    = map_dbl(stats, pluck, 4),
      q4    = map_dbl(stats, pluck, 5),
      max   = max(!!sym(feat))
    ) %>%
    arrange(med, q3, q4, max) %>%
    pull(grps) %>%
    rev() %>%
    as.character()
  
  df_in <- df_in %>%
    mutate(!!sym(grps) := fct_relevel(!!sym(grps), lvls))
  
  # Create boxplots
  bx <- df_in %>%
    ggplot(aes(!!sym(grps), !!sym(feat), fill = !!sym(grps), color = !!sym(grps))) + 
    geom_boxplot(outlier.size = 0.3, alpha = 0.4) +
    scale_fill_manual(values = box_clrs) +
    scale_color_manual(values = box_clrs) +
    labs(x = "cluster") +
    djvdj_theme() +
    theme(
      plot.margin     = margin(r = 15, b = 15, l = 15),
      legend.position = "none",
      axis.title.x    = element_blank(),
      axis.title.y    = element_blank()
      # axis.title.x    = element_text(size = txt_pt1),
    )
  
  # Create final figure
  res <- plot_grid(
    u, bx,
    ncol        = 1,
    rel_heights = c(1, 0.5),
    align       = "v",
    axis        = "rl"
  )
  
  res
}


# Marker figures ----

#' Run gprofiler
#' 
#' @param gene_list List of input genes.
#' @param genome Genome to use for identifying GO terms.
#' @param gmt_id GMT ID to use instead of genome name.
#' @param p_max Maximum adjusted p-value for GO term to be included in results.
#' @param GO_size Minimum term size for GO term.
#' @param intrsct_size Minimum number of genes that intersect GO term.
#' @param order_query Run an ordered query where more significant genes are at
#' the top of the list.
#' @param dbases GO databases to query.
#' @param file_path File path to write results.
#' @param overwrite Overwrite existing file.
#' @param ... Additional arguments to pass to gprofiler2.
#' @return tibble containing GO terms.
#' @export
run_gprofiler <- function(gene_list, genome = NULL, gmt_id = NULL, p_max = 0.05, GO_size = 10, intrsct_size = 10, 
                          order_query = FALSE, dbases = c("GO:BP", "GO:MF", "KEGG"), file_path = NULL,
                          overwrite = FALSE, ...) {
  
  # Check for empty gene list
  if (is_empty(gene_list)) {
    return(as_tibble(NULL))
  }
  
  # Check arguments
  if (is.null(genome) && is.null(gmt_id)) {
    stop("Must specifiy genome or gmt_id.")
  }
  
  if (!is.null(file_path) && file.exists(file_path) && !overwrite) {
    warning("Loading existing file: ", file_path, ". Set overwrite = TRUE to rerun.")
    
    return(read_tsv(file_path))
  }
  
  # Use gmt id
  if (!is.null(gmt_id)) {
    genome <- gmt_id
    dbases <- NULL
  }
  
  # Run gProfileR
  res <- gene_list %>%
    gost(
      organism       = genome,
      sources        = dbases,
      domain_scope   = "annotated",
      significant    = TRUE,
      user_threshold = p_max,
      ordered_query  = order_query,
      ...
    )
  
  # Format and sort output data.frame
  res <- as_tibble(res$result)
  
  if (!is_empty(res)) {
    
    if (!is.null(GO_size)) {
      res <- res %>%
        filter(term_size > GO_size)
    }
    
    if (!is.null(intrsct_size)) {
      res <- res %>%
        filter(intersection_size > intrsct_size)
    }
    
    res <- res %>%
      arrange(source, p_value)
  }
  
  # Write table
  if (!is.null(file_path) && nrow(res) > 0) {
    res %>%
      select(-parents) %>%
      write_tsv(file_path)    
  }
  
  res
}

#' Find markers with presto
#' 
#' @param sobj_in Seurat object.
#' @param grp_column meta.data column containing labels to use for grouping cells.
#' @param exclude_grp Name of cell group to exclude from comparison.
#' @param groups_use Vector of cell groups to use for comparison.
#' @param filt_group Filter results to only include markers for groups specified by filt_group.
#' @param fc_range Vector of length 2 containing minimum and maximum log fold change for markers.
#' @param auc_min Minimum AUC for markers.
#' @param p_max Maximum adjusted p-value for filtering markers.
#' @param pct_in_min Minimum percentage of cells within the group that express each marker.
#' @param pct_out_max Maximum percentage of cells outside the group that express each marker.
#' @param file_path File path to write results.
#' @param overwrite Overwrite existing file specified by file_path.
#' @param ... Additional arguments to pass to wilcoxauc.
#' @return tibble
#' @export
find_markers <- function(sobj_in, grp_column = NULL, exclude_grp = NULL,
                         groups_use = NULL, filt_group = NULL,
                         fc_range = c(0.25, Inf), auc_min = 0.5, p_max = 0.05,
                         pct_in_min = 50, pct_out_max = 100, file_path = NULL,
                         overwrite = FALSE, ...) {
  
  if (!is.null(file_path) && file.exists(file_path) && !overwrite) {
    warning("Loading existing file: ", file_path, ". Set overwrite = TRUE to rerun.")
    
    return(read_tsv(file_path))
  }
  
  if (!is.null(exclude_grp) && is.null(groups_use)) {
    sobj_in <- sobj_in %>%
      subset(subset = !!sym(grp_column) != exclude_grp)
  }
  
  res <- sobj_in %>%
    wilcoxauc(
      group_by   = grp_column,
      groups_use = groups_use,
      ...
    ) %>%
    as_tibble() %>%
    filter(
      padj    < p_max,
      logFC   > fc_range[1],
      logFC   < fc_range[2],
      auc     > auc_min,
      pct_in  > pct_in_min,
      pct_out < pct_out_max
    ) %>%
    arrange(desc(logFC)) %>%
    rename(gene = feature)
  
  if (!is.null(filt_group)) {
    res <- res %>%
      filter(group %in% filt_group)
  }
  
  # Write table
  if (!is.null(file_path) && nrow(res) > 0) {
    write_tsv(res, file_path)
  }
  
  res
}

#' Run FindConservedMarkers for all cell identities
#' 
#' @param sobj_in Seurat object.
#' @param ident_1 Identity class to define markers for.
#' @param idnet_2 Second identity class for comparison. If NULL use all other
#' cells for comparison.
#' @param grp_var Grouping variable for finding markers.
#' @param p_max Maximum p-value for filtering marker genes.
#' @param min_fc Minimum log2 fold change for filtering marker genes.
#' @param fc_range Vector of length 2 containing minimum and maximum log2 fold
#' change for markers. Use this argument to select positive or negative
#' markers.
#' @param file_path File path to write results.
#' @param filt_regex Remove genes from results that match the provided regular
#' expression
#' @param overwrite Overwrite existing file specified by file_path.
#' @return tibble containing marker genes.
#' @export
find_conserved_markers <- function(sobj_in, ident_1 = NULL, ident_2 = NULL,
                                   grp_var, p_max = 0.05, fc_range = c(log2(1.25), Inf),
                                   file_path = NULL, filt_regex = NULL,
                                   overwrite = FALSE) {
  
  if (is.null(ident_1) && !is.null(ident_2)) {
    stop("Must specify ident_1 when ident_2 is given.")
  }
  
  if (!is.null(file_path) && file.exists(file_path) && !overwrite) {
    warning("Loading existing file: ", file_path, ". Set overwrite = TRUE to rerun.")
    
    return(read_tsv(file_path))
  }
  
  if (length(fc_range) != 2) {
    stop("fc_range must be a vector of length 2")
  }
  
  if (!fc_range[1] < fc_range[2]) {
    stop("The first value of fc_range must be smaller than the second value")
  }
  
  if (is.null(ident_1)) {
    ident_1 <- unique(Idents(sobj_in))
  }
  
  # Set FC filtering parameters
  only_pos <- FALSE
  fc_lim   <- 0
  
  if (all(fc_range > 0)) {
    only_pos <- TRUE
    fc_lim   <- fc_range[1]
    
  } else if (all(fc_range < 0)) {
    fc_lim <- abs(fc_range[2])
  }
  
  # Find markers
  res <- ident_1 %>%
    map_dfr(~ {
      mks <- sobj_in %>%
        FindConservedMarkers(
          ident.1         = .x,
          ident.2         = ident_2,
          grouping.var    = grp_var,
          logfc.threshold = fc_lim,
          only.pos        = only_pos
        ) %>%
        rownames_to_column("gene")
      
      # Filter and format output table
      fc_clmns <- grep("_avg_log2FC$", colnames(mks), value = TRUE)
      
      if (nrow(mks) > 0) {
        mks <- mks %>%
          rowwise() %>%
          mutate(
            avg_log2FC = mean(!!!syms(fc_clmns)),
            ident_1    = .x,
            ident_2    = ident_2
          ) %>%
          filter(
            across(
              all_of(fc_clmns),
              ~ .x > fc_range[1] & .x < fc_range[2]
            ),
            max_pval < p_max
          ) %>%
          ungroup()
      }
      
      mks
    })
  
  # Remove genes that match filt_regex
  if (!is.null(filt_regex)) {
    res <- res %>%
      filter(!grepl(filt_regex, gene))
  }
  
  # Write table
  if (!is.null(file_path) && nrow(res) > 0) {
    write_tsv(res, file_path)
  }
  
  res
}

#' Find marker genes and create summary figure
#' 
#' @param sobj_in Seurat object.
#' @param marker_fn Function to use for finding markers.
#' @param n_genes Number of marker genes to plot.
#' @param box_x Group to use for x-axis of boxplots.
#' @param box_nrow Number of rows for arranging boxplots.
#' @param box_ncol Number of columns for arranging boxplots.
#' @param box_aspect_ratio Aspect ratio to set for boxplots.
#' @param box_cols Colors for boxplots.
#' @param go_genome Genome to use for finding GO terms.
#' @param go_cols Colors for GO plots.
#' @param rel_h Relative height of boxplot and GO term panels.
#' @param plot_vln Should violin plots be generated instead of boxplots?
#' @param show_go Should panels showing GO terms be included?
#' @param show_x_labels Should x-axis labels be shown?
#' @param file_dir Directory to save output files.
#' @param file_prefix Prefix to add to output files.
#' @param overwrite Should gene lists be overwritten if they already exist?
#' @param gene_col Column in marker table containing gene names
#' @param ... Arguments to pass to marker_fn.
#' @return ggplot object
#' @export
create_marker_fig <- function(sobj_in, marker_fn, n_genes = 40, box_x,
                              box_nrow = 4, box_ncol = 10, box_cols = NULL,
                              go_genome = "mmusculus", n_terms = 10,
                              go_cols = NULL, rel_h = c(2, 0.5),
                              box_aspect_ratio = 1.3, plot_vln = FALSE,
                              show_go = FALSE, show_x_labels = FALSE,
                              file_dir = NULL, file_prefix = NULL,
                              overwrite = FALSE, gene_col = "gene", ...) {
  
  # Markers
  m_file <- file_dir
  
  if (!is.null(m_file)) {
    m_file <- str_c(file_prefix, "markers.tsv", sep = "_")
    m_file <- here(file_dir, m_file)
  }
  
  marks <- sobj_in %>%
    marker_fn(
      ...,
      file_path = m_file,
      overwrite = overwrite
    )
  
  N_MARKERS <<- nrow(marks)
  
  top_marks <- marks %>%
    head(n_genes) %>%
    pull(gene_col)
  
  # Boxplots
  boxes <- sobj_in %>%
    create_feat_boxes(
      feats         = top_marks,
      grp_column    = box_x,
      box_colors    = box_cols,
      panels_n_row  = box_nrow,
      panels_n_col  = box_ncol,
      plot_vln      = plot_vln,
      show_x_labels = show_x_labels,
      aspect.ratio  = box_aspect_ratio
    ) +
    theme(plot.margin = unit(c(0.2, 0.2, 1, 0.2), "cm"))
  
  # GO terms
  go_file <- file_dir
  
  if (!is.null(go_file)) {
    go_file <- str_c(file_prefix, "GO.tsv", sep = "_")
    go_file <- here(file_dir, go_file)
  }
  
  go_df <- marks %>%
    pull(gene_col) %>%
    run_gprofiler(
      genome      = go_genome,
      order_query = TRUE,
      file_path   = go_file,
      overwrite   = overwrite
    )
  
  N_GO_TERMS <<- nrow(go_df)
  
  # Return boxplots if show_go = FALSE
  if (!show_go) {
    return(boxes)
  }
  
  # Create dot plots
  go_dot <- go_df %>%
    create_bubbles(
      n_terms     = n_terms,
      plot_colors = go_cols
    ) +
    theme(legend.title = element_text(size = txt_pt2))
  
  # Create final figure
  res <- plot_grid(
    boxes, go_dot,
    rel_heights = rel_h,
    ncol        = 1
  )
  
  res
}

#' Create summary figures for treatment group
#' 
#' @param so_in Seurat object.
#' @param treat Treatment group to filter for. The provided label should be
#' present in the 'treatment' meta.data column. If NULL, all cells will be
#' included.
#' @param grp_clmn meta.data column containing groups. Markers will be
#' identified for each group in grp_clmn.
#' @param rep_clmn meta.data column containing replicate IDs. This is used to 
#' identify markers shared across replicates.
#' @param p_max Maximum p-value for filtering marker genes.
#' @param fc_range Vector of length 2 containing minimum and maximum log2 fold
#' change for markers. Use this argument to select positive or negative
#' markers.
#' @param file_dir Directory to save output files.
#' @param file_prefix Prefix to add to output files.
#' @param box_cols Character vector containing colors to use for boxplots.
#' @param n_genes Number of genes to plot.
#' @param box_nrow Number of rows for arranging boxplots.
#' @param box_ncol Number of columns for arranging boxplots.
#' @param box_aspect_ratio Aspect ratio to set for boxplots.
#' @param rel_height Relative height of boxplot and GO term panels.
#' @param footer Text to include after plot. This can be used to add spacing
#' between the next section of the document.
#' @param ... Additional arguments to pass to create_marker_fig().
#' @return ggplot object
#' @export
create_treatment_marker_fig <- function(so_in, treat = NULL, grp_clmn, rep_clmn = "rep", p_max = 0.05,
                                        fc_range = c(log2(1.25), Inf), file_dir = NULL, file_prefix = NULL,
                                        box_cols, n_genes = 48, box_nrow = 8, box_ncol = 6, header = "####",
                                        footer = "\n\n<br>\n\n", ...) {
  
  # Input data
  so_m <- so_in
  
  if (!is.null(treat)) {
    so_m <- so_in %>%
      subset(treatment == treat)
  }
  
  Idents(so_m) <- so_m %>%
    FetchData(grp_clmn)
  
  # Cell types to test
  grps <- so_m %>%
    get_cell_types(grp_clmn, "sample")
  
  # Create figure panels
  for (grp in grps) {
    cat("\n", header, " ", grp, "\n", sep = "")
    
    prfx <- str_c(c(file_prefix, grp), collapse = "_")
    prfx <- str_replace_all(prfx, " ", "_")
    prfx <- str_remove_all(prfx, "[\\)\\(]")
    
    marker_fig <- so_m %>%
      create_marker_fig(
        marker_fn   = find_conserved_markers,
        ident_1     = grp,
        grp_var     = rep_clmn,
        p_max       = p_max,
        fc_range    = fc_range,
        box_x       = grp_clmn,
        box_cols    = box_cols,
        n_genes     = n_genes,
        box_ncol    = box_ncol,
        box_nrow    = box_nrow,
        file_dir    = file_dir,
        file_prefix = prfx,
        ...
      )
    
    abs_fc <- fc_range[is.finite(fc_range)]
    abs_fc <- abs(2 ^ abs_fc)
    
    cat(
      N_MARKERS, "marker genes with a p-value <", p_max,
      "and an absolute fold change >", abs_fc, "were identified.",
      N_GO_TERMS, "GO terms were identified."
    )
    
    print(marker_fig)
    cat(footer)
  }
}

#' Create summary figures comparing treatments within each group
#' 
#' e.g. use this to compare mock vs CHIKV for each cell type or cluster
#' 
#' @param so_in Seurat object.
#' @param treat Treatment group to filter for. The provided label should be
#' present in the 'treatment' meta.data column.
#' @param grp_clmn meta.data column containing groups. Markers will be
#' identified by comparing treatments for each group in grp_clmn.
#' @param rep_clmn meta.data column containing replicate IDs. This is used to 
#' identify markers shared across replicates.
#' @param p_max Maximum p-value for filtering marker genes.
#' @param fc_range Vector of length 2 containing minimum and maximum log2 fold
#' change for markers. Use this argument to select positive or negative
#' markers.
#' @param file_dir Directory to save output files.
#' @param file_prefix Prefix to add to output files.
#' @param box_cols Character vector containing colors to use for boxplots.
#' @param n_genes Number of genes to plot.
#' @param box_nrow Number of rows for arranging boxplots.
#' @param box_ncol Number of columns for arranging boxplots.
#' @param box_aspect_ratio Aspect ratio to set for boxplots.
#' @param rel_height Relative height of boxplot and GO term panels.
#' @param footer Text to include after plot. This can be used to add spacing
#' between the next section of the document.
#' @param ... Additional arguments to pass to create_marker_fig().
#' @return ggplot object
#' @export
create_treatment_vs_marker_fig <- function(so_in, treat, grp_clmn, rep_clmn = "rep", p_max = 0.05,
                                           fc_range = c(log2(1.25), Inf), file_dir = NULL, file_prefix = NULL,
                                           box_cols, n_genes = 80, box_nrow = 8, box_ncol = 10,
                                           footer = "\n\n<br>\n\n", ...) {
  
  # Input data
  Idents(so_in) <- so_in %>%
    FetchData("treatment")
  
  treat_nm <- levels(unique(Idents(so_in)))
  treat_nm <- c(treat, treat_nm[treat_nm != treat])
  treat_nm <- str_c(treat_nm, collapse = "-")
  
  # Cell types to test
  grps <- so_in %>%
    get_cell_types(grp_clmn, "sample")
  
  # Create figure panels
  for (grp in grps) {
    cat("\n##### ", grp, "\n", sep = "")
    
    prfx <- str_c(file_prefix, treat_nm, grp, sep = "_")
    prfx <- str_replace_all(prfx, " ", "_")
    prfx <- str_remove_all(prfx, "[\\)\\(]")
    
    marker_fig <- so_in %>%
      subset(subset = !!sym(grp_clmn) == grp) %>%
      
      create_marker_fig(
        marker_fn   = find_conserved_markers,
        ident_1     = treat,
        grp_var     = rep_clmn,
        p_max       = p_max,
        fc_range    = fc_range,
        box_x       = "sample",
        box_cols    = box_cols,
        n_genes     = n_genes,
        box_nrow    = box_nrow,
        box_ncol    = box_ncol,
        file_dir    = file_dir,
        file_prefix = prfx,
        ...
      )
    
    abs_fc <- fc_range[is.finite(fc_range)]
    abs_fc <- abs(2 ^ abs_fc)
    
    cat(
      N_MARKERS, "marker genes with a p-value <", p_max,
      "and an absolute fold change >", abs_fc, "were identified.",
      N_GO_TERMS, "GO terms were identified."
    )
    
    print(marker_fig)
    cat(footer)
  }
}







