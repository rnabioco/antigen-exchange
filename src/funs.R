
# Utility ----

.check_mtime <- function(upstream_file, new_file) {
  make_file <- !file.exists(new_file)
  
  if (!make_file) {
    up_mtime <- file.info(upstream_file)$mtime
    mtime    <- file.info(new_file)$mtime
    
    if (is.na(up_mtime)) cli_abort("{.file {upstream_file}} does not exist")
    
    make_file <- up_mtime > mtime
  }
  
  make_file
}

#' Save figure as TIFF
#' 
#' @param plt Plot to save
#' @param prefix File prefix to use for naming
#' @param w Width
#' @param h Height
#' @param mtplyr Fraction to increase final image size by, this is helpful
#' since for TIFFs the font size appears slightly larger when saved 
#' @param dir Directory path to save file
#' @param dev Device to use for saving image
#' @export
.save_fig <- function(plt, filename, w = 13, h = 12, mtplyr = 0,
                      dir = here("results/figures/images"), dev = NULL) {
  
  if (!is.null(dev)) filename <- str_c(filename, ".", dev)
  
  plt %>%
    ggsave(
      filename = here(dir, filename),
      width    = w + (w * mtplyr),
      height   = h + (h * mtplyr),
      dpi      = 600,
      bg       = "white"
    )
}

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

#' Load slimmed down Seurat object
#' 
#' @param file File path to qs object
#' @param diet Return slimmed down object
load_obj <- function(file, diet = TRUE) {
  obj <- qread(file)
  
  if (diet) {
    obj <- obj %>%
      DietSeurat(
        assays = "RNA",
        counts = FALSE
      )
  }
  
  gc()
  
  obj
}

#' Export counts and meta.data tables
#' 
#' @param sobj_in Seurat object
#' @param assays Assays to include in counts matrix
#' @param feat_type Feature type for specified assay, provide a feature type
#' for each assay
#' @param gene_prefix Prefix to add to gene names for counts matrix
#' @param columns meta.data columns to export
#' @param out_dir Output directory
#' @param file_prefix Prefix to add to output files
#' @return Counts and meta.data tables
#' @export
export_matrices <- function(sobj_in, assays = "RNA", feat_type = "Gene Expression",
                            gene_prefix = "", columns, out_dir, file_prefix = "") {
  
  if (length(gene_prefix) == 1) {
    gene_prefix <- rep(gene_prefix, length(assays))
  }
  
  # Format count matrices
  assays <- purrr::set_names(feat_type, assays)
  
  mat <- assays %>%
    purrr::imap(~ Seurat::GetAssayData(sobj_in, assay = .y, "counts"))
  
  gene_names <- mat %>%
    purrr::map(rownames) %>%
    c(recursive = TRUE, use.names = FALSE)
  
  feat_type <- mat %>%
    purrr::imap(~ rep(assays[[.y]], nrow(.x))) %>%
    c(recursive = TRUE, use.names = FALSE)
  
  mat <- mat %>%
    purrr::reduce(rbind)
  
  barcodes <- colnames(mat)
  
  # Create HDF5 file
  counts_out <- file.path(out_dir, str_c(file_prefix, "count_matrix.h5"))
  
  h5file <- H5File$new(counts_out, mode = "w")
  
  h5file$create_group("matrix")
  h5file$create_group("matrix/features")
  
  h5file[["matrix/barcodes"]]              <- barcodes
  h5file[["matrix/features/name"]]         <- gene_names
  h5file[["matrix/features/feature_type"]] <- feat_type
  
  # Save matrix data (in sparse format)
  # * convert to COO format (i, j, x)
  triplet <- summary(mat)
  
  h5file[["matrix/data"]]    <- triplet$x                       # Non-zero values
  h5file[["matrix/indices"]] <- triplet$i - 1                   # Adjust indices for 0-based indexing
  h5file[["matrix/indptr"]]  <- c(0, cumsum(table(triplet$j)))  # Compressed column pointer
  h5file[["matrix/shape"]]   <- c(nrow(mat), ncol(mat))         # Matrix dimensions
  
  h5file$close()
  
  # Write meta.data table
  meta_out <- file.path(out_dir, str_c(file_prefix, "metadata.tsv.gz"))
  
  sobj_in@meta.data %>%
    tibble::as_tibble(rownames = "cell_id") %>%
    dplyr::select(any_of(columns)) %>%
    readr::write_tsv(meta_out)
}

#' Run hypergeometric test
#' 
#' Arguments match those used for dhyper()
#' 
#' @param x number of white balls drawn
#' @param k number of total balls drawn
#' @param m number of white balls in urn
#' @param n number of black balls in urn
#' @param alt alternative hypothesis, 'greater' tests whether more white balls
#' were drawn than expected
#' @export
.calc_fisher <- function(x, k, m, n, alt = "two.sided") {
  tot <- m + n
  k   <- k - x
  m   <- m - x
  n   <- n - k
  
  # Example contingency table
  # the sum of the matrix should equal the total number of cells
  # 23  244  | 267
  # 51  3235 | 3286
  #
  # 74  3479 | 3553
  
  mat <- c(x, k, m, n) %>%
    matrix(nrow = 2)
  
  if (sum(mat) != tot) {
    stop(
      "To create contingency table, the following must be TRUE: ",
      "x + (k - x) + (m - x) + (n - k + x) == m + n"
    )
  }
  
  res <- mat %>%
    fisher.test(alternative = alt)
  
  res$p.value
}

#' Format p values for labels
#' 
#' modified from djvdj
.format_pvalue <- function(p, digits = 1, cutoffs = NULL, show_decimal = 0.1) {
  
  if (p == 0) {
    p <- str_c("italic(p) < 1*x*10^-16")
    
    return(p)
  }
  
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
  
  clsts <- so_in
  
  if (is(so_in, "Seurat")) {
    clsts <- so_in %>%
      FetchData(unique(c(feats, clst_col, type_col)))
  }
  
  num_feats <- clsts %>%
    keep(is.numeric) %>%
    colnames()
  
  num_feats <- feats[feats %in% num_feats]
  chr_feats <- feats[!feats %in% num_feats]
  
  clsts <- clsts %>%
    group_by(!!sym(clst_col)) %>%
    summarize(
      across(all_of(num_feats), summary_fn),
      across(all_of(chr_feats), unique),
      .groups = "drop"
    )
  
  n_clsts <- nrow(clsts)
  
  clsts <- clsts %>%
    filter({{filt}}) %>%
    pull(clst_col) %>%
    as.character()
  
  if (n_distinct(so_in[[clst_col]]) != n_clsts) {
    warning("multiple values of one of the `feats` of type character are present for some clusters")
  }
  
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

#' Integrate Seurat objects using harmony
#' 
#' @param obj Input object
#' @param group_vars Variable(s) to using for grouping cells
#' @noRd
integrate_sobjs <- function(so_in, group_vars = "mouse", dims = 1:40,
                            resolution = c(1, 3, 5, 10),
                            umap.method = "umap-learn") {
  
  rm_clmns <- c(str_c("UMAP_", 1:2), str_c("hUMAP_", 1:2))
  
  res <- so_in %>%
    mutate_meta(dplyr::select, -any_of(rm_clmns)) %>%
    FindVariableFeatures(nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA() %>%
    RunUMAP(dims = dims, umap.method = umap.method) %>%
    AddMetaData(FetchData(., c("UMAP_1", "UMAP_2")))
  
  res <- res %>%
    RunHarmony(
      group.by.vars = group_vars,
      dims = dims,
      plot_convergence = TRUE
    ) %>%
    RunUMAP(
      dims = dims,
      reduction      = "harmony",
      reduction.name = "humap",
      reduction.key  = "hUMAP_"
    ) %>%
    FindNeighbors(reduction = "harmony", dims = dims) %>%
    FindClusters(resolution = resolution) %>%
    AddMetaData(FetchData(., c("hUMAP_1", "hUMAP_2")))
  
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

#' Calculate antigen score by subtracting B/T cell signal
#' 
#' @param obj Input object
#' @param clmns Columns containing antigen counts
#' @param control_types Cell types to use for calculating control signal
#' @param type_clmn Column containing cell types, this should contain the types
#' provided to control_types
#' @param norm_factors Named vector where names are sample names found in
#' `sample_clmn` and values are normalization factors for each sample,
#' e.g. total counts expressed in millions. If NULL, signal is not transformed
#' prior to correcting for background.
#' @param sample_clmn Column containing sample names that match those provided
#' by `norm_factors`
#' @param q Quantile to use as cutoff for control signal
#' @param operation Function to use for normalizing by control signal, this 
#' should take two arguments, the first is the signal to normalize,
#' the second is the control signal.
#' @param slot Slot in object to fetch data, default is counts.
#' @param final_trans Final transformation to perform after normalizing and
#' correcting background
#' @param suffix Suffix to add to new columns
#' @noRd
calc_ag_score <- function(obj, clmns, control_types = c("B cells", "T cells"),
                          type_clmn = "cell_type", sample_clmn = "directory",
                          norm_factors = NULL, q = 0.75,
                          operation = `-`, slot = "counts",
                          final_trans = log1p, suffix = "_score"
) {
  
  # Fetch clmns
  assay             <- DefaultAssay(obj)
  DefaultAssay(obj) <- "ADT"
  
  obj <- obj %>%
    AddMetaData(FetchData(., clmns, slot = slot))
  
  DefaultAssay(obj) <- assay
  
  # Normalize using norm_factors
  # * if norm_factors not provided, values are not transformed
  norm_clmns <- clmns
  
  if (!is.null(norm_factors)) {
    norm_clmns <- str_c(clmns, "_norm")
    
    obj <- obj %>%  
      mutate_meta(~ {
        .x %>%
          group_by(!!sym(sample_clmn)) %>%
          mutate(
            across(all_of(unname(clmns)), ~ {
              fctr <- norm_factors[!!sym(sample_clmn)]
              .x / fctr
            },
            .names = "{.col}_norm"
            )
          ) %>%
          ungroup()
      })
  }
  
  # Calculate background signal
  if (!is.null(control_types)) {
    obj <- obj %>%
      mutate_meta(~ {
        .x %>%
          group_by(!!sym(sample_clmn)) %>%
          mutate(across(
            all_of(norm_clmns),
            ~ {
              missing <- control_types[!control_types %in% !!sym(type_clmn)]
              
              if (!is_empty(missing)) {
                cli_warn("{missing} w{?as/ere} not found in {type_clmn}.")
              }
              
              idx <- !!sym(type_clmn) %in% control_types  # rows to use
              res <- quantile(.x[idx], q)                 # set cutoff
              
              res
            },
            .names = str_c("{.col}_threshold")
          )) %>%
          ungroup()
      })
    
    # Subtract background signal
    obj <- obj %>%
      mutate_meta(~ {
        .x %>%
          group_by(!!sym(sample_clmn)) %>%
          mutate(across(
            all_of(norm_clmns),
            ~ {
              clmn   <- str_c(cur_column(), "_threshold")
              thresh <- obj@meta.data[cur_group_rows(), clmn]
              thresh <- unique(thresh)
              
              if (length(thresh) > 1) {
                cli_abort("Threshold must be a single value for each group")
              }
              
              res <- .x - thresh
              
              res[res < 0] <- 0
              
              final_trans(res)
            },
            .names = str_c("{.col}", suffix)
          )) %>%
          ungroup()
      })
    
    norm_clmns <- str_c(norm_clmns, suffix)
  }
  
  # Rename final columns
  final_clmns <- set_names(
    norm_clmns,
    names(clmns) %||% str_c(clmns, suffix)
  )
  
  res <- obj %>%
    mutate_meta(~ dplyr::rename(.x, !!!final_clmns))
  
  res
}

#' Format and normalize Ag meta.data columns
#' 
#' @param so_in Input object
#' @param norm_factors Named vector containing normalization factors for
#' `calc_ag_score()`
#' @param ag_key List of named vectors containing labels for each Ag tag for
#' each experiment. List names should be experiment name, vector names should
#' be new label for the Ag tag.
#' @noRd
format_ag_data <- function(so_in, norm_factors, ag_key, ...) {
  
  # Assign labels for Ag tags (ABPS2, ABPS4)
  # * the identity of each tag differs depending on the experiment
  res <- so_in %>%
    mutate_meta(~ {
      dat <- .x
      
      c("3wk", "6wk") %>%
        walk(~ {
          clmn <- str_c("Ag_", .x)
          
          dat <<- dat %>%
            group_by(exp) %>%
            mutate(
              !!sym(clmn) := ag_key[[as.character(cur_group())]][.x],
              
              !!sym(clmn) := case_when(
                !!sym(clmn) == "ABPS2" ~ ABPS2,
                !!sym(clmn) == "ABPS4" ~ ABPS4
              )
            ) %>%
            ungroup()
        })
      
      dat
    })
  
  # Calculate Ag scores
  ag_clmns <- c("ovalbumin", "Ag_3wk", "Ag_6wk")
  
  res <- res %>%
    mutate_meta(mutate, across(all_of(ag_clmns), ~ replace_na(.x, 0))) %>%
    calc_ag_score(
      clmns        = ag_clmns,
      type_clmn    = "cell_type",
      sample_clmn  = "directory",
      norm_factors = norm_factors,
      ...
    )
  
  # Format Ag score columns
  res <- res %>%
    mutate_meta(
      mutate,
      Ag_score = case_when(
        exp   == "exp-0" ~ ovalbumin_score,
        mouse == "3wk" ~ Ag_3wk_score,
        mouse == "6wk" ~ Ag_6wk_score,
        
        # for 6wk-3wk use highest score
        mouse == "6wk-3wk" & Ag_3wk_score >= Ag_6wk_score ~ Ag_3wk_score,
        mouse == "6wk-3wk" & Ag_6wk_score >= Ag_3wk_score ~ Ag_6wk_score
      ),
      
      # these columns are used for antigen-exchange plots
      Ag_score_3 = ifelse(vaccination == "dual", Ag_3wk_score, Ag_score),
      Ag_score_6 = ifelse(vaccination == "dual", Ag_6wk_score, Ag_score),
      
      tm_3 = ifelse(vaccination == "dual", 21, tm),
      tm_6 = ifelse(vaccination == "dual", 42, tm)
    )
}

#' Identify GO terms
#' 
#' By default all expressed genes are used as background for the cell type
#' This will include results for all terms with any overlap regardless of
#' significance, we need all terms for plotting enrichment scores
#' 
#' Terms with no overlap will not be included in results
#' 
#' @param genes Named list of differentially expressed genes for each cell type,
#' names should correspond to cell type labels found in type_clmn
#' @param so_in Seurat object so use for determining background gene set, if
#' `NULL` all genes are used
#' @param type_clmn Column in Seurat object containing cell types, these should
#' match the cell types provided by genes
#' @param max_term_size Maximum term size to include in results
#' @param min_term_size Minimum term size to include in results
#' @param db GO database to use, specify 'ALL' to use all three
#' @param simplify_terms Should GO terms be simplified to removed redundant
#' results, this has a significant effect on performance
#' @param n_bkgd Number of top expressed genes to include for the background
#' gene set
#' @param exclude_genes Regular expression to use for removing genes before
#' performing GO analysis, by default ribo protein and mito genes are excluded
#' @param org_db Bioconductor organism annotation database, should be class
#' OrgDb
#' @param sim_data Similarity data to use for simplifying GO terms, generated
#' using [clusterProfiler::godata()]
#' @param file Path to output file for saving results, do not include extension
get_go <- function(genes, so_in, cell_types = NULL, type_clmn = "subtype",
                   max_term_size = 750, min_term_size = 10, db = "BP",
                   simplify_terms = TRUE, n_bkgd = Inf,
                   exclude_genes = "^(Rp[sl]|mt-)", org_db = org.Mm.eg.db,
                   sim_data = go_sim_data, file = NULL) {
  
  # Check for save GO file
  # this function will save a tsv and a GO object
  obj_file <- str_c(file, ".qs")
  tbl_file <- str_c(file, ".tsv")
  
  if (!is.null(file) && file.exists(obj_file)) {
    cli::cli_alert("Loading file {.file {obj_file}}")
    
    return(qread(obj_file))
  }
  
  # Set cell types to use for background genes
  if (is.null(cell_types)) {
    cell_types <- names(genes)
    
  } else if (length(cell_types) != length(genes)) {
    cli_abort("The length of `cell_types` must match the length of `genes`")
  }
  
  # Identify background gene set for each cell type
  # this includes all genes with >0 counts for any cell
  bkgd <- genes %>%
    map2(cell_types, ~ {
      bkgd <- NULL
      
      # so_in is required to determine background gene list
      # if genes is not named, all cells in so_in will be used
      if (!is.null(so_in)) {
        if (is.character(.y) && (.y %in% so_in[[type_clmn]][[1]])) {
          bkgd <- subset(so_in, !!sym(type_clmn) == .y)
          
        } else {
          cli_warn("{.y} not in {type_clmn}, all genes used for background.")
          
          bkgd <- so_in
        }
        
        bkgd <- bkgd@assays$RNA@data %>%
          rowMeans() %>%
          sort(decreasing = TRUE) %>%
          head(n_bkgd)
        
        bkgd <- names(bkgd[bkgd > 0])
        
      } else {
        cli::cli_warn(
          "Background gene set can only be determined if so_in is provided"
        )
      }
      
      bkgd
    })
  
  # Identify GO terms for each gene list
  # * exclude specified genes
  # * the geneRatio column may not include some input genes if they are not
  #   included in the "GOALL" gene universe used by clusterProfiler
  go <- genes %>%
    imap(~ {
      .x <- .x[!grepl(exclude_genes, .x)]
      
      g <- .x %>%
        enrichGO(
          keyType      = "SYMBOL",
          OrgDb        = org_db,
          universe     = bkgd[[.y]],
          maxGSSize    = max_term_size,
          minGSSize    = min_term_size,
          pvalueCutoff = 1.1,
          qvalueCutoff = 1.1,
          ont          = db
        )
      
      # Add background genes for each term to object
      # * this column is formatted in the same way as the geneID column
      # * these genes could be used to calculate overall fold enrichment for
      #   each cell type for GO clusters, i.e. determine total number of unique
      #   background genes overlapping all terms in the cluster
      bkgd_gns <- g@result$ID %>%
        clusterProfiler::bitr(         # used to fetch gene symbols for GO terms
          fromType = "GOALL", toType = "SYMBOL",
          OrgDb = org.Mm.eg.db
        ) %>%
        filter(SYMBOL %in% bkgd[[.y]]) %>%
        split(.$GOALL)
      
      bkgd_gns <- bkgd_gns %>%
        map_chr(~ {
          .x$SYMBOL %>%
            unique() %>%
            str_c(collapse = "/")
        })
      
      g %>%
        mutate(bgID = bkgd_gns[ID])
    })
  
  # Merge GO objects
  if (!is.null(names(go))) {
    go <- merge_result(go)
    
    go@fun <- "enrichGO"
    
  } else if (length(go) == 1) {
    go <- go[[1]]
  }
  
  # Simplify terms
  # * this collapses terms that are very similar
  # * this is very slow
  if (simplify_terms && !is.list(go)) {
    go <- go %>%
      clusterProfiler::simplify(semData = sim_data)
  }
  
  # Calculate enrichment scores
  go <- go %>%
    mutate(
      n_ovlp       = as.numeric(str_extract(GeneRatio, "^[0-9]+")),
      tot_genes    = as.numeric(str_extract(GeneRatio, "[0-9]+$")),
      n_bg_ovlp    = as.numeric(str_extract(BgRatio,   "^[0-9]+")),
      tot_bg_genes = as.numeric(str_extract(BgRatio,   "[0-9]+$")),
      enrichment   = (n_ovlp / tot_genes) / (n_bg_ovlp / tot_bg_genes)
    )
  
  # Save results
  if (!is.null(file)) {
    qsave(go, obj_file)
    
    go@compareClusterResult %>%
      write_tsv(tbl_file)
  }
  
  go
}

# Plotting ----

#' Create boxplots comparing RF module scores for predicted Ag classes
#' 
#' Used for figure 2
create_rf_module_boxes <- function(df_in, module_clmns, module_title, clrs,
                                   p_method = wilcox.test, p_alt = "two.sided",
                                   p_hjust = c(0.7, 0.5), p_vjust = c(1.2, 2.4),
                                   n_label = FALSE, pred_clmn = "ml_pred_1_grp",
                                   ...
) {
  
  # Format data for boxplots  
  typs <- str_remove(module_clmns, "_[^_]+$")
  
  dat <- df_in %>%
    pivot_longer(all_of(module_clmns)) %>%
    mutate(module_type = str_remove(name, "_[^_]+$")) %>%
    filter(subtype == module_type) %>%
    mutate(
      subtype           = fct_relevel(subtype, typs),
      !!sym(pred_clmn) := fct_relevel(!!sym(pred_clmn), names(clrs))  # must set as factor to
    )                                                # calculate correlation
  
  # Calculate Spearman correlation
  p_dat <- dat %>%
    group_by(subtype, tm) %>%
    summarize(corr = list(tidy(cor.test(
      as.numeric(!!sym(pred_clmn)), value,
      method      = "spearman",
      continuity  = TRUE,
      exact       = FALSE,
      alternative = p_alt
    ))), .groups = "drop") %>%
    unnest(corr) %>%
    mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
    rowwise() %>%
    mutate(
      p_adj = .format_pvalue(p_adj, show_decimal = Inf),
      p_adj = str_c("italic(p) == ", p_adj),
      r_lab = str_c("italic(r[s]) == ", round(estimate, 2))
    ) %>%
    ungroup()
  
  # Format n labels
  # the order of n labels in n_dat determines order on plot
  # need to make sure n_dat is ordered based on factor level which is set
  # upstream
  n_dat <- dat %>%
    group_by(tm, subtype, !!sym(pred_clmn)) %>%
    summarize(n = n_distinct(.cell_id), .groups = "drop") %>%
    mutate(n_lab = label_comma()(n)) %>%
    arrange(!!sym(pred_clmn))
  
  # Create boxplots
  y_exp <- c(0.05, 0.25)
  
  if (n_label) y_exp[1] <- 0.1
  
  res <- dat %>%
    ggplot(aes(as.character(tm), value, fill = !!sym(pred_clmn))) +
    geom_boxplot(
      outlier.size = 0.1,
      position     = position_dodge2(preserve = "single"),
      key_glyph    = draw_key_point,
      ...
    ) +
    
    # r label
    geom_text(
      aes(y = Inf, fill = NULL, label = r_lab),
      data = p_dat,
      size = 10 / .pt,
      hjust = p_hjust[1],
      vjust = p_vjust[1],
      parse = TRUE
    ) +
    
    # p label
    geom_text(
      aes(y = Inf, fill = NULL, label = p_adj),
      data  = p_dat,
      size  = 10 / .pt,
      hjust = p_hjust[2],
      vjust = p_vjust[2],
      parse = TRUE
    ) +
    
    scale_fill_manual(values = clrs) +
    scale_y_continuous(expand = expansion(y_exp)) +
    facet_wrap(~ subtype, nrow = 1, scales = "free_y") +
    labs(x = "days post immunization", y = module_title) +
    guides(fill = guide_legend(override.aes = list(shape = 22, size = 5))) +
    base_theme +
    theme(
      aspect.ratio    = 0.85,
      plot.title      = element_text(hjust = 0.5),
      legend.position = "bottom",
      legend.title    = element_blank(),
      plot.margin     = margin(t = 0, r = 5, b = 5, l = 5)
    )
  
  # n label
  if (n_label) {
    res <- res +
      geom_text(
        aes(y = -Inf, fill = NULL, label = n_lab),
        data     = n_dat,
        position = position_dodge2(preserve = "single", width = 1),
        size     = 10 / .pt,
        vjust    = -0.25
      )
  }
  
  res
}

#' Create UMAPs showing module scores for each timepoint
#' 
#' Used for figures 2 and 5
create_rf_module_umaps <- function(df_in, typ, module_clmn, module_title,
                                   ag_class_clmn, pt_size = 1, facet_var = "tm",
                                   equal_scales = TRUE, outline = FALSE,
                                   circle_cells = TRUE, show_scale_title = FALSE,
                                   ag_clrs = c(
                                     "lightblue", "lightblue",
                                     "white", "#D7301F"
                                   ),
                                   class_clrs = c(
                                     `Ag-high` = "#D7301F",
                                     `Ag-low`  = "#56B4E9",
                                     other     = "white"
                                   ),
                                   lab_fn = as_labeller(function(x) str_c("day ", x)),
                                   umap_labels = TRUE) {
  
  # Format plot data
  if (length(pt_size) == 1) pt_size <- rep(pt_size, 2)
  
  dat <- df_in %>%
    mutate(
      !!sym(ag_class_clmn) := ifelse(
        subtype == typ,
        as.character(!!sym(ag_class_clmn)),
        as.character(NA)
      ),
      !!sym(ag_class_clmn) := fct_relevel(!!sym(ag_class_clmn), names(class_clrs)),
      !!sym(facet_var) := as.character(!!sym(facet_var))
    )
  
  # Set coordinates for drawing circle
  circ_dat <- dat %>%
    filter(subtype == typ) %>%
    group_by(!!sym(facet_var)) %>%
    summarize(across(c(hUMAP_1, hUMAP_2), mean), .groups = "drop")
  
  # Plot module scores
  grps         <- sort(unique(df_in[[facet_var]]))
  facet_frm    <- as.formula(str_c("~ ", facet_var))
  scale_limits <- guides <- NULL
  
  if (equal_scales) {
    scale_limits <- c(min(dat[[module_clmn]]), max(dat[[module_clmn]]))
    guides       <- "collect"
  }
  
  ag_mod_u <- grps %>%
    imap(~ {
      u <- dat %>%
        filter(!!sym(facet_var) == .x) %>%
        create_sig_umap(
          dat_col      = module_clmn,
          grp_col      = facet_var,
          clrs         = ag_clrs,
          pt_size      = pt_size[1],
          pt_stroke    = 0.7,
          scale_limits = scale_limits
        ) +
        facet_wrap(facet_frm, labeller = lab_fn) +
        theme(legend.text = element_text(size = 10))
      
      if (show_scale_title) {
        u <- u +
          guides(fill = guide_colorbar(
            ticks = FALSE, title.position = "left", title = module_title
          )) +
          theme(
            legend.title = element_text(size = ttl_pt2, hjust = 0.5),
            legend.key.height = unit(6, "pt"),
            legend.key.width  = unit(30, "pt")
          )
      }
      
      if (circle_cells) {
        u <- u +
          geom_circle(
            aes(x0 = hUMAP_1, y0 = hUMAP_2, fill = NULL, color = NULL, r = 2.5),
            data = filter(circ_dat, !!sym(facet_var) == .x),
            linetype = 2
          )
      }
      
      if (umap_labels) {
        u <- u +
          scale_y_continuous(expand = expansion(c(0.15, 0.05)))
        
        if (.y == 1) {
          u <- .add_umap_labels(u)
        } else {
          u <- .add_umap_labels(u, include = "x")
        }
      }
      
      u
    }) %>%
    wrap_plots(nrow = 1, guides = guides) &
    theme(legend.position = "bottom")
  
  if (!show_scale_title) {
    ag_mod_u <- ag_mod_u +
      plot_annotation(
        title = module_title,
        theme = theme(
          plot.title = element_text(hjust = 0.5, size = 18)
        )
      )
  }
  
  # Plot Ag-high cells
  ag_class_u <- grps %>%
    imap(~ {
      plt <- dat %>%
        filter(!!sym(facet_var) == .x) %>%
        create_grp_umap(
          dat_col  = ag_class_clmn,
          grp_col  = facet_var,
          clrs     = class_clrs,
          lvls     = names(class_clrs),
          size     = pt_size[2],
          na_color = "grey85"
        ) +
        facet_wrap(facet_frm, labeller = lab_fn) +
        theme(
          legend.position = c(0.5, -0.1),
          legend.direction = "horizontal"
        )
      
      if (umap_labels) {
        plt <- plt +
          scale_y_continuous(expand = expansion(c(0.15, 0.05)))
        
        if (.y == 1) {
          plt <- .add_umap_labels(plt)
        } else {
          plt <- .add_umap_labels(plt, include = "x")
        }
      }
      
      plt
    }) %>%
    wrap_plots(nrow = 1)
  
  # Return list of UMAPs
  res <- list(ag_mod_u, ag_class_u)
}

#' Create boxplots comparing single and dual vaccinated mice
#' 
#' Used for figures 3 and 4
create_exchange_boxes <- function(df_in, tm = NULL, tm_clmn, ag_clmn,
                                  typs = NULL, clrs, labs = waiver(),
                                  facet_vars = "subtype",
                                  p_method = wilcox.test, p_vars = NULL,
                                  p_hjust = 0.5, p_vjust = 0.5, p_size = 10,
                                  p_cutoff = Inf, include_n = TRUE) {
  
  # Set plot colors
  ln_clrs   <- clrs
  ln_clrs[] <- "black"
  lgnd_clr  <- "black"
  
  if (!is.null(tm)) {
    ln_clrs[names(ln_clrs) != tm] <- "grey65"
    clrs[names(clrs) != tm]       <- "white"
    lgnd_clr                      <- clrs[[tm]]
  }
  
  # Format plot data
  x_facet <- dplyr::last(facet_vars)
  
  dat <- df_in %>%
    mutate(
      !!sym(x_facet) := fct_reorder(
        !!sym(x_facet), !!sym(ag_clmn), mean,
        .desc = TRUE
      ),
      vaccination = fct_relevel(vaccination, vacs)
    )
  
  if (!is.null(typs)) {
    dat <- dat %>%
      filter(!!sym(x_facet) %in% typs)
  }
  
  # Data to draw line connecting timepoints
  ln_dat <- dat %>%
    filter(vaccination == "single") %>%
    group_by(!!!syms(c(tm_clmn, facet_vars))) %>%
    summarize(!!sym(ag_clmn) := median(!!sym(ag_clmn)), .groups = "drop")
  
  # p-values for vaccination
  y_max <- max(pull(dat, ag_clmn))
  
  p_dat <- dat %>%
    group_by(!!!syms(c(tm_clmn, facet_vars))) %>%
    filter(all(vacs %in% vaccination)) %>%
    summarize(
      pval = p_method(
        (!!sym(ag_clmn))[vaccination == vacs[1]],
        (!!sym(ag_clmn))[vaccination == vacs[2]],
        alternative = "less"
      )$p.value,
      !!sym(vacs[1]) := sum(vaccination == vacs[1]),
      !!sym(vacs[2]) := sum(vaccination == vacs[2]),
      y  = max(!!sym(ag_clmn)),
      y  = y + (y_max * 0.15),
      .groups = "drop"
    )
  
  # Adjust p-values separately for each replicate
  if (!is.null(p_vars)) {
    p_dat <- p_dat %>%
      group_by(!!!syms(p_vars))
  }
  
  p_dat <- p_dat %>%
    mutate(padj = p.adjust(pval, method = "BH")) %>%
    rowwise() %>%
    mutate(
      plab = .format_pvalue(padj),
      plab = str_c("italic(p) == ", plab)
    ) %>%
    ungroup() %>%
    filter(padj < p_cutoff)
  
  # Set facet variables
  if (length(facet_vars) == 1) facet_vars <- c("", facet_vars)
  
  facet_vars <- str_c(facet_vars, collapse = "~")
  
  # Create boxplots
  res <- dat %>%
    ggplot(aes(!!sym(tm_clmn), !!sym(ag_clmn))) +
    
    geom_line(
      data = ln_dat,
      color = "grey85", linewidth = 0.5, linetype = 2
    ) +
    geom_boxplot(
      aes(
        fill  = as.character(!!sym(tm_clmn)),
        color = as.character(!!sym(tm_clmn)),
        alpha = vaccination
      ),
      outlier.size = 0.2, width = 5, key_glyph = draw_key_point,
      position = position_dodge2(preserve = "single")
    ) +
    geom_text(
      aes(y = y, label = plab, alpha = NULL),
      show.legend = FALSE,
      data  = p_dat,
      parse = TRUE,
      size  = p_size / .pt,
      hjust = p_hjust,
      vjust = p_vjust
    ) +
    facet_grid(as.formula(facet_vars)) +
    scale_alpha_manual(
      values = c(single = 0.35, dual = 1),
      labels = labs
    ) +
    scale_fill_manual(values = clrs) +
    scale_color_manual(values = ln_clrs) +
    guides(
      fill = "none", color = "none",
      alpha = guide_legend(override.aes = list(
        shape = 15, size = 4, color = lgnd_clr
      ))
    ) +
    labs(x = "days post immunization", y = "Ag-score") +
    base_theme +
    theme(aspect.ratio = 0.7)
  
  if (include_n) {
    n_dat <- p_dat %>%
      pivot_longer(all_of(vacs), names_to = "vaccination", values_to = "n") %>%
      mutate(
        vaccination = fct_relevel(vaccination, vacs),
        n           = label_comma()(n),
        hjust       = if_else(vaccination == vacs[1], 1, 0)
      )
    
    res <- res +
      scale_y_continuous(expand = expansion(c(0.17, 0.05))) +
      geom_text(
        aes(y = -Inf, label = n, alpha = NULL, hjust = hjust),
        show.legend = FALSE,
        data  = n_dat,
        size  = p_size / .pt,
        vjust = -0.3,
        position = position_dodge2(width = 3)
      )
  }
  
  res
}

#' Create heatmaps comparing single and dual vaccinated mice
#' 
#' Used for figures 3 and 4
create_exchange_heats <- function(df_in, ag_clmn, ttl, clr,
                                  # label_fn = ~ str_replace(.x, " day", "d")
                                  label_fn = ~ str_remove(.x, "^Group ")
) {
  
  # Format data for heatmaps  
  dat <- df_in %>%
    mutate(subtype = fct_reorder(subtype, !!sym(ag_clmn), mean))
  
  # Create heatmaps
  res <- dat %>%
    ggplot(aes(mouse, subtype, fill = !!sym(ag_clmn))) +
    geom_tile() +
    labs(x = "immunization group") +
    scale_fill_gradientn(colours = c("lightblue", "white", clr)) +
    scale_x_discrete(labels = label_fn) +
    guides(fill = guide_colorbar(
      ticks = FALSE,
      title = ttl,
      barheight = unit(150, "pt"), barwidth = unit(10, "pt")
    )) +
    base_theme +
    theme(
      aspect.ratio = n_distinct(dat$subtype) / n_distinct(dat$mouse),
      plot.margin  = margin(0, 15, 0, 15),
      axis.title.y = element_blank(),
      axis.text.y  = element_text(size = ttl_pt2)
    )
  
  res
}

#' Create scatter plots comparing Ag-scores for dual vaccinated mice
#' 
#' Used for figure 3 
create_exchange_scatter <- function(df_in,
                                    x = c(Ag_3wk_score = "Ag-score (21 day)"),
                                    y = c(Ag_6wk_score = "Ag-score (42 day)"),
                                    clrs = lec_clrs) {
  
  # Set x y variables
  x_var <- names(x)
  y_var <- names(y)
  
  # Calculate correlation
  lab_dat <- df_in %>%
    group_by(subtype, exp) %>%
    summarize(
      cor_res = list(cor.test(!!sym(x_var), !!sym(y_var), method = "spearman")),
      cor     = map_dbl(cor_res, ~ .x$estimate[[1]]),
      p_val   = map_dbl(cor_res, ~ .x$p.value),
      n       = n_distinct(.cell_id),
      n       = label_comma()(n),
      .groups = "drop"
    ) %>%
    mutate(p_adj = p.adjust(p_val, method = "BH")) %>%
    rowwise() %>%
    mutate(
      p_lab = .format_pvalue(p_adj),
      r_lab = str_c("italic(r[s]) == ", round(cor, digits = 2)),
      n_lab = str_c("n = ", n)
    )
  
  # Create scatter plot
  res <- df_in %>%
    ggplot(aes(!!sym(x_var), !!sym(y_var), color = subtype)) +
    geom_point(size = 1.5) +
    geom_smooth(
      method    = "lm",
      linetype  = 2,
      linewidth = 0.5,
      color     = "black",
      se        = FALSE,
      formula   = "y ~ x"
    ) +
    geom_text(
      aes(x = -Inf, y = Inf, label = r_lab),
      data  = lab_dat,
      hjust = -0.1,
      vjust = 1.2,
      color = "black",
      size  = txt_pt2 / .pt,
      parse = TRUE
    ) +
    geom_text(
      aes(x = -Inf, y = Inf, label = n_lab),
      data  = lab_dat,
      hjust = -0.1,
      vjust = 2.7,
      color = "black",
      size  = txt_pt2 / .pt
    ) +
    facet_wrap(~ subtype, nrow = 1) +
    
    scale_color_manual(values = clrs) +
    labs(x = unname(x), y = unname(y)) +
    base_theme +
    theme(
      aspect.ratio = 1,
      legend.position = "none"
    )
  
  res
}

#' Create boxplots comparing RF module scores and Ag-classes 
#' 
#' Used for figure 3
create_module_boxes <- function(df_in, ag_class_clmn = "Ag_class_dual",
                                module_clmns, module_title,
                                clrs = lec_clrs, lvls, p_alt = "two.sided",
                                p_method = "correlation",
                                p_x = 1,
                                p_hjust = c(0, 0),
                                p_vjust = c(1.2, 1.9),
                                ...) {
  
  # Format plotting data
  # cell types filtered based on module_clmns names
  # module score is only plotted for the corresponding cell type
  # group all low cells together under single mouse, otherwise three boxplots
  # are shown
  typs <- str_remove(module_clmns, "_[^_]+$")
  
  dat <- df_in %>%
    pivot_longer(all_of(module_clmns)) %>%
    mutate(module_type = str_remove(name, "_[^_]+$")) %>%
    filter(subtype == module_type) %>%
    mutate(
      subtype = fct_relevel(subtype, typs),
      !!sym(ag_class_clmn) := fct_relevel(!!sym(ag_class_clmn), lvls)
    ) %>%
    arrange(!!sym(ag_class_clmn)) %>%
    mutate(
      mouse_class = fct_inorder(str_c(mouse, "-", subtype, "-", !!sym(ag_class_clmn)))
    )
  
  n_grps <- n_distinct(dat[[ag_class_clmn]])
  
  # Calculate p-values
  if (!identical(p_method, "correlation")) {
    p_dat <- dat %>%
      split(.$subtype) %>%
      map(~ {
        pairwise.wilcox.test(
          .x$value, .x[[ag_class_clmn]],
          alternative     = p_alt,
          p.adjust.method = "none"  # adjust below using all p-values
        ) %>%
          broom::tidy() %>%
          mutate(
            ymax = max(.x$value),
            ymin = min(.x$value)
          )
      }) %>%
      bind_rows(.id = "subtype") %>%
      filter(group1 == "double-high") %>%
      mutate(
        p_adj   = p.adjust(p.value, method = "BH"),
        ln_x    = match(group2, lvls),
        ln_xend = match(group1, lvls),
        y_fctr  = ((ymax - ymin) * 0.1),
        ln_y    = ymax + (y_fctr * ln_x),
        y       = ymax + (y_fctr * (ln_x + 1)),
        x       = ln_x + 0.5,
        ln_x    = ln_x + 0.1,
        ln_xend = ln_xend - 0.1
      ) %>%
      rowwise() %>%
      mutate(p_adj = .format_pvalue(p_adj)) %>%
      ungroup()
    
    # Calculate Spearman correlation
  } else {
    p_dat <- dat %>%
      group_by(subtype) %>%
      summarize(
        corr = list(tidy(cor.test(
          as.numeric(!!sym(ag_class_clmn)), value,
          method      = "spearman",
          continuity  = TRUE,
          exact       = FALSE,
          alternative = p_alt
        ))),
        .groups = "drop"
      ) %>%
      unnest(corr) %>%
      mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
      rowwise() %>%
      mutate(
        p_adj = .format_pvalue(p_adj, show_decimal = Inf),
        p_adj = str_c("italic(p) == ", p_adj),
        r_lab = str_c("italic(r[s]) == ", round(estimate, 2))
      ) %>%
      ungroup()
  }
  
  # Format n labels for x axis
  n_dat <- dat %>%
    group_by(!!sym(ag_class_clmn), mouse_class, subtype) %>%
    summarize(n = n(), .groups = "drop") %>%
    mutate(
      n = str_c(!!sym(ag_class_clmn), "\nn = ", scales::label_comma()(n))
    )
  
  n_lab <- set_names(n_dat$n, n_dat$mouse_class)
  
  # Set p label position
  # by default center on plot
  p_x <- p_x %||% ((n_grps / 2) + 0.5)
  
  if (length(p_x) == 1)     p_x     <- rep(p_x, 2)
  if (length(p_hjust) == 1) p_hjust <- rep(p_hjust, 2)
  if (length(p_vjust) == 1) p_vjust <- rep(p_vjust, 2)
  
  # Create boxplots
  x_idx   <- 1:n_grps
  x_hjust <- (x_idx - min(x_idx)) / (max(x_idx) - min(x_idx))
  # x_hjust <- rev(x_hjust)
  
  res <- dat %>%
    ggplot(aes(mouse_class, value, fill = subtype)) +
    geom_boxplot(outlier.size = 0.1, ...) +
    
    facet_wrap(~ subtype, scales = "free") +
    scale_fill_manual(values = clrs) +
    scale_x_discrete(labels = n_lab) +
    scale_y_continuous(expand = expansion(c(0.05, 0.07))) +
    labs(y = module_title) +
    base_theme +
    theme(
      aspect.ratio    = 1.5,
      legend.position = "none",
      axis.text.x     = element_text(angle = 90, hjust = 1, vjust = x_hjust),
      # axis.text.x     = element_text(angle = 90, hjust = 1),
      axis.title.x    = element_blank()
    )
  
  if (!identical(p_method, "correlation")) {
    res <- res +
      geom_segment(
        aes(x = ln_x, xend = ln_xend, y = ln_y, yend = ln_y),
        data = p_dat,
        color = ln_col
      ) +
      
      geom_text(
        aes(x = x, y = y, label = p_adj),
        # direction = "x", seed = 42, force = 0.1, force_pull = 0, box.padding = 0.5,
        data = p_dat,
        parse = TRUE,
        size = 9 / .pt
      )
    
  } else {
    res <- res +
      geom_text(
        aes(x = !!(p_x[1]), y = Inf, fill = NULL, label = r_lab),
        size = 10 / .pt,
        data = p_dat,
        hjust = p_hjust[1],
        vjust = p_vjust[1],
        parse = TRUE
      ) +
      geom_text(
        aes(x = !!(p_x[2]), y = Inf, fill = NULL, label = p_adj),
        size = 10 / .pt,
        data = p_dat,
        hjust = p_hjust[2],
        vjust = p_vjust[2],
        parse = TRUE
      )
  }
  
  res
}

#' Create scatter plots comparing RF module scores and Ag-scores
#'
#' Used for figure 4
create_module_scatter <- function(df_in, ag_clmns, module_clmns, module_title,
                                  ag_class_clmn = "Ag_class_dual", clrs,
                                  label_x = -Inf, label_hjust = -0.1,
                                  label_start = -0.5, label_pad = 0.15) {
  
  # Format plotting data
  # * cell types are filtered based on the module_clmns names
  typs <- str_remove(module_clmns, "_[^_]+$")
  
  dat <- df_in %>%
    pivot_longer(
      all_of(ag_clmns),
      names_to = "ag_clmn", values_to = "Ag-score"
    ) %>%
    pivot_longer(
      all_of(module_clmns),
      names_to = "module_name", values_to = "module_score"
    ) %>%
    mutate(module_type = str_remove(module_name, "_[^_]+$")) %>%
    filter(subtype == module_type) %>%
    mutate(
      subtype               = fct_relevel(subtype, typs),
      !!sym(ag_class_clmn) := fct_relevel(!!sym(ag_class_clmn), rev(names(clrs)))
    )
  
  # Format r label data
  lab_dat <- dat %>%
    group_by(ag_clmn) %>%
    mutate(
      min_score = min(`Ag-score`),
      max_score = max(`Ag-score`)
    ) %>%
    .calc_cor_pvals(
      format_data = FALSE,
      keep_columns = c("min_score", "max_score")
    ) %>%
    rowwise() %>%
    mutate(
      r_lab = str_c("italic(r[s]) == ", round(s_cor, digits = 2)),
      p_lab = str_c("italic(p) == ", .format_pvalue(s_p_adj)),
      n_lab = str_c("n == ", n)
    )
  
  # Set label coordinates
  lab_dat <- lab_dat %>%
    mutate(
      y_fctr = (max_score - min_score) * label_pad,
      r_y    = max_score + (y_fctr * (label_start + 2)),
      p_y    = max_score + (y_fctr * (label_start + 1)),
      n_y    = max_score + (y_fctr * label_start)
    )
  
  # Labelling function
  lab_fn <- function(x) {
    modify <- x %in% names(ag_labs_2)
    
    x[modify] <- ag_labs_2[x[modify]]
    
    x
  }
  
  # Create scatter plots
  res <- dat %>%
    ggplot(aes(module_score, `Ag-score`, color = !!sym(ag_class_clmn))) +
    geom_point(size = 0.7) +
    geom_smooth(
      method = "lm", formula = y ~ x,
      color = "black", linetype = 2, linewidth = 1, alpha = 0.1
    ) +
    facet_grid(
      ag_clmn ~ subtype,
      scales = "free", switch = "y",
      labeller = as_labeller(lab_fn)
    ) +
    scale_color_manual(values = clrs) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    labs(x = module_title) +
    base_theme +
    theme(
      aspect.ratio    = 1,
      strip.placement = "outside",
      legend.position = "bottom",
      legend.title    = element_blank(),
      axis.title.y    = element_blank()
    )
  
  # Add labels
  lab_args <- list(
    y = str_c(c("r", "p", "n"), "_y"),
    label = str_c(c("r", "p", "n"), "_lab")
  )
  
  lab_args %>%
    pwalk(~ {
      res <<- res +
        geom_text(
          aes(!!label_x, !!sym(.x), color = NULL, label = !!sym(.y)),
          data  = lab_dat,
          hjust = label_hjust,
          parse = TRUE,
          show.legend = FALSE
        )
    })
  
  res
}

.calc_cor_pvals <- function(df_in, ag_clmns, module_clmns, format_data = TRUE,
                            keep_columns = NULL) {
  
  if (format_data) {
    df_in <- df_in %>%
      pivot_longer(
        all_of(ag_clmns),
        names_to = "ag_clmn", values_to = "Ag-score"
      ) %>%
      pivot_longer(
        all_of(module_clmns),
        names_to = "module_name", values_to = "module_score"
      ) %>%
      mutate(module_type = str_remove(module_name, "_[^_]+$")) %>%
      filter(subtype == module_type)
  }
  
  .calc_cor_p <- function(df_in, x, y, prefix = "s", method = "spearman") {
    df_in %>%
      mutate(
        res = list(cor.test(!!sym(x), !!sym(y), method = method)),
        cor = map_dbl(res, ~ .x$estimate[[1]]),
        p   = map_dbl(res, ~ .x$p.value)
      ) %>%
      select(-res) %>%
      rename_with(.cols = all_of(c("cor", "p")), ~ str_c(prefix, "_", .x))
  }
  
  keep_clmns <- c(
    "subtype", "module_name", "ag_clmn",
    "s_cor", "s_p", "p_cor", "p_p", "n"
  )
  
  res <- df_in %>%
    group_by(subtype, module_name, ag_clmn) %>%
    .calc_cor_p("module_score", "Ag-score", prefix = "s", method = "spearman") %>%
    .calc_cor_p("module_score", "Ag-score", prefix = "p", method = "pearson") %>%
    mutate(n = n_distinct(.cell_id)) %>%
    distinct(!!!syms(keep_clmns), .keep_all = TRUE) %>%
    ungroup() %>%
    select(all_of(c(keep_clmns, keep_columns))) %>%
    mutate(across(c(s_p, p_p), ~ p.adjust(.x, method = "BH"), .names = "{.col}_adj"))
  
  res
}

#' Create boxplots comparing Ag-score for each LEC and DC subset
#' 
#' Used for figure 1
create_ag_boxes <- function(df_in, clrs, facet_vars = "subtype",
                            include_n = TRUE) {
  
  # Line data
  # * make sure tm is numeric
  grp_vars <- unique(c(facet_vars, "tm"))
  
  df_in <- df_in %>%
    mutate(tm = as.numeric(as.character(tm)))
  
  ln_dat <- df_in %>%
    group_by(!!!syms(grp_vars)) %>%
    summarize(Ag_score = median(Ag_score), .groups = "drop")
  
  # Create boxplots
  x_facet <- dplyr::last(facet_vars)
  
  if (length(facet_vars) == 1) facet_vars <- c("", facet_vars)
  
  facet_vars <- str_c(facet_vars, collapse = "~")
  
  res <- df_in %>%
    ggplot(aes(tm, Ag_score, fill = !!sym(x_facet))) +
    geom_line(
      data = ln_dat,
      color = "grey85", linewidth = 0.5, linetype = 2
    ) +
    geom_boxplot(
      aes(group = mouse),
      outlier.size = 0.2,
      width = 3,
      key_glyph = draw_key_point
    ) +
    
    guides(
      fill     = guide_legend(override.aes = list(size = 4, shape = 22), order = 1),
      linetype = guide_legend(override.aes = list(color = "black"), title = NULL)
    ) +
    facet_grid(as.formula(facet_vars)) +
    
    scale_fill_manual(values = clrs) +
    labs(x = "days post immunization", y = "Ag-score") +
    base_theme +
    theme(
      aspect.ratio    = 0.9,
      legend.position = "none",
      plot.margin     = margin(0, 15, 0, 15)
    )
  
  # Add N values to plots
  if (include_n) {
    n_dat <- df_in %>%
      group_by(!!!syms(grp_vars)) %>%
      summarize(n = label_comma()(n()), .groups = "drop")
    
    res <- res +
      geom_text_repel(
        aes(y = Inf, label = n),
        data  = n_dat,
        seed  = 42,
        force = 0.005,
        size  = 8 / .pt,
        vjust = 1.2
      ) +
      scale_y_continuous(expand = expansion(c(0.05, 0.15))) +
      theme(plot.margin = margin(15, 15, 30, 15))
  }
  
  res
}

#' Create bargraphs summarizing the fraction of Ag-competent cells
#' 
#' Used for figure 5
create_chikv_class_bars <- function(df_in, x = "treatment", clrs = treat_clrs,
                                    rep_clmn = "rep", ttl, n_p = NULL,
                                    pred_clmn = "ml_pred_1_grp") {
  
  # Set colors
  grps <- unique(df_in[[x]])
  clrs <- clrs[names(clrs) %in% grps]
  
  # Format data for bargraphs
  df_in <- df_in %>%
    ungroup() %>%
    mutate(
      !!sym(x) := fct_relevel(!!sym(x), names(clrs)),
      x_clmn    = str_c(subtype, "-", !!sym(x))
    ) %>%
    arrange(!!sym(x)) %>%
    mutate(x_clmn = fct_inorder(x_clmn))
  
  dat <- df_in %>%
    group_by(!!!syms(c("x_clmn", x, "subtype", rep_clmn))) %>%
    summarize(
      frac_comp = sum(!!sym(pred_clmn) == "Ag-competent") / n(),
      .groups = "drop"
    )
  
  # Perform fisher exact test
  p_dat <- df_in %>%
    group_by(x_clmn, subtype, !!sym(x)) %>%
    mutate(n = n()) %>%
    group_by(x_clmn, subtype, !!sym(x), n) %>%
    summarize(
      n_comp = sum(!!sym(pred_clmn) == "Ag-competent"),
      .groups = "drop"
    ) %>%
    group_by(subtype) %>%
    mutate(
      frac_comp = n_comp / n,
      tot_comp  = sum(n_comp),
      tot       = sum(n)
    ) %>%
    ungroup()
  
  p_vals <- p_dat %>%
    rowwise() %>%
    mutate(p = .calc_fisher(
      x = n_comp,
      k = n,
      m = tot_comp,
      n = tot - tot_comp
    )) %>%
    ungroup() %>%
    distinct(subtype, p) %>%    # there should be a single p-value per cell type
    mutate(
      p_adj = p.adjust(
        p,
        method = "BH",
        n = n_p %||% length(p)  # specifying n_p is useful when only plotting a 
      )                         # single panel and showing the other panels in
    ) %>%                       # other figures
    rowwise() %>%
    mutate(p_adj = str_c("italic(p) == ", .format_pvalue(p_adj))) %>%
    ungroup()
  
  if (nrow(p_vals) > n_distinct(p_vals$subtype)) {
    cli_abort("Malformed p-value data, check calculations")
  }
  
  # Format x-axis labels
  x_labs <- p_dat %>%
    mutate(
      n = str_c("n = ", label_comma()(n)),
      n = str_c(!!sym(x), "\n", n)
    )
  
  x_labs <- set_names(x_labs$n, x_labs$x_clmn)
  
  # Create bargraphs
  plt_aes <- aes(
    x    = x_clmn,
    y    = frac_comp,
    fill = !!sym(x)
  )
  
  if (!is.null(rep_clmn)) plt_aes$alpha <- sym(rep_clmn)
  
  res <- dat %>%
    ggplot(plt_aes) +
    geom_col(position = position_dodge2(padding = 0.15)) +
    
    geom_text(
      aes(1.5, Inf, fill = NULL, alpha = NULL, label = p_adj),
      data  = p_vals,
      parse = TRUE,
      vjust = 1.3,
      size  = 12 / .pt
    ) +
    
    facet_wrap(~ subtype, scales = "free_x") +
    scale_fill_manual(values = clrs) +
    scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1), expand = expansion(c(0.05, 0.2))) +
    scale_x_discrete(labels = x_labs) +
    guides(alpha = "none") +
    labs(y = ttl) +
    base_theme +
    theme(
      aspect.ratio    = 1.5,
      legend.position = "none",
      strip.clip      = "off",
      axis.title.x    = element_blank(),
      axis.text.x     = element_text(size = ttl_pt2)
    )
  
  if (!is.null(rep_clmn)) {
    al_vals <- rep(1, n_distinct(df_in[[rep_clmn]]))
    
    res <- res +
      scale_alpha_manual(values = al_vals)
  }
  
  res
}

#' Create boxplots comparing Ag module scores for mock vs CHIKV
#' 
#' Used for figure 5
create_chikv_module_boxes <- function(df_in, x = "treatment", module_clmns,
                                      clrs = treat_clrs, ttl, n_p = NULL, ...) {
  
  # Format colors
  grps <- unique(df_in[[x]])
  clrs <- clrs[names(clrs) %in% grps]
  
  # Format data for plotting
  if (length(grps) > 2) {
    cli_warn("More than two groups present in `x`, using the first two")
  }
  
  dat <- df_in %>%
    pivot_longer(all_of(module_clmns)) %>%
    mutate(
      module_subtype = str_remove(name, "_[^_]+$"),
      module_type    = str_remove(name, "^.+_"),
      module_type    = str_c("Ag-", module_type),
      !!sym(x)      := fct_relevel(!!sym(x), names(clrs)),
      x_clmn         = str_c(subtype, "-", !!sym(x))
    ) %>%
    filter(subtype == module_subtype) %>%
    arrange(!!sym(x)) %>%
    mutate(x_clmn = fct_inorder(x_clmn))
  
  # Calculate p-values
  p_dat <- dat %>%
    group_by(subtype, module_type) %>%
    summarize(
      p = wilcox.test(
        value[!!sym(x) == grps[1]],
        value[!!sym(x) == grps[2]]
      )$p.value,
      .groups = "drop"
    ) %>%
    mutate(
      p_adj = p.adjust(
        p,
        method = "BH",
        n = n_p %||% length(p)  # specifying n_p is useful when only plotting a 
      )                         # single panel and showing the other panels in
    ) %>%                       # other figures
    rowwise() %>%
    mutate(p_adj = str_c("italic(p) == ", .format_pvalue(p_adj))) %>%
    ungroup()
  
  # Format x-axis labels
  x_labs <- dat %>%
    group_by(!!sym(x), x_clmn) %>%
    summarize(n = n_distinct(.cell_id), .groups = "drop") %>%
    mutate(
      n = str_c("n = ", label_comma()(n)),
      n = str_c(!!sym(x), "\n", n)
    )
  
  x_labs <- set_names(x_labs$n, x_labs$x_clmn)
  
  # Create boxplots
  res <- dat %>%
    ggplot(aes(x_clmn, value, fill = !!sym(x), alpha = rep)) +
    geom_boxplot(outlier.size = 0.5, ...) +
    
    geom_text(
      aes(1.5, Inf, fill = NULL, alpha = NULL, label = p_adj),
      data  = p_dat,
      parse = TRUE,
      vjust = 1.3,
      size  = 12 / .pt
    ) +
    
    facet_grid(module_type ~ subtype, scales = "free") +
    scale_fill_manual(values = clrs) +
    scale_alpha_manual(values = rep(1, n_distinct(df_in$rep))) +
    scale_y_continuous(expand = expansion(c(0.05, 0.1))) +
    scale_x_discrete(labels = x_labs) +
    guides(alpha = "none") +
    labs(y = ttl) +
    base_theme +
    theme(
      aspect.ratio    = 1.5,
      legend.position = "none",
      axis.title.x    = element_blank(),
      axis.text.x     = element_text(size = ttl_pt2)
    )
  
  res
}

create_gn_plots <- function(so_in, p_data,
                            x = "ml_pred_1_grp", plt_clrs, x_lvls,
                            n_gns = 10, top_gns = NULL,
                            p_test = wilcox.test, draw_line = TRUE, pt_size = 1,
                            sort = TRUE) {
  
  # Set genes to plot
  gns <- p_data %>%
    distinct(gene, class) %>%
    mutate(top = gene %in% top_gns) %>%
    group_by(class) %>%
    mutate(rank = row_number()) %>%
    arrange(desc(top), rank) %>%
    slice(1:n_gns) %>%
    ungroup() %>%
    arrange(rank)
  
  gns <- set_names(gns$class, gns$gene)
  
  p_data <- p_data %>%
    filter(gene %in% names(gns))
  
  # Format data
  n_mice <- n_distinct(so_in$mouse)
  
  cell_type <- unique(so_in$subtype)
  
  if (length(cell_type) > 1) cli_abort("Object contains more than one cell type")
  
  plt_dat <- so_in %>%
    FetchData(c("mouse", "tm", "subtype", x, names(gns))) %>%
    as_tibble(rownames = ".cell_id") %>%
    pivot_longer(all_of(names(gns)), names_to = "gene") %>%
    left_join(p_data, by = c("tm", x, "subtype", "gene")) %>%
    
    mutate(
      gene      = fct_relevel(gene, names(gns)),
      tm        = str_c("day ", tm),
      x         = str_c(!!sym(x), "-", tm),
      !!sym(x) := fct_relevel(!!sym(x), x_lvls),
    ) %>%
    arrange(!!sym(x)) %>%
    mutate(x = fct_inorder(x))
  
  plt_dat <- plt_dat %>%
    group_by(gene, mouse, tm, subtype, !!sym(x), x) %>%
    summarize(
      n       = n_distinct(.cell_id),
      med     = median(value),
      q1      = quantile(value, 0.25),
      q3      = quantile(value, 0.75),
      p_adj   = unique(p_adj),
      .groups = "drop"
    ) %>%
    group_by(gene) %>%
    mutate(range = max(q3) - min(q1)) %>%
    ungroup() %>%
    mutate(
      class = gns[gene],
      p_y = q3 + (range * 0.2)
    )
  
  # Format x-axis labels
  x_lab <- plt_dat %>%
    mutate(
      xlab = scales::label_comma()(n),
      xlab = str_c(!!sym(x), " (n = ", xlab, ")")
    ) %>%
    distinct(xlab, x)
  
  x_lab <- set_names(x_lab$xlab, x_lab$x)
  
  # Create boxplots
  plt_dat <- split(plt_dat, plt_dat$class)
  
  res <- plt_dat %>%
    imap(~ {
      cls <- .y
      
      plt <- .x %>%
        ggplot(aes(x, med, color = !!sym(x))) +
        geom_segment(
          aes(x = x, xend = x, y = q1, yend = q3),
          linewidth = pt_size, color = ln_col
        ) +
        geom_point(aes(size = !!sym(x)))
      
      if (draw_line) {
        plt <- plt +
          geom_smooth(
            aes(x = as.numeric(as.factor(!!sym(x)))),
            method = "lm", linewidth = 0.25, linetype = 2,
            se = FALSE, formula = y ~ x
          )
      }
      
      plt <- plt +
        geom_point(
          aes(y = p_y),
          data = ~ filter(.x, class == cls, p_adj < 0.05),
          color = "black", shape = 6, stroke = 1, size = pt_size
        ) +
        facet_grid(gene ~ tm, scales = "free", switch = "y") +
        scale_size_manual(values = c(pt_size, pt_size * 1.25, pt_size * 1.25)) +
        scale_color_manual(values = plt_clrs) +
        scale_y_continuous(
          breaks = ~ c(ceiling(min(.x)), floor(max(.x))),
          sec.axis = dup_axis(name = str_c(cell_type, " Ag-", .y))
        ) +
        scale_x_discrete(labels = x_lab) +
        base_theme +
        theme(
          aspect.ratio       = 1,
          plot.margin        = margin(15, 15, 15, 15),
          legend.position    = "none",
          panel.border       = element_rect(color = "grey75"),
          strip.placement    = "outside",
          strip.text.y.left  = element_text(face = "italic"),
          axis.ticks.y.right = element_blank(),
          axis.title         = element_blank(),
          axis.title.y.right = element_text(angle = 270, size = ttl_pt2 * 1.2),
          axis.text.y.right  = element_blank(),
          axis.text.x        = element_text(angle = 45, hjust = 1)
        )
      
      # Adjust theme
      if (.y == dplyr::first(names(plt_dat))) {
        plt <- plt +
          theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
          )
        
      } else {
        plt <- plt +
          theme(
            strip.text.x = element_blank(),
            plot.margin  = margin(20, 15, 15, 15)
          )
      }
      
      plt
    })
  
  res
}

create_sig_umap <- function(so_in, dat_col, grp_col = "mouse", clrs,
                            pt_size = 0.15, pt_stroke = 0.6,
                            scale_limits = NULL, outline = TRUE) {
  res <- so_in %>%
    plot_scatter(
      dat_col,
      "hUMAP_1", "hUMAP_2",
      group_col    = grp_col,
      outline      = outline,
      plot_colors  = clrs,
      panel_nrow   = 1,
      size         = pt_size,
      stroke       = pt_stroke,
      label_params = list(size = ttl_pt2)
    ) +
    guides(fill = guide_colorbar(
      ticks = FALSE, title.position = "top",
      barheight = unit(6, "pt"), barwidth = unit(120, "pt"), 
    )) +
    base_theme +
    theme(
      plot.margin     = margin(),
      plot.title      = element_text(size = 18),
      legend.margin   = margin(),
      legend.position = "bottom",
      legend.title    = element_blank(),
      axis.title      = element_blank(),
      axis.ticks      = element_blank(),
      axis.text       = element_blank(),
      aspect.ratio    = 0.9
    )
  
  if (!is.null(scale_limits)) {
    res <- res +
      scale_fill_gradientn(colours = clrs, limits = scale_limits)
  }
  
  res
}

create_grp_umap <- function(so_in, dat_col, grp_col = "mouse", clrs,
                            lvls = c("Ag-competent", "Ag-low", "Ag-high", "other"),
                            ...) {
  
  # Format n labels  
  n_labs <- so_in %>%
    filter(!is.na(!!sym(dat_col))) %>%
    group_by(!!sym(dat_col)) %>%
    summarize(n = n()) %>%
    mutate(
      n_lab = str_c(!!sym(dat_col), "\nn = ", n)
    )
  
  n_labs <- set_names(n_labs$n_lab, n_labs[[dat_col]])
  
  lvls <- lvls[lvls %in% names(n_labs)]
  
  # Create UMAPs
  # * add colors separately since ggplot2 3.5.1 automatically includes `NA` in
  #   plot legend
  res <- so_in %>%
    plot_scatter(
      dat_col,
      "hUMAP_1", "hUMAP_2",
      group_col    = grp_col,
      plot_lvls    = lvls,
      panel_nrow   = 1,
      label_params = list(size = ttl_pt2),
      ...
    ) +
    guides(fill = guide_legend(
      label.position = "bottom", reverse = TRUE,
      override.aes = list(size = 4)
    )) +
    
    scale_color_manual(
      values   = clrs,
      limits   = lvls,
      breaks   = lvls,
      labels   = n_labs,
      na.value = "grey80"
    ) +
    
    base_theme +
    theme(
      plot.margin     = margin(),
      plot.title      = element_text(size = ttl_pt2),
      legend.margin   = margin(),
      legend.position = "bottom",
      legend.title    = element_blank(),
      axis.title      = element_blank(),
      axis.ticks      = element_blank(),
      axis.text       = element_blank(),
      aspect.ratio    = 0.9
    )
  
  res
}

#' Format n values for plot labels
.format_n_label <- function(n) {
  ifelse(
    n >= 1000,
    str_c(round(n / 1000, 1), "k"),
    as.character(n)
  )
}

.create_umap <- function(dat, dat_clmn = "subtype", clrs, add_axis_title = TRUE,
                         ...) {
  
  # Format data for n labels
  n_dat <- dat %>%
    group_by(!!sym(dat_clmn)) %>%
    summarize(
      n       = n_distinct(.cell_id),
      hUMAP_1 = mean(hUMAP_1),
      hUMAP_2 = mean(hUMAP_2),
      .groups = "drop"
    ) %>%
    mutate(
      !!sym(dat_clmn) := fct_reorder(!!sym(dat_clmn), n, .desc = TRUE),
      !!sym(dat_clmn) := fct_relevel(!!sym(dat_clmn), "unassigned", after = Inf)
    ) %>%
    arrange(!!sym(dat_clmn)) %>%
    mutate(
      rank = row_number(!!sym(dat_clmn)),
      rank = if_else(!!sym(dat_clmn) == "unassigned", NA, rank),
      lgnd = if_else(
        !!sym(dat_clmn) != "unassigned",
        str_c(rank, "-", !!sym(dat_clmn)),
        !!sym(dat_clmn)
      ),
      lgnd = str_c(lgnd, "\n   (", label_comma()(n), ")")
    )
  
  # Create UMAP
  plt <- dat %>%
    plot_scatter(
      dat_clmn, "hUMAP_1", "hUMAP_2",
      size         = 0.1,
      label_params = list(size = ttl_pt1),
      panel_nrow   = 2,
      plot_lvls    = levels(n_dat[[dat_clmn]])
    ) +
    geom_text_repel(
      aes(label = rank),
      data     = na.omit(n_dat),
      color    = "black",
      fontface = "bold",
      size     = txt_pt2 / .pt,
      seed     = 42,
      force    = 0.001,
      min.segment.length = Inf,
      ...
    ) +
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 4), reverse = TRUE)) +
    scale_color_manual(
      values = clrs,
      labels = set_names(n_dat$lgnd, n_dat[[dat_clmn]])
    ) +
    scale_y_continuous(expand = expansion(c(0.1, 0.05))) +
    umap_theme +
    theme(
      aspect.ratio    = 0.9,
      legend.position = "right",
      legend.title    = element_blank(),
      legend.text     = element_text(size = txt_pt2)
    )
  
  # Set coords for axis lines and labels
  if (add_axis_title) {
    plt <- plt %>%
      .add_umap_labels(
        arrows     = TRUE,
        arrow_frac = 0.3
      )
  }
  
  plt
}

.add_umap_labels <- function(plt, labels = NULL, include = c("x", "y"),
                             hjust = 0.5, pad = margin(0, 5, 0, 5, "pt"),
                             arrows = FALSE, arrow_frac = 0.25) {
  
  x <- as.character(plt$mapping$x)[2]
  y <- as.character(plt$mapping$y)[2]
  
  if (is.null(labels)) {
    labels <- c(x, y) %>%
      map_chr(~ {
        .x %>%
          str_remove_all(".+(?=UMAP)") %>%
          str_replace_all("[_-]", " ")
      })
  }
  
  if (!arrows) {
    plt <- plt +
      scale_y_continuous(expand = expansion(c(0.15, 0.05))) +
      labs(x = labels[1], y = labels[2])
    
    if ("x" %in% include) {
      plt <- plt +
        theme(axis.title.x = element_text(hjust = hjust)) +
        theme(
          axis.title.x = element_textbox(
            fill    = "white",
            margin  = margin(t = -8),
            padding = pad
          )
        )
    }
    
    if ("y" %in% include) {
      plt <- plt +
        theme(axis.title.y = element_text(hjust = hjust)) +
        theme(
          axis.title.y = element_textbox(
            orientation = "left-rotated",
            fill    = "white",
            margin  = margin(b = -10.8),
            padding = pad
          )
        )
    }
    
    return(plt)
  }
  
  # Add titles with arrows
  dat   <- plt$data
  rat   <- plt$theme$aspect.ratio %||% 1
  x_min <- min(dat[[x]]) - (diff(range(dat[[x]])) * 0.05)
  y_min <- min(dat[[y]]) - (diff(range(dat[[y]])) * 0.05)
  x_q1  <- min(dat[[x]]) + (diff(range(dat[[x]])) * arrow_frac * rat)
  y_q1  <- min(dat[[y]]) + (diff(range(dat[[y]])) * arrow_frac)
  
  ln_dat <- tibble(
    x     = rep(x_min, 2),
    y     = rep(y_min, 2),
    xend  = c(x_q1, x_min),
    yend  = c(y_min, y_q1),
    angle = c(0, 90),
    vjust = c(1.2, -0.2),
    label = labels
  )
  
  plt <- plt +
    geom_segment(
      aes(x = x, xend = xend, y = y, yend = yend),
      data = ln_dat,
      color = "black",
      linewidth = ln_pt,
      arrow = arrow(length = unit(5, "pt"))
    ) +
    geom_text(
      aes(x = xend, y = yend, label = label, vjust = vjust, angle = angle),
      data = ln_dat,
      color = "black",
      hjust = 1.2
    )
  
  plt
}

#' Plot image as ggplot2 object
#'
#' @param img Raster array, e.g. output from readTIFF
#' @param t,r,b,l Plot margins  
.plot_img <- function(img, t = 0, r = 0, b = 0, l = 0) {
  
  img <- rasterGrob(img, interpolate = TRUE)
  
  # Adjust plot limits if necessary
  img_lims <- tibble(x = c(0, 1), y = c(0, 1))
  
  # Create an empty plot for placing the image, adjusted to fill the area
  res <- ggplot(img_lims, aes(x, y)) + 
    annotation_custom(
      grob = img,
      xmin = -Inf, xmax = Inf,
      ymin = -Inf, ymax = Inf
    ) +
    xlim(c(-Inf, Inf)) +
    ylim(c(-Inf, Inf)) +
    theme_void() +
    theme(plot.margin = margin(t, r, b, l))
  
  res
}
