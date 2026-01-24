#' @title Data Input/Output Functions for STRIDER
#' @description Functions for reading and writing spatial transcriptomics data
#' @name data_io
NULL

#' Read 10X Genomics HDF5 File
#'
#' @description Read gene expression data from 10X Genomics HDF5 format.
#'   Supports both Cell Ranger v2 and v3 formats.
#'
#' @param filename Path to the HDF5 file
#' @param use.names Use gene names instead of IDs (default: TRUE)
#' @param unique.features Make feature names unique (default: TRUE)
#'
#' @return A list containing:
#'   \itemize{
#'     \item{matrix}{Sparse dgCMatrix of counts}
#'     \item{features}{Gene names/IDs}
#'     \item{barcodes}{Cell/spot barcodes}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data <- read_10x_h5("filtered_feature_bc_matrix.h5")
#' }
read_10x_h5 <- function(filename, use.names = TRUE, unique.features = TRUE) {
  # Check file exists
  if (!file.exists(filename)) {
    rlang::abort(paste0("File not found: ", filename))
  }
  
  # Try hdf5r first, then rhdf5
  if (requireNamespace("hdf5r", quietly = TRUE)) {
    return(.read_10x_h5_hdf5r(filename, use.names, unique.features))
  } else if (requireNamespace("rhdf5", quietly = TRUE)) {
    return(.read_10x_h5_rhdf5(filename, use.names, unique.features))
  } else {
    rlang::abort(
      "Either 'hdf5r' or 'rhdf5' package is required to read HDF5 files.\n",
      "Install with: install.packages('hdf5r') or BiocManager::install('rhdf5')"
    )
  }
}

#' @keywords internal
.read_10x_h5_hdf5r <- function(filename, use.names, unique.features) {
  h5file <- hdf5r::H5File$new(filename, mode = "r")
  on.exit(h5file$close_all(), add = TRUE)
  
  # Check for matrix group
  if (h5file$exists("matrix")) {
    # Cell Ranger v3 format
    group <- h5file[["matrix"]]
    
    # Read matrix data
    data <- group[["data"]][]
    indices <- group[["indices"]][]
    indptr <- group[["indptr"]][]
    shape <- group[["shape"]][]
    
    # Read features
    if (group$exists("features")) {
      feature_group <- group[["features"]]
      if (use.names && feature_group$exists("name")) {
        features <- feature_group[["name"]][]
      } else {
        features <- feature_group[["id"]][]
      }
    } else {
      features <- group[["gene_names"]][]
    }
    
    # Read barcodes
    barcodes <- group[["barcodes"]][]
    
  } else {
    rlang::abort("Unrecognized HDF5 format. Expected 'matrix' group.")
  }
  
  # Convert bytes to character if needed
  if (is.raw(features[1])) {
    features <- rawToChar(features)
  }
  if (inherits(features, "H5R")) {
    features <- as.character(features)
  }
  features <- as.character(features)
  barcodes <- as.character(barcodes)
  
  # Create sparse matrix (0-indexed to 1-indexed)
  mat <- Matrix::sparseMatrix(
    i = indices + 1L,
    p = indptr,
    x = as.numeric(data),
    dims = rev(shape),
    dimnames = list(features, barcodes)
  )
  
  # Make features unique if requested
  if (unique.features) {
    rownames(mat) <- make.unique(rownames(mat))
    features <- rownames(mat)
  }
  
  list(
    matrix = mat,
    features = features,
    barcodes = barcodes
  )
}

#' @keywords internal
.read_10x_h5_rhdf5 <- function(filename, use.names, unique.features) {
  # Read using rhdf5
  data <- rhdf5::h5read(filename, "matrix/data")
  indices <- rhdf5::h5read(filename, "matrix/indices")
  indptr <- rhdf5::h5read(filename, "matrix/indptr")
  shape <- rhdf5::h5read(filename, "matrix/shape")
  barcodes <- rhdf5::h5read(filename, "matrix/barcodes")
  
  # Try to read features
  tryCatch({
    if (use.names) {
      features <- rhdf5::h5read(filename, "matrix/features/name")
    } else {
      features <- rhdf5::h5read(filename, "matrix/features/id")
    }
  }, error = function(e) {
    features <<- rhdf5::h5read(filename, "matrix/gene_names")
  })
  
  features <- as.character(features)
  barcodes <- as.character(barcodes)
  
  # Create sparse matrix
  mat <- Matrix::sparseMatrix(
    i = indices + 1L,
    p = indptr,
    x = as.numeric(data),
    dims = rev(shape),
    dimnames = list(features, barcodes)
  )
  
  if (unique.features) {
    rownames(mat) <- make.unique(rownames(mat))
    features <- rownames(mat)
  }
  
  list(
    matrix = mat,
    features = features,
    barcodes = barcodes
  )
}

#' Write 10X Genomics HDF5 File
#'
#' @description Write gene expression data to 10X Genomics HDF5 format.
#'
#' @param matrix Sparse matrix (genes x cells/spots)
#' @param features Gene names
#' @param barcodes Cell/spot barcodes
#' @param filename Output file path
#' @param overwrite Overwrite existing file (default: FALSE)
#'
#' @return Invisible NULL
#' @export
#'
#' @examples
#' \dontrun{
#' write_10x_h5(mat, genes, cells, "output.h5")
#' }
write_10x_h5 <- function(matrix, features, barcodes, filename, overwrite = FALSE) {
  if (file.exists(filename) && !overwrite) {
    rlang::abort(paste0("File already exists: ", filename, 
                        ". Use overwrite = TRUE to replace."))
  }
  
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    rlang::abort("Package 'hdf5r' is required. Install with: install.packages('hdf5r')")
  }
  
  # Convert to dgCMatrix
  mat <- methods::as(matrix, "dgCMatrix")
  
  # Create HDF5 file
  h5file <- hdf5r::H5File$new(filename, mode = "w")
  on.exit(h5file$close_all(), add = TRUE)
  
  # Create matrix group
  mat_group <- h5file$create_group("matrix")
  
  # Write matrix data
  mat_group[["data"]] <- mat@x
  mat_group[["indices"]] <- mat@i
  mat_group[["indptr"]] <- mat@p
  mat_group[["shape"]] <- rev(dim(mat))
  mat_group[["barcodes"]] <- barcodes
  
  # Create features group
  feat_group <- mat_group$create_group("features")
  feat_group[["id"]] <- features
  feat_group[["name"]] <- features
  feat_group[["feature_type"]] <- rep("Gene Expression", length(features))
  
  invisible(NULL)
}

#' Read Count Matrix from Text File
#'
#' @description Read gene expression count matrix from tab-separated text file.
#'   First column should contain gene names, first row should contain cell/spot IDs.
#'
#' @param filename Path to the text file
#' @param sep Field separator (default: tab)
#' @param header Whether file has header (default: TRUE)
#' @param row.names Column containing row names (default: 1)
#' @param sparse Convert to sparse matrix (default: TRUE)
#'
#' @return A list containing:
#'   \itemize{
#'     \item{matrix}{Count matrix (sparse or dense)}
#'     \item{features}{Gene names}
#'     \item{barcodes}{Cell/spot barcodes}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data <- read_count_matrix("gene_count.txt")
#' }
read_count_matrix <- function(filename, sep = "\t", header = TRUE, 
                               row.names = 1, sparse = TRUE) {
  if (!file.exists(filename)) {
    rlang::abort(paste0("File not found: ", filename))
  }
  
  vmessage("Reading count matrix: ", filename)
  
  # Use data.table for fast reading
  dt <- data.table::fread(
    filename, 
    sep = sep, 
    header = header,
    data.table = TRUE,
    showProgress = FALSE
  )
  
  # Extract row names (genes)
  if (is.numeric(row.names)) {
    features <- as.character(dt[[row.names]])
    dt <- dt[, -row.names, with = FALSE]
  } else if (is.character(row.names)) {
    features <- as.character(dt[[row.names]])
    dt[, (row.names) := NULL]
  } else {
    features <- paste0("Gene_", seq_len(nrow(dt)))
  }
  
  # Get column names (barcodes)
  barcodes <- colnames(dt)
  
  # Convert to matrix
  mat <- as.matrix(dt)
  rownames(mat) <- features
  colnames(mat) <- barcodes
  storage.mode(mat) <- "numeric"
  
  # Convert to sparse if requested
  if (sparse) {
    mat <- as_sparse(mat)
  }
  
  vmessage("  Loaded ", nrow(mat), " genes x ", ncol(mat), " cells/spots")
  
  list(
    matrix = mat,
    features = features,
    barcodes = barcodes
  )
}

#' Read Spatial Coordinates
#'
#' @description Read spatial coordinates from tab-separated text file.
#'   First column should contain spot IDs, followed by coordinates.
#'
#' @param filename Path to the text file
#' @param sep Field separator (default: tab)
#' @param header Whether file has header (default: TRUE)
#' @param spot.col Column containing spot IDs (default: 1)
#' @param coord.cols Columns containing coordinates (default: c(2, 3))
#' @param sample.col Column containing sample IDs (optional)
#'
#' @return A data.frame with columns: spot, x, y (and optionally sample)
#' @export
#'
#' @examples
#' \dontrun{
#' coords <- read_spatial_coords("location.txt")
#' }
read_spatial_coords <- function(filename, sep = "\t", header = TRUE,
                                 spot.col = 1, coord.cols = c(2, 3),
                                 sample.col = NULL) {
  if (!file.exists(filename)) {
    rlang::abort(paste0("File not found: ", filename))
  }
  
  vmessage("Reading spatial coordinates: ", filename)
  
  df <- data.table::fread(
    filename,
    sep = sep,
    header = header,
    data.table = FALSE
  )
  
  # Extract spot IDs
  spots <- as.character(df[[spot.col]])
  
  # Extract coordinates
  x <- as.numeric(df[[coord.cols[1]]])
  y <- as.numeric(df[[coord.cols[2]]])
  
  result <- data.frame(
    spot = spots,
    x = x,
    y = y,
    row.names = spots,
    stringsAsFactors = FALSE
  )
  
  # Add sample column if specified
  if (!is.null(sample.col)) {
    result$sample <- as.character(df[[sample.col]])
  }
  
  vmessage("  Loaded coordinates for ", nrow(result), " spots")
  
  result
}

#' Read Cell Type Annotations
#'
#' @description Read cell type annotations from tab-separated text file.
#'   File should have two columns: cell ID and cell type.
#'
#' @param filename Path to the annotation file
#' @param sep Field separator (default: tab)
#' @param header Whether file has header (default: FALSE)
#'
#' @return A named vector mapping cell IDs to cell types
#' @export
#'
#' @examples
#' \dontrun{
#' celltypes <- read_celltype_annotations("celltype.txt")
#' }
read_celltype_annotations <- function(filename, sep = "\t", header = FALSE) {
  if (!file.exists(filename)) {
    rlang::abort(paste0("File not found: ", filename))
  }
  
  vmessage("Reading cell type annotations: ", filename)
  
  df <- data.table::fread(
    filename,
    sep = sep,
    header = header,
    data.table = FALSE
  )
  
  if (ncol(df) < 2) {
    rlang::abort("Annotation file must have at least 2 columns: cell ID and cell type")
  }
  
  celltypes <- stats::setNames(as.character(df[[2]]), as.character(df[[1]]))
  
  vmessage("  Loaded annotations for ", length(celltypes), " cells")
  vmessage("  Cell types: ", paste(unique(celltypes), collapse = ", "))
  
  celltypes
}

#' Read Gene List
#'
#' @description Read a list of genes from a text file (one gene per line).
#'
#' @param filename Path to the gene list file
#'
#' @return A character vector of gene names
#' @export
#'
#' @examples
#' \dontrun{
#' genes <- read_gene_list("markers.txt")
#' }
read_gene_list <- function(filename) {
  if (!file.exists(filename)) {
    rlang::abort(paste0("File not found: ", filename))
  }
  
  genes <- readLines(filename, warn = FALSE)
  genes <- trimws(genes)
  genes <- genes[nchar(genes) > 0]  # Remove empty lines
  
  vmessage("  Loaded ", length(genes), " genes")
  
  genes
}

#' Write Results to File
#'
#' @description Write STRIDER results to tab-separated text files.
#'
#' @param x Data to write (matrix or data.frame)
#' @param filename Output file path
#' @param row.names Include row names (default: TRUE)
#' @param col.names Include column names (default: TRUE)
#'
#' @return Invisible NULL
#' @export
write_results <- function(x, filename, row.names = TRUE, col.names = TRUE) {
  vmessage("Writing results to: ", filename)
  
  # Create directory if needed
  dir_path <- dirname(filename)
  if (!dir.exists(dir_path) && dir_path != ".") {
    dir.create(dir_path, recursive = TRUE)
  }
  
  # Convert sparse to dense if needed
  if (is_sparse(x)) {
    x <- as.matrix(x)
  }
  
  data.table::fwrite(
    data.table::as.data.table(x, keep.rownames = row.names),
    file = filename,
    sep = "\t",
    col.names = col.names,
    row.names = FALSE
  )
  
  invisible(NULL)
}
