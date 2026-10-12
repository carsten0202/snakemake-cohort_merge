king_file <- snakemake@input[["king_matrix"]]
pcrelate_file <- snakemake@input[["pcrelate_data"]]
output_file <- snakemake@output[["report"]]
minimum_kinship <- snakemake@params[["minimum_kinship"]]

king <- as.matrix(
    read.table(
        king_file,
        header = TRUE,
        row.names = 1,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
)

if (nrow(king) != ncol(king)) {
    stop("KING kinship matrix is not square")
}
if (anyDuplicated(rownames(king)) || anyDuplicated(colnames(king))) {
    stop("KING kinship matrix contains duplicate sample IDs")
}
if (!setequal(rownames(king), colnames(king))) {
    stop("KING kinship matrix row and column sample IDs differ")
}

king <- king[, rownames(king), drop = FALSE]
suppressWarnings(storage.mode(king) <- "double")
if (anyNA(king)) {
    stop("KING kinship matrix contains missing or non-numeric values")
}
if (!isTRUE(all.equal(king, t(king), tolerance = 1e-8, check.attributes = FALSE))) {
    stop("KING kinship matrix is not symmetric")
}

pcrelate_environment <- new.env(parent = baseenv())
loaded_objects <- load(pcrelate_file, envir = pcrelate_environment)
if (!("mypcrelate" %in% loaded_objects)) {
    stop("PC-Relate data does not contain an object named mypcrelate")
}

mypcrelate <- get("mypcrelate", envir = pcrelate_environment, inherits = FALSE)
if (is.null(mypcrelate$kinBtwn)) {
    stop("mypcrelate does not contain kinBtwn pair estimates")
}

pcrelate <- as.data.frame(mypcrelate$kinBtwn, stringsAsFactors = FALSE)
required_columns <- c("ID1", "ID2", "kin")
missing_columns <- setdiff(required_columns, names(pcrelate))
if (length(missing_columns) > 0) {
    stop(
        paste0(
            "PC-Relate kinBtwn is missing required columns: ",
            paste(missing_columns, collapse = ", ")
        )
    )
}

pcrelate$ID1 <- as.character(pcrelate$ID1)
pcrelate$ID2 <- as.character(pcrelate$ID2)
suppressWarnings(pcrelate$kin <- as.numeric(pcrelate$kin))

sample_ids <- rownames(king)
first_index <- match(pcrelate$ID1, sample_ids)
second_index <- match(pcrelate$ID2, sample_ids)
if (anyNA(first_index) || anyNA(second_index)) {
    unknown_ids <- unique(
        c(pcrelate$ID1[is.na(first_index)], pcrelate$ID2[is.na(second_index)])
    )
    stop(
        paste0(
            "PC-Relate contains sample IDs absent from the KING matrix: ",
            paste(unknown_ids, collapse = ", ")
        )
    )
}
if (any(first_index == second_index)) {
    stop("PC-Relate kinBtwn contains self-pairs")
}

pair_first <- pmin(first_index, second_index)
pair_second <- pmax(first_index, second_index)
pcrelate_keys <- paste(pair_first, pair_second, sep = ":")
if (anyDuplicated(pcrelate_keys)) {
    stop("PC-Relate kinBtwn contains duplicate sample pairs")
}

king_indices <- which(
    upper.tri(king) & king >= minimum_kinship,
    arr.ind = TRUE
)
king_keys <- if (nrow(king_indices) > 0) {
    paste(king_indices[, "row"], king_indices[, "col"], sep = ":")
} else {
    character()
}
pcrelate_keys_above_threshold <- pcrelate_keys[
    !is.na(pcrelate$kin) & pcrelate$kin >= minimum_kinship
]
reported_keys <- unique(c(king_keys, pcrelate_keys_above_threshold))

relationship_class <- function(kinship) {
    ifelse(
        is.na(kinship),
        NA_character_,
        ifelse(
            kinship >= 0.354,
            "duplicate_or_mz_twin",
            ifelse(
                kinship >= 0.177,
                "first_degree",
                ifelse(
                    kinship >= 0.0884,
                    "second_degree",
                    ifelse(
                        kinship >= 0.0442,
                        "third_degree",
                        "distant_or_unrelated"
                    )
                )
            )
        )
    )
}

if (length(reported_keys) == 0) {
    report <- data.frame(
        sample_1 = character(),
        sample_2 = character(),
        king_kinship = character(),
        king_relationship = character(),
        pcrelate_kinship = character(),
        pcrelate_relationship = character()
    )
} else {
    pair_indices <- do.call(
        rbind,
        strsplit(reported_keys, ":", fixed = TRUE)
    )
    storage.mode(pair_indices) <- "integer"

    king_kinship <- king[pair_indices]
    pcrelate_match <- match(reported_keys, pcrelate_keys)
    pcrelate_kinship <- pcrelate$kin[pcrelate_match]
    strength <- pmax(king_kinship, pcrelate_kinship, na.rm = TRUE)
    ordering <- order(
        -strength,
        pair_indices[, 1],
        pair_indices[, 2]
    )

    format_kinship <- function(values) {
        ifelse(is.na(values), NA_character_, sprintf("%.6g", values))
    }

    report <- data.frame(
        sample_1 = sample_ids[pair_indices[, 1]],
        sample_2 = sample_ids[pair_indices[, 2]],
        king_kinship = format_kinship(king_kinship),
        king_relationship = relationship_class(king_kinship),
        pcrelate_kinship = format_kinship(pcrelate_kinship),
        pcrelate_relationship = relationship_class(pcrelate_kinship),
        stringsAsFactors = FALSE
    )[ordering, ]
}

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write.table(
    report,
    file = output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)
