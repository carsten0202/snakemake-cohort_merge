log_connection <- file(snakemake@log[[1]], open = "wt")
sink(log_connection)
sink(log_connection, type = "message")

manifest <- read.delim(
    snakemake@input[["manifest"]], colClasses = "character", check.names = FALSE,
    na.strings = character()
)
sample_ids <- manifest$merged_sample_id
environment <- new.env(parent = baseenv())
loaded <- load(snakemake@input[["pcrelate_data"]], envir = environment)
if (!("mypcrelate" %in% loaded)) {
    stop("PC-Relate data does not contain an object named mypcrelate")
}
mypcrelate <- get("mypcrelate", envir = environment, inherits = FALSE)
if (is.null(mypcrelate$kinBtwn)) {
    stop("mypcrelate does not contain kinBtwn pair estimates")
}
pairs <- as.data.frame(mypcrelate$kinBtwn, stringsAsFactors = FALSE)
if (!all(c("ID1", "ID2", "kin") %in% names(pairs))) {
    stop("PC-Relate kinBtwn is missing required columns ID1, ID2, kin")
}
first <- match(as.character(pairs$ID1), sample_ids)
second <- match(as.character(pairs$ID2), sample_ids)
if (anyNA(first) || anyNA(second)) {
    stop("PC-Relate contains sample IDs absent from the sample manifest")
}
if (any(first == second)) {
    stop("PC-Relate kinBtwn contains self-pairs")
}
keys <- paste(pmin(first, second), pmax(first, second), sep = ":")
if (anyDuplicated(keys)) {
    stop("PC-Relate kinBtwn contains duplicate sample pairs")
}
kinship <- suppressWarnings(as.numeric(as.character(pairs$kin)))
selected <- is.finite(kinship) & kinship >= snakemake@params[["minimum_kinship"]]
report <- data.frame(
    sample_1 = sample_ids[pmin(first, second)[selected]],
    sample_2 = sample_ids[pmax(first, second)[selected]],
    pcrelate_kinship = kinship[selected],
    stringsAsFactors = FALSE
)
report <- report[order(report$sample_1, report$sample_2), ]
write.table(
    report, file = snakemake@output[["pairs"]], sep = "\t",
    quote = FALSE, row.names = FALSE
)
cat(nrow(report), "qualifying duplicate/twin pairs\n")
sink(type = "message")
sink()
close(log_connection)
