#' Deciphers SNV Signatures
#' 
#' This function loads an SNV data file and runs the mutation signature deciphering
#' process on it. It returns the output of nnls_exposures.
#' @param vcf_file Path to an VCF file, suitable for import_vcfs()
#' @param snv_file Path to an SNV file, suitable for import_snvs()
#' @param genome ID of the genome being used (default: 'hg19')
#' @param fractions If TRUE, mutation signature exposures are normalized to sum to 1 (default: 'FALSE')
#' @param silent If TRUE, does not print exposures output (default: 'FALSE')
#' @export

run_snv <- function(vcf_file=NULL, snv_file=NULL, genome='hg19', fractions=FALSE, iterations = 1000, silent=FALSE) {
    if (is.null(vcf_file) && is.null(snv_file)) {
        stop('run_snv requires at least one of vcf_file or snv_file')
    } else if (! is.null(vcf_file) && !is.null(snv_file)) {
        warning('Both vcf_file and snv_file provided. snv_file will be ignored.')
        vr = import_vcf(vcf_file, genome)
    } else if (is.null(snv_file)) {
        vr = import_vcf(vcf_file, genome)
    } else if (is.null(vcf_file)) {
        vr = import_snvs(snv_file, genome)
    }

    vr = sanitize_vranges(vr)

    catalog <- mutation_catalog(vr)
    ref <- get_reference_signatures()
    exposures <- nnls_exposures(catalog, ref, fraction=fractions, montecarlo=TRUE, iterations = iterations)

    if (!silent) {
        print(exposures)
    }

}
