#!/usr/bin/env Rscript

' Usage: test_multiple.R -l LOHFILES -o OUTPUT [ -p PREVIOUS -r GENOME -n NCORES -s SUCCESS -f FAILURE --log LOGPATH --debug ]

Options:
    -l --loh LOHFILES       Path to file containing the paths (one per line) to loh segs files from APOLLOH
    -o --output OUTPUT      Output TSV file destination
    -p --previous PREVIOUS  Path to previously computed results to update
    -r --genome GENOME      Genome version (default: hg19)
    -n --ncores NCORES      Number of cores to run with
    -s --success SUCCESS    Path to file which will be created if script executes successfully
    -f --failure FAILURE    Path to file which will be created if script does NOT execute successfully
    --log LOGPATH           Path to log file
    --debug                 Runs in debug mode
' -> doc

library(docopt)
library(futile.logger)
args <- docopt(doc)

loh_files <- readLines(args[['loh']])
loh_files <- loh_files[file.exists(loh_files)]

if (is.null(args[['genome']])) {
    args[['genome']] = 'hg19'
}

library(snow)

if (args[['debug']]) {
    library(devtools)
    load_all('scripts/hrdtools')
    n <- readLines('scripts/hrdtools/NAMESPACE')
    for (p in unique(gsub('.*?\\((.*?)[,)].*', '\\1', n[grepl('^import', n)]))) library(p, character.only = TRUE)
} else {
    suppressMessages(library(hrdtools))
}

tryCatch({
    test_multiple(loh_files, args[['output']], args[['genome']],
              number_of_cores = args[['ncores']], previous_output = args[['previous']],
              log_path = args[['log']])
    writeLines(c(sprintf('Successfully computed HRD scores at %s', args[['output']])), args[['success']])
}, error = function(e) {
    writeLines(c(sprintf('Failed to calculate new HRD scores for %s. Error: %s', args[['output']], e)), args[['failure']])
})
