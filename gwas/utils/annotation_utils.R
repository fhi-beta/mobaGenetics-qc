##
#
# This script contains functions for the query of Phenoscanner and Ensembl. 
# 
##


# Ensembl API

server <- "https://rest.ensembl.org"
server37 <- "https://grch37.rest.ensembl.org"


# List of variants that make VEP crash

vepBlackList <- c("rs138993410")


# Default cache folders

tempFolder <- "~/.gwas_pipeline"
phenoScannerCacheFolder <- file.path(tempFolder, "phenoscanner_cache")
ensemblCacheFolder <- file.path(tempFolder, "ensembl_cache")

if (!dir.exists(tempFolder)) {
  dir.create(tempFolder)
}
if (!dir.exists(phenoScannerCacheFolder)) {
  dir.create(phenoScannerCacheFolder)
}
if (!dir.exists(ensemblCacheFolder)) {
  dir.create(ensemblCacheFolder)
}


#' Fixes custom ids that crash the API.
#' 
#' @param variantId the id of the variant
#' 
#' @return the cleaned id of the variant
cleanRsId <- function(
    variantIds
) {
    
    # check for rs ids with a suffix
    
    suffixI <- !is.na(variantIds) & is.character(variantIds) & startsWith(variantIds, "rs") & str_detect(variantIds, '_')
    index <- ifelse(!is.na(variantIds) & is.character(variantIds), str_locate(variantIds, pattern = "_")[, 1], 0)
    
    cleanIds <- ifelse(!is.na(variantIds), variantIds, "")
    cleanIds <- ifelse(is.character(variantIds), variantIds, as.character(variantIds))
    cleanIds <- ifelse(suffixI, str_sub(variantIds, end = index - 1), variantIds)
    
    # Remove special characters 
    
    cleanIds <- str_replace_all(cleanIds, "/", "_")
    cleanIds <- str_replace_all(cleanIds, ":", "_")
    
    return(cleanIds)
    
}



#' Writes generic files from a phenoscanner query.
#' 
#' @param variantId the id of the variant to query
#' @param altAllele the alternative allele
#' @param phenoscannerFolder the folder where to save the phenoscanner results
#' 
#' @return returns a label for the table or "No Results"
getPhenoscannerDocs <- function(
    variantId,
    altAllele = NULL,
    phenoscannerFolder
) {
  
  if (!is.na(variantId) & is.character(variantId)) {
    
    variantId <- cleanRsId(variantId)
    
    resultsDF <- NULL
    
    noResultsTextFile <- file.path(phenoScannerCacheFolder, "noResults.txt")
    phenoscannerCacheFile <- file.path(phenoScannerCacheFolder, paste0(variantId, ".gz"))
    phenoscannerTextFile <- file.path(phenoscannerFolder, paste0(variantId, ".gz"))
    
    if (file.exists(phenoscannerCacheFile) && file.exists(phenoscannerTextFile)) {
        
        resultsDF <- read.table(
            file = gzfile(phenoscannerCacheFile),
            header = T,
            stringsAsFactors = F,
            sep = "\t",
            comment.char = ""
        )
        write.table(
            x = resultsDF,
            file = gzfile(phenoscannerTextFile),
            col.names = T,
            row.names = F,
            quote = T,
            sep = "\t"
        )
        
    } else {
        
        if (file.exists(noResultsTextFile)) {
            
            noResultsVariants <- readLines(noResultsTextFile)
            
        } else {
            
            noResultsVariants <- c()
            
        }
        
        if (! variantId %in% noResultsVariants) {
            
            phenoScannerResult <- phenoscanner(
                snpquery = variantId,
                proxies = "EUR",
                pvalue = 5e-8
            )
            
            resultsDF <- phenoScannerResult$results
            
            if (nrow(resultsDF > 0)) {
                
                write.table(
                    x = resultsDF,
                    file = gzfile(phenoscannerCacheFile),
                    col.names = T,
                    row.names = F,
                    quote = T,
                    sep = "\t"
                )
                write.table(
                    x = resultsDF,
                    file = gzfile(phenoscannerTextFile),
                    col.names = T,
                    row.names = F,
                    quote = T,
                    sep = "\t"
                )
            } else {
                
                write(
                    x = variantId,
                    file = noResultsTextFile,
                    append = T
                )
                
            }
        }
    }
    
    if (!is.null(resultsDF) && nrow(resultsDF) > 0) {
        
        if (!is.null(altAllele)) {
            
            resultsDF$beta <- as.numeric(resultsDF$beta)
            resultsDF$beta <- ifelse(resultsDF$ref_a2 == altAllele, resultsDF$beta, -resultsDF$beta)
            
        }
        
        phenoScannerFile <- file.path(phenoscannerFolder, paste0(variantId, ".md"))
        
        write(x = paste0("# Phenoscanner - ", variantId, "\n"), file = phenoScannerFile, append = F)
        
        write(x = paste0("Phenoscanner executed on ", format(Sys.time(), "%d.%m.%Y"), " querying ", variantId, "\n"), file = phenoScannerFile, append = T)
        
        write(x = paste0("| SNP | proxy | r² (EUR) | Trait | Beta | se | p |"), file = phenoScannerFile, append = T)
        write(x = paste0("| --- | ----- | -------- | ----- | ---- | -- | - |"), file = phenoScannerFile, append = T)
        
        for (i in 1:nrow(resultsDF)) {
            
            write(x = paste0("| ", resultsDF$snp[i], " | ", resultsDF$rsid[i], " | ", resultsDF$r2[i], " | ", resultsDF$trait[i], " | ", resultsDF$beta[i], " | ", resultsDF$se[i], " | ", resultsDF$p[i], " |"), file = phenoScannerFile, append = T)
            
        }
        
        write(x = "\n", file = phenoScannerFile, append = T)
        write(x = paste0("[[Download]](", paste0(variantId, ".gz"), ")\n"), file = phenoScannerFile, append = T)
        
        uniqueTraits <- unique(resultsDF$trait)
        
        if (length(uniqueTraits) <= 3) {
            
            traitsLabel <- paste(sort(uniqueTraits), collapse = ", ")
            
        } else {
            
            traitsLabel <- "[...]"
            
        }
        
        return(traitsLabel)
        
    } else {
        
        return("No Results")
        
    }
  } else {
    
    print(paste0("Variant ID not supported: ", variantId))
    
  }
}

#' Writes generic files from a vep query.
#' 
#' @param variantId the id of the variant to query
#' @param chr the chromosome of the variant to query
#' @param bp the pb of the variant to query
#' @param ensemblFolder the folder where to save the VEP results
#' @param maxDistance the distance to use around variants of interest for gene annotation
#' 
#' @return returns a label for the table or "No Results"
getEnsemblDocs <- function(
    variantId,
    chr,
    bp,
    ensemblFolder,
    maxDistance = 100000
) {
    
  if (!is.na(variantId) & is.character(variantId)) {
  
    variantId <- cleanRsId(variantId)
    
    noResultsVepFile <- file.path(ensemblCacheFolder, "noResults_vep.txt")
    vepCacheFile <- file.path(ensemblCacheFolder, paste0(variantId, "_vep.json.gz"))
    vepResultsFile <- file.path(ensemblFolder, paste0(variantId, "_vep.json.gz"))
    
    vepJson <- NULL
    
    if (file.exists(vepCacheFile)) {
        
        vepJson <- fromJSON(file = gzfile(vepCacheFile))
        jsonString <- toJSON(vepJson)
        write(
            x = jsonString, 
            file = gzfile(vepResultsFile)
        )
        
    } else if (startsWith(x= variantId, prefix = "rs")) {
        
        if (file.exists(noResultsVepFile)) {
            
            noResultsVariants <- readLines(noResultsVepFile)
            
        } else {
            
            noResultsVariants <- c()
            
        }
        
        if (! variantId %in% noResultsVariants && !variantId %in% vepBlackList) {
            
            ext <- paste0("/vep/human/id/", variantId, "?")
            
            for (try in 1:100) {
                
                tryCatch(
                    {
                        r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
                        
                        responseCode <- status_code(r)
                        
                        if (responseCode == 200) {
                            
                            vepJson <- fromJSON(toJSON(content(r)))
                            
                            break()
                            
                        } else {
                            
                            headers <- headers(r)
                            
                            remaining <- headers[["x-ratelimit-remaining"]]
                            reset <- headers[["x-ratelimit-reset"]]
                            
                            if (!is.null(remaining) && !is.null(reset) && remaining <= 1) {
                                
                                Sys.sleep(reset + 1)
                                
                            } else {
                                
                                Sys.sleep(try)
                                
                            }
                        }
                    },
                    error = function(e) {
                        
                        print(e)
                        vepJson <- NULL
                        
                    }
                )
            }
            
            if (is.null(vepJson)) {
                
                return("Error")
                
            }
            
            if (length(vepJson) > 0) {
                
                jsonString <- toJSON(vepJson)
                write(
                    x = jsonString, 
                    file = gzfile(vepCacheFile)
                )
                write(
                    x = jsonString, 
                    file = gzfile(vepResultsFile)
                )
                
            } else {
                
                write(
                    x = variantId,
                    file = noResultsVepFile,
                    append = T
                )
                
            }
        }
    }
    
    
    noResultsGeneFile <- file.path(ensemblCacheFolder, "noResults_gene.txt")
    geneCacheFile <- file.path(ensemblCacheFolder, paste0(variantId, "_gene.json.gz"))
    geneResultsFile <- file.path(ensemblFolder, paste0(variantId, "_gene.json.gz"))
    
    geneJson <- NULL
    
    if (file.exists(geneCacheFile)) {
        
        geneJson <- fromJSON(file = gzfile(geneCacheFile))
        jsonString <- toJSON(geneJson)
        write(
            x = jsonString, 
            file = gzfile(geneResultsFile)
        )
        
    } else {
        
        if (file.exists(noResultsGeneFile)) {
            
            noResultsVariants <- readLines(noResultsGeneFile)
            
        } else {
            
            noResultsVariants <- c()
            
        }
        
        if (! variantId %in% noResultsVariants) {
            
            if (chr == 23) {
                
                chr <- 'X'
                
            }
            
            minBp <- max(bp - maxDistance, 0)
            maxBp <- bp + maxDistance
            
            region <- paste0(chr, ":", minBp, "-", maxBp)
            
            ext <- paste0("/overlap/region/human/", region, "?feature=gene")
            
            for (try in 1:100) {
                
                tryCatch(
                    {
                        
                        r <- GET(paste(server37, ext, sep = ""), content_type("application/json"))
                        
                        responseCode <- status_code(r)
                        
                        if (responseCode == 200) {
                            
                            geneJson <- fromJSON(toJSON(content(r)))
                            
                            break()
                            
                        } else {
                            
                            headers <- headers(r)
                            
                            remaining <- headers[["x-ratelimit-remaining"]]
                            reset <- headers[["x-ratelimit-reset"]]
                            
                            if (!is.null(remaining) && !is.null(reset) && remaining <= 1) {
                                
                                Sys.sleep(reset + 1)
                                
                            } else {
                                
                                Sys.sleep(try)
                                
                            }
                        }
                    },
                    error = function(e) {
                        
                        print(e)
                        vepJson <- NULL
                        
                    }
                )
            }
            
            if (is.null(geneJson)) {
                
                return("Error")
                
            }
            
            if (length(geneJson) > 0) {
                
                jsonString <- toJSON(geneJson)
                write(
                    x = jsonString, 
                    file = gzfile(geneResultsFile)
                )
                write(
                    x = jsonString, 
                    file = gzfile(geneCacheFile)
                )
                
            } else {
                
                write(
                    x = variantId,
                    file = noResultsGeneFile,
                    append = T
                )
                
            }
        }
    }
    
    if (length(vepJson) > 0 || length(geneJson) > 0) {
        
        genesLabel <- "[...]"
        genesDF <- NULL
        
        ensemblFile <- file.path(ensemblFolder, paste0(variantId, ".md"))
        
        write(x = paste0("# ", variantId, "\n"), file = ensemblFile, append = F)
        
        bpMin <- max(bp - maxDistance, 0)
        bpMax <- bp + maxDistance
        
        if (startsWith(variantId, prefix = "rs")) {
            
            write(x = paste0("[[GWAS Catalog]](https://www.ebi.ac.uk/gwas/variants/", variantId, ")[[UCSC]](https://genome.ucsc.edu/cgi-bin/hgTracks?position=chr", chr, ":", bpMin, "-", bpMax, "&addHighlight=hg19.chr", chr, "%3A123065528%2D123066028%23fcfcac&hgFind.matches=", variantId, "&db=hg19)[[Ensembl]](https://grch37.ensembl.org/Homo_sapiens/Variation/Explore?r=", chr, ":", bp, "-", bp, ";v=", variantId, ";vdb=variation)"), file = ensemblFile, append = T)
            
        } else {
            
            write(x = paste0("[[Ensembl]](https://grch37.ensembl.org/Homo_sapiens/Location/View?r=", chr, "%3A", bpMin, "-", bpMax, ")"), file = ensemblFile, append = T)
            
        }
        
        if (length(geneJson) > 0) {
            
            write(x = paste0("## Closest genes\n"), file = ensemblFile, append = T)
            
            write(x = paste0("Genes withing 100 kbp\n"), file = ensemblFile, append = T)
            
            genesDF <- data.frame(
                direction = character(length(geneJson)),
                distance = numeric(length(geneJson)),
                name = character(length(geneJson)),
                strand = character(length(geneJson)),
                start = numeric(length(geneJson)),
                end = numeric(length(geneJson)),
                biotype = character(length(geneJson)),
                desciption = character(length(geneJson)),
                id = character(length(geneJson)),
                stringsAsFactors = F
            )
            
            for (geneI in 1:length(geneJson)) {
                
                geneDetails <- geneJson[[geneI]]
                
                if (geneDetails$start <= bp && geneDetails$end >= bp) {
                    
                    distance <- 0
                    direction <- "Intragenic"
                    
                } else if (geneDetails$start > bp) {
                    
                    distance <- geneDetails$start - bp
                    direction <- "Downstream"
                    
                } else if (geneDetails$end < bp) {
                    
                    distance <- bp - geneDetails$end
                    direction <- "Upstream"
                    
                }
                
                description <- ""
                
                if (!is.null(geneDetails$description) && length(geneDetails$description) > 0) {
                    
                    description <- geneDetails$description
                    
                }
                
                name <- ""
                
                if (!is.null(geneDetails$external_name) && length(geneDetails$external_name) > 0) {
                    
                    name <- geneDetails$external_name
                    
                }
                
                biotype <- ""
                
                if (!is.null(geneDetails$biotype) && length(geneDetails$biotype) > 0) {
                    
                    biotype <- geneDetails$biotype
                    
                }
                
                geneId <- ""
                
                if (!is.null(geneDetails$gene_id) && length(geneDetails$gene_id) > 0) {
                    
                    geneId <- geneDetails$gene_id
                    
                }
                
                genesDF$direction[geneI] <- direction
                genesDF$distance[geneI] <- distance
                genesDF$name[geneI] <- name
                genesDF$strand[geneI] <- geneDetails$strand
                genesDF$start[geneI] <- geneDetails$start
                genesDF$end[geneI] <- geneDetails$end
                genesDF$biotype[geneI] <- biotype
                genesDF$desciption[geneI] <- description
                genesDF$id[geneI] <- geneId
                
            }
            
            genesDF %>%
                filter(
                    name != "" | id != "" | description != "" | biotype != ""
                ) %>% 
                mutate(
                    signedDistance = ifelse(direction == "Upstream", -distance, distance)
                ) %>% 
                arrange(
                    signedDistance
                ) -> genesDF
            
            
            write(x = paste0("| Direction | Distance | Name | Strand | Start | End | Biotype | Description | ID |"), file = ensemblFile, append = T)
            write(x = paste0("| --------- | -------- | ---- | ------ | ----- | --- | ------- | ----------- | -- |"), file = ensemblFile, append = T)
            
            for (i in 1:nrow(genesDF)) {
                
                write(x = paste0("| ", genesDF$direction[i], " | ", genesDF$distance[i], " | ", genesDF$name[i], " | ", genesDF$strand[i], " | ", genesDF$start[i], " | ", genesDF$end[i], " | ", genesDF$biotype[i], " | ", genesDF$desciption[i], " | ", genesDF$id[i], " |"), file = ensemblFile, append = T)
                
            }
            
            write(x = "\n", file = ensemblFile, append = T)
            
            genesDF %>% 
                filter(
                    name != ""
                ) -> genesDF
            
            if (nrow(genesDF) > 0) {
                
                i <- which.min(genesDF$distance)
                genesLabel <- genesDF$name[i]
                
            }
        }
        
        if (length(vepJson) > 0) {
            
            if ("transcript_consequences" %in% names(vepJson)) {
                
                consequences <- vepJson[["transcript_consequences"]]
                
                if (!is.null(consequences)) {
                    
                    write(x = paste0("## Transcript Consequences\n"), file = ensemblFile, append = T)
                    
                    consequencesDF <- data.frame(
                        symbol = character(length(consequences)),
                        biotype = character(length(consequences)),
                        consequence = character(length(consequences)),
                        impact = character(length(consequences)),
                        geneId = character(length(consequences)),
                        transcriptId = character(length(consequences)),
                        stringsAsFactors = F
                    )
                    
                    for (consequenceI in 1:length(consequences)) {
                        
                        consequence <- consequences[[consequenceI]]
                        
                        consequencesDF$symbol[consequenceI] <- ifelse(is.null(consequence[["gene_symbol"]]), "", paste(consequence[["gene_symbol"]], collapse = ","))
                        consequencesDF$biotype[consequenceI] <- ifelse(is.null(consequence[["biotype"]]), "", paste(consequence[["biotype"]], collapse = ","))
                        consequencesDF$impact[consequenceI] <- ifelse(is.null(consequence[["impact"]]), "", paste(consequence[["impact"]], collapse = ","))
                        consequencesDF$consequence[consequenceI] <- ifelse(is.null(consequence[["consequence_terms"]]), "", paste(consequence[["consequence_terms"]], collapse = ","))
                        consequencesDF$geneId[consequenceI] <- ifelse(is.null(consequence[["gene_id"]]), "", paste(consequence[["gene_id"]], collapse = ","))
                        consequencesDF$transcriptId[consequenceI] <- ifelse(is.null(consequence[["transcript_id"]]), "", paste(consequence[["transcript_id"]], collapse = ","))
                        
                    }
                    
                    consequencesDF %>% 
                        mutate(
                            impact = tolower(impact),
                            impactFactor = factor(impact, levels = c("high", "moderate", "low", "modifier"))
                        ) %>%
                        arrange(
                            impactFactor, consequence, biotype, symbol
                        ) %>% 
                        distinct() -> consequencesDF
                    
                    
                    write(x = paste0("| Symbol | Biotype | Consequence | Gene ID | Transcript ID |"), file = ensemblFile, append = T)
                    write(x = paste0("| ------ | ------- | ----------- | ------- | ------------- |"), file = ensemblFile, append = T)
                    
                    for (i in 1:nrow(consequencesDF)) {
                        
                        write(x = paste0("| ", consequencesDF$biotype[i], " | ", consequencesDF$symbol[i], " | ", consequencesDF$consequence[i], " | ", consequencesDF$geneId[i], " | ", consequencesDF$transcriptId[i], " |"), file = ensemblFile, append = T)
                        
                    }
                    
                    write(x = "\n", file = ensemblFile, append = T)
                    
                    for (impact in c("high", "moderate", "low", "modifier")) {
                        
                        uniqueGenes <- unique(consequencesDF$symbol[consequencesDF$impact == impact])
                        
                        if (length(uniqueGenes) > 0) {
                            
                            if (!is.null(genesDF) && sum(uniqueGenes %in% genesDF$name) > 0) {
                                
                                bestGenesDF <- genesDF %>%
                                    filter(
                                        name %in% uniqueGenes
                                    )
                                
                                i <- which.min(bestGenesDF$distance)
                                genesLabel <- bestGenesDF$name[i]
                                
                            } else if (length(uniqueGenes) <= 5) {
                                
                                genesLabel <- paste(sort(uniqueGenes), collapse = ", ")
                                
                            }
                            
                            break()
                            
                        }
                    }
                }
            }
            
            
            if ("regulatory_feature_consequences" %in% names(vepJson)) {
                
                consequences <- vepJson[["regulatory_feature_consequences"]]
                
                if (!is.null(consequences)) {
                    
                    write(x = paste0("## Regulatory Features Consequences\n"), file = ensemblFile, append = T)
                    
                    consequencesDF <- data.frame(
                        biotype = character(length(consequences)),
                        consequence = character(length(consequences)),
                        featureId = character(length(consequences)),
                        stringsAsFactors = F
                    )
                    
                    for (consequenceI in 1:length(consequences)) {
                        
                        consequence <- consequences[[consequenceI]]
                        
                        consequencesDF$biotype[consequenceI] <- ifelse(is.null(consequence[["biotype"]]), "", paste(consequence[["biotype"]], collapse = ","))
                        consequencesDF$consequence[consequenceI] <- ifelse(is.null(consequence[["consequence_terms"]]), "", paste(consequence[["consequence_terms"]], collapse = ","))
                        consequencesDF$featureId[consequenceI] <- ifelse(is.null(consequence[["regulatory_feature_id"]]), "", paste(consequence[["regulatory_feature_id"]], collapse = ","))
                        
                    }
                    
                    consequencesDF %>% 
                        distinct() -> consequencesDF
                    
                    
                    write(x = paste0("| Biotype | Consequence | Regulatory Feature ID |"), file = ensemblFile, append = T)
                    write(x = paste0("| ------- |  ---------- | --------------------- |"), file = ensemblFile, append = T)
                    
                    for (i in 1:nrow(consequencesDF)) {
                        
                        write(x = paste0("| ", consequencesDF$biotype[i], " | ", consequencesDF$consequence[i], " | ", consequencesDF$featureId[i], " |"), file = ensemblFile, append = T)
                        
                    }
                    
                    write(x = "\n", file = ensemblFile, append = T)
                    
                }
            }
        }
        
        if (length(vepJson) > 0 || length(geneJson) > 0) {
            
            download <- "Download: "
            
            if (length(vepJson) > 0) {
                
                download <- paste0(download, "[[VEP]](", variantId, "_vep.json.gz)")
                
            }
            
            if (length(geneJson) > 0) {
                
                download <- paste0(download, "[[Genes]](", variantId, "_gene.json.gz)")
                
            }
            
            write(x = paste0(download, "\n"), file = ensemblFile, append = T)
            
        }
        
        return(genesLabel)
        
    } else {
        
        return("No Results")
        
    }
  } else {
    
    print(paste0("Variant ID not supported: ", variantId))
    return("Variant ID not supported")
    
  }
}
