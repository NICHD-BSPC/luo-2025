# Number of chimeric fragments

# number of total reads and mapped reads
total_mapped <- function(fn, fragment_type){
  x <- read.csv(fn, sep='t', header=FALSE)
  return(list(total = x[1,1],
              mapped = x[2,1]))
}

# number of fragments
fragments <- function(fn, fragment_type){
  file.i <- readLines(fn)
  if (fragment_type == 'total') {
    return(length(file.i))
  } else {
    return(length(grep(fragment_type, file.i)))
  }
}


# Interaction plots

print_svg <- function(sp, pxsize=400) {
  lk <- paste0('../data/rnaseq_rilseq/', sp, '/', sp, '_circos_plot.svg')
  cat(paste0('<img width="', pxsize, 'px" src="', lk, '"></img>'))
  cat(paste0('\n\n- [', sp, '_circos_plot.svg](', lk, ')\n\n'))
}



# get_RNAnames: matching RNA names between significant fragments and all chimera fragemnts


convert_to_Granges <- function(df, read_len=0){
  # from 4 columns df (chrom, start, end, strand)
  # return a Grange object
  names(df) <- c('RNA_chr', 'RNA_first_start', 'RNA_last_start', 'RNA_strand')
  return(GRanges(seqnames=as.character(df$RNA_chr),
                 range=IRanges(start=as.numeric(df$RNA_first_start),
                               end=as.numeric(df$RNA_last_start) + read_len),
                 strand=as.character(df$RNA_strand))
  )
}


get_hit <- function(sub.chi.df, sub.fragment.df) {
  # convert the dfs of significant chimeras and all chimeras (fragments)
  # to GRanges so can find their overlap
  # return the matrix of overlaps
  chimera_i <- convert_to_Granges(sub.chi.df)
  fragments_i <- convert_to_Granges(sub.fragment.df)
  hit_i <- as.matrix(findOverlaps(fragments_i, chimera_i))
  return(hit_i)
  }


getChimericFragmentsByGRanges <- function(fragment.df, chi.df, read_len=0) {
  # find the overlaps (get_hit) in fragments for the significant chimeras,
  # for RNA1 and RNA2 of interactions
  # return the target fragments (the unique set between the 2 hit matrices)
  #
  # RNA one in fragment should match to RNA one in chimeras.
  # Only start position is used to find overlap
  hit_one <- get_hit(chi.df[RNA1_cols],
                     fragment.df[RNA1_cols_fg])

  # RNA two in fragment should match to RNA two in chimeras
  hit_two <- get_hit(chi.df[RNA2_cols],
                         fragment.df[RNA2_cols_fg])

  # Find the unique set of same rows in two hit matrices
  hit_rows <- data.table(rbind(hit_one, hit_two), key="queryHits,subjectHits")
  hit_rows <- unique(hit_rows[duplicated(hit_rows)])

  target_fragments <- fragment.df[as.numeric(hit_rows$queryHits),]

  return (target_fragments )
}

getOrderedGeneName <- function(annotation) {
  # Order RNA names based on their genomic coordinates (1-based on plus
  # strand only)
  annotation <- sort(annotation);
  RNA_names <- paste0(annotation$gene_id, collapse=":");
  return (RNA_names)
}

findGeneNamesByOverlaps <- function(RNA_gr, annotation, min_ratio = 0.7) {
  # given a Granges of target fragments, and a Granges of annotation,
  # get the genes associated with the target fragments by overlap with
  # the annotation
  # special cases when several genes, depending on min_ratio, which is the ratio threshold of
  # a number of fragments in a gene to "win" and return only this gene
  # otherwise returns all overlapping genes
  overlaps <- as.matrix(findOverlaps(RNA_gr, annotation))

  if(nrow(overlaps) == 0) {
    stop("the annotation includes IGRs, there should not be non-overlap, but no overlap was found")
  }

  hit_count <- as.data.frame(table(overlaps[,'subjectHits']))
  id_index  <- as.numeric(as.character(hit_count[,'Var1']))
  hit_freq  <- as.numeric(as.character(hit_count[,'Freq']))

  if(length(id_index) == 1) {
    # If only 1 gene overlaps, add it to RNA_name
    RNA_name <- annotation$gene_id[id_index];
  } else {
  # if more than one gene/IGR overlaps with RNA_gr, compare the ratio
    if(min_ratio == 0) {
      # no min_ratio, append all genes names
      RNA_name <- paste0(annotation$gene_id[id_index], collapse=":")
    } else {
    # min_ratio specified, only use the ones > min ratio
      max_count <- max(hit_freq)
      if((max_count/sum(hit_freq)) >= min_ratio) {
        # if one gene is over the ratio threshold, use this one only
        # NOTE: not sure I understand this part
        row_index <- which(hit_freq == max_count)
        id_row <- id_index[row_index]
        RNA_name <- annotation$gene_id[id_row]
      } else {
        # if no gene satisfies the min_ratio, order them
        RNA_name <- getOrderedGeneName(annotation[id_index])
      }
    }
  }
  return (RNA_name)
}



# funcs_unified
getChimeraHits <- function(unified.df, chi.df, read_len=0) {
  # find the overlaps (get_hit) in fragments for the significant chimeras,
  hit_one <- get_hit(unified.df[RNA1_cols],
                     chi.df[RNA1_cols])

  # RNA two in fragment should match to RNA two in chimeras
  hit_two <- get_hit(unified.df[RNA2_cols],
                     chi.df[RNA2_cols])

  # Find the unique set of same rows in two hit matrices
  hit_rows <- data.table(rbind(hit_one, hit_two), key="queryHits,subjectHits")
  hit_rows <- as.matrix(unique(hit_rows[duplicated(hit_rows)]))

  hit_status <- cbind(1:dim(unified_df)[1],
                      rep(0, dim(unified_df)[1]))
  hit_status[match(hit_rows[,2], hit_status[,1]), 2] <- hit_rows[,1]

  return (hit_status)
  }




# sig.summary add annotation
get_ann_col <- function(df.i.col, ann.col, ann.df = ann_df) {
  # df.i.col: df.i[['RNA1.EcoCyc.ID']]
  # ann.col : 'type'
  # returns the corresponding column in annotation to df.i.col
  # split by ':' if necessary, when several genes are in df.i.col
  df.i.col <- strsplit(df.i.col, ':')
  df.i.col <- lapply(df.i.col, function(x) {
                                 if (length(x) == 1){ ann.df[x, ann.col] %>% as.character()}
                                 else {tmp <-  unlist(lapply(x, function(y){ann.df[y, ann.col]}));
                                       return(paste(tmp,collapse=":"))}
                                 }) %>%
                                 unlist()
  return(df.i.col)
}


# save the results in excel file Summary_RILSeq_S-chimeras
exported_excel <- function(sig.chimera, unified_fns, file='Summary_RILSeq_S-chimeras.xlsx'){
  wb <- openxlsx::createWorkbook()

  for (gp in names(unified_fns)) {
    for (sp in unified_fns[[gp]]) {
      # sheetnames can't be > 31 char
      sheetname <- gsub('Total_', '', sp)
      sheetname <- gsub('Repeat', '', sheetname)
      sheetname <- gsub('Stationary', 'Stat', sheetname)
      openxlsx::addWorksheet(wb, sheetname)
      openxlsx::writeData(wb, sheetname, x=sig.chimeras[[sp]] %>% select(any_of(target_cols)), withFilter=TRUE)
      }
    sheetname <- gsub('Total_', '', gp)
    sheetname <- gsub('Repeat', '', sheetname)
    sheetname <- gsub('Stationary', 'Stat', sheetname)
    openxlsx::addWorksheet(wb, sheetname)
    openxlsx::writeData(wb, sheetname, x=sig.chimeras[[gp]][common_target_cols], withFilter=TRUE)

  }
  openxlsx::saveWorkbook(wb, file=file, overwrite=TRUE)
}


# save the results in excel file Summary_interactions_with_IP
exported_excel2 <- function(raw.counts, sp.tabs, ip.tab, file='Summary_interactions_with_IP.xlsx'){
  wb <- openxlsx::createWorkbook()

  # add the raw and normalized counts
  for (gp in names(raw.counts)) {
    openxlsx::addWorksheet(wb, tabs.names[[gp]])
    openxlsx::writeData(wb, tabs.names[[gp]], x=raw.counts[[gp]], withFilter=FALSE)
  }

  # add the IP results
  openxlsx::addWorksheet(wb, 'IP_only')
  openxlsx::writeData(wb, 'IP_only', x=ip.tab, withFilter=FALSE)

  # add individual sample group tabs
  for (sp in names(sp.tabs)) {
    openxlsx::addWorksheet(wb, sp)
    openxlsx::writeData(wb, sp, x=sp.tabs[[sp]], withFilter=FALSE)
  }

  openxlsx::saveWorkbook(wb, file=file, overwrite=TRUE)
}


