#################################################################
#STAT Functions
#################################################################



#'Title: Coverage Statistics
#' @description 
#'This function calculate coverage statistics for samples. You can filter by a specific trailing character.<br>
#'OSS: expected SAMPLE and MEAN_COVERAGE columns
#'@param df Data frame containing sample coverage data from metrics (metrics/all.wgs_metrics.txt)
#'@param charTerm usually A or B for norm and tumor
#'@export
covStats = function(df, charTerm){
    df %>% ungroup %>% filter(str_ends(SAMPLE, charTerm)) %>% 
    summarize(
		median_value = median(MEAN_COVERAGE),
		sd_value = sd(MEAN_COVERAGE),
		min_value = min(MEAN_COVERAGE),
		max_value = max(MEAN_COVERAGE)) %>% round
}

#'Title: Purity Statistics
#'@description
#'Calculate purity statistics <br>
#'Not take in account 0.3 and NA values <br>
#'OSS: expected Purity column
#'@param df Data frame containing sample statistics of data from facets (facets/cncf/all.summary.txt)
#'@param charTerm usually A or B for norm and tumor or ORG for organoids.
#'@export
purStats = function(df,charTerm=""){
	if (charTerm == "") {
	# Filter by the character term and remove NA and 0.3 values
		df %>% ungroup %>%
		filter(!is.na(Purity) & Purity!=0.3) %>% 
		summarize(
			median_value = median(Purity),
			sd_value = sd(Purity),
			min_value = min(Purity),
			max_value = max(Purity))  %>% round(digits=2)
  	}else {
	# Filter by the character term and remove NA and 0.3 values
		df %>% ungroup %>% filter(str_ends(tumor, charTerm)) %>%
		filter(!is.na(Purity) & Purity!=0.3) %>% 
		summarize(
			median_value = median(Purity),
			sd_value = sd(Purity),
			min_value = min(Purity),
			max_value = max(Purity))  %>% round(digits=2)
   	}
}

#'Title: FGA (Fraction of the genome altered)
#' @description
#' This function calculates the fraction of the genome altered (FGA) based on cncf data from facets. <br>
#' OSS: expected TUMOR_NORMAL, start, end and cnlr.median columns
#'@param df Data frame containing cncf info from facets (facets/all.cncf.txt)
#'@param threshold Value to define a segment altered or not.
#'It is the absolute value of cnlr.median. (default 0.2)
#'@export
FGA = function(df, threshold = 0.2){
	df %>% ungroup %>% group_by(TUMOR_NORMAL)  %>% 
	mutate(
		length = end-start, 
		altered = abs(cnlr.median) > threshold)  %>% 
	summarise(
		FGA = (sum(length[altered]) / sum(length))*100)
}

#'Title: Whole Genome Douplication (WGD)
#' @description
#'This function calculates the whole genome duplication (WGD) status based on ploidy values from facets.
#'OSS: expected tumor and Ploidy columns
#'@param facets Data frame containing facets data (facets/all.summary.txt)
#' @param threshold Value to define a sample as WGD or diploid (default 2.6)
#'@export
WGD = function(facets, threshold=2.6){
	facets %>% ungroup %>%
	mutate(WGD = ifelse(Ploidy>=threshold,"WGD","diploid")) %>% 
	select(tumor, WGD)  %>% 
	rename(TUMOR_SAMPLE = tumor)
} 

#'Title: Tumor Mutational Burden (TMB)
#' @description
#'This function calculates the TMB based on the number of mutations per sample.
#'OSS :It assumes that the input data frame has a column named TUMOR_SAMPLE.
#'@param df Data frame containing mutations per sample (mutect2/strelka2 etc.)
#'@export
TMB = function(df){
   df %>% 
	count(TUMOR_SAMPLE, name = "mut_count") %>% 
	mutate(TMB = round(mut_count / 3000,2) ) %>% 
	select(-mut_count)
}


#'Title: Chek TERT promoter mutations
#' @description
#'This function shows TERT promoter mutations C228T and C250T.
#'These mutations are common in gliomas and other cancers.
#'It filters the input data frame for mutations in the TERT promoter region at specific positions
#'@param df Data frame containing mutations per sample (mutect2/strelka2 etc.)
#'@export
checkTERT = function(df){
    df %>% filter(CHROM	== "chr5", POS == 1295113 | POS == 1295135)
}

#################################################################
#ONCOPLOT Functions
#################################################################


#' Rename Annotation
#' @description
#' Rename the annotation string based on SnpEff data.
#'
#' This function splits the annotation string (potentially multiple annotations joined by '&'),
#' looks up their ranks in the provided SnpEff dataframe, and chooses the highest-priority
#' (lowest rank) annotation. It then recodes some specific labels to 'Other' for consistency.
#'
#' @param annotation Character string to be renamed. Typically a concatenation of variant effect terms.
#' @param snpeff Data frame containing SnpEff effect ontology and ranks. Default is the bundled snpeff.xlsx
#'        from the package's `ref` folder. Need to be passed as argument!
#'
#' @return A character string with the renamed annotation.
#'
#' @export
rename_anno = function(annotation, snpeff){
    annotation_list <- unlist(strsplit(annotation, split="&")[[1]])
    if(length(annotation_list) > 1){
        val = min(subset(snpeff, Effect_seq_ontology %in% annotation_list)$Broad_classes_rank,na.rm=TRUE)
        new_anno = ifelse(val!=Inf,
                          subset(snpeff, Effect_seq_ontology %in% annotation_list) %>%
                              filter(Broad_classes_rank == val) %>%
                              pull(Broad_classes),
                          "Not_interesting")
    } else {
        new_anno = ifelse(annotation_list!=".",
                          subset(snpeff, Effect_seq_ontology %in% annotation_list) %>% pull(Broad_classes),
                          "Not_interesting")
        if(is.na(new_anno)){
            new_anno = "Not_interesting"
        }
    }
    if(new_anno %in% c('Interaction_locus','Start_stop_lost','SV','Promoter')){
        new_anno = "Other"
    }
    return(new_anno)
}


#' Make Mutation Landscape
#' @description
#' Create a mutation landscape matrix from a dataframe of mutations.Samples are cols, genes are rows.
#'
#' This function constructs a mutation landscape matrix where each cell contains the combined annotations
#' for a specific gene in a specific sample. It filters the mutations based on the provided landmark genes,
#' and combines annotations for each sample and gene. The resulting matrix is cleaned to remove unnecessary
#' annotations and ensure proper formatting.
#' @param df  Data frame containing mutations with columns TUMOR_SAMPLE, ANN....GENE, and annotation. Filtere removing not interesting mutations.
#' @param landmark Character vector of landmark genes to filter the mutations.
#' @return A data frame representing the mutation landscape, where rows are genes and columns are samples
#'
#' @export
make_mut_landscape = function(df, landmark){

	# Specify row names and column names
	rownames <- unique(df$TUMOR_SAMPLE)
	colnames <- landmark

	# Create a data frame with all elements as NA
	landscape <- data.frame(matrix("", nrow = length(rownames), ncol = length(colnames),
							dimnames = list(rownames, colnames)))
	for(driver in colnames){
		gene = df  %>% 
			filter(str_detect(ANN....GENE, paste0("\\b",driver,"\\b")))

		# Create the combined annotation for each sample
		new_anno <- gene %>% 
			group_by(TUMOR_SAMPLE) %>% 
			summarize(Combined_Anno = paste(unique(annotation), collapse = ";"))
		
		# Loop through each sample (by row name in landscape)
		for(sample in rownames(landscape)){
			# Check if the sample exists in the new_anno data frame
			matching_anno <- new_anno %>% filter(TUMOR_SAMPLE == sample)
			
			# If there is a match, assign the Combined_Anno value, otherwise assign ""
			landscape[sample, driver] <- ifelse(nrow(matching_anno) > 0,
												matching_anno %>% pull(Combined_Anno),
												"")
		}
	}
	variant = t(landscape)
	variant_filtered <- as.data.frame(
	apply(variant, c(1, 2), function(x) {
		if (is.character(x)) {
		x <- gsub("Not_interesting", "", x) # Remove "notInteresting"
		x <- gsub(";;", ";", x)           # Clean up double semicolons
		x <- gsub("^;|;$", "", x)         # Remove leading/trailing semicolons
		}
		x
	})
	)
	return(variant_filtered)
}

#'Title: Draw Oncoplot
#'@description 
#'Create annotated complex heatmap (oncoplot)
#'@param mat_filtered matrix containing relevant variants
#'@param anno_df annotation data frame containing sample annotations
#'@param title Title of the oncoplot
#'@param rm_empty_rows Logical, if TRUE remove empty rows (default FALSE)
#'@export
oncoplot = function(mat_filtered, anno_df, title,rm_empty_rows = FALSE, keep_gene_order = FALSE, related_samples=FALSE){

    # Define the original column order
    #original_column_order <- colnames(mat_filtered)

    # Custom alteration functions
    alter_fun <- list(
        background = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#CCCCCC", col = NA)),
        Nonsense = function(x, y, w, h) grid.polygon(
                unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w), 
                unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h),
                gp = gpar(fill = col["Nonsense"], col = NA)),
        Missense = function(x, y, w, h) grid.polygon(
                unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
                unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
                gp = gpar(fill = col["Missense"], col = NA)),
        Splice_site = function(x, y, w, h) grid.segments(x - w*0.5, y - h*0.5, x + w*0.5, y + h*0.5, gp = gpar(lwd = 3 ,col = col["Splice_site"])),
        Frameshift = function(x, y, w, h) grid.points(x, y, pch = 10,gp = gpar(col = col["Frameshift"])),
        Other = function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = NA, lwd = 2)),
        Inframe = function(x, y, w, h) grid.segments(x - w*0.5, y + h*0.5, x + w*0.5, y - h*0.5, gp = gpar(lwd = 3, col = col["Inframe"])))

    anno_list <- list()
    col_list <- list()
    
	if(!keep_gene_order){
		row_counts <- apply(mat_filtered != "", 1, sum)  # Count non-empty elements in each row
    	ranked_rows <- rownames(mat_filtered)[order(-row_counts)]  # Order row names by descending counts
	}else{
		ranked_rows <- rownames(mat_filtered)
	}
	# save header names for top annotation
	top_anno = colnames(anno_df)

    if("Histotype" %in% top_anno){
        histotype_vector = anno_df %>% pull(Histotype)
        histotype_vector <- factor(histotype_vector, levels = c("pancreato-biliary", "intestinal", "adenosquamous", "mixed"))
        col_list$Histotype <- c(
            "pancreato-biliary" = "#66C2A5",
            "intestinal" = "#FC8D62",
            "adenosquamous" = "#8DA0CB",
            "mixed" = "#E78AC3"
        )
    }
    if("WGD" %in% top_anno){
        wgd_vector = anno_df %>% pull(WGD)
        wgd_vector <- factor(wgd_vector, levels = c("diploid", "WGD"))
        col_list$WGD <- c("diploid" = "#B3B3B3", "WGD" = "#E31A1C")
    }
    if("MSI" %in% top_anno){
        msi_vector = anno_df %>% pull(MSI)
        msi_vector <- factor(msi_vector)
        col_list$MSI <- colorRamp2(c(0, 20), c("white", "red"))
    }
    if("POLE_POLD" %in% top_anno){
        polepold_vector <- anno_df %>% pull("POLE_POLD")
        polepold_vector <- factor(polepold_vector, levels = c(TRUE, FALSE))  # Set levels explicitly
        col_list$POLE_POLD <- c("TRUE" = "red", "FALSE" = "white")
        anno_df$POLE_POLD <- polepold_vector  # Ensure it's added back to anno_df for HeatmapAnnotation
    }
    if("HRD" %in% top_anno){
        lst_vector = anno_df %>% pull(HRD)
        lst_vector <- factor(lst_vector, levels = c("low", "high"))
        col_list$HRD <- c("low" = "white", "high" = "red")
    }
    if("FGA" %in% top_anno){
        fga_vector = cut(
        anno_df %>% pull(FGA),
        breaks = c(-Inf, 10, 20, Inf),
        labels = c("<10", "10–20", ">20")
        )
        fga_vector <- factor(fga_vector, levels = c("<10", "10–20", ">20"))

        # Add to col_list with appropriate colors
        col_list$FGA <- c(
            "<10" = "white",
            "10–20" = "orange",
            ">20" = "red"
        )

        # Replace column in anno_df so it is included in HeatmapAnnotation
        anno_df$FGA <- fga_vector
    }

	# Create a named vector for TMB counts
	column_counts <- anno_df %>% pull(TMB)

    top_annotation <- HeatmapAnnotation(
        TMB = anno_barplot(column_counts, border = FALSE, gp = gpar(fill = "black")),
        df = anno_df  %>% select(-c(TMB,TUMOR_SAMPLE)),
        col = col_list,
        annotation_name_side = "left"
    )

        
    # Get the full Dark2 palette
    full_palette <- brewer.pal(8, "Dark2")  # "Dark2" has up to 8 colors

    # Select specific colors by their indices (e.g., 1, 3, 5, 6, 8)
    selected_colors <- full_palette[c(1, 8, 3, 4, 6)]
        
    # Assign colors to mutation types
    col <- structure(selected_colors, names = c('Nonsense', 'Missense', 'Splice_site', 'Inframe', 'Frameshift'))

	if(related_samples){
		group_labels <- gsub("ORG|UD","",colnames(mat_filtered))
		group_breaks <- unique(group_labels)
	}else{
		group_labels <- anno_df %>% pull(TUMOR_SAMPLE)
	}
    # Plot
    p = oncoPrint(
        mat_filtered,
        alter_fun = alter_fun,
        alter_fun_is_vectorized = FALSE,
        col = col,
        column_title = paste0("oncoprint with ",title),
        show_column_names = TRUE, 
        column_names_gp = gpar(fontsize = 8),
		column_order = colnames(mat_filtered),
        row_order = ranked_rows,
        remove_empty_rows = rm_empty_rows,
        top_annotation = top_annotation,
		column_split = group_labels
        )

    return(p)
}