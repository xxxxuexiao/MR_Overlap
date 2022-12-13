library(tidyverse)
library(data.table)
setDTthreads(20)

# if data.frame/data.table read it in directly, otherwise read from path
read_sumstats <- function(path)
{
	if(is.data.frame(path)){ return(path) } else{ return(fread(path, fill = TRUE)) }
}

# @gwasFilePaths: a list with analyzed gwas or gwas paths
# @traits: phenotype strings for related gwas
merge_sumstats <- function(gwasFilePaths, traits, merged_path=NULL)
{	
	if(is.null(merged_path)){
		gwasList = list()
		for (i in 1:length(gwasFilePaths))
		{
			
			gwas = read_sumstats(gwasFilePaths[[i]])
			if(!is.na(traits[i])){
				pheno = traits[i]
				## the first 5 columns are c("SNP", "CHR", "BP", "A1", "A2"), and rename other columns with trait name is it's not na	
				cols = as.character(lapply(colnames(gwas)[-(1:5)], function(x) {paste0(x, "_", pheno)}))
				colnames(gwas) <- c(colnames(gwas)[1:5], cols)
			}
			gwasList = c(gwasList, list(gwas))
		}
		merged = try(gwasList %>% reduce(inner_join, by = c("SNP", "CHR", "BP", "A1", "A2")))
		## if no "CHR", "BP", "A1", "A2" in gwasList, merge data by SNP
		if(class(merged)[1] == "try-error"){merged = gwasList %>% reduce(inner_join, by = "SNP")}	
	} else{
		merged = read_sumstats(merged_path)
	}
	
	return(merged)
}




