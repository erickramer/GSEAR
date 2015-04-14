gsea <- function(rnk, 
                 gmt, 
                 gmt.path=NULL, 
                 nperm=1e3, 
                 output.dir=NULL,
                 output.name=NULL){
  
  gsea.jar = "~/Projects/GSEA/src/gsea/gsea2-2.0.13.jar"
  
  output.dir = if(is.null(output.dir)) tempdir() else output.dir
  
  gmt.file = if(!is.null(gmt.path)) paste(gmt.path, gmt, sep="") else
    paste("~/Projects/GSEA/src/gsea/gmt/", gmt, sep="")
  
  time = Sys.time()
  time = as.character(time)
  time = gsub("-", "_", time, fixed=T)
  time = gsub(" ", "_", time, fixed=T)
  time = gsub(":", "_", time, fixed=T)
  
  rnk.file = paste(output.dir, "/scores_", time, ".rnk", sep="")
  
  output.name = if(is.null(output.name)) paste("result_", time, sep="") else 
    output.name
  
  if(!file.exists(gmt.file)) stop("gmt file not found")
  
  write.table(data.frame(gene=names(rnk),
                         rnk=rnk),
              quote=F,
              sep="\t",
              row.names=F,
              col.names=F,
              file=rnk.file)
  
  cmd = paste("java -Xmx2g -cp", gsea.jar, "xtools.gsea.GseaPreranked",
              "-gmx", gmt.file,
              "-rnk", rnk.file,
              "-rpt_label", output.name,
              "-nperm", nperm,
              "-out", output.dir,
              "-collapse false -set_max 5000 -scoring_scheme classic")
  system(cmd)
  
  # read in findings
  
  result.dir = grep(output.name, list.files(output.dir), value=T) 
  result.id = strsplit(result.dir, ".", fixed=T)[[1]][3]
  
  if(grepl("error", result.dir)) stop("GSEA hit an error")
  
  pos = read.delim(paste(output.dir, "/", 
                         result.dir, "/",
                         "gsea_report_for_na_pos_", result.id, ".xls",
                         sep=""))
  
  neg = read.delim(paste(output.dir, "/", 
                         result.dir, "/",
                         "gsea_report_for_na_neg_", result.id, ".xls",
                         sep=""))
  
  rbind(pos, neg)
}