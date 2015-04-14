parse.gmt <- function(filename){
  f = file(filename)
  lines = readLines(f)
  info = strsplit(lines, "\t", fixed=T)
  
  names(info) = sapply(info, function(x) x[[1]])
  sets = lapply(info, function(x) x[3:length(x)])
  return(sets)
}


library("igraph")
library("plyr")

load("graph.Rdata")
comm = multilevel.community(g)

comm.genes = V(g)$description[comm$membership == 2581]
genes = unique(V(g)$description)

sets = parse.gmt("c2.cgp.v4.0.symbols.gmt")
sets.matrix = llply(sets, function(x) genes %in% x, .progress="text")

comm.vector = genes %in% comm.genes 

f.tables = llply(sets.matrix, function(x) table(x, comm.vector),
                 .progress="text")
names(f.tables) = names(sets)

good.tables = laply(f.tables, function(x) length(x) == 4)
f.tables = f.tables[good.tables]

res = llply(f.tables, fisher.test, .progress="text")

p.values = laply(res, function(x) x$p.value) 
df = data.frame(set=names(res), p=p.values)

df = df[order(df$p), ]