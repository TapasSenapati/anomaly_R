require(igraph)
args <- commandArgs(trailingOnly = TRUE)
setwd("~/Desktop/R_code/")

#dirname <- args[1]
dirname <- "enron"

total_graphs = length(list.files(paste("./",dirname,"/",sep = "")))   #total number of graphs in time series

no_vertices = 0 #total no of vertex in graph
vertex_list = vector("list",total_graphs) #list of vertices in each graph in timeseries

no_edges = list()  #list of total edges in each graph in timeseries
graph = vector("list",total_graphs)    #list of graphs      

#create all time series graphs
for(i in 1:total_graphs)
{
  filename = paste("./",dirname,"/",i-1,sep = "") 
  vertex_list[[i]] <- read.table(filename, header=T, quote="\"")
  header <- names(vertex_list[[i]])
  no_vertices <- as.numeric(substring(header[1],2))
  no_edges[i] <- as.numeric(substring(header[2],2))
  graph[[i]] <- graph.data.frame(vertex_list[[i]], directed=F)
}


ged = list()
#calculate graph edit distance
for(i in 1:(total_graphs-1))
{
  Vg = no_vertices
  Vh = no_vertices
  Vg_Vh = length(intersect(V(graph[[i]])$name, V(graph[[i+1]])$name))
  Eg = no_edges[i]
  Eh = no_edges[i+1]
  temp1 <- c(paste(vertex_list[[i]][,1], vertex_list[[i]][,2]), paste(vertex_list[[i]][,2], vertex_list[[i]][,1]))
  temp2 <- c(paste(vertex_list[[i+1]][,1], vertex_list[[i+1]][,2]), paste(vertex_list[[i+1]][,2], vertex_list[[i+1]][,1]))
  Eg_Eh = length(intersect(temp1,temp2))/2
  ged[i] = Vg+Vh-2*Vg_Vh+as.numeric(Eg)+as.numeric(Eh)-2*Eg_Eh
}

ged_x = seq(from = 1, to = total_graphs-1, by = 1)
median_ged = median(as.numeric(ged))
sd_ged = sd(as.numeric(ged))
upper_thres = median_ged + 2*sd_ged
lower_thres = median_ged - 2*sd_ged
plot(ged_x,ged,type="l",xlab="days", ylab="edit distance")
abline(h=upper_thres,col="red",lty=2)
abline(h=lower_thres,col="red",lty=2)
title(main="Graph Edit distance", col.main="blue")