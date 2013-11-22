require(igraph)
require(compiler)

args <- commandArgs(trailingOnly = TRUE)
setwd("C:\Users\Tapas\Documents\GitHub\anomaly_R")

#dirname <- args[1]
dirname <- "enron"

total_graphs = length(list.files(paste("./",dirname,"/",sep = "")))   #total number of graphs in time series

no_vertices = list()   #list of total edges in each graph in timeseries
vertex_list = vector("list",total_graphs) #list of vertices in each graph in timeseries

no_edges = list()  #list of total edges in each graph in timeseries
graph = vector("list",total_graphs)    #list of graphs      

#create all time series graphs
for(i in 1:total_graphs){
  filename = paste("./",dirname,"/",i-1,sep = "") 
  vertex_list[[i]] <- read.table(filename, header=T, quote="\"")
  header <- names(vertex_list[[i]])
  no_vertices[i] <- as.numeric(substring(header[1],2))
  no_edges[i] <- as.numeric(substring(header[2],2))
  graph[[i]] <- graph.data.frame(vertex_list[[i]], directed=F)
}

##########################################################################################################################
###########################                     GRAPH EDIT DISTANCE                              #########################
##########################################################################################################################

calculate_ged <- function(x,y) {
  #cat(sprintf("x y: %d %d", x,y))
  Vg = no_vertices[x]
  Vh = no_vertices[y]
  Vg_Vh = no_vertices[x]
  Eg = no_edges[x]
  Eh = no_edges[y]
  temp1 <- c(paste(vertex_list[[x]][,1], vertex_list[[x]][,2]), paste(vertex_list[[x]][,2], vertex_list[[x]][,1]))
  temp2 <- c(paste(vertex_list[[y]][,1], vertex_list[[y]][,2]), paste(vertex_list[[y]][,2], vertex_list[[y]][,1]))
  Eg_Eh = length(intersect(temp1,temp2))/2
  result <- as.numeric(Vg)+as.numeric(Vh)-2*as.numeric(Vg_Vh)+as.numeric(Eg)+as.numeric(Eh)-2*Eg_Eh
  return(result)
}

ged = list()
#calculate graph edit distance
for(i in 1:(total_graphs-1))
  ged[i] = calculate_ged(i,(i+1))

ged_x = seq(from = 1, to = total_graphs-1, by = 1)
med_ged = median(as.numeric(ged))
sd_ged = sd(as.numeric(ged))
upper_thres = med_ged + 2*sd_ged
lower_thres = med_ged - 2*sd_ged

#write to output file
sink("./output/ged_outfile.txt")
cat(upper_thres)
cat(" ")
cat(lower_thres)
for(i in 1:(total_graphs-1)){
  cat("\n")
  cat(i)
  cat(" ")
  cat(paste(ged[i]))
}  
sink()

#png("./plots/ged.#png")
dev.off()

plot(ged_x,ged,type="l",xlab="days", ylab="edit distance")
abline(h=upper_thres,col="red",lty=2)
abline(h=lower_thres,col="red",lty=2)
title(main="Graph Edit distance", col.main="blue")


##########################################################################################################################
###########################               MEDIAN GRAPH EDIT DISTANCE                              ########################
##########################################################################################################################

#return median graph with windows size 5
calculate_g_bar <- function(x) {
  val1 = calculate_ged(x,(x+1)) + calculate_ged(x,(x+2)) + calculate_ged(x,(x+3)) + calculate_ged(x,(x+4))
  val2 = calculate_ged((x+1),x) + calculate_ged((x+1),(x+2)) + calculate_ged((x+1),(x+3)) + calculate_ged((x+1),(x+4))
  val3 = calculate_ged((x+2),x) + calculate_ged((x+2),(x+1)) + calculate_ged((x+2),(x+3)) + calculate_ged((x+2),(x+4))
  val4 = calculate_ged((x+3),x) + calculate_ged((x+3),(x+1)) + calculate_ged((x+3),(x+2)) + calculate_ged((x+3),(x+4))
  val5 = calculate_ged((x+4),x) + calculate_ged((x+4),(x+1)) + calculate_ged((x+4),(x+2)) + calculate_ged((x+4),(x+3))
  new_list = c(val1,val2,val3,val4,val5)
  index <- which.min(new_list)
  if(index==1)
    return(x)
  else if(index==2)
    return((x+1))
  else if(index==3)
    return((x+2))
  else if(index==4)
    return((x+3))
  else
    return((x+4))
}

set_median_graph = list()

for(i in 1:(total_graphs-5)){
  set_median_graph[i] = calculate_g_bar(i)
}

median_ged = list()
#calculate median graph edit distance
for(i in 1:(total_graphs-5))
  median_ged[i] = calculate_ged(as.numeric(set_median_graph[i]),(i+5))

median_ged_x = seq(from = 5, to = total_graphs-1, by = 1)
med_median_ged = median(as.numeric(median_ged))
sd_median_ged = sd(as.numeric(median_ged))
upper_thres = med_median_ged + 2*sd_median_ged 
lower_thres = med_median_ged - 2*sd_median_ged 

#write to output file
sink("./output/median_ged_outfile.txt")
cat(upper_thres)
cat(" ")
cat(lower_thres)
for(i in 1:(total_graphs-5)){
  cat("\n")
  cat(i)
  cat(" ")
  cat(paste(median_ged[i]))
}  
sink()

#png("./plots/median_ged.#png")
dev.off()
plot(median_ged_x,median_ged,type="l",xlab="days", ylab="median graph edit distance")
abline(h=upper_thres,col="red",lty=2)
abline(h=lower_thres,col="red",lty=2)
title(main="Median Graph Edit distance", col.main="blue")

##########################################################################################################################
###########################                     ENTROPY                                          #########################
##########################################################################################################################

entropy<-function(i)
{
    normweight = 1/no_edges[[i]];
    currGraphEntropy = -(no_edges[[i]])*(normweight - log(normweight));
}

entropy_distance = c();

for(i in 1:(total_graphs-1))
  entropy_distance[i] = entropy(i) - entropy(i+1)

entropy_distance[is.nan(entropy_distance)] = 0

median_entropy_x = seq(from = 1, to = (total_graphs-1), by = 1)
med_entropy_distance = median(entropy_distance)
sd_entropy_distance = sd(entropy_distance)
upper_thres = med_entropy_distance + 2*sd_entropy_distance
lower_thres = med_entropy_distance - 2*sd_entropy_distance

#write to output file
sink("./output/entropy_outfile.txt")
cat(upper_thres)
cat(" ")
cat(lower_thres)
for(i in 1:(total_graphs)){
  cat("\n")
  cat(i)
  cat(" ")
  cat(paste(entropy_distance[i]))
}  
sink()

#png("./plots/entropy.#png")
dev.off()
plot(ged_x,ged,type="l",xlab="days", ylab="Entropy distance")
abline(h=upper_thres,col="red",lty=2)
abline(h=lower_thres,col="red",lty=2)
title(main="Entropy distance", col.main="blue")


##########################################################################################################################
###########################                     SPECTRAL DISTANCE                                #########################
##########################################################################################################################
spectralDistance<-function()
{
  prevGraphEigenValues = c();
  currGraphEigenValues = rep(0,5);
  distance = rep(0,total_graphs);
  k=3;#consider k largest eigen values
  for(i in 1:total_graphs)
  {	
    g_laplacian = graph.laplacian(graph[[i]]);
    if(nrow(g_laplacian) > 0)#sometimes the graph laplacian is zero length?
    {
      if(dirname == "enron")
      {
        eigvalues = eigen(g_laplacian, only.values=TRUE)$values;
        currGraphEigenValues = eigvalues[1:k];
      }
      else
      {
        func <- function(x, extra=NULL) { as.vector(g_laplacian %*% x) };
        currGraphEigenValues = arpack(func, options=list(n=vcount(graph[[i]]), nev=k, 
                                                         ncv=8, which="LM", maxiter=200))$values
      }
      if (length(prevGraphEigenValues) == k)#ignore the first one
      {
        sumdiffsq  = 0;
        sumcurrsq = 0;
        sumprevsq = 0;
        minsq = 0;
        for(j in 1:k)
        {
          sumdiffsq = sumdiffsq + (abs(currGraphEigenValues[j]) - abs(prevGraphEigenValues[j]))^2;
          sumcurrsq = sumcurrsq + abs(currGraphEigenValues[j])^2;
          sumprevsq = sumprevsq + abs(prevGraphEigenValues[j])^2;
        }
        minsq = min(sumcurrsq , sumprevsq);
        #cat("Minsq: ", minsq);
        distance[i] = (sumdiffsq/minsq)^0.5;
        #print(distance[i]);
      }
      prevGraphEigenValues = currGraphEigenValues;
      currGraphEigenValues = rep(0,5);
      #print(i)
    }
  }
  return(distance)
}


spectral_distance <- spectralDistance()
spectral_distance[is.na(spectral_distance)] <- 0

spectral_distance_x = seq(from = 1, to = (total_graphs), by = 1)
med_spectral_distance = median(spectral_distance)
sd_spectral_distance= sd(spectral_distance)
upper_thres = med_spectral_distance + 2*sd_spectral_distance
lower_thres = med_spectral_distance - 2*sd_spectral_distance

#write to output file
sink("./output/spectral_distance_outfile.txt")
cat(upper_thres)
cat(" ")
cat(lower_thres)
for(i in 1:(total_graphs)){
  cat("\n")
  cat(i)
  cat(" ")
  cat(paste(spectral_distance[i]))
}  
sink()

#png("./plots/entropy.#png")
dev.off()
plot(spectral_distance_x,spectral_distance,type="l",xlab="days", ylab="Spectral distance")
abline(h=upper_thres,col="red",lty=2)
abline(h=lower_thres,col="red",lty=2)
title(main="Spectral distance", col.main="blue")
