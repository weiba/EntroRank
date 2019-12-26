load('data/new_InfluenceGraph.Rdata')
load('data/Prostate.Rdata')###prostate/breast/lung
compartment=read.csv('data/Gen_compartment_together.csv')
main_function=function(influenceGraph)
{
  inf=influenceGraph[intersect(colnames(patMutMatrix),row.names(influenceGraph)),intersect(colnames(patMutMatrix),row.names(influenceGraph))]
  tt=normalize_grph_freq(inf,patMutMatrix)
  #tt=normalize_grph_freq(ECC_tissue_matrix,gene_freq)
  ent=compute_entrophy(tt$total_graph,compartment,tt$frequency)
  #IF=calculate_IF(ent$subgroup,ent$H_G)
  
  final_result=combine_sort_new(ent$subgroup,compartment,inf)
  final_result
}
normalize_grph_freq=function(total_graph,patMatrix)
{
#   fre_nor=(as.numeric(frequency[,2])-min(as.numeric(frequency[,2])))/(max(as.numeric(frequency[,2]))-min(as.numeric(frequency[,2])))
#   frequency[,2]=fre_nor
  gen_fr=colSums(patMatrix)
  gene_frequence=gen_fr/nrow(patMatrix)
  gene_frequence=cbind(names(gene_frequence),gene_frequence)
  total_graph=(total_graph-min(total_graph))/(max(total_graph)-min(total_graph))
  list(total_graph=total_graph,frequency=gene_frequence)
}
compute_entrophy=function(total_graph,seperate,frequency)
{
  com=unique(seperate[,2])
  group=list()
  for(i in 1:length(com))#####divide matrix into 11 compartment subgroups
  {
    com_i=seperate[which(seperate[,2]==com[i]),1]
    group_i=intersect(colnames(total_graph),com_i)
    group[[i]]=cbind(Gene=group_i,compartment=as.character(com[i]))
  }
  H_G=c()
 for(j in 1:length(group))###calculate the entrophy for each sub-group
 {
   sub_group=total_graph[group[[j]][,1],group[[j]][,1]]
   row.names(frequency)=frequency[,1]
   i_entropy=c()
   
   ##
   for(k in 1:nrow(sub_group))
   {
     
     n=length(which(sub_group[k,]!=0))####number of the neighbor nodes
     
    #######################below is the add_result
#      wij=sum(sub_group[k,])
#      yj=sum(as.numeric(frequency[colnames(sub_group)[which(sub_group[k,]!=0)],2]))
#       
#     Wi=sum(wij+yj)/n+as.numeric(frequency[which(frequency[,1]==row.names(sub_group)[k]),2])
 ##################   
    #######below is multi_result
#      wij=sub_group[k,which(sub_group[k,]!=0)]
#      
#       yj=frequency[colnames(sub_group)[which(sub_group[k,]!=0)],]
#      if(is.null(dim(yj))==TRUE)
#      {
#        Wi=sum(wij*as.numeric(yj[2]))/n*as.numeric(frequency[row.names(sub_group)[k],2])
#      }
#      else
#      {
#        Wi=sum(wij*as.numeric(yj[,2]))/n*as.numeric(frequency[row.names(sub_group)[k],2])
#      }
   #########################below is the cross_entropy of variation frequency
     wij=sub_group[k,which(sub_group[k,]!=0)]
     
     yj=frequency[colnames(sub_group)[which(sub_group[k,]!=0)],]
     if(is.null(dim(yj))==TRUE)
     {
       #Wi=sum(wij*abs(log(as.numeric(frequency[row.names(sub_group)[k],2])/as.numeric(yj[2]))))*as.numeric(frequency[row.names(sub_group)[k],2])
       #Wi=sum(wij*abs(as.numeric(frequency[row.names(sub_group)[k],2])*log(as.numeric(yj[2]))))*as.numeric(frequency[row.names(sub_group)[k],2])
       Wi=sum(wij*log(as.numeric(frequency[row.names(sub_group)[k],2])/as.numeric(yj[2]))*as.numeric(frequency[row.names(sub_group)[k],2]))
#        yi=as.numeric(frequency[row.names(sub_group)[k],2])
#        log_value=log(yi/as.numeric(yj[2]))
#        Wi=sum(wij*1/exp(abs(yi*log_value)))
     }
     else
     {
       #Wi=sum(wij*abs(log(as.numeric(frequency[row.names(sub_group)[k],2])/as.numeric(yj[,2]))))*as.numeric(frequency[row.names(sub_group)[k],2])
       #Wi=sum(wij*abs(as.numeric(frequency[row.names(sub_group)[k],2])*log(as.numeric(yj[,2]))))*as.numeric(frequency[row.names(sub_group)[k],2])
       Wi=sum(wij*log(as.numeric(frequency[row.names(sub_group)[k],2])/as.numeric(yj[,2]))*as.numeric(frequency[row.names(sub_group)[k],2]))
#        yi=as.numeric(frequency[row.names(sub_group)[k],2])
#        log_value=log(yi/as.numeric(yj[,2]))
#        Wi=sum(wij*1/exp(abs(yi*log_value)))
      
     }
#browser()
     ##########################
     
     Pi=n/(length(which(sub_group!=0))/2)
     #Pi=sum(sub_group[k,])/sum(sub_group[which(sub_group!=0)])
     Ei=Wi*Pi*log2(1/Pi)
     
     #Ei=-1*Wi*Pi*log2(Pi)
     if(n==0||Pi==0)
     {
       Ei=0
     }
     i_entropy=rbind(i_entropy,entropy=Ei)
     
   }
   group[[j]]=cbind(group[[j]],entrophy=i_entropy)
    H_k=sum(as.numeric(group[[j]][,3]))
    H_G=rbind(H_G,cbind(as.character(group[[j]][1,2]),H_k))
 }
 list(subgroup=group,H_G=H_G)
}

calculate_IF=function(subgraph,H_G)
{
  final_result=c()
  for(k in 1:length(subgraph))
    {
    IF=c()
     for(i in 1:nrow(subgraph[[k]]))
    {
      Ei=as.numeric(subgraph[[k]][i,3])
      EN_i=as.numeric(H_G[k,2])-Ei
      if(log2(EN_i/Ei)!=0&log2(EN_i/Ei)!=-Inf)
      {
        IF_i=EN_i/(log2(EN_i/Ei))
      }
      else
      {
        IF_i=0
      }
      IF=rbind(IF,cbind(subgraph[[k]][i,1],IF_i))
     }
    final_result=rbind(final_result,IF)
  }
  final_result
}
combine_sort=function(finalresult)
{finalresult1=c()
  for(k in 1:length(finalresult))
  {finalresult1=rbind(finalresult1,finalresult[[k]])
  }
  
comp=compartment[which(duplicated(compartment[,2])==FALSE),2]
  com=c()
  for(i in 1:length(comp))
  {
    com=rbind(com,cbind(as.character(comp[i]),length(which(compartment[,2]==comp[i]))))
  }
  com[,2]=as.numeric(com[,2])/max(as.numeric(com[,2]))
 
  for(j in 1:nrow(finalresult1))
  {
    finalresult1[j,3]=as.numeric(com[which(finalresult1[j,2]==com[,1]),2])*as.numeric(finalresult1[j,3])
  }
  
  
  #####
  gene=unique(finalresult1[,1])
  total_gene=c()
  for(i in 1:length(gene))
  {
    #G_v=max(as.numeric(finalresult1[which(finalresult1[,1]==gene[i]),3]))###the maximize value as the final result
    G_v=sum(as.numeric(finalresult1[which(finalresult1[,1]==gene[i]),3]))###sum up the entropy of each gene.
    total_gene=rbind(total_gene,cbind(gene[i],G_v))
  }
  total_gene1=total_gene[order(as.numeric(total_gene[,2]),decreasing=T),]
  
  
  
  return(total_gene1)
}

combine_sort_new=function(finalresult,compartment,ECC_tissue_matrix)
{
  finalresult1=c()
  for(i in 1:length(finalresult))
  {
    finalresult1=rbind(finalresult1,finalresult[[i]])
  }
  #######weight of compartment size
#   comp=compartment[which(duplicated(compartment[,2])==FALSE),2]
#   com=c()
#   for(i in 1:length(comp))
#   {
#     com=rbind(com,cbind(as.character(comp[i]),length(which(compartment[,2]==comp[i]))))
#   }
#   com[,2]=as.numeric(com[,2])/max(as.numeric(com[,2]))

  ######below is weight of the compartment edges
  comp=compartment[which(duplicated(compartment[,2])==FALSE),2]
  com=c()
  
  for(i in 1:length(comp))
  {
    gen_nam=finalresult1[which(comp[i]==finalresult1[,2]),1]
    com_graph=ECC_tissue_matrix[gen_nam,gen_nam]
    com_value=length(which(com_graph!=0))/length(which(ECC_tissue_matrix!=0))
    com=rbind(com,cbind(as.character(comp[i]),com_value))
 }
  #########
    for(j in 1:nrow(finalresult1))
    {
      finalresult1[j,3]=as.numeric(com[which(finalresult1[j,2]==com[,1]),2])*as.numeric(finalresult1[j,3])
    }
  gene=unique(finalresult1[,1])
  total_gene=c()
  for(i in 1:length(gene))
  {
    G_v=max(as.numeric(finalresult1[which(finalresult1[,1]==gene[i]),3]))###the maximize value as the final result
    #G_v=sum(as.numeric(finalresult1[which(finalresult[,1]==gene[i]),3]))###sum up the entropy of each gene.
    total_gene=rbind(total_gene,cbind(gene[i],G_v))
  }
  total_gene1=total_gene[order(as.numeric(total_gene[,2]),decreasing=T),]
  
  return(total_gene1)
}

