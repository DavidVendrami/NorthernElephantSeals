########################################
## Get genetic loads

# parse files:
for i in rep*
do
cd $i
for j in gen_*.txt
do
sed -n '/Mutations:/,/Individuals:/p' $j | grep -v 'Mutations:' | grep -v 'Individuals' > ${j%.txt}_parsed.txt
done
cd ../
done

##############################

# Remove neutral mutations if present
# for i in rep*
# do
# cd $i
# for j in gen_*parsed.txt
# do
# grep -v 'm1' $j > ${j%.txt}_parsed_noN.txt
# done
# cd ../
# done


################################

# Genetic load, realized genetic load, masked genetic load, absolute numbers

for i in rep*
do
cd $i
ls -1 *_parsed.txt > file_list.txt
cd ../
done

###############################

for i in rep_*
do
cd $i
Rscript /prj/furseal-genome/David/NES/SLiM/SLiM_3.7.1/WF/Parse_simulations_output_RAD.R # Edit the file to point at the right Ne table
Rscript /prj/furseal-genome/David/NES/SLiM/SLiM_3.7.1/WF/Drift_Load_RAD.R # Edit the file to point at the right Ne table
cd ../
done

# Parse_simulations_output.R
# R
# 
# files<-read.table("file_list.txt",h=F)
# if (dim(files)[1]==25){
# 
# Nes<-read.table("/prj/furseal-genome/David/NES/SLiM/SLiM_3.7.1/WF/Ne_RAD_table.txt",h=T)
# GL<-c() # Genetic load
# RGL<-c() # Realized genetic load
# MGL<-c() # Masked genetic load
# LD<-c() # Number of mildly deleterious mutations
# HD<-c() # Number of highly deleterious mutations
# fLD<-c() # Number of fixed mildly deleterious mutations
# fHD<-c() # Number of fixed highly deleterious mutations
# 
# 
# for (i in 1:length(files[,1])){
# 
# data<-read.table(paste(files[i,1]),h=F)
# Ne<-Nes[i,2]
# colnames(data)<-c("id1","id2","mut","pos","s","h","pop","bla","N")
# 
# #data$q<-data$N/(2*Ne)
# #data$p<-1-data$q
# #data$gl<-2*data$h*data$s*data$p*data$q+data$s*(data$q)^2
# #GL[i]<-sum(data$gl)
# 
# af<-data$N/(2*Ne)
# data$q<-af
# gl<-data$s*data$q
# data$gl<-gl
# GL[i]<-sum(data$gl)
# 
# 
# rgl_p1<-(data$q^2)*data$s
# rgl_p2<-data$q*(1-data$q)*data$h*data$s
# srgl_p1<-sum(rgl_p1)
# srgl_p2<-2*sum(rgl_p2)
# RGL[i]<-srgl_p1+srgl_p2
# 
# MGL[i]<-GL[i]-RGL[i]
# 
# LD[i]<-length(which(data$mut=="m2"))
# 
# HD[i]<-length(which(data$mut=="m3"))
# 
# fLD[i]<-length(which(data$mut=="m2"&data$N==(Ne*2)))
# 
# fHD[i]<-length(which(data$mut=="m3"&data$N==(Ne*2)))
# }
# 
# write.table(GL,"GL.txt",quote=F,col.names=F,row.names=F,sep="\t")
# write.table(RGL,"RGL.txt",quote=F,col.names=F,row.names=F,sep="\t")
# write.table(MGL,"MGL.txt",quote=F,col.names=F,row.names=F,sep="\t")
# write.table(LD,"LD.txt",quote=F,col.names=F,row.names=F,sep="\t")
# write.table(HD,"HD.txt",quote=F,col.names=F,row.names=F,sep="\t")
# write.table(fLD,"fLD.txt",quote=F,col.names=F,row.names=F,sep="\t")
# write.table(fHD,"fHD.txt",quote=F,col.names=F,row.names=F,sep="\t")
# 
# } else {
# 
# Nes<-read.table("/prj/furseal-genome/David/NES/SLiM/SLiM_3.7.1/WF/Ne_RAD_table.txt",h=T)
# GL<-c() # Genetic load
# RGL<-c() # Realized genetic load
# MGL<-c() # Masked genetic load
# LD<-c() # Number of mildly deleterious mutations
# HD<-c() # Number of highly deleterious mutations
# fLD<-c() # Number of fixed mildly deleterious mutations
# fHD<-c() # Number of fixed highly deleterious mutations
# 
# 
# for (i in 1:length(files[,1])){
# 
# data<-read.table(paste(files[i,1]),h=F)
# Ne<-Nes[i,2]
# colnames(data)<-c("id1","id2","mut","pos","s","h","pop","bla","N")
# 
# #data$q<-data$N/(2*Ne)
# #data$p<-1-data$q
# #data$gl<-2*data$h*data$s*data$p*data$q+data$s*(data$q)^2
# #GL[i]<-sum(data$gl)
# 
# af<-data$N/(2*Ne)
# data$q<-af
# gl<-data$s*data$q
# data$gl<-gl
# GL[i]<-sum(data$gl)
# 
# 
# rgl_p1<-(data$q^2)*data$s
# rgl_p2<-data$q*(1-data$q)*data$h*data$s
# srgl_p1<-sum(rgl_p1)
# srgl_p2<-2*sum(rgl_p2)
# RGL[i]<-srgl_p1+srgl_p2
# 
# MGL[i]<-GL[i]-RGL[i]
# 
# LD[i]<-length(which(data$mut=="m2"))
# 
# HD[i]<-length(which(data$mut=="m3"))
# 
# fLD[i]<-length(which(data$mut=="m2"&data$N==(Ne*2)))
# 
# fHD[i]<-length(which(data$mut=="m3"&data$N==(Ne*2)))
# }
# 
# GL[(length(files[,1])+1):25]<-NA
# RGL[(length(files[,1])+1):25]<-NA
# MGL[(length(files[,1])+1):25]<-NA
# LD[(length(files[,1])+1):25]<-NA
# HD[(length(files[,1])+1):25]<-NA
# fLD[(length(files[,1])+1):25]<-NA
# fHD[(length(files[,1])+1):25]<-NA
# 
# write.table(GL,"GL.txt",quote=F,col.names=F,row.names=F,sep="\t")
# write.table(RGL,"RGL.txt",quote=F,col.names=F,row.names=F,sep="\t")
# write.table(MGL,"MGL.txt",quote=F,col.names=F,row.names=F,sep="\t")
# write.table(LD,"LD.txt",quote=F,col.names=F,row.names=F,sep="\t")
# write.table(HD,"HD.txt",quote=F,col.names=F,row.names=F,sep="\t")
# write.table(fLD,"fLD.txt",quote=F,col.names=F,row.names=F,sep="\t")
# write.table(fHD,"fHD.txt",quote=F,col.names=F,row.names=F,sep="\t")}

# Drift_Load.R
# R
# 
# files<-read.table("file_list.txt",h=F)
# if (dim(files)[1]==25){
# 
# Nes<-read.table("/prj/furseal-genome/David/NES/SLiM/SLiM_3.7.1/WF/Ne_RAD_table.txt",h=T)
# DL<-c() # drift load
# 
# data<-read.table(paste(files[1,1]),h=F)
# Ne<-Nes[1,2]
# colnames(data)<-c("id1","ID","mut","pos","s","h","pop","GEN","N")
# 
# drift<-data[data$N==(2*Ne),]
# ids<-drift$ID
# DL[1]<-sum(drift$s)
# 
# for (i in 2:length(files[,1])){
# 
# data<-read.table(paste(files[i,1]),h=F)
# Ne<-Nes[i,2]
# colnames(data)<-c("id1","ID","mut","pos","s","h","pop","GEN","N")
# 
# drift<-data[data$N==(2*Ne),]
# #ind<-which(drift$ID%in%ids==FALSE)
# #ids<-sort(c(ids,drift$ID[ind]))
# #drift<-drift[ind,]
# DL[i]<-sum(drift$s)
# }
# 
# #write.table(DL,"DL.txt",quote=F,col.names=F,row.names=F,sep="\t")
# 
# } else {
# 
# Nes<-read.table("/prj/furseal-genome/David/NES/SLiM/WF/Ne_table.txt",h=T)
# DL<-c() # drift load
# 
# data<-read.table(paste(files[1,1]),h=F)
# Ne<-Nes[1,2]
# colnames(data)<-c("id1","ID","mut","pos","s","h","pop","GEN","N")
# 
# drift<-data[data$N==(2*Ne),]
# ids<-drift$ID
# DL[1]<-sum(drift$s)
# 
# for (i in 1:length(files[,1])){
# 
# data<-read.table(paste(files[i,1]),h=F)
# Ne<-Nes[i,2]
# colnames(data)<-c("id1","id2","mut","pos","s","h","pop","bla","N")
# 
# drift<-data[data$N==(2*Ne),]
# #ind<-which(drift$ID%in%ids==FALSE)
# #ids<-sort(c(ids,drift$ID[ind]))
# #drift<-drift[ind,]
# DL[i]<-sum(drift$s)
# }
# 
# DL[(length(files[,1])+1):25]<-NA
# 
# #write.table(DL,"DL.txt",quote=F,col.names=F,row.names=F,sep="\t")
# }

####################################################

# Merge outputs from different simulations
paste rep_*/GL.txt > GL_all.txt
paste rep_*/RGL.txt > RGL_all.txt
paste rep_*/MGL.txt > MGL_all.txt
paste rep_*/DL.txt > DL_all.txt

##################################################

