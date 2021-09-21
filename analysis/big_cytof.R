library(mclust)
library(umap)
library(MASS)
library(sp)

getSampleID = function(x){
  relevant = strsplit(x,"-")[[1]][1]
  rlist = strsplit(relevant," ")[[1]]
  btype = rlist[1]
  sno = sprintf("%02d",as.numeric(rlist[2]))
  return(paste(btype,sno,sep="_"))
}

# https://stackoverflow.com/questions/30542128/circular-shifting-arrays-in-r-by-distance-n
shifter <- function(x, n = 1) {
  if (n == 0) x else c(tail(x, -n), head(x, n))
}

makepols = function(pol,plot=FALSE){
 a = which.min(pol$x)
 b = which.max(pol$x)
 first = min(a,b)
 second = max( a,b)
 N = length(pol$x)
 shiftind = shifter(1:N,first)
 xvals = pol$x[shiftind]
 yvals = pol$y[shiftind]
 switchpoint = second - first
 upperx = c(xvals[N-1],xvals[1:switchpoint],xvals[switchpoint],xvals[N-1])
 uppery = c(99999999,yvals[1:switchpoint],yvals[switchpoint],99999999)
  
 lowerx = c(xvals[switchpoint],xvals[switchpoint:N],xvals[N],xvals[N])
 lowery = c(-999999999,yvals[switchpoint:N],yvals[N],-999999999)

 if(plot){
  plot(xvals,yvals,type="b",col="black")
  points(upperx,uppery,type="l",col="red",lwd=2)
  points(lowerx,lowery,type="l",col="blue",lwd=2)
 }
 res = list(up=list(x=upperx,y=uppery),down=list(x=lowerx,y=lowery))
 return(res)
}

froots=c("01_Osteocalcin","02_Osteocalcin_VDAC_mask_n","03_Osteocal_VDAC_mask_without_n","04_VDAC_mask_osteocalcin_witho")
froots = c("02_Osteocalcin_VDAC_mask_n")

metad=read.delim("CyToF_Analysis_meta.txt",sep="\t",stringsAsFactors=FALSE)
rownames(metad) = metad$Sample

for(froot in froots){
print(froot)
fname = paste(froot,".txt",sep="")
dir.create(froot)
pngname = file.path(froot,paste(froot,"%03d.png"))
pdfname = paste(froot,"_stripped.pdf",sep="")

sink(paste(froot,"CombinationsQuantified.txt",sep="_"))

pdf(pdfname,height=8.27, width=11.69, useDingbats=FALSE)
#png(pngname,height=827*2, width=1169*2, pointsize=12*(827/480)*2)
 amyv = read.delim("AmySignallingVDAC.txt",header=FALSE,stringsAsFactors=FALSE)

 dat = read.delim(fname,sep="\t",stringsAsFactors=FALSE,na.strings="N/A")
 dat$Sample = sapply(dat$Item.Name,getSampleID)
 samps = sort(unique(dat$Sample))

 # Drop strong outlier ROIs.  Still need to get to bottom of how these occur.
 dropROI = c("Hip 6-1-6","Hip 6-1-7","Hip 5-1-8","Hip 7-1-1","Hip 2-2-2","Hip 3-3-1","Hip 3-3-2","Hip 3-3-3","Hip 3-3-4","Hip 3-3-5","Hip 3-3-6","Hip 3-3-7","Hip 3-3-8","Hip 3-3-9")
 dat = dat[!dat$Item.Name%in%dropROI,]

 op=par(mfrow=c(1,2))
 plot(density(log(dat$Mean.VDAC)),ylim=c(0,1.5),main="log(VDAC) distributions",lwd=2)
 points(density(log(amyv$V1)),type="l",col="red",lwd=2)
 legend("topright",legend=c("DanH","AmyV"),col=c("black","red"),lwd=2)

 plot(density(dat$Mean.VDAC),ylim=c(0,1.0),main="VDAC distributions",lwd=2)
 points(density(amyv$V1),type="l",col="red",lwd=2)
 legend("topright",legend=c("DanH","AmyV"),col=c("black","red"),lwd=2)
 par(op)

 mitochan = "Mean.VDAC"
 chans = colnames(dat)[grepl("Mean.",colnames(dat))]
 chans = chans[chans!=mitochan]
 #chans = paste("Mean.",c("SDHA","NDUFB8","GRIM19","COX4","MTCO1"),sep="")
 
 res = data.frame(matrix(0.0, ncol = length(chans), nrow = length(samps),dimnames=list(samps,chans)))
 res$Age = metad[samps,"Age"]
 res$Sex = metad[samps,"Sex"] 

 res_sdha = res

 res_up = res
 res_down = res

 metad = metad[samps,]
 metad = metad[order(metad$Age),]

 ctrls = c("Paediatric_01","Femur_01")
 #samps = samps[!samps%in%ctrls]
 samps = metad$Sample

 for (s in samps){
 print(paste(froot,s))
 sdat = dat[dat$Sample==s,]
 sdat$Mean.TOM22 = as.numeric(sdat$Mean.TOM22)
 sdat$Item.Name = sub(".*? ", "", sdat$Item.Name)

 mlabs = paste(s,metad[s,"Age"],metad[s,"Sex"])
 mlab = paste(mlabs,"\nN =",length(unique(sdat$ID)))

 cdat = dat[dat$Sample%in%ctrls,]
 cdat$Mean.TOM22 = as.numeric(cdat$Mean.TOM22)

 mat = as.matrix(sdat[,c(chans[!chans%in%c("Mean.DNA1","Mean.DNA2","Mean.TOM22")],mitochan)])
 colnames(mat)=gsub("Mean.","",colnames(mat))

 #mb = Mclust(mat,1:2)
 #sdat$cluster = mb$classification
 #sdat$strength = apply(mb$z,1,max)
 #sdat$ccol = ifelse(sdat$cluster==1,rgb(1,0,0,0.15),rgb(0,0,1,0.15)) 
 
 #colvec = rainbow(mb$G)
 #sdat$ccol = colvec[sdat$cluster]

 #op = par(mfrow=c(3,4),mar=c(4,4,1,1))
 op = par(mfrow=c(2,3),mar=c(4,4,1,1))
 #for(chan in chans){
 #  xaxrng = c(0.001,max(as.numeric(dat[,mitochan]),na.rm=TRUE))
 #  yaxrng = c(0.001,max(as.numeric(dat[,chan]),na.rm=TRUE))
 #  plot(as.numeric(sdat[,mitochan]),as.numeric(sdat[,chan]),xlim=xaxrng, ylim=yaxrng,xlab=mitochan,ylab=chan,pch=16,col=rgb(0,0,0,0.05),log="",cex=0.5)
 #}
 #plot.new()
 #text(0.5,0.5,paste(s,"(linear)",sep="\n"), cex = 2.0,col="red")
 for(chan in chans[chans!="Mean.TOM22"]){
   stripchart(sdat[[chan]]~sdat$Item.Name,vertical=TRUE,ylab=chan,method="jitter",jitter=0.1,pch=16,col=rgb(0,0,0,0.15),las=2,main=mlabs,cex.axis=0.75)
 }
 plot.new()
 plot.new()
 for(chan in chans[chans!="Mean.TOM22"]){
   lchan = log(sdat[[chan]])
   itemname = sdat$Item.Name
   stripchart(lchan[is.finite(lchan)]~itemname[is.finite(lchan)],vertical=TRUE,ylab=paste("log(",chan,")",sep=""),method="jitter",jitter=0.1,pch=16,col=rgb(0,0,0,0.15),las=2,main=mlabs,cex.axis=0.75)
 }
 plot.new()
 plot.new()
 par(op)

 figcomps = c("(CI)","(CI)","(CII)","(CIV)","(CIV)","ATP")
 figlabs = c("A","B","C","D","E","F")
 figchans = c("Mean.NDUFB8","Mean.GRIM19","Mean.SDHA","Mean.COX4","Mean.MTCO1","Mean.OSCP")
 names(figcomps)=figchans
 names(figlabs)=figchans


 patres=list()

 op = par(mfrow=c(2,4),mar=c(4,4,1,1))
 for(mitochan in c("Mean.VDAC")){
 #for(chan in chans[chans!="Mean.TOM22"]){
 for(chan in figchans){
   xaxrng = c(0.001,max(as.numeric(dat[,mitochan]),na.rm=TRUE))
   yaxrng = c(0.001,max(as.numeric(dat[,chan]),na.rm=TRUE))

   xpat = log(as.numeric(sdat[,mitochan]))
   ypat = log(as.numeric(sdat[,chan]))

   xctrl = log(as.numeric(cdat[,mitochan]))
   yctrl = log(as.numeric(cdat[,chan]))
   
   dmat = log(cdat[,c(mitochan,chan)])
   dmat = dmat[is.finite(rowSums(dmat)),]
   # https://stackoverflow.com/questions/16225530/contours-of-percentiles-on-level-plot
   x = dmat[,mitochan]
   y = dmat[,chan]
   dens = kde2d(x, y, n=200); ## estimate the z counts

   prob = c(0.95, 0.5)
   dx = diff(dens$x[1:2])
   dy = diff(dens$y[1:2])
   sz = sort(dens$z)
   c1 = cumsum(sz) * dx * dy
   levels = sapply(prob, function(x) {
    approx(c1, sz, xout = 1 - x)$y
   })

   clines = contourLines(dens,levels=list(levels[1]))
   pip = function(x) return(point.in.polygon(xpat,ypat,x$x,x$y))
   inners = lapply(clines,pip)
   best_inner = which.max(as.numeric(lapply(inners,sum)))

   u_d = makepols(clines[[best_inner]])
   over_exp = pip(u_d$up)
   under_exp = pip(u_d$down)
   patres[[paste(chan,"over",sep="_")]] = over_exp
   patres[[paste(chan,"under",sep="_")]] = under_exp

   inner = inners[best_inner][[1]][is.finite(xpat)&is.finite(ypat)]==1
   prop = sum(!inner)/length(inner)
   propU = sum(over_exp)/length(over_exp)
   propD = sum(under_exp)/length(under_exp)
   if(mitochan=="Mean.VDAC"){
     res[s,chan] = prop
     res_up[s,chan] = propU
     res_down[s,chan] = propD
   }else{res_sdha[s,chan] = prop}
   mlabp = paste(signif(100*(prop),2),"% different, ",signif(100*(propU),2),"% over & ",signif(100*(propD),2),"% under", sep="")

   plot(xctrl,yctrl,xlim=log(xaxrng), ylim=log(yaxrng), main=mlabp,
      xlab=paste("log(",mitochan,")",sep=""),ylab=paste("log(",chan,")",sep=""),pch=16,col=rgb(0,0,0,0.05),
      cex=0.5,cex.axis=1.25,cex.lab=1.4,type="n",cex.main=0.9)

   #points(xpat,ypat,pch=16,col=sdat$ccol,cex=0.5)
   points(xpat[is.finite(xpat)&is.finite(ypat)],ypat[is.finite(xpat)&is.finite(ypat)],pch=16,col=ifelse(inner,rgb(0,0,0,0.05),rgb(1,0,0,0.1)),cex=0.5)

   #lapply(clines, lines, lwd=4, col="green")
   contour(dens, levels=levels, labels=prob, add=T,lwd=2)
   #mtext(side=3, line=-2, text=figlabs[chan], adj=0.02, outer=T,cex=2)
   mtext(side=3, line=-1.75, text=figcomps[chan], adj=0.02, outer=F,cex=1.5)

   
 }
 plot.new()
 plot.new()
 text(0.5,0.5,mlab, cex = 2.5,col="black")
 if(FALSE){
 cm = log(cdat[,mitochan])
 pm = log(sdat[,mitochan])
 mdc = density(cm)
 mdp = density(pm)
 dmax = max(c(mdc$y,mdp$y))
 dxlim = range(c(cm[is.finite(cm)],pm[is.finite(pm)]))
 plot(mdc,xlim=dxlim,ylim=c(0,dmax),main=paste(s,"(log)",sep=" "),xlab=paste("log(",mitochan,")",sep=""))
 points(mdp,type="l",col="red")
 
 um = umap(mat[,colnames(mat)[!colnames(mat)%in%c("DNA.1","DNA.2")]])
 plot(um$layout,pch=16,col=sdat$ccol,xlab="UMAP 1",ylab="UMAP 2")
 }
 }
 par(op)

 patdf = do.call(cbind.data.frame,patres)
 N = dim(patdf)[1]
 colnames(patdf)=gsub("Mean.","",colnames(patdf))
 fchans = gsub("Mean.","",figchans)
 same_complex = list(c("NDUFB8","GRIM19"),c("COX4","MTCO1"))
 directs = list(c("over","over"),c("under","under"),c("over","under"))
 
 for(fchan in fchans){
   A = paste(fchan,"over",sep="_")
   B = paste(fchan,"under",sep="_")
   cat(paste(s,A,B),sep="\n") 
   print(sum(patdf[[A]]&patdf[[B]])/N) # Should be impossible to be classified as over and under (should be zero)
   for (x in c(A,B)){
     cat(paste(s,x),sep="\n")
     cat(sum(patdf[[x]])/N,sep="\n") # Should be as reported in .pdf
   }
 }
 # Should be consistency between proteins from same complex? (should be greater than zero)
 allpairs = combn(fchans,2)

 for (i in 1:ncol(allpairs)){
  a = allpairs[1,i]
  b = allpairs[2,i]
  if(sum(unlist(lapply(same_complex,identical,c(a,b))))) cat("SAME COMPLEX",sep="\n")
  for (direct in directs){ 
    A = paste(a,direct[1],sep="_")
    B = paste(b,direct[2],sep="_")
    cat(paste(s,A,B),sep="\n") 
    cat(sum(patdf[[A]]&patdf[[B]])/N,sep="\n") # If directs same should be high, else should be low
  }
  cat("",sep="\n")
 }

 #sd = sdat[,c(mitochan,chans)]
 #range(sd,na.rm=TRUE)

 #xcoords = 1:dim(sd)[2]
 #xlabs = c(mitochan,chans)

 #plot(NULL,xlim=c(1,length(chans)+1),ylim=c(0.001,100),xlab="",ylab="Raw expression",log="y",xaxt="n",main=paste("Slide",slideno,"Profiles"))
 #axis(1, at = xcoords,labels=xlabs,las=2)
 #for(i in 1:dim(sd)[1]){
 # points(xcoords,sd[i,],type="b",col=rgb(0,0,0,0.1))
 #}

 #sd = sdat[,chans]/sdat[,mitochan]

 #xcoords = 1:dim(sd)[2]
 #xlabs = chans

 #plot(NULL,xlim=c(1,length(chans)+1),ylim=c(0.001,100),xlab="",ylab=paste("Ratio to",mitochan),log="y",xaxt="n",main=paste("Slide",slideno,"Profiles"))
 #axis(1, at = xcoords,labels=xlabs,las=2) 
 #for(i in 1:dim(sd)[1]){
 # points(xcoords,sd[i,],type="b",col=rgb(0,0,0,0.1))
 # points(xcoords,sd[i,],col=rgb(1,0,0,0.1),pch=16)
 #}
 }
par(op)
op = par(mfrow=c(2,3))
for(resdf in list(res,res_up,res_down)){
if(resdf == res){rlab = "Different from controls"}
if(resdf == res_up){rlab = "Above controls"}
if(resdf == res_down){rlab = "Below controls"}
for(ch in chans){
  mlab = paste(rlab,gsub("Mean.","",ch),sep="\n")
  mlab = gsub("OSCP","ATP Synthase",mlab)
  
  plot(resdf$Age,100*resdf[[ch]],pch=16,cex=2,col=ifelse(resdf$Sex=="f",rgb(1,0,0,0.5),rgb(0,0,1,0.5)),xlab = "Age (y)", ylab="Cells different to controls (%)"
  ,ylim=c(0,30),xlim=c(0,100),main=mlab,cex.lab=1.45,cex.axis=1.45,cex.main=2.0)
  abline(h=5,lty=2,lwd=3)
}
plot.new()
legend("topleft",legend=c("Male","Female"),pch=16,col=c(rgb(0,0,1,0.5),rgb(1,0,0,0.5)),cex=2)
#if(identical(resdf,res)) {text(0,0,"x-axis: VDAC1",cex=2,pos=4)}else{text(0,0,"x-axis: SDHA",cex=2,pos=4)}
}
par(op)

dev.off()
sink()
}


pdf("Fig4.pdf",height=8.27, width=11.69, useDingbats=FALSE)
resdf = res_down
chans = c("Mean.NDUFB8","Mean.GRIM19","Mean.SDHA","Mean.MTCO1","Mean.COX4","Mean.OSCP")
op = par(mfrow=c(2,3),mar=c(5,5,2,1))
for(ch in chans){
  mlab = gsub("Mean.","",ch)
  mlab = gsub("OSCP","ATP Synthase",mlab)
  
  plot(resdf$Age,100*resdf[[ch]],pch=16,cex=2,col=ifelse(resdf$Sex=="f",rgb(1,0,0,0.5),rgb(0,0,1,0.5)),xlab = "Age (y)", ylab="Cells below controls (%)"
  ,ylim=c(0,17.5),xlim=c(0,100),main=mlab,cex.lab=2.45,cex.axis=1.45,cex.main=2.0)
  abline(h=5,lty=2,lwd=3)
}
#plot.new()
#legend("topleft",legend=c("Male","Female"),pch=16,col=c(rgb(0,0,1,0.5),rgb(1,0,0,0.5)),cex=2)
#if(identical(resdf,res)) {text(0,0,"x-axis: VDAC1",cex=2,pos=4)}else{text(0,0,"x-axis: SDHA",cex=2,pos=4)}
par(op)
dev.off()




