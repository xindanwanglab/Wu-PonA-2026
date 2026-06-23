r=read.table("BWX6205_LB_rep1-10k_overlap.hm_corrected")
for (k in c(0.01)){
#Define lower and upper cut-off
zlower = 0
zupper = k

#Function to trim the matrix down to min and max cut-off points
image_mod = function(m, zlower, zupper){
  norm_m = m
  for(i in 1:nrow(m) ){
    for (j in 1:ncol(m) ){
      if(m[i,j] >= zupper){
        norm_m[i,j] = zupper
      }
      if(m[i,j] <= zlower){
        norm_m[i,j] = zlower
      }
      if(m[i,j] < zupper & m[i,j] > zlower){
        norm_m[i,j] = m[i,j]
      }
    }
  }
  return(norm_m)
}
#######################################
#White blue purple black linear scheme
colvec= rgb(
  rbind(
    c(255, 255, 255)/255 ,
    c(162, 192, 222)/255 ,
    c(140, 137, 187)/255 ,
    c(140, 87, 167)/255 ,
    c(140, 45, 143)/255 ,
    c(120, 20, 120)/255 ,
    c(90, 15, 90)/255 ,
    c(60, 10, 60)/255 ,
    c(30, 5, 30)/255 ,
    c(0, 0, 0)/255
  )
)
m = as.matrix(r)
#Fill in the diagonal and +1/-1 subdiagonal with max cut-off point
diag(m) = zupper
diag(m[-nrow(m), -1]) = zupper
diag(m[-1, -nrow(m)]) = zupper
norm_m = m

par(pty = "s")
image(seq(0, 422, 1), seq(0, 422,1),image_mod(norm_m,zlower, zupper), col=colorRampPalette(colvec)(256),xlim=c(150, 250), ylim=c(150, 250), xlab="genome position (kb)", ylab="genome position (kb)", useRaster=TRUE,axes=FALSE, xaxs="i")
axis(1, at=c(0, 100, 150, 200,250, 300, 400), labels=c(0,1000,1500, 2000,2500, 3000,4000))
axis(2, at=c(0, 100, 150, 200,250, 300, 400), labels=c(0,1000,1500, 2000,2500, 3000,4000))
axis(3, at=c(0, 100, 150, 200,250, 300, 400), labels=c(0,1000,1500, 2000,2500, 3000,4000))
axis(4, at=c(0, 100, 150, 200,250, 300, 400), labels=c(0,1000,1500, 2000,2500, 3000,4000))
box()

#title(ylab="Reads Per Million", xlab="genome position")

dev.copy2pdf(file= paste ("BWX6205_LB_rep1-10k_terCentered_zoom", k, ".pdf", sep=""), useDingbats=FALSE)
}