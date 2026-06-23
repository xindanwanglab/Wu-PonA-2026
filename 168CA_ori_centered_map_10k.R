r=read.table("Errington_BWX6205_WT168_rep2-10k_overlap.hm_corrected")
for (k in c(0.005, 0.007, 0.01, 0.02)){
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

#Generate the matrix with appropriate cutoff value
n= image_mod(norm_m, zlower, zupper)
# if 10kb bin, then the dimension of n is 404x404

#generate the 2by2 then 4by4 matrix
g = cbind(rbind(n,n), rbind(n,n))

#now need to cut this matrix to put ori in the middle of the map
#bin1: 0-10, bin2:10-20, bin3:20-30
starting_bin = 201
h = g[(starting_bin-1): (starting_bin + nrow(m)-2), (starting_bin-1): (starting_bin + nrow(m)-2)]

par(pty = "s")
image(seq(0, 422, 1), seq(0, 422,1),h, col=colorRampPalette(colvec)(256), zlim=c(zlower, zupper),xlab="genome position (kb)", ylab="genome position (kb)", useRaster=TRUE,axes=FALSE, xaxs="i")
axis(1, at=c(0, 100, 200, 222.1, 322.1, 422.1), labels=c(2000, 3000, 4000,0, 1000, 1999))
axis(2, at=c(0, 100, 200, 222.1, 322.1, 422.1), labels=c(2000, 3000, 4000,0, 1000, 1999))
axis(3, at=c(0, 100, 200, 222.1, 322.1, 422.1), labels=c(2000, 3000, 4000,0, 1000, 1999))
axis(4, at=c(0, 100, 200, 222.1, 322.1, 422.1), labels=c(2000, 3000, 4000,0, 1000, 1999))

box()

#title(ylab="Reads Per Million", xlab="genome position")

dev.copy2pdf(file= paste ("Test_", k, ".pdf", sep=""), useDingbats=FALSE)
}