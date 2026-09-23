# Input is csv from CLC genomics - in the reads for the strain, right click on the graph and hit "export as CSV"
#  Need to remove commas from exported CSV file
input1 = read.csv ("CA168_1_exp_both_clean.csv", header=T)
# Normalization 
input1_reorder = c(input1$Value[2000001:4215606], input1$Value[1:1994394])
input1_reorder_last = c(input1$Value[1994395:2000000])
# Change the 39 to be the number of bp per read that you get. 
input1_reorder_norm = input1_reorder/  ((sum(input1_reorder, input1_reorder_last)/150)  / 1e6)
input1_reorder_last_norm = input1_reorder_last/  ((sum(input1_reorder, input1_reorder_last)/150)  / 1e6)
input1_bin = colSums(matrix(input1_reorder_norm,nrow=10000))/10000
input1_bin_last = sum(input1_reorder_last_norm)/5605
input1_bin_whole = c(input1_bin, input1_bin_last)

input2 = read.csv ("d_ponA_s1exp_both_clean.csv", header=T)
# Normalization 
input2_reorder = c(input2$Value[2000001:4215606], input2$Value[1:1994394])
input2_reorder_last = c(input2$Value[1994395:2000000])
# Change the 39 to be the number of bp per read that you get. 
input2_reorder_norm = input2_reorder/  ((sum(input2_reorder, input2_reorder_last)/150)  / 1e6)
input2_reorder_last_norm = input2_reorder_last/  ((sum(input2_reorder, input1_reorder_last)/150)  / 1e6)
input2_bin = colSums(matrix(input2_reorder_norm,nrow=10000))/10000
input2_bin_last = sum(input2_reorder_last_norm)/5605
input2_bin_whole = c(input2_bin, input2_bin_last)

# Plotting the MFA plot
plot(seq(0, 421, 1), input1_bin_whole, type="l", lwd=1.5, axes=FALSE, ylab="Reads per Million", xlab="genome position", col="black", ylim = c(0, 100))
lines(seq(0, 421, 1), input2_bin_whole, ylim=c(0, 100), type="l", lwd=1.5, axes=FALSE, xlab="", ylab="", col="blue")
axis(1, at=c(0, 100, 200, 222, 322, 422), labels=c(2000, 3000, 4000, 0, 1000, 2000))
axis(2)
# abline(v = 415, col = "red")
box()

# saves file
dev.copy2pdf(file="WT1 to d_ponA_s1exp_both_consistentYlim.pdf", width = 5, height = 5)
