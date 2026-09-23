# Input is csv from CLC genomics - in the reads for the strain, right click on the graph and hit "export as CSV"
#  Need to remove commas from exported CSV file
input1 = read.csv ("332382_dnaBtsdelponAcat_1_clean.csv", header=T)
# Normalization 
input1_reorder = c(input1$Value[1:4210000])
input1_reorder_last = c(input1$Value[4210001:4215606])
# Change the 39 to be the number of bp per read that you get. 
input1_reorder_norm = input1_reorder/  ((sum(input1_reorder, input1_reorder_last)/39)  / 1e6)
input1_reorder_last_norm = input1_reorder_last/  ((sum(input1_reorder, input1_reorder_last)/39)  / 1e6)
input1_bin = colSums(matrix(input1_reorder_norm,nrow=10000))/10000
input1_bin_last = sum(input1_reorder_last_norm)/5605
input1_bin_whole = c(input1_bin, input1_bin_last)

#

# Plotting the MFA plot
plot(seq(0, 421, 1), input1_bin_whole, type="l", lwd=1.5, axes=FALSE, ylim=c(6,14), xlab="", ylab="Reads per Million", col="black")
axis(1, at=c(0, 100, 200, 300, 400), labels=c(0, 1000, 2000, 3000, 4000))
axis(2)


box()

title(ylab="", xlab="genome position")

# saves file
dev.copy2pdf(file="332382_dnaBtsdelponAcat_1_10k.pdf", width = 10, height = 5)








