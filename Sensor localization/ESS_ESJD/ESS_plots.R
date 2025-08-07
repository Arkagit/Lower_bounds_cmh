########## Multivariate ESS
load("RAB_ess_data.Rdata")

ess_ram = matrix(0, ncol = length(samp_size), nrow = repet)
ess_mh = matrix(0, ncol = length(samp_size), nrow = repet)
ess_bark = matrix(0, ncol = length(samp_size), nrow = repet)

for(i in 1:repet){
	ess_ram[i, ] = as.numeric(alpha_list[[i]][[2]])
	ess_mh[i, ] = as.numeric(alpha_list[[i]][[6]])
	ess_bark[i, ] = as.numeric(alpha_list[[i]][[10]])
}

colMeans(ess_ram); colMeans(ess_mh); colMeans(ess_bark)


################ Componentwise ESS

ess_comp_ram = matrix(0, ncol = 4, nrow = repet)
ess_comp_mh = matrix(0, ncol = 4, nrow = repet)
ess_comp_bark = matrix(0, ncol = 4, nrow = repet)

for(i in 1:repet){
	ess_comp_ram[i, ] = as.numeric(alpha_list[[i]][[3]])
	ess_comp_mh[i, ] = as.numeric(alpha_list[[i]][[7]])
	ess_comp_bark[i, ] = as.numeric(alpha_list[[i]][[11]])
}

colMeans(ess_comp_ram); colMeans(ess_comp_mh); colMeans(ess_comp_bark)


################ Multivariate ESJD
esjd_ram = matrix(0, ncol = 9, nrow = repet)
esjd_mh = matrix(0, ncol = 9, nrow = repet)
esjd_bark = matrix(0, ncol = 9, nrow = repet)

for(i in 1:repet){
	esjd_ram[i, ] = as.numeric(alpha_list[[i]][[4]])
	esjd_mh[i, ] = as.numeric(alpha_list[[i]][[8]])
	esjd_bark[i, ] = as.numeric(alpha_list[[i]][[12]])
}

colMeans(esjd_ram); colMeans(esjd_mh); colMeans(esjd_bark)

################ Componentwise ESJD

Esjd_comp_ram = matrix(0, nrow = repet, ncol = 4)
Esjd_comp_mh = matrix(0, nrow = repet, ncol = 4)
Esjd_comp_bark = matrix(0, nrow = repet, ncol = 4)

for(i in 1:4){
	Esjd_comp_ram[,i] = rowMeans(esjd_ram[ ,c(2*i-1, 2*i)])
	Esjd_comp_mh[,i] = rowMeans(esjd_mh[ ,c(2*i-1, 2*i)])
	Esjd_comp_bark[,i] = rowMeans(esjd_bark[ ,c(2*i-1, 2*i)])
}

Esjd_comp_ram; Esjd_comp_mh; Esjd_comp_bark



set.seed(42)
group <- rep(c("CRAM", "CMH", "CRAB"), each = 100)
values1 <- c(Esjd_comp_ram[,1], Esjd_comp_mh[,1], Esjd_comp_bark[,1])
values2 <- c(Esjd_comp_ram[,2], Esjd_comp_mh[,2], Esjd_comp_bark[,2])
values3 <- c(Esjd_comp_ram[,3], Esjd_comp_mh[,3], Esjd_comp_bark[,3])
values4 <- c(Esjd_comp_ram[,4], Esjd_comp_mh[,4], Esjd_comp_bark[,4])

# Boxplot
pdf("Component_esjd.pdf", width = 20, height = 20)
par(mfrow = c(2,2), mar = c(6, 8, 6, 4))  

boxplot(values1 ~ group,
        col = c("skyblue", "tomato", "lightgreen"),
        main = "1st Location",
        xlab = "CARB algorithms",
        ylab = "Mean ESJD",
        cex.lab = 3,
        cex.main = 3,
        cex.axis = 2)

boxplot(values2 ~ group,
        col = c("skyblue", "tomato", "lightgreen"),
        main = "2nd Location",
        xlab = "CARB algorithms",
        ylab = "Mean ESJD",
        cex.lab = 3,
        cex.main = 3,
        cex.axis = 2)

boxplot(values3 ~ group,
        col = c("skyblue", "tomato", "lightgreen"),
        main = "3rd Location",
        xlab = "CARB algorithms",
        ylab = "Mean ESJD",
        cex.lab = 3,
        cex.main = 3,
        cex.axis = 2)

boxplot(values4 ~ group,
        col = c("skyblue", "tomato", "lightgreen"),
        main = "4th Location",
        xlab = "CARB algorithms",
        ylab = "Mean ESJD",
        cex.lab = 3,
        cex.main = 3,
        cex.axis = 2)

dev.off()
