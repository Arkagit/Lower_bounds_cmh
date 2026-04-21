########## Multivariate
load("RAB_ess_data.Rdata")


################ Multivariate ESJD
esjd_ram_1 = matrix(0, ncol = 9, nrow = repet)
esjd_ram_2 = matrix(0, ncol = 9, nrow = repet)
esjd_ram_3 = matrix(0, ncol = 9, nrow = repet)

esjd_mh_1 = matrix(0, ncol = 9, nrow = repet)
esjd_mh_2 = matrix(0, ncol = 9, nrow = repet)
esjd_mh_3 = matrix(0, ncol = 9, nrow = repet)

esjd_bark_1 = matrix(0, ncol = 9, nrow = repet)
esjd_bark_2 = matrix(0, ncol = 9, nrow = repet)
esjd_bark_3 = matrix(0, ncol = 9, nrow = repet)

for(i in 1:repet){
	esjd_ram_1[i, ] = as.numeric(alpha_list[[i]][[1]][[4]])
        esjd_ram_2[i, ] = as.numeric(alpha_list[[i]][[2]][[4]])
        esjd_ram_3[i, ] = as.numeric(alpha_list[[i]][[3]][[4]])

	esjd_mh_1[i, ] = as.numeric(alpha_list[[i]][[1]][[8]])
        esjd_mh_2[i, ] = as.numeric(alpha_list[[i]][[2]][[8]])
        esjd_mh_3[i, ] = as.numeric(alpha_list[[i]][[3]][[8]])

	esjd_bark_1[i, ] = as.numeric(alpha_list[[i]][[1]][[12]])
        esjd_bark_2[i, ] = as.numeric(alpha_list[[i]][[2]][[12]])
        esjd_bark_3[i, ] = as.numeric(alpha_list[[i]][[3]][[12]])
}

#colMeans(esjd_ram_1); colMeans(esjd_mh); colMeans(esjd_bark)

Table = matrix(0, nrow = 3, ncol = 4*dim(j.scale)[1])

for(i in 1:4*dim(j.scale)[1]){
        dum1 = c(sum(colMeans(esjd_ram_1)[1:2]), sum(colMeans(esjd_ram_1)[3:4]), sum(colMeans(esjd_ram_1)[5:6]),
                sum(colMeans(esjd_ram_1)[7:8]))
        dum2 = c(sum(colMeans(esjd_ram_2)[1:2]), sum(colMeans(esjd_ram_2)[3:4]), sum(colMeans(esjd_ram_2)[5:6]),
                sum(colMeans(esjd_ram_2)[7:8]))
        dum3 = c(sum(colMeans(esjd_ram_3)[1:2]), sum(colMeans(esjd_ram_3)[3:4]), sum(colMeans(esjd_ram_3)[5:6]),
                sum(colMeans(esjd_ram_3)[7:8]))
        Table[1,] = c(dum1, dum2, dum3)

        dum1 = c(sum(colMeans(esjd_mh_1)[1:2]), sum(colMeans(esjd_mh_1)[3:4]), sum(colMeans(esjd_mh_1)[5:6]),
                sum(colMeans(esjd_mh_1)[7:8]))
        dum2 = c(sum(colMeans(esjd_mh_2)[1:2]), sum(colMeans(esjd_mh_2)[3:4]), sum(colMeans(esjd_mh_2)[5:6]),
                sum(colMeans(esjd_mh_2)[7:8]))
        dum3 = c(sum(colMeans(esjd_mh_3)[1:2]), sum(colMeans(esjd_mh_3)[3:4]), sum(colMeans(esjd_mh_3)[5:6]),
                sum(colMeans(esjd_mh_3)[7:8]))
        Table[2,] = c(dum1, dum2, dum3)

        dum1 = c(sum(colMeans(esjd_bark_1)[1:2]), sum(colMeans(esjd_bark_1)[3:4]), sum(colMeans(esjd_bark_1)[5:6]),
                sum(colMeans(esjd_bark_1)[7:8]))
        dum2 = c(sum(colMeans(esjd_bark_2)[1:2]), sum(colMeans(esjd_bark_2)[3:4]), sum(colMeans(esjd_bark_2)[5:6]),
                sum(colMeans(esjd_bark_2)[7:8]))
        dum3 = c(sum(colMeans(esjd_bark_3)[1:2]), sum(colMeans(esjd_bark_3)[3:4]), sum(colMeans(esjd_bark_3)[5:6]),
                sum(colMeans(esjd_bark_3)[7:8]))
        Table[3,] = c(dum1, dum2, dum3)
}


################ Componentwise ESJD

Esjd_comp_ram_1 = matrix(0, nrow = repet, ncol = 4)
Esjd_comp_ram_2 = matrix(0, nrow = repet, ncol = 4)
Esjd_comp_ram_3 = matrix(0, nrow = repet, ncol = 4)

Esjd_comp_mh_1 = matrix(0, nrow = repet, ncol = 4)
Esjd_comp_mh_2 = matrix(0, nrow = repet, ncol = 4)
Esjd_comp_mh_3 = matrix(0, nrow = repet, ncol = 4)

Esjd_comp_bark_1 = matrix(0, nrow = repet, ncol = 4)
Esjd_comp_bark_2 = matrix(0, nrow = repet, ncol = 4)
Esjd_comp_bark_3 = matrix(0, nrow = repet, ncol = 4)

for(i in 1:4){
	Esjd_comp_ram_1[,i] = rowMeans(esjd_ram_1[ ,c(2*i-1, 2*i)])
        Esjd_comp_ram_2[,i] = rowMeans(esjd_ram_2[ ,c(2*i-1, 2*i)])
        Esjd_comp_ram_3[,i] = rowMeans(esjd_ram_3[ ,c(2*i-1, 2*i)])

	Esjd_comp_mh_1[,i] = rowMeans(esjd_mh_1[ ,c(2*i-1, 2*i)])
        Esjd_comp_mh_2[,i] = rowMeans(esjd_mh_2[ ,c(2*i-1, 2*i)])
        Esjd_comp_mh_3[,i] = rowMeans(esjd_mh_3[ ,c(2*i-1, 2*i)])

	Esjd_comp_bark_1[,i] = rowMeans(esjd_bark_1[ ,c(2*i-1, 2*i)])
        Esjd_comp_bark_2[,i] = rowMeans(esjd_bark_2[ ,c(2*i-1, 2*i)])
        Esjd_comp_bark_3[,i] = rowMeans(esjd_bark_3[ ,c(2*i-1, 2*i)])
}

#Esjd_comp_ram; Esjd_comp_mh; Esjd_comp_bark

######CRAM

set.seed(42)
group <- rep(c("h = 0.01", "h = 1", "h = 10"), each = repet)
values1 <- c(Esjd_comp_ram_1[,1], Esjd_comp_ram_2[,1], Esjd_comp_ram_3[,1])
values2 <- c(Esjd_comp_ram_1[,2], Esjd_comp_ram_2[,2], Esjd_comp_ram_3[,2])
values3 <- c(Esjd_comp_ram_1[,3], Esjd_comp_ram_2[,3], Esjd_comp_ram_3[,3])
values4 <- c(Esjd_comp_ram_1[,4], Esjd_comp_ram_2[,4], Esjd_comp_ram_3[,4])

# Boxplot
pdf("Component_esjd_ram.pdf", width = 20, height = 20)
par(mfrow = c(2, 2), mar = c(6, 8, 6, 4))  

boxplot(values1 ~ group,
        col = c("skyblue", "tomato", "lightgreen"),
        main = "1st Location",
        xlab = "CRAM algorithms",
        ylab = "Mean ESJD",
        cex.lab = 3,
        cex.main = 3,
        cex.axis = 2)

boxplot(values2 ~ group,
        col = c("skyblue", "tomato", "lightgreen"),
        main = "2nd Location",
        xlab = "CRAM algorithms",
        ylab = "Mean ESJD",
        cex.lab = 3,
        cex.main = 3,
        cex.axis = 2)

boxplot(values3 ~ group,
        col = c("skyblue", "tomato", "lightgreen"),
        main = "3rd Location",
        xlab = "CRAM algorithms",
        ylab = "Mean ESJD",
        cex.lab = 3,
        cex.main = 3,
        cex.axis = 2)

boxplot(values4 ~ group,
        col = c("skyblue", "tomato", "lightgreen"),
        main = "4th Location",
        xlab = "CRAM algorithms",
        ylab = "Mean ESJD",
        cex.lab = 3,
        cex.main = 3,
        cex.axis = 2)

dev.off()


###### CMH

group <- rep(c("h = 0.01", "h = 1", "h = 10"), each = repet)
values1 <- c(Esjd_comp_mh_1[,1], Esjd_comp_mh_2[,1], Esjd_comp_mh_3[,1])
values2 <- c(Esjd_comp_mh_1[,2], Esjd_comp_mh_2[,2], Esjd_comp_mh_3[,2])
values3 <- c(Esjd_comp_mh_1[,3], Esjd_comp_mh_2[,3], Esjd_comp_mh_3[,3])
values4 <- c(Esjd_comp_mh_1[,4], Esjd_comp_mh_2[,4], Esjd_comp_mh_3[,4])

# Boxplot
pdf("Component_esjd_mh.pdf", width = 20, height = 20)
par(mfrow = c(2, 2), mar = c(6, 8, 6, 4))  

boxplot(values1 ~ group,
        col = c("skyblue", "tomato", "lightgreen"),
        main = "1st Location",
        xlab = "CMH algorithms",
        ylab = "Mean ESJD",
        cex.lab = 3,
        cex.main = 3,
        cex.axis = 2)

boxplot(values2 ~ group,
        col = c("skyblue", "tomato", "lightgreen"),
        main = "2nd Location",
        xlab = "CMH algorithms",
        ylab = "Mean ESJD",
        cex.lab = 3,
        cex.main = 3,
        cex.axis = 2)

boxplot(values3 ~ group,
        col = c("skyblue", "tomato", "lightgreen"),
        main = "3rd Location",
        xlab = "CMH algorithms",
        ylab = "Mean ESJD",
        cex.lab = 3,
        cex.main = 3,
        cex.axis = 2)

boxplot(values4 ~ group,
        col = c("skyblue", "tomato", "lightgreen"),
        main = "4th Location",
        xlab = "CMH algorithms",
        ylab = "Mean ESJD",
        cex.lab = 3,
        cex.main = 3,
        cex.axis = 2)

dev.off()


###### CBARK

group <- rep(c("h = 0.01", "h = 1", "h = 10"), each = repet)
values1 <- c(Esjd_comp_bark_1[,1], Esjd_comp_bark_2[,1], Esjd_comp_bark_3[,1])
values2 <- c(Esjd_comp_bark_1[,2], Esjd_comp_bark_2[,2], Esjd_comp_bark_3[,2])
values3 <- c(Esjd_comp_bark_1[,3], Esjd_comp_bark_2[,3], Esjd_comp_bark_3[,3])
values4 <- c(Esjd_comp_bark_1[,4], Esjd_comp_bark_2[,4], Esjd_comp_bark_3[,4])

# Boxplot
pdf("Component_esjd_bark.pdf", width = 20, height = 20)
par(mfrow = c(2, 2), mar = c(6, 8, 6, 4))  

boxplot(values1 ~ group,
        col = c("skyblue", "tomato", "lightgreen"),
        main = "1st Location",
        xlab = "CBARK algorithms",
        ylab = "Mean ESJD",
        cex.lab = 3,
        cex.main = 3,
        cex.axis = 2)

boxplot(values2 ~ group,
        col = c("skyblue", "tomato", "lightgreen"),
        main = "2nd Location",
        xlab = "CBARK algorithms",
        ylab = "Mean ESJD",
        cex.lab = 3,
        cex.main = 3,
        cex.axis = 2)

boxplot(values3 ~ group,
        col = c("skyblue", "tomato", "lightgreen"),
        main = "3rd Location",
        xlab = "CBARK algorithms",
        ylab = "Mean ESJD",
        cex.lab = 3,
        cex.main = 3,
        cex.axis = 2)

boxplot(values4 ~ group,
        col = c("skyblue", "tomato", "lightgreen"),
        main = "4th Location",
        xlab = "CBARK algorithms",
        ylab = "Mean ESJD",
        cex.lab = 3,
        cex.main = 3,
        cex.axis = 2)

dev.off()



