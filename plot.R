
setwd('C:\\Users\\Xiao\\Dropbox\\projects\\range_size_dist\\results\\final_results\\')
taxa = c('amphibians', 'mammals')
continents = c('global', 'Africa', 'Eurasia', 'North America', 'South America')
cols = rainbow(length(continents))

# 1. S vs partial S 
png('plots\\partialS.png', width = 480, height = 720)
par(mfcol = c(3, 2))
for (taxon in taxa){
    emp_dir = paste('partial_S\\', taxon, '.csv', sep = '')
    emp_partialS = read.csv(emp_dir, header = F)
    plot(1, type="n", xlim=c(1, 4), ylim=c(0, 1), 
         xlab = 'Range size quartile', ylab = 'r, empirical', cex.axis = 1.5,
         cex.lab = 1.5, main = taxon, cex.main = 1.5)
    for (i in 1:length(continents)){
        continent = continents[i]
        cont_r = emp_partialS[emp_partialS[, 1] == continent, 2:5]
        lines(4:1, cont_r, col = cols[i], lwd = 2, type = 'b', pch = 19)
    }
    if (taxon == taxa[1])
        legend('topleft', continents, col = cols, lwd = 2)
    
    scatter_dir = paste('partial_S\\', taxon, '_sim_scattered.csv', sep = '')
    scatter_partialS = read.csv(scatter_dir, header = F)
    plot(1, type="n", xlim=c(1, 4), ylim=c(0, 1), 
         xlab = 'Range size quartile', ylab = 'r, simulated (scattered)', 
         cex.axis = 1.5, cex.lab = 1.5, main = '')
    for (i in 1:length(continents)){
        continent = continents[i]
        cont_r = scatter_partialS[scatter_partialS[, 1] == continent, 2:5]
        quant = apply(cont_r, 2, function(x) quantile(x, c(0.025, 0.975)))
        polygon(c(4:1,1:4),c(quant[1, ],rev(quant[2, ])),
                col= adjustcolor(cols[i],alpha.f=0.3), border = NA)
    }
    
    contiguous_dir = paste('partial_S\\', taxon, '_sim_contiguous.csv', sep = '')
    contiguous_partialS = read.csv(contiguous_dir, header = F)
    plot(1, type="n", xlim=c(1, 4), ylim=c(0, 1), 
         xlab = 'Range size quartile', ylab = 'r, simulated (contiguous)', 
         cex.axis = 1.5, cex.lab = 1.5, main = '')
    for (i in 1:length(continents)){
        continent = continents[i]
        cont_r = contiguous_partialS[contiguous_partialS[, 1] == continent, 2:5]
        quant = apply(cont_r, 2, function(x) quantile(x, c(0.025, 0.975)))
        polygon(c(4:1,1:4),c(quant[1, ],rev(quant[2, ])),
                col= adjustcolor(cols[i],alpha.f=0.3), border = NA)
    }
}
dev.off()

# 2. S vs S(q)
png('plots\\sq_corr.png', width = 960, height = 480)
par(mfrow = c(1, 2))
for (taxon in taxa){
    sq_corr_dir = paste('sq\\', taxon, '_corr.csv', sep = '')
    sq_corr = read.csv(sq_corr_dir, header = F)
    plot(1, type="n", xlim=c(-1, 1), ylim=c(-0.1, 1), 
         xlab = 'q', ylab = 'Correlation, S(q) vs S', 
         cex.axis = 1.5, cex.lab = 1.5, main = taxon, cex.main = 1.5)
    for (i in 1:length(continents)){
       corr_cont = sq_corr[sq_corr[, 1] == continents[i], ]
       lines(corr_cont[, 2], corr_cont[, 3], col = cols[i], lwd = 2, 
             type = 'b', pch = 19)
    }
    if (taxon == taxa[1])
        legend('bottomleft', continents, col = cols, lwd = 2)
}
dev.off()

# 3. Range size distributions
png('plots\\range_dist.png', width = 960, height = 480)
par(mfrow = c(1, 2))
for (taxon in taxa){
    range_dir = paste('sp_dist\\', taxon, '_range_size.csv', sep = '')
    range_taxon = read.csv(range_dir)
    plot(1, type="n", xlim=c(10, log(max(range_taxon[, 2:ncol(range_taxon)]))), 
         ylim=c(0, 0.2), xlab = 'Log(range size)', ylab = 'Density', 
         cex.axis = 1.5, cex.lab = 1.5, main = taxon, cex.main = 1.5, 
         log = 'x')
    for (i in 1:length(continents)){
        range_cont = range_taxon[, i + 1]
        range_cont = range_cont[range_cont > 0]
        lines(density(log(range_cont)), lwd = 2, col = cols[i])
    }
    if (taxon == taxa[1])
        legend('topleft', continents, col = cols, lwd = 2)
}
dev.off()

# 4. Model S(q) with environmental variables
# By variable
png('plots\\sq_lm.png',  units="in", width=5, height=12.5, res = 300)
par(mfcol = c(5, 2))
for (taxon in taxa){
    lm_dir = paste('sq\\', taxon, '_sq_lm.csv', sep = '')
    lm_taxon = read.csv(lm_dir)
    for (i in 1:5){
        dat_r2 = lm_taxon[, c(1, 2, i + 6)]
        plot(1, type="n", xlim=c(-1, 1), ylim=c(0, 2 * max(dat_r2[, 3])), 
             xlab = 'q', ylab = names(lm_taxon)[i + 6], 
             cex.axis = 1.5, cex.lab = 1.5, main = taxon, cex.main = 1.5)
        for (j in 1:length(continents)){
            continent = continents[j]
            dat_r2_cont = dat_r2[dat_r2$continent == continent, ]
            lines(dat_r2_cont$q, dat_r2_cont[, 3], col = cols[j], 
                  lwd = 2, type = 'l')
        }
        if ((taxon == taxa[1]) & (i == 1))
            legend('topleft', continents, col = cols, cex = 0.7, lwd = 2)
    }
}
dev.off()

# By continent
png('plots\\sq_lm_cont.png',  units="in", width=5, height=12.5, res = 300)
par(mfcol = c(5, 2))
for (taxon in taxa){
    lm_dir = paste('sq\\', taxon, '_sq_lm.csv', sep = '')
    lm_taxon = read.csv(lm_dir)
    for (i in 1:5){
        continent = continents[i]
        dat_r2 = lm_taxon[lm_taxon$continent == continent, ]
        plot(1, type="n", xlim=c(-1, 1), ylim=c(0, 1.5 * max(dat_r2[, 7:11])), 
             xlab = 'q', ylab = 'R-squared', main = continent,
             cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
        for (j in 7:11){
            dat_r2_col = dat_r2[, j]
            lines(dat_r2$q, dat_r2_col, col = cols[j - 6], lwd = 2, type = 'l')
        }
        if ((taxon == taxa[1]) & (i == 1))
            legend('topleft', names(lm_taxon)[7:11], col = cols, cex = 0.7, lwd = 2)
    }
}
dev.off()
