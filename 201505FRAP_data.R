# 201505 Ark FRAP data analysis with R

# read file list
constructs = c('N725', 'N725sim', 'N865', 'N865sim', 'N865simhis')

base_path = '~/Documents/Projects/RING_SUMO/Raw_Data/'
all_files = list(
  N725 = c(Sys.glob(paste0(base_path, '20150518FRAP_ArkNsim_data/*N725_A4*.txt')),
           Sys.glob(paste0(base_path, '20150407U2OS_EGFP-ArkN_FRAP_summary/*N725*.txt'))),
  N725sim = Sys.glob(paste0(base_path, '20150518FRAP_ArkNsim_data/*N725sim_A5*.txt')),
  N865 = c(Sys.glob(paste0(base_path, '20150518FRAP_ArkNsim_data/*N865_A1*.txt')),
           Sys.glob(paste0(base_path, '20150407U2OS_EGFP-ArkN_FRAP_summary/*N865*.txt'))),
  N865sim = Sys.glob(paste0(base_path, '20150518FRAP_ArkNsim_data/*N865sim_A2*.txt')),
  N865simhis = Sys.glob(paste0(base_path, '20150518FRAP_ArkNsim_data/*N865simhis_A3*.txt'))
  )

add_labels = sapply(constructs, function(i) cbind(all_files[[i]], i))
all_data = as.data.frame(Reduce(rbind, add_labels))
names(all_data) = c('file_name', 'construct')
all_data[, 1] = as.vector(all_data[, 1])

# FRAP model formula
# using basic non-linear least square method nls
library(minpack.lm)
fit = function(f){df = read.delim(f); 
                  start = which(df[, 'time_corrected'][1:10] == 0);
                  end = dim(df)[1] - 2;
                  x = df[start:end, 'time_corrected']; 
                  y = df[start:end, 'region1_.'];
                  Y0 = y[1];
                  frap_model = y ~ plateau - (plateau - Y0) * (2.0 ^ (- x / t_half)) + decay * x;
                  nls_fit = nlsLM(frap_model, start=list(plateau=0.5, t_half=20, decay=0));
                  mobile_fraction = (as.numeric(coef(nls_fit)['plateau']) - Y0) / (1 - Y0);
                  return(c(Y0=Y0, coef(nls_fit)['t_half'], 
                           mobile_fraction=mobile_fraction))}

# data collection and plotting
res = cbind(all_data, as.data.frame(t(sapply(all_data$file_name, fit))))

plot.new()
par(mfrow=c(1, 1), bty='n', mgp=c(2.4, 1, 0))
boxplot(res$t_half ~ res$construct, ylim=c(0, 160), 
        ylab='half-recovery time (sec)', outline=FALSE, boxwex=0.6)
stripchart(res$t_half ~ res$construct, ylim=c(0, 160), 
           vertical=TRUE, pch=16, col='blue',
           method='jitter', jitter=0.12, cex=1, add=TRUE)

plot.new()
boxplot(res$mobile_fraction ~ res$construct, ylim=c(0, 1), ylab='mobile fraction',
        outline=FALSE, boxwex=0.6)
stripchart(res$mobile_fraction ~ res$construct, ylim=c(0, 8),  
           vertical=TRUE, pch=16, col='blue',
           method='jitter', jitter=0.12, cex=1, add=TRUE)

# mtext("FRAP statistics", side = 3, line = -2, outer = TRUE)

# pairwise t tests
pairwise.t.test(res$t_half, res$construct, pool.sd=FALSE, paired=FALSE, p.adjust='none')
pairwise.t.test(res$mobile_fraction, res$construct, pool.sd=FALSE, paired=FALSE, p.adjust='none')


