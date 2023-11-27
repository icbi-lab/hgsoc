library(ggplot)
data = read.csv('COXPH_CONT_CLIN_CONT_MyeloidCheckpoints1.txt',sep = '\t',header = FALSE)
data = within(data, data1 <- data.frame(do.call('rbind', strsplit(as.character(data$V6),'-',fixed=TRUE))))

label = data$V1
mean = data$V4
lower = data$data1$X1
upper = data$data1$X2

df <- data.frame(label, mean, lower, upper)
df$label <- factor(df$label, levels=rev(df$label))


ggplot(data=df, aes(x=as.factor(df$label), y=as.numeric(df$mean), ymin=as.numeric(df$lower), ymax=as.numeric(df$upper))) +
  geom_pointrange()  +
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Genes") + ylab("HR (95% CI)") +
  theme_bw()  # use a white background
