# Load required libraries
library(survival)
library(sigmoid)
library(cutpointr)

# Define a function to dichotomize a numeric vector 'x' based on quantiles 'pr'.
# Values above the quantile are set to 1, others to 2.
hiloquant <- function(x, pr) {
  qu <- as.vector(quantile(x, pr))
  a <- which(x > qu)
  b <- which(x <= qu)
  x[a] <- 1
  x[b] <- 2
  return(x)
}

# Another function to dichotomize data, possibly for validation or comparison.
hiloquant2 <- function(x, pr) {
  qu <- as.vector(quantile(x, pr))
  a <- which(x > qu)
  b <- which(x <= qu)
  x[a] <- 1
  x[b] <- 2
  return(x)
}

# Function to dichotomize a numeric vector 'x' using a predefined threshold 'pr'.
hiloquant3 <- function(x, pr) {
  a <- which(x >= pr)
  b <- which(x < pr)
  x[a] <- 1
  x[b] <- 2
  return(x)
}

# Function to split a continuous variable 'x' into two groups based on a cutoff 'cf'.
SplitBRCAness <- function(x, cf) {
  qu <- as.vector((x))
  a <- which(x > cf)
  b <- which(x <= cf)
  x[a] <- 1
  x[b] <- 2
  return(x)
}


# Load the dataset from a text file with headers, using comma as separator and period as decimal.
D <- read.table("DataMatrix_VulnerabilityScore_and_expression_MUI.txt", header=TRUE, sep=",", dec=".", check.names=FALSE)

# Convert overall survival from years to months by multiplying by 12.
SURV <- D$OS * 12
SURV1 <- SURV
# Extract the time of death or last follow-up event indicator from the dataset.
EVENT1 <- D$TOD

# Perform univariate Cox proportional hazards modeling for the genes with the indices 226, 28295, 28292.
for (i in c(226, 28295, 28292)){
  # Select the expression level for the ith gene.
  EXPR1 <- D[[i]]
  # Identify non-missing data points for survival time, event, and gene expression.
  inde1 <- which(!is.na(EXPR1) & !is.na(EVENT1) & !is.na(SURV1)) 
  # Run Cox regression using the survival and event data against gene expression.
  s <- coxph(Surv(as.numeric(SURV1[inde1]), as.numeric(EVENT1[inde1])) ~ EXPR1[inde1])
  # Extract and format the hazard ratios and confidence intervals from the model summary.
  HR1 <- sprintf("%.2f", round(summary(s)$conf.int[1], 2))
  HR2 <- sprintf("%.2f", round(summary(s)$conf.int[2], 2))
  HRL <- summary(s)$conf.int[3]
  HRU <- summary(s)$conf.int[4]
  # Extract and format the p-value from the Cox regression summary.
  P <- sprintf("%.2E", summary(s)$coefficients[5])
  # Create a confidence interval string.
  CI <- paste(sprintf("%.2f", round(HRL, 2)), "-", sprintf("%.2f", round(HRU, 2)), sep = "")
  # Prepare a text line containing the formatted results.
  text1 <- paste(names(D)[i], "OS", length(inde1), HR1, HR2, CI, P, sep = "\t")
  # Write the results line to a text file without quotes and with tab separation.
  write.table(text1, file="COXPH_CONT.txt", append=TRUE, quote=FALSE, sep="\t", eol="\n", na="NA", dec=".", row.names=FALSE, col.names=FALSE)
} 


# The loop below performs similar operations to the first loop but focuses on finding the optimal cut-point for gene expression
# that maximizes the rank statistic for survival difference (the best split for dichotomizing the gene expression data).
for (i in c(226,28295,28292)){
  # Again, extract expression data for the current gene.
  EXPR1 <- D[[i]]
  # Ensure data completeness for the variables of interest.
  inde1 <- which(!is.na(EXPR1) & !is.na(EVENT1) & !is.na(SURV1))
  # Initialize a vector to store rank statistics for each potential cut-point.
  rankstat <- NULL
  # Test each percentile from 20% to 80% as a potential cut-point for dichotomizing expression data.
  for (inda in 20:80) {
    qr <- inda/100
    exeo <- hiloquant2(EXPR1[inde1], qr)
    rankstat[inda] <- survdiff(Surv(SURV1[inde1], EVENT1[inde1]) ~ factor(exeo))[[5]]
  }
  # Identify the percentile (cut-point) that gives the maximum rank statistic.
  OP <- which(rankstat == max(rankstat, na.rm=TRUE))[1]
  op <- OP/100
  # Dichotomize the expression data at the optimal percentile.
  ex <- hiloquant2(EXPR1[inde1], op)
  # Perform Cox regression using the dichotomized data.
  s <- coxph(Surv(as.numeric(SURV1[inde1]), as.numeric(EVENT1[inde1])) ~ factor(ex))
  # Extract and format the results as before.
  HR1 <- sprintf("%.2f", round(summary(s)$conf.int[1], 2))
  HR2 <- sprintf("%.2f", round(summary(s)$conf.int[2], 2))
  HRL <- summary(s)$conf.int[3]
  HRU <- summary(s)$conf.int[4]
  PV <- summary(s)$coefficients[5]
  P <- sprintf("%.2E", PV)
  # Adjust the p-value for further analysis or reporting.
  PA <- PV * (-1.63) * (1 + 2.35 * log(PV))
  CI <- paste(sprintf("%.2f", round(HRL, 2)), "-", sprintf("%.2f", round(HRU, 2)), sep="")
  # Concatenate results into a single string for output.
  text1 <- paste(names(D)[i], "OS", length(which(ex == 1)), length(which(ex == 2)), HR1, HR2, CI, P, PA, sep="\t")
  # Write the results to a different text file, specifically for dichotomized data.
  write.table(text1, file="COXPH_DICHOTOM_MAXRANKSTAT.txt", append=TRUE, quote=FALSE, sep="\t", eol="\n", na="NA", dec=".", row.names=FALSE, col.names=FALSE)
}


# The next block of code introduces additional clinical variables into the Cox regression model performing multivariate Cox regression model.
for (i in c(226,28295,28292)){
  # Convert gene expression data to numeric and extract additional clinical variables.
  EXPR1 <- as.numeric(D[[i]])
  EXPR2 <- D$FIGO_einfach  # FIGO stage simplified.
  EXPR3 <- D$Alter  # Patient age.
  # Filter for complete case data across all variables of interest.
  inde1 <- which((!is.na(EXPR1) & !is.na(EVENT1) & !is.na(SURV1)) & !is.na(EXPR2) & !is.na(EXPR3))
  # Run a multivariate Cox regression with the gene expression and clinical variables.
  s <- coxph(Surv(as.numeric(SURV1[inde1]), as.numeric(EVENT1[inde1])) ~ EXPR1[inde1] + EXPR2[inde1] + EXPR3[inde1])
  # Print summary of the Cox model to the console for inspection.
  print(summary(s))
  # Extract and format the results as before.
  HR1 <- sprintf("%.2f", round(summary(s)$conf.int[1], 2))
  HR2 <- sprintf("%.2f", round(summary(s)$conf.int[2], 2))
  HRL <- summary(s)$conf.int[3]
  HRU <- summary(s)$conf.int[4]
  P <- sprintf("%.2E", summary(s)$coefficients[5])
  CI <- paste(sprintf("%.2f", round(HRL, 2)), "-", sprintf("%.2f", round(HRU, 2)), sep="")
  # Prepare the results string.
  text1 <- paste(names(D)[i], "OS", length(inde1), HR1, HR2, CI, P, sep="\t")
  # Write the results to a file, with a filename that includes the gene name.
  write.table(text1, file=paste0("COXPH_CONT_CLIN_CONT", colnames(D)[i], ".txt"), append=TRUE, quote=FALSE, sep="\t", eol="\n", na="NA", dec=".", row.names=FALSE, col.names=FALSE)
}



#************************************
# Set up a PDF device to save the plots, specifying the file location and dimensions.
pdf(file="/data/projects/2020/OvarianCancerHH/OV_Project/Results/MUI/OS_opt_new.pdf", width=12, height=5, useDingbats = FALSE)

# Set graphical parameters for the plots: layout, margins, line width, color, etc.
par(mfrow=c(1,2), tcl=0.5, cex=1.1, lwd=1.5, lend="square", ljoin="mitre", las=1, mgp=c(2,0.3,0), mar=c(6,6,3,2), oma=c(1,1,1,1))

# Define colors for plotting different groups or categories in the survival plots.
colred <- c("#DD0000") # Red color for one group
coldark <- c("#0000DD") # Blue color for another group
colextra <- c("#FFFFFF") # White color, potentially for background or other elements
#************************************

# Set a maximum time point for plotting survival on the x-axis.
mx1 <- 140

# Loop through the specified columns (gene indices) for survival analysis.
for (i in c(28298,28295)){
  # Extract the expression data for the current gene.
  EXPR1 <- D[[i]]
  # Subset the data to remove NA values for expression, event, and survival.
  inde1 <- which(!is.na(EXPR1) & !is.na(EVENT1) & !is.na(SURV1))
  # Subset the survival and event data to the same cases.
  SURV <- SURV1[inde1]
  EVENT <- EVENT1[inde1]
  
  # Loop to create plots for different methods of dichotomizing expression data.
  for (jj in 1:3) {
    
    # For the first method, dichotomize data at the median expression value.
    if (jj == 1) {
      ex <- hiloquant(EXPR1[inde1], 0.5)
    }
    
    # For the second method, use an optimal cutpoint based on Youden's index.
    if (jj == 2) {
      op <- cutpointr(EXPR1[inde1], EVENT1[inde1], pos_class = 1, neg_class = 0, method = maximize_metric, metric = youden)
      ex <- hiloquant3(EXPR1[inde1], op[[2]])
    }
    
    # For the third method, find the cutpoint that maximizes the rank statistic.
    if (jj == 3) {
      rankstat <- NULL
      for (inda in 20:80) {
        qr <- inda / 100
        exeo <- hiloquant2(EXPR1[inde1], qr)
        rankstat[inda] <- survdiff(Surv(SURV1[inde1], EVENT1[inde1]) ~ factor(exeo))[[5]]
      }
      OP <- which(rankstat == max(rankstat, na.rm=TRUE))[1]
      op <- OP / 100
      ex <- hiloquant2(EXPR1[inde1], op)
    }
    
    # Fit survival curves based on the dichotomized data.
    DIS.bygene <- survfit(Surv(SURV, EVENT) ~ ex)
    # Perform a log-rank test to compare the survival curves.
    DIS.bygene.logrank <- survdiff(Surv(SURV, EVENT) ~ ex) 
    
    # Calculate the p-value from the log-rank test.
    PSURV1 <- (1 - pchisq(DIS.bygene.logrank[[5]], df=1))
    # Adjust the p-value using a transformation (purpose-specific and requires context).
    PA <- PSURV1 * (-1.63) * (1 + 2.35 * log(PSURV1))
    
    # Calculate the hazard ratio (HR) from the log-rank test.
    HR <- round(((DIS.bygene.logrank$obs[1]/DIS.bygene.logrank$exp[1])/(DIS.bygene.logrank$obs[2]/DIS.bygene.logrank$exp[2])), 2)
    # Extract the log hazard ratio and its variance.
    lHR1 <- (DIS.bygene.logrank$obs[1] - DIS.bygene.logrank$exp[1]) / (DIS.bygene.logrank$var[1, 1])
    HR1 <- exp(1)^lHR1
    V <- abs(DIS.bygene.logrank$var[1, 1])
    # Calculate the upper and lower confidence limits for the hazard ratio.
    upHR1 <- exp(1)^(lHR1 + 1.96 / sqrt(V))
    loHR1 <- exp(1)^(lHR1 - 1.96 / sqrt(V))
    CI <- paste(sprintf("%.2f", round(loHR1, 2)), sprintf("%.2f", round(upHR1, 2)), sep = "-")
    # Extract the number of subjects at risk in each group from the log-rank test.
    n1 <- as.vector(DIS.bygene.logrank[[1]][1])
    m1 <- as.vector(DIS.bygene.logrank[[1]][2])
    # Plot the survival curves with appropriate colors and line widths.
    plot(DIS.bygene, col = c(colred, coldark), lwd = 1.5, xlab = "Time (months)", ylab = "Overall survival (probability)", xlim = c(0, mx1), cex.lab = 1.1, cex.axis = 1.1, main = "", axes = FALSE, mark.time = TRUE)
    # Add axis labels and a box around the plot.
    axis(side = 2, las = 1, cex.axis = 1.1, tck = -0.015, lwd = 1.5, mgp = c(6, 0.4, 0))
    axis(side = 1, cex.axis = 1.1, tck = -0.015, lwd = 1.5, mgp = c(3, 0.4, 0))
    box(lwd = 1.5)
    # Add annotations to the plot for the p-value and hazard ratio.
    # ...
    # Add annotations for the number of subjects at risk at different time points.
    # ...
  }
}

# Close the PDF device to save the plots to the file.
dev.off()


# Survival analysis where BRCAness is dichotomized based on prediction probability.
#************************************
# Set up the PDF device to save the plots, specifying the file location and dimensions.
pdf(file="OS_BRCAness.pdf", width=12, height=5, useDingbats = FALSE)

# Set the graphical parameters for the layout of the plot area, including the number of rows and columns of plots,
# text label sizes, line widths, and margin sizes.
par(mfrow=c(1,2), tcl=0.5, cex=1.1, lwd=1.5, lend="square", ljoin="mitre", las=1, mgp=c(2,0.3,0), mar=c(6,6,3,2), oma=c(1,1,1,1))

# Define colors to be used in the plots for different groups.
colred <- c("#DD0000") # Red color for high group
coldark <- c("#0000DD") # Blue color for low group
colextra <- c("#FFFFFF") # White color for extra group
#************************************

# Set the maximum x-axis limit for the survival plots.
mx1 <- 140

# Loop over a specific column index (233) to perform survival analysis.
for (i in c(233)){
  # Extract the expression data for the current gene or index.
  EXPR1 <- D[[i]]
  # Identify the subset of data without missing values.
  inde1 <- which(!is.na(EXPR1) & !is.na(EVENT1) & !is.na(SURV1))
  # Subset the survival and event data for the identified cases.
  SURV <- SURV1[inde1]
  EVENT <- EVENT1[inde1] 
  # Dichotomize the BRCAness based on a predefined cutoff probability.
  ex <- SplitBRCAness(EXPR1[inde1], 0.5266)

  # Calculate the survival curves for the dichotomized groups.
  DIS.bygene <- survfit(Surv(SURV, EVENT) ~ ex)
  # Perform a log-rank test to compare the survival curves.
  DIS.bygene.logrank <- survdiff(Surv(SURV, EVENT) ~ ex)			
  
  # Calculate the p-value from the log-rank test.
  PSURV1 <- (1 - pchisq(DIS.bygene.logrank[[5]], df=1))
  # Apply an adjustment to the p-value (the logic behind the adjustment formula should be validated).
  PA <- PSURV1 * (-1.63) * (1 + 2.35 * log(PSURV1))
  
  # Calculate the hazard ratio (HR) and its confidence interval (CI) from the log-rank test.
  HR <- round(((DIS.bygene.logrank$obs[1] / DIS.bygene.logrank$exp[1]) / (DIS.bygene.logrank$obs[2] / DIS.bygene.logrank$exp[2])), 2)
  lHR1 <- (DIS.bygene.logrank$obs[1] - DIS.bygene.logrank$exp[1]) / (DIS.bygene.logrank$var[1, 1])
  HR1 <- exp(1) ^ lHR1
  V <- abs(DIS.bygene.logrank$var[1, 1])
  upHR1 <- exp(1) ^ (lHR1 + 1.96 / sqrt(V))
  loHR1 <- exp(1) ^ (lHR1 - 1.96 / sqrt(V))
  CI <- paste(sprintf("%.2f", round(loHR1, 2)), sprintf("%.2f", round(upHR1, 2)), sep="-")
  n1 <- as.vector(DIS.bygene.logrank[[1]][1])
  m1 <- as.vector(DIS.bygene.logrank[[1]][2])
  
  # Plot the survival curves using the defined colors and line widths.
  plot(DIS.bygene, col=c(colred, coldark), lwd=1.5, xlab="Time (months)", ylab="Overall survival (probability)", xlim=c(0, mx1), cex.lab=1.1, cex.axis=1.1, main="", axes=FALSE, mark.time=TRUE) 				
  # Add axis labels and a box around the plot.
  axis(side=2, las=1, cex.axis=1.1, tck=-0.015, lwd=1.5, mgp=c(6,0.4,0))
  axis(side=1, cex.axis=1.1, tck=-0.015, lwd=1.5, mgp=c(3,0.4,0))
  box(lwd=1.5)
  
  # Add annotations for the p-value and hazard ratio on the plot.
  TXD <- paste("logrank p=", sprintf("%.3f", PSURV1), sep="")
  text(40, 0.9, TXD, col="black", cex=1.1, pos=4)
  text(40, 0.8, paste("HR=", sprintf("%.2f", HR), " (", CI, ")", sep=""), col="black", cex=1.1, pos=4)
  
  # Add annotations for the number of subjects at risk at different time points below the plot.
  text(-7, -0.24, "No. at risk", xpd=TRUE, col="black", cex=1.1, pos=1)
  text(-15, -0.33, "High:", xpd=TRUE, col=colred, cex=1.1, pos=1)
  text(-16, -0.42, "Low:", xpd=TRUE, col=coldark, cex=1.1, pos=1)
  text(mx1 / 2, 1.20, colnames(D)[i], xpd=TRUE, col="black", cex=1.1, pos=1)
  text(0, -0.33, n1, xpd=TRUE, col=colred, cex=1.1, pos=1)
  text(0, -0.42, m1, xpd=TRUE, col=coldark, cex=1.1, pos=1)
  
  # Annotate the number of subjects still at risk at specified time intervals.
  for (j in timeseq) {
    nr1 <- n1 - length(which(OSA < j))
    nr2 <- m1 - length(which(OSB < j))
    text(j, -0.33, nr1, xpd=TRUE, col=colred, cex=1.1, pos=1)
    text(j, -0.42, nr2, xpd=TRUE, col=coldark, cex=1.1, pos=1)
  }
}

# Close the PDF device, saving the plot to the specified file.
dev.off()
