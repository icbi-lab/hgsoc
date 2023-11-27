# Load necessary libraries for survival analysis and cutpoint analysis
library(survival)
library(sigmoid)
library(cutpointr)

# Define a function to dichotomize a vector based on quantiles.
hiloquant <- function(x, pr) {
  # Calculate the quantile value for the given probability
  qu <- as.vector(quantile(x, pr))
  # Identify which values are above the quantile
  a <- which(x > qu)
  # Identify which values are at or below the quantile
  b <- which(x <= qu)
  # Assign a value of 1 to values above the quantile
  x[a] <- 1
  # Assign a value of 2 to values at or below the quantile
  x[b] <- 2
  return(x)
}

# Define a similar function for dichotomization as hiloquant.
# This might be used for validation or to ensure reproducibility.
hiloquant2 <- function(x, pr) {
  # Calculate the quantile value for the given probability
  qu <- as.vector(quantile(x, pr))
  # Identify which values are above the quantile
  a <- which(x > qu)
  # Identify which values are at or below the quantile
  b <- which(x <= qu)
  # Assign a value of 1 to values above the quantile
  x[a] <- 1
  # Assign a value of 2 to values at or below the quantile
  x[b] <- 2
  return(x)
}

# Define a function to dichotomize a vector based on a specific threshold.
hiloquant3 <- function(x, pr) {
  # Identify which values are above or equal to the threshold
  a <- which(x >= pr)
  # Identify which values are below the threshold
  b <- which(x < pr)
  # Assign a value of 1 to values above or equal to the threshold
  x[a] <- 1
  # Assign a value of 2 to values below the threshold
  x[b] <- 2
  return(x)
}

# Define a function to dichotomize BRCAness based on a cutoff value.
SplitBRCAness <- function(x, cf) {
  # The input 'x' is assumed to be a vector of BRCAness scores.
  # Convert the input vector to a numeric vector
  qu <- as.vector((x))
  # Identify which scores are above the cutoff
  a <- which(x > cf)
  # Identify which scores are at or below the cutoff
  b <- which(x <= cf)
  # Assign a value of 1 to scores above the cutoff
  x[a] <- 1
  # Assign a value of 2 to scores at or below the cutoff
  x[b] <- 2
  return(x)
}


# Read the dataset from a specified file path into a dataframe 'D'.
# The dataset is expected to contain vulnerability scores and gene expressions.
D <- read.table("DataMatrix_VulnerabilityScore_and_expression.csv", header=TRUE, sep=",", dec=".", check.names=FALSE)

# Extract survival time in months and event occurrence from 'D'
SURV <- D$OS_MONTHS
SURV1 <- SURV
EVENT1 <- D$OS

# Univariate Cox regression analysis for continous  gene expression data
for (i in c(73,89,28186,28184)){
  # Extract the expression levels for the current gene by column index
  EXPR1 <- D[[i]]
  
  # Find indices of complete cases where expression levels, event indicator, and survival times are not missing
  inde1 <- which(!is.na(EXPR1) & !is.na(EVENT1) & !is.na(SURV1)) 
  
  # Perform Cox proportional hazards regression for the non-missing cases
  # This models the hazard of the event as a function of gene expression
  s <- coxph(Surv(as.numeric(SURV1[inde1]), as.numeric(EVENT1[inde1])) ~ EXPR1[inde1])
  
  # Extract the lower and upper bounds of the hazard ratio (HR) confidence interval and round them to 2 decimal places
  HR1 <- sprintf("%.2f", round(summary(s)$conf.int[1], 2))
  HR2 <- sprintf("%.2f", round(summary(s)$conf.int[2], 2))
  
  # Extract the precise lower and upper bounds of the HR confidence interval without rounding
  HRL <- summary(s)$conf.int[3]
  HRU <- summary(s)$conf.int[4]
  
  # Format the p-value of the model to scientific notation
  P <- sprintf("%.2E", summary(s)$coefficients[5])
  
  # Concatenate the rounded confidence interval bounds into a string
  CI <- paste(sprintf("%.2f", round(HRL, 2)), "-", sprintf("%.2f", round(HRU, 2)), sep="")
  
  # Create a string that combines the gene name, overall survival (OS), number of cases, HR, and p-value
  text1 <- paste(names(D)[i], "OS", length(inde1), HR1, HR2, CI, P, sep="\t")
  
  # Write the results to a text file, appending to it if it already exists. This creates a tab-separated values file without quotes.
  write.table(text1, file="COXPH_CONT.txt", append=TRUE, quote=FALSE, sep="\t", eol="\n", na="NA", dec=".", row.names=FALSE, col.names=FALSE)
}


# Univariate Cox regression analysis for dichotomized gene expression data
for (i in c(73,89,28186,28184)) {
  # Extract the expression levels for the current gene by its column index
  EXPR1 <- D[[i]]
  
  # Identify the subset of data that has non-missing values for gene expression, survival, and event data
  inde1 <- which(!is.na(EXPR1) & !is.na(EVENT1) & !is.na(SURV1))
  
  # Dichotomize the gene expression data based on the median value
  # using the hiloquant2 function, where values above the median are set to 1 and below or equal are set to 2
  ex <- hiloquant2(EXPR1[inde1], 0.5)
  
  # Fit a Cox proportional hazards model using the dichotomized gene expression data as a factor
  s <- coxph(Surv(as.numeric(SURV1[inde1]), as.numeric(EVENT1[inde1])) ~ factor(ex))
  
  # Extract the estimated hazard ratios (HR) and their confidence intervals (CI) from the Cox model
  # and format these values to two decimal places
  HR1 <- sprintf("%.2f", round(summary(s)$conf.int[1], 2))
  HR2 <- sprintf("%.2f", round(summary(s)$conf.int[2], 2))
  HRL <- summary(s)$conf.int[3] # Lower limit of the confidence interval
  HRU <- summary(s)$conf.int[4] # Upper limit of the confidence interval
  
  # Convert the p-value of the model to scientific notation for readability
  P <- sprintf("%.2E", summary(s)$coefficients[5])
  
  # Construct the confidence interval string by concatenating the lower and upper limits
  CI <- paste(sprintf("%.2f", round(HRL, 2)), "-", sprintf("%.2f", round(HRU, 2)), sep="")
  
  # Create a summary string containing the gene name, sample size for each group, HR, CI, and p-value
  # This string is tab-separated to prepare for output to a text file
  text1 <- paste(names(D)[i], "OS", length(which(ex==1)), length(which(ex==2)), HR1, HR2, CI, P, sep="\t")
  
  # Write the summary string to a text file. If the file already exists, the new data is appended.
  # The output format is specified without quotes, with tab separators, and with a newline at the end of each gene's results.
  write.table(text1, file="COXPH_DICHOTOM.txt", 
              append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE, col.names = FALSE)
}

#Univariate Cox regression analysis for dichotomized gene expression data with optimization of dichotomization cutoff
# Loop over selected gene indices for advanced survival analysis
for (i in c(73,89,28186,28184)){
   # Extract the expression data for the ith gene
   EXPR1 <- D[[i]]
   # Identify the subset of data without missing values for expression, event, and survival
   inde1 <- which(!is.na(EXPR1) & !is.na(EVENT1) & !is.na(SURV1))
   
   # Initialize an empty vector to store rank statistics for different quantiles
   rankstat <- NULL
   
   # Loop through quantiles from 20% to 80% to find the optimal cutpoint for dichotomizing gene expression
   for (inda in 20:80) {
     qr <- inda/100 # Convert the quantile to a proportion
     # Dichotomize gene expression data at the current quantile
     exeo <- hiloquant2(EXPR1[inde1], qr)
     # Perform survival difference analysis and store the rank statistic for the current quantile
     rankstat[inda] <- survdiff(Surv(SURV1[inde1], EVENT1[inde1]) ~ factor(exeo))[[5]] 
   }
   
   # Find the quantile with the maximum rank statistic, which will be used as the cutpoint
   OP <- which(rankstat == max(rankstat, na.rm=TRUE))[1]
   op <- OP/100 # Convert the index to a proportion
   # Dichotomize the expression data at the optimal cutpoint
   ex <- hiloquant2(EXPR1[inde1], op)
   
   # Perform Cox proportional hazards regression using the dichotomized data
   s <- coxph(Surv(as.numeric(SURV1[inde1]), as.numeric(EVENT1[inde1])) ~ factor(ex))
   
   # Extract hazard ratios (HR) and confidence intervals (CI) from the model and format them
   HR1 <- sprintf("%.2f", round(summary(s)$conf.int[1], 2))
   HR2 <- sprintf("%.2f", round(summary(s)$conf.int[2], 2))
   HRL <- summary(s)$conf.int[3]
   HRU <- summary(s)$conf.int[4]
   # Format the p-value from the model for output
   PV <- summary(s)$coefficients[5]
   P <- sprintf("%.2E", PV)
   # Adjust the p-value using a transformation (reason for adjustment should be clear from the study design)
   PA <- PV * (-1.63) * (1 + 2.35 * log(PV))
   # Concatenate the confidence interval limits into a string
   CI <- paste(sprintf("%.2f", round(HRL, 2)), "-", sprintf("%.2f", round(HRU, 2)), sep = "")
   
   # Create a formatted string with gene name, number of observations, HRs, CI, p-value, and adjusted p-value
   text1 <- paste(names(D)[i], "OS", length(which(ex == 1)), length(which(ex == 2)), HR1, HR2, CI, P, PA, sep = "\t")
   
   # Write the results to a text file, appending if it already exists, with specified formatting
   write.table(text1, file="COXPH_DICHOTOM_MAXRANKSTAT.txt", append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
}


#Multivariate Cox regression analysis adding clinical paramets as covariates for continous gene expression data
for (i in c(1733,23093,17409,14145,14163,16627,13978,14696,19145,16377,6886,12038,18312,3328,24910,1023,20665)){
  
  # Subset the data to only include patients with available tumor residual disease information
  D1 <- D[!is.na(D$TUMOR_RESIDUAL_DISEASE),]
  # Extract the expression level for the current gene
  EXPR1 = D1[[i]]
  
  # Extract survival time and event information from the subsetted data
  SURV2 = D1$OS_MONTHS
  EVENT2 = D1$OS

  # Recode clinical stage from string to numeric format for Cox regression
  # Each stage is converted to a numeric level: IIA, IIB, IIC -> 2; IIIA, IIIB, IIIC -> 3; IV -> 4
  D1$CLINICAL_STAGE[D1$CLINICAL_STAGE == 'IIA'] <- 2
  D1$CLINICAL_STAGE[D1$CLINICAL_STAGE == 'IIB'] <- 2
  D1$CLINICAL_STAGE[D1$CLINICAL_STAGE == 'IIC'] <- 2
  D1$CLINICAL_STAGE[D1$CLINICAL_STAGE == 'IIIA'] <- 3
  D1$CLINICAL_STAGE[D1$CLINICAL_STAGE == 'IIIB'] <- 3
  D1$CLINICAL_STAGE[D1$CLINICAL_STAGE == 'IIIC'] <- 3
  D1$CLINICAL_STAGE[D1$CLINICAL_STAGE == 'IV'] <- 4

  # Extract BRCAness status and relevel the factor so 'noBRCAness' is the reference category
  EXPR2 <- D1$BRCAness
  EXPR2 = relevel(as.factor(EXPR2), 'noBRCAness')

  # Create a binary variable for tumor residual disease, with 'NO MACROSCOPIC DISEASE' as 0, otherwise 1
  EXPR3 <- c()
  for (x in D1$TUMOR_RESIDUAL_DISEASE){
    if (x == "NO MACROSCOPIC DISEASE"){
      EXPR3 <- c(EXPR3, 0)
    } else {
      EXPR3 <- c(EXPR3, 1)
    }
  }
  
  # Extract patient age
  EXPR4 <- D1$AGE
  # Convert the clinical stage to a factor for inclusion in the model
  EXPR5 <- as.factor(D1$CLINICAL_STAGE)
  
  # Find indices of cases with complete data across all clinical and gene expression variables
  inde1 <- which((!is.na(EXPR1) & !is.na(EVENT2) & !is.na(SURV2)) & !is.na(EXPR2) & !is.na(EXPR3) & !is.na(EXPR4) & !is.na(EXPR5))
  
  # Perform Cox regression with the gene expression and clinical variables
  s <- coxph(Surv(as.numeric(SURV2[inde1]), as.numeric(EVENT2[inde1])) ~ EXPR1[inde1] + EXPR2[inde1])
  # Print the name of the gene and the summary of the Cox regression model to the console
  print(names(D1)[i])
  print(summary(s))
  
  # Extract and format the hazard ratio (HR) and confidence intervals (CI) from the Cox model summary
  HR1 <- sprintf("%.2f", round(summary(s)$conf.int[1], 2))
  HR2 <- sprintf("%.2f", round(summary(s)$conf.int[2], 2))
  HRL <- summary(s)$conf.int[1, 3]
  HRU <- summary(s)$conf.int[1, 4]
  
  # Convert the p-value of the model to scientific notation
  P <- sprintf("%.2E", summary(s)$coefficients[1, 5])
  # Construct the confidence interval string
  CI <- paste(sprintf("%.2f", round(HRL, 2)), "-", sprintf("%.2f", round(HRU, 2)), sep = "")
  
  # Create a string of the results for output
  text1 <- paste(names(D1)[i], "OS", length(inde1), HR1, HR2, CI, P, sep = "\t")
  
  # Write the string to a text file, appending if it already exists
  write.table(text1, file=paste0("COXPH_CONT_CLIN_CONT.txt"), append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
}

## Add Kaplan Meier curves for several parameters with different dichotomisation methods 
#************************************
# Set up a PDF device for plotting survival curves. This is where the plots will be saved.
pdf(file="OS_Kaplan_Meier.pdf", width=12, height=5, useDingbats = FALSE)

# Configure the plotting area to have 1 row and 2 columns of plots,
# with specific margins and text sizes for the axes and labels.
par(mfrow=c(1,2), tcl=0.5, cex=1.1, lwd=1.5, lend="square", ljoin="mitre", las=1, mgp=c(2,0.3,0), mar=c(6,6,3,2), oma=c(1,1,1,1))

# Define colors for plotting to visually distinguish different groups in the survival curves.
colred <- c("#DD0000") # Red color for one group
coldark <- c("#0000DD") # Blue color for another group
colextra <- c("#FFFFFF") # White color, not used in this particular block

# Maximum x-axis limit for the survival plot, presumably to standardize the plot range.
mx1 <- 140

# This line is trying to access a column 'ratio_CYT_C1QA' in the data frame 'D', but it doesn't do anything with it in this block of code.
# D$ratio_CYT_C1QA

# Loop through a series of gene indices to perform survival analysis
for (i in c(73,89,28186,28184)){
  # Extract expression data for the current gene.
  EXPR1 <- D[[i]]
  # Identify complete cases where expression data, event data, and survival data are not missing.
  inde1 <- which(!is.na(EXPR1) & !is.na(EVENT1) & !is.na(SURV1))
  # Subset survival time and event indicator for the complete cases.
  SURV <- SURV1[inde1]
  EVENT <- EVENT1[inde1]

  # Loop through three methods of dichotomizing expression data for survival analysis.
  for (jj in 1:3) {
    
    # First method: Dichotomize data at the median.
    if (jj == 1) {
      ex <- hiloquant(EXPR1[inde1], 0.5)
    }
    
    # Second method: Use the cutpointr package to find an optimal cutpoint based on maximizing the Youden's index.
    if (jj == 2) {
      op <- cutpointr(EXPR1[inde1], EVENT1[inde1], pos_class = 1, neg_class = 0, method = maximize_metric, metric = youden)
      ex <- hiloquant3(EXPR1[inde1], op$optimal_cutpoint)
    }
    
    # Third method: Find the cutpoint that maximizes the rank statistic for survival difference.
    if (jj == 3) {
      rankstat <- NULL
      for (inda in 20:80) {
        qr <- inda/100
        exeo <- hiloquant2(EXPR1[inde1], qr)
        rankstat[inda] <- survdiff(Surv(SURV1[inde1], EVENT1[inde1]) ~ factor(exeo))$chisq
      }
      OP <- which(rankstat == max(rankstat, na.rm = TRUE))[1]
      op <- OP / 100
      ex <- hiloquant2(EXPR1[inde1], op)
    }
    
    # Fit a survival curve for the dichotomized groups.
    DIS.bygene <- survfit(Surv(SURV, EVENT) ~ ex)
    # Perform a log-rank test to compare the survival curves between groups.
    DIS.bygene.logrank <- survdiff(Surv(SURV, EVENT) ~ ex)
    
    # Calculate the p-value from the log-rank test.
    PSURV1 <- (1 - pchisq(DIS.bygene.logrank$chisq, df = 1))
    # Apply a formula to adjust the p-value (likely study-specific and needs context for interpretation).
    PA <- PSURV1 * (-1.63) * (1 + 2.35 * log(PSURV1))
    
    # Calculate the hazard ratio (HR) from the survival curves.
    HR <- round(((DIS.bygene.logrank$obs[1] / DIS.bygene.logrank$exp[1]) / (DIS.bygene.logrank$obs[2] / DIS.bygene.logrank$exp[2])), 2)
    # Calculate the log hazard ratio and derive its confidence interval.
    lHR1 <- (DIS.bygene.logrank$obs[1] - DIS.bygene.logrank$exp[1]) / (DIS.bygene.logrank$var[1, 1])
    HR1 <- exp(1)^lHR1
    V <- abs(DIS.bygene.logrank$var[1, 1])
    upHR1 <- exp(1)^(lHR1 + 1.96 / sqrt(V))
    loHR1 <- exp(1)^(lHR1 - 1.96 / sqrt(V))
    # Concatenate the confidence interval limits into a string.
    CI <- paste(sprintf("%.2f", round(loHR1, 2)), sprintf("%.2f", round(upHR1, 2)), sep = "-")
    
    # Extract the number of individuals at risk in each group from the log-rank test.
    n1 <- as.vector(DIS.bygene.logrank$n[1])
    m1 <- as.vector(DIS.bygene.logrank$n[2])
    
    # Plot the survival curves on the configured PDF device.
    plot(DIS.bygene, col = c(colred, coldark), lwd = 1.5, xlab = "Time (months)", ylab = "Overall survival (probability)", xlim = c(0, mx1), cex.lab = 1.1, cex.axis = 1.1, main = "", axes = FALSE, mark.time = TRUE)
    # Add axis labels and a box around the plot.
    axis(side = 2, las = 1, cex.axis = 1.1, tck = -0.015, lwd = 1.5, mgp = c(6, 0.4, 0))
    axis(side = 1, cex.axis = 1.1, tck = -0.015, lwd = 1.5, mgp = c(3, 0.4, 0))
    box(lwd = 1.5)

    # Annotate the plot with the p-value, adjusted p-value, and hazard ratio.
    if (PSURV1<0.001) {
          pid<-sprintf("%.1E", PSURV1)
          pid<-sub("E","x10",pid)
          pid1<-strsplit(pid, "-")
          
          TXD1<-bquote("logrank p="*.(pid1[[1]][1])^-.(as.numeric(pid1[[1]][2]))*" ("*.(sprintf ("%.3f",PA))*")")
          TXD<-bquote("logrank p="*.(pid1[[1]][1])^-.(as.numeric(pid1[[1]][2])))
          
        } else {
          TXD<-paste("logrank p=",sprintf ("%.3f",PSURV1),sep="")
          TXD1<-paste("logrank p=",sprintf ("%.3f",PSURV1), " (",sprintf ("%.3f",PA),")",sep="")
        }
        if (jj==3) {
          text(40,0.9,TXD1, col="black", cex=1.1, pos=4)
        } else {
          text(40,0.9,TXD, col="black", cex=1.1, pos=4)
          
        }
        text(40,0.8,paste("HR=",sprintf("%.2f",HR)," (",CI,")",sep=""), col="black", cex=1.1, pos=4)
        
    # Annotate the plot with the number of individuals at risk at various time points.
    OSA<-SURV[which(ex==1)]
    OSB<-SURV[which(ex==2)]
    timeseq<-seq(20,mx1,by=20)
    text(-7,-0.24,"No. at risk",xpd=TRUE,col="black",cex=1.1,pos=1)
    text(-15,-0.33,"High:",xpd=TRUE,col=colred,cex=1.1,pos=1)
    text(-16,-0.42,"Low:",xpd=TRUE,col=coldark,cex=1.1,pos=1)
    text(mx1/2,1.20,colnames(D)[i],xpd=TRUE,col="black",cex=1.1,pos=1)
    text(0,-0.33,n1,xpd=TRUE,col=colred,cex=1.1,pos=1)
    text(0,-0.42,m1,xpd=TRUE,col=coldark,cex=1.1,pos=1)
    for (j in timeseq) {
      nr1<-n1-length(which(OSA<j))
      nr2<-m1-length(which(OSB<j))
      text(j,-0.33,nr1,xpd=TRUE,col=colred,cex=1.1,pos=1)
      text(j,-0.42,nr2,xpd=TRUE,col=coldark,cex=1.1,pos=1)	
    }
  }
}

# Close the PDF device, saving the plot.
dev.off() 

# Survival analysis Kaplan Meier curves where BRCAness is dichotomized based on prediction probability.
#************************************
# Set up a PDF device for plotting survival curves and saving the output to a file.
pdf(file="OS_BRCAness.pdf", width=12, height=5, useDingbats = FALSE)

# Configure the plotting parameters for a 1x2 panel layout.
par(mfrow=c(1,2), tcl=0.5, cex=1.1, lwd=1.5, lend="square", ljoin="mitre", las=1, mgp=c(2,0.3,0), mar=c(6,6,3,2), oma=c(1,1,1,1))

# Define colors for plotting different groups in the survival curves.
colred <- c("#DD0000") # Red for one group
coldark <- c("#0000DD") # Blue for another group
colextra <- c("#FFFFFF") # White, unused in this snippet

#************************************

# Set the maximum x-axis limit for the survival plots.
mx1 <- 140

# Loop through the specified column indices in the data frame to perform survival analysis.
for (i in c(125)){
  # Extract the expression levels for the current gene index.
  EXPR1 <- D[[i]]
  # Select the subset of data that is complete (non-missing) for expression, event, and survival.
  inde1 <- which(!is.na(EXPR1) & !is.na(EVENT1) & !is.na(SURV1))
  # Subset the survival time and event data for the complete cases.
  SURV <- SURV1[inde1]
  EVENT <- EVENT1[inde1]
  
  # Apply the SplitBRCAness function to dichotomize the gene expression data based on a cutoff.
  ex <- SplitBRCAness(EXPR1[inde1], 0.526)
  
  # Fit a survival curve using the Kaplan-Meier estimator for the dichotomized expression data.
  DIS.bygene <- survfit(Surv(SURV, EVENT) ~ ex)
  # Perform a log-rank test to assess the difference in survival between the two groups.
  DIS.bygene.logrank <- survdiff(Surv(SURV, EVENT) ~ ex)
  
  # Calculate the p-value from the log-rank test.
  PSURV1 <- (1 - pchisq(DIS.bygene.logrank[[5]], df = 1))
  # Perform an adjustment on the p-value (the formula for adjustment should be provided by study design or domain knowledge).
  PA <- PSURV1 * (-1.63) * (1 + 2.35 * log(PSURV1))
  
  # Calculate the hazard ratio (HR) for the groups and derive the confidence intervals (CI).
  HR <- round(((DIS.bygene.logrank$obs[1] / DIS.bygene.logrank$exp[1]) / (DIS.bygene.logrank$obs[2] / DIS.bygene.logrank$exp[2])), 2)
  lHR1 <- (DIS.bygene.logrank$obs[1] - DIS.bygene.logrank$exp[1]) / (DIS.bygene.logrank$var[1, 1])
  HR1 <- exp(lHR1)
  V <- abs(DIS.bygene.logrank$var[1, 1])
  upHR1 <- exp(lHR1 + 1.96 / sqrt(V))
  loHR1 <- exp(lHR1 - 1.96 / sqrt(V))
  CI <- paste(sprintf("%.2f", round(loHR1, 2)), sprintf("%.2f", round(upHR1, 2)), sep = "-")
  
  # Extract the number of individuals at risk in each group from the log-rank test.
  n1 <- as.vector(DIS.bygene.logrank$n[1])
  m1 <- as.vector(DIS.bygene.logrank$n[2])
  
  # Plot the survival curves with the defined colors and add labels and annotations.
  plot(DIS.bygene, col = c(colred, coldark), lwd = 1.5, xlab = "Time (months)", ylab = "Overall survival (probability)", xlim = c(0, mx1), cex.lab = 1.1, cex.axis = 1.1, main = "", axes = FALSE, mark.time = TRUE)
  axis(side = 2, las = 1, cex.axis = 1.1, tck = -0.015, lwd = 1.5, mgp = c(6, 0.4, 0))
  axis(side = 1, cex.axis = 1.1, tck = -0.015, lwd = 1.5, mgp = c(3, 0.4, 0))
  box(lwd = 1.5)
  
  # Annotate the plot with the log-rank p-value and hazard ratio.
  TXD <- paste("logrank p=", sprintf("%.3f", PSURV1), sep = "")
  text(40, 0.9, TXD, col = "black", cex = 1.1, pos = 4)
  text(40, 0.8, paste("HR=", sprintf("%.2f", HR), " (", CI, ")", sep = ""), col = "black", cex = 1.1, pos = 4)
  
  # Annotate the plot with the number of individuals at risk at various time points below the plot.
  OSA <- SURV[which(ex == 1)]
  OSB <- SURV[which(ex == 2)]
  timeseq <- seq(20, mx1, by = 20)
  text(-7, -0.24, "No. at risk", xpd = TRUE, col = "black", cex = 1.1, pos = 1)
  text(-15, -0.33, "High:", xpd = TRUE, col = colred, cex = 1.1, pos = 1)
  text(-16, -0.42, "Low:", xpd = TRUE, col = coldark, cex = 1.1, pos = 1)
  text(mx1 / 2, 1.20, colnames(D)[i], xpd = TRUE, col = "black", cex = 1.1, pos = 1)
  text(0, -0.33, n1, xpd = TRUE, col = colred, cex = 1.1, pos = 1)
  text(0, -0.42, m1, xpd = TRUE, col = coldark, cex = 1.1, pos = 1)
  for (j in timeseq) {
    nr1 <- n1 - length(which(OSA < j))
    nr2 <- m1 - length(which(OSB < j))
    text(j, -0.33, nr1, xpd = TRUE, col = colred, cex = 1.1, pos = 1)
    text(j, -0.42, nr2, xpd = TRUE, col = coldark, cex = 1.1, pos = 1)  
  }
}

# Close the PDF device, completing the plot and saving the file.
dev.off() 
