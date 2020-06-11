#Downloading the file gene_expression.tsv from the Github repository.
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part1_files/gene_expression.tsv",
              destfile = "try.tsv")

#Read the file gene_expression.tsv Gene accession numbers considered as row names
a <- read.table("try.tsv", header = TRUE, row.names = 1)

#Displaying the values of the first six genes
head(a)

#Finding the mean of given data
a$Mean <- rowMeans(a)

#Displaying the values for first six genes
head(a)

#Listing 10 the genes with highest mean expression
list <- a[order(-a$Mean),]
head(list,10)

#Number of genes with a mean value less than 10.
nrow( subset(a, a$Mean<10))

#histogram plot of the means
a$Mean <- as.matrix(a)
range(a$Mean)
hist(a$Mean)
hist(as.matrix(a$Mean),10, xlab = "Mean", breaks = 50, col = "blue", xlim = c(0,50000))

#download file growth_data.csv
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part1_files/growth_data.csv",
              destfile = "try2.csv")

#Read the file growth_data.csv.
b <- read.table("try2.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE)
colnames(b)

#Subset
b[1:50,]
northeast <- b[1:50, ]
#Calculate mean and standard deviation
mean(northeast$Circumf_2004_cm)
sd(northeast$Circumf_2004_cm)
mean(northeast$Circumf_2019_cm)
sd(northeast$Circumf_2019_cm)

#subset
b[51:100,]
southwest <- b[51:100, ]

#Calculate mean and standard deviation at the beginning and end of study period
mean(southwest$Circumf_2004_cm)
sd(southwest$Circumf_2004_cm)
mean(southwest$Circumf_2019_cm)
sd(southwest$Circumf_2019_cm)

#Displaying box plot
boxplot(southwest$Circumf_2004_cm,southwest$Circumf_2019_cm,northeast$Circumf_2004_cm,northeast$Circumf_2019_cm,
        names = c("SW2004","SW2019","NE2004","NE2019"),ylab="Cirumference (cm)",
        main="Tree Circumference Growth at Two Plantation Sites")

#Mean tree growth over 10 years in both sites
GrowthSW <- (southwest$Circumf_2019_cm-southwest$Circumf_2009_cm)     
GrowthNE <- (northeast$Circumf_2019_cm-northeast$Circumf_2009_cm)
mean(GrowthSW)
mean(GrowthNE)
head(b)

#Findings for t.test
res <- t.test(GrowthSW,GrowthNE, var.equal = FALSE)
res

#Findings for wilcox.test
res <- wilcox.test(GrowthSW,GrowthNE)
res
