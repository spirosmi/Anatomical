##############################################################################
FAMD CODE
##############################################################################
#Modified from Kassambara, Alboukadel. Practical guide to principal component methods in R: PCA, M (CA), FAMD, MFA, HCPC, factoextra. Vol. 2. Sthda, 2017.

library(FactoMineR) #FDMA
library(factoextra) #FDMA
library(missMDA) #impute data

set.seed(1234)
getwd()
setwd("")

#Data
data<-read.csv("data.csv",head=T)
rownames(data) <- c("ANAT-001", "ANAT-002", "ANAT-003", "ANAT-004", "ANAT-005", "ANAT-006", "ANAT-007", "ANAT-008", "ANAT-009", "ANAT-010", "ANAT-011", "ANAT-012", "ANAT-013", "ANAT-014", "ANAT-019", "ANAT-020")
df <- data[1:16,2:65]
head(df[, 1:7], 4)
str(df)
df[4:64] <- lapply(df[4:64], as.numeric)
str(df)

### with missing values
data.impute <- imputeFAMD(df, ncp=5) 

#FAMD
res.famd <- FAMD(df,tab.disj=data.impute$tab.disj, ncp=64) 
print(res.famd)

#Eigenvalues
eig.val<- get_eigenvalue(res.famd)
head(eig.val)
eig.val
fviz_screeplot(res.famd)

res.famd$var

#FAMD results for variables
var <- get_famd_var(res.famd)

# Coordinates of variables
head(var$coord)

# Cos2: quality of representation on the factor map
head(var$cos2)

# Contributions to the  dimensions
head(var$contrib)

res.famd$quanti.var

M<-as.data.frame(res.famd$var$contrib)

# Plot of variables
plot<- fviz_famd_var(res.famd, col.var="cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, title="Variable Associations")
plot

# Contribution to the first dimension
plot1<-fviz_contrib(res.famd, "var", axes = 1, top=35, color="#00AFBB", fill="#00AFBB", title="Contribution of Variables to Dimension 1")
plot1

# Contribution to the second dimension
plot2<-fviz_contrib(res.famd, "var", axes = 2, top=25, color="#E7B800", fill="#E7B800", title="Contribution of Variables to Dimension 2")
plot2

#Extract the results for quantitative variables
quanti.var <- get_famd_var(res.famd, "quanti.var")
quanti.var
fviz_famd_var(res.famd, "quanti.var", repel = TRUE,
              col.var = "black")

plot3<-fviz_famd_var(res.famd, "quanti.var", col.var = "contrib", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              repel = TRUE, title="Quantitative Variables Correlations")
plot3

# Color by cos2 values: quality on the factor map
fviz_famd_var(res.famd, "quanti.var", col.var = "cos2",
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
              repel = TRUE)

#Extract the results for qualitative variables
quali.var <- get_famd_var(res.famd, "quali.var")
quali.var 
fviz_famd_var(res.famd, "quali.var", col.var = "contrib", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)

#Graph of Individuals
ind <- get_famd_ind(res.famd)
ind
plot4<-fviz_famd_ind(res.famd, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, title="Individuals & Qualitative Biplot")
plot4 

