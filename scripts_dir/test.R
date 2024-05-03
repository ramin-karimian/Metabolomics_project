print("Hiiiiii")
# library(pairwiseAdonis)
# install.packages("vegan")
library(vegan)

args = commandArgs(trailingOnly=TRUE)

# c1 = 'A'
# t1 = '0'
# c2 = 'A'
# t2 = '90'

c1 = args[1]
t1 = args[2]
c2 = args[3]
t2 = args[4]
df = read.csv(sprintf("../outputs/pca_%s%s-%s%s/temp/df.csv",c1,t1,c2,t2 ),row.names = 1)
# print(df)

ad <- adonis2(df[,1:(length(df)-2)]~individual * time,df[,(length(df)-1):(length(df))], permutations = 999, method="euclidean") # bray
#ad <- adonis2(df[,1:(length(df)-2)]~ individual * time,list(df[,(length(df))]), permutations = 999, method="euclidean")
print(ad)


