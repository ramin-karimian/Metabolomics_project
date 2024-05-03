# library(pairwiseAdonis)
library(vegan)

args = commandArgs(trailingOnly=TRUE)

c1 = args[1]
t1 = args[2]
c2 = args[3]
t2 = args[4]
col1 = args[5]
col2 = args[6]
group = args[7]


print(sprintf("R data is: ../outputs/pca_%s%s-%s%s_%s/temp/df.csv",c1,t1,c2,t2, group ))
df = read.csv(sprintf("../outputs/pca_%s%s-%s%s_%s/temp/df.csv",c1,t1,c2,t2, group ),row.names = 1)
# print(df)
if (col2 =='time') {
  ad <- adonis2(df[,1:(length(df)-2)]~ individual * time,df[,(length(df)-1):(length(df))], permutations = 999, method="euclidean") # bray

} else if (col2 =='condition') {
    ad <- adonis2(df[,1:(length(df)-2)]~ individual * condition,df[,(length(df)-1):(length(df))], permutations = 999, method="euclidean") # bray
}

#ad <- adonis2(df[,1:(length(df)-2)]~ individual * time,list(df[,(length(df))]), permutations = 999, method="euclidean")
print(ad)

write.csv(ad,sprintf('../outputs/pca_%s%s-%s%s_%s/Rpermanova_results.csv',c1,t1,c2,t2,group))

