# PCA-Clustering-Berry-Dataset


* Dataset : 
Berry South Africa titled, "Blood Trancriptional Profiles of TB in South Africa."

It records the expression profiling by array. It aims to compare transcriptional response to Tuberculosis in regions
of different incidence or prevalence.

* Methods : 
COnducted Principal Component Analysis and Heirarchial Clustering on Berry Data set. 



* MISCELLANEOUS ANALYSIS :
In order to restore the proper ordering of the columns of a given matrix, X, there was analysis done. The top & bottom
30 rows seemed to have a gradient pattern. By sorting the rows according to their mean, we hoped to get an
order. The heatmap after sorting the matrix by their mean indicated a clear image of the picture and were reversed
photos of each other.

An alternative way of the attempt would be using the Principal Component Analysis. From the graph plotted
with the various Principal components (fit graph), we can observe that Principal Component 1 contributes maximum
variance. 

First we play with whole data set, the solution seems to be in a mess. Then try out with only the top 30
rows, the solution seems to be better. Finally we treat 30 rows from both the top and bottom, we get much better
solution.
