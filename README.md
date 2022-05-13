# BF528-project-5
Repo containing code written for the final project of the BF528 course.  
Roles performed: Project 3- Programmer and Biologist.  
The scripts ```featureCounts.qsub``` and ```create_counts.R``` are for creating the counts matrices. ```deseq.R``` was used for performing differential gene expression analysis. ```heatmap.R``` was used to plot the heatmaps for the Biologist role.  
How to run-
1. To run the ```featureCounts.qsub``` you need to give two command line inputs- the first for the name of the output file and the second for the name of the input file.  
2. To run ```create_counts.R``` you need to import the results from the featureCounts.  
3. To run ```deseq.R``` you need the counts produced by the ```create_counts.R``` script.  It also requires the associated metadata file.
4. To run ```heatmaps.R``` you need the normalized counts matrixes produced by the ```deseq.R``` script. 
