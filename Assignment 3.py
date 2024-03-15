# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 18:40:55 2024

@author: venkata pradeep kumar athota
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# a. Load the Data
df1 = pd.read_excel(r"C:\Users\pradeepchowdary\anaconda\Lib\programming\Gene_Expression_Data.xlsx")
df2 = pd.read_csv(r"C:\Users\pradeepchowdary\anaconda\Lib\programming\Gene_Information.csv")
df3 = pd.read_csv(r"C:\Users\pradeepchowdary\anaconda\Lib\programming\Sample_Information.tsv", sep = '\t')


# b. Change Sample Names
# Create a mapping from sample names to phenotypes
mapping = df3['group'].to_dict()
# Rename columns in the gene expression DataFrame
df1_renamed = df1.rename(columns=mapping)


# c. Split Merged Data and Split the DataFrame into 'tumor' and 'normal'
tumor_columns = [col for col in df1_renamed.columns if 'tumor' in col or col == 'Probe_ID']
tumor_df = df1_renamed[tumor_columns]
normal_columns = [ col for col in df1_renamed.columns if 'normal' in col or col == 'Probe_ID']
normal_df = df1_renamed[normal_columns]

# d.Compute the average expression for all genes from the 2 data sets from part d
tumor_avg_expression = tumor_df.mean(axis=1,numeric_only=True)
normal_avg_expression = normal_df.mean(axis=1,numeric_only=True)

# e.Determine the fold change for each Probe between the two groups ((Tumour – Control) / Control)
fold_change = (tumor_avg_expression - normal_avg_expression) / normal_avg_expression

#Use the data from part e and “Gene_Information.csv” to identify all genes fold change magnitude (absolute value) was greater than 5
fold_change_genes = pd.DataFrame({'Probe_ID': df2['Probe_ID'], 'Fold_Change': fold_change})
fold_change_genes = fold_change_genes[abs(fold_change_genes['Fold_Change']) > 5]

#Add a column to the result of part f to include if the gene was higher expressed in “Normal” or “Tumor” samples
fold_change_genes['Higher_Expression'] = np.where(fold_change_genes['Fold_Change'] > 0, 'Tumor', 'Normal')

# 2nd question
#merging the df2 and fold change
fold_change_genes = pd.merge(fold_change_genes, df2[['Probe_ID', 'Chromosome']], on='Probe_ID', how='left')

# a.Perform exploratory data analysis on the genes from part 1g
plt.scatter(range(len(fold_change_genes['Fold_Change'])), fold_change_genes['Fold_Change'])
plt.title('Scatter Plot of Fold Change Values')
plt.xlabel('Fold Change')
plt.show()

# b. Histogram of DEGs by Chromosome
plt.figure(figsize=(10, 6))
plt.hist(fold_change_genes['Chromosome'].astype(str),color = 'yellow')
plt.title('Distribution of DEGs by Chromosome')
plt.xlabel('Chromosome')
plt.ylabel('Frequency')
plt.show()

# c. Histogram of DEGs by Chromosome Segregated by Sample Type
normal_genes = fold_change_genes[fold_change_genes['Higher_Expression'] == 'Normal']
tumor_genes = fold_change_genes[fold_change_genes['Higher_Expression'] == 'Tumor']
plt.figure(figsize=(12, 6))
plt.hist(normal_genes['Chromosome'].astype(str), color='blue', label='Normal', bins=20)
plt.hist(tumor_genes['Chromosome'].astype(str), color='red', label='Tumor', bins=20)
plt.title('Distribution of DEGs by Chromosome Segregated by Sample Type')
plt.xlabel('Chromosome')
plt.ylabel('Frequency')
plt.legend(title='Higher Expression')
plt.show()

# d. Bar Chart of Upregulated and Downregulated DEGs
upregulated_percentage = (fold_change_genes[fold_change_genes['Higher_Expression'] == 'Tumor'].shape[0] / fold_change_genes.shape[0]) * 100
downregulated_percentage = 100 - upregulated_percentage
plt.bar(['Upregulated', 'Downregulated'], [upregulated_percentage, downregulated_percentage], color=['black', 'orange'])
plt.xlabel('Gene Expression')
plt.ylabel('Percentage')
plt.title('Percentage of DEGs Upregulated and Downregulated in Tumor Samples')
plt.show()

# e. Heatmap of Gene Expression by Sample
# drop the probe_id column
dfx = df1_renamed.drop(['Probe_ID'], axis =1)
plt.figure(figsize=(12, 8))
sns.heatmap(dfx.corr(), cmap='viridis', annot=True)
plt.title('Heatmap of Gene Expression by Sample')
plt.xlabel('Sample')
plt.ylabel('Gene')
plt.show()

# f. Clustermap of Gene Expression by Sample
plt.figure(figsize=(12, 8))
sns.clustermap(dfx, cmap='viridis')
plt.title('Clustermap of Gene Expression by Sample')
plt.xlabel('Sample')
plt.ylabel('Gene')
plt.show()

# g .write a few sentences of analysis
#scatter plot of fold change indicates the significant variation in gene expression between tumor and normal samples exhibiting foldcahnges greater than 5.
# Histogram of DEG by chromosome indicates the non uniform distribution of DEG's across the chromosome implicating the tumors genomic regions
# the higher percentage of upregulated genes in tumor samples compared to down regulated suggests shifts towards gene expression levels in tumor samples.
# the heat map and cluster map provide the understanding of the correaltion and clustering  patterns of gene expression. this analysis comprehensively emphasize the importance of understanding of gene expression in biology.
