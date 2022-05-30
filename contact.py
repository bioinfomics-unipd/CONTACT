#!/usr/bin/env python
# coding: utf-8




import pandas as pd
import numpy as np
import sys
import seaborn as sns
import matplotlib.pyplot as plt






# check command line, if not specified use 5 as minimum number of genomes

if len(sys.argv)!=6:
    if len(sys.argv)==5 and type(sys.argv[-1])==str:
        sys.argv.append(5)
    else:
        print('Syntax: contact.py checkM.txt classification.txt ncbi_taxonomy.txt tax_level min_number_genomes')
        exit()





# load NCBI database, with no header
ncbi_data=pd.read_table(sys.argv[3],low_memory=False, header=None)





# convert column 6 and 7 (genome size and GC content) into float type, replace special character '-' with nan if present

for col in list(ncbi_data.columns)[6:8]:
    if ncbi_data[col].dtype!='float64':
        try:
            ncbi_data[col].astype(float)
        except Exception:
            ncbi_data[col]=ncbi_data[col].replace('-',np.nan)
        ncbi_data[col]=ncbi_data[col].astype(float)





# load checkM output table, header is line 0 
checkM=pd.read_table(sys.argv[1], header=0, engine='python').set_index('Bin Id')





# load conversion table with classification from gtdb to ncbi db
classification=pd.read_table(sys.argv[2], header=0).set_index('user_genome')





# extract MAGs from classification table assigned at least at user-defined taxonomic level

l=[]
for i in range(len(classification)):
    val=False
    if sys.argv[4]=='family':
        val='f__;' in classification['NCBI classification'][i]
        if val==False:
            l.append(i)
    elif sys.argv[4]=='genus':
        val='g__;' in classification['NCBI classification'][i]
        if val==False:
            l.append(i)
    elif sys.argv[4]=='species':
        val=classification['NCBI classification'][i].endswith('s__')
        if val==False:
            l.append(i)
        

tax_assigned=classification.iloc[l,1]





# create dictionary with each MAG as key and taxonomy as value (user-defined level)

dict_taxonomy={}
for i in range(len((tax_assigned))):
    if sys.argv[4]=='family':
        dict_taxonomy[tax_assigned.index[i]]=tax_assigned[i].split(';')[4]
    elif sys.argv[4]=='genus':
        dict_taxonomy[tax_assigned.index[i]]=tax_assigned[i].split(';')[5]
    elif sys.argv[4]=='species':
        dict_taxonomy[tax_assigned.index[i]]=tax_assigned[i].split(';')[6]





# define function to extract from NCBI db organisms matching taxonomy of each MAG (considering user_defined taxonomic level)

def extract_org_name(org_name):
    mask=ncbi_data[23].str.contains(org_name)
    return ncbi_data[mask]

# define function to calculate statistics for genome size and GC content (mean and standard deviation) and number of extracted organisms, return a tuple of 5 elements 

def stats(extracted):
    numb_genomes=extracted.shape[0]
    avg_genome_size=np.mean(extracted[6])
    avg_GC_content=np.mean(extracted[7])
    std_genome_size=np.std(extracted[6])
    std_GC_content=np.std(extracted[7])
    return (numb_genomes,avg_genome_size,avg_GC_content,std_genome_size,std_GC_content)





# create dictionary of tuples: each MAG is a key associated to a tuple calculated by stat function
# statistics are calculated only if the number of extracted organisms from db is higher than the user-defined minimum number (default = 5)

dict_stat={}
new='.'
for k in dict_taxonomy.keys():
    tax=dict_taxonomy.get(k)
    name=tax.split('_')[2]
    if name!=new and name!='':
        ex=extract_org_name(name)
        if ex.shape[0]>=int(sys.argv[5]):
            s=stats(ex)
            dict_stat[tax]=s
            new=name





# prepare a new dataframe 
output_table=pd.DataFrame(columns=['ID','taxonomy','genome_size_STD_from_avg','expected_missing_exceeding_portion','completeness','contamination','genome_size','expected_genome_size','avg_genome_size_NCBI','std_genome_size_NCBI','GC_content','GC_content_STD_from_avg','avg_GC_content_NCBI','std_GC_content_NCBI','number_of_genomes','size_diff','GC_diff','genome_size_STD_from_avg_>1','GC_content_STD_from_avg_>1']).set_index('ID')





# fill the columns of the dataframe with taxonomy, GC content, genome size, completeness, contamination for each MAG
# new columns are created containing number of genomes, average and standard deviation for genome size and GC content from stat function

for k in dict_taxonomy.keys():
    tax=dict_taxonomy.get(k)
    if tax in dict_stat.keys():
        output_table.loc[k,'taxonomy']=tax
        output_table.loc[k,'GC_content']=checkM.loc[k,'GC']
        output_table.loc[k,'genome_size']=checkM.loc[k,'Genome size (bp)']/1000000
        output_table.loc[k,'completeness']=checkM.loc[k,'Completeness']
        output_table.loc[k,'contamination']=checkM.loc[k,'Contamination'] 
        output_table.loc[k,'avg_genome_size_NCBI']=dict_stat[tax][1]
        output_table.loc[k,'std_genome_size_NCBI']=dict_stat[tax][3]
        output_table.loc[k,'avg_GC_content_NCBI']=dict_stat[tax][2]
        output_table.loc[k,'std_GC_content_NCBI']=dict_stat[tax][4]
        output_table.loc[k,'number_of_genomes']=dict_stat[tax][0]





# other columns are the the difference in genome size and GC content and the number of STDs each MAG is far from the average
# boolean columns to easily detect those that are more distant than 1 std
# the expected genome size is calculated by subtracting the contaminated portion from the genome size and dividing by the completeness
# the expected missing or exceeding portion is calculated subtracting the expected genome size by the average from NCBI

output_table['size_diff']=(output_table['genome_size']>=output_table['avg_genome_size_NCBI']).replace(to_replace=[True,False],value=['larger','smaller'])
output_table['GC_diff']=(output_table['GC_content']>=output_table['avg_GC_content_NCBI']).replace(to_replace=[True,False],value=['higher','lower'])
output_table['genome_size_STD_from_avg']=(output_table['genome_size']-output_table['avg_genome_size_NCBI'])/output_table['std_genome_size_NCBI']
output_table['GC_content_STD_from_avg']=(output_table['GC_content']-output_table['avg_GC_content_NCBI'])/output_table['std_GC_content_NCBI']
output_table['genome_size_STD_from_avg_>1']=abs(output_table['genome_size_STD_from_avg'])>=1
output_table['GC_content_STD_from_avg_>1']=abs(output_table['GC_content_STD_from_avg'])>=1
output_table['expected_genome_size']=(output_table['genome_size']-(output_table['genome_size']*(output_table['contamination']/100)))/(output_table['completeness']/100)
output_table['expected_missing_exceeding_portion']=output_table['expected_genome_size']-output_table['avg_genome_size_NCBI']





# output table is created in the working directory in csv format named with the taxonomic level and the minimum number of genomes considered when performing analyses
output_table.to_csv('output_table_'+str(sys.argv[4])+'_'+str(sys.argv[5])+'.csv')





# a scatterplot is created in the working directory showing the completeness of each genome as function of the number of std its size is distant from the average
# points are colored based on the contamination level
# the file is saved as a pdf and named with the taxonomic level and the minimum number of genomes used

plt.figure()
sns.scatterplot(x='completeness',y='genome_size_STD_from_avg',data=output_table, hue='contamination').axhline(y=0)
plt.savefig('scatterplot_'+str(sys.argv[4])+'_'+str(sys.argv[5])+'.pdf')





# another scatterplot is showing the completeness as function of the expected missing or exceeding genomic portion

plt.figure()
sns.scatterplot(x='expected_missing_exceeding_portion',y='completeness',data=output_table).axvline(x=0, color='r')
plt.savefig('scatterplot2_'+str(sys.argv[4])+'_'+str(sys.argv[5])+'.pdf')







