import pandas as pd
import numpy as np

#type root directory that the Union CpG ade and ah folders are in
# ex. directory = 'C:\\Users\\tug39776\\OneDrive - Temple University\\wgbs\\'
directory = 'C:\\Users\\tug39776\\OneDrive - Temple University\\wgbs\\'

#enter gene, chromosome number in quotations, and start/end position for promoter
gene = "TNFRSF19_02_"
chrom_num = "13"
start_pos = 23577883
end_pos = 23579883

#set file names and check output at runtime that the files are correct
file_ade = directory+"ade\\chr"+chrom_num+".txt"
file_ah = directory+"ah\\chr"+chrom_num+".txt"
print(file_ade)
print(file_ah)

#read in text file for the chromosome into 2 dataframes
data_ade = pd.read_table(file_ade, header=None, names=["chromosome","position","CpGs","Methylated"])
data_ah = pd.read_table(file_ah, header=None, names=["chromosome","position","CpGs","Methylated"])

#trim the data for just the selected promoter region
trim_ade = data_ade[(data_ade["position"] >= start_pos) & (data_ade["position"] <= end_pos)]
trim_ah = data_ah[(data_ah["position"] >= start_pos) & (data_ah["position"] <= end_pos)]

#set the output files and print the trimmed dataset to 'gene_name_adh/ah.csv' in the root direction
#header : chromosome position num_CpGs num_methylated
file_ade_out = directory+gene+"_ade.csv"
file_ah_out = directory+gene+"_ah.csv"

#write raw data files
trim_ah.to_csv(file_ah_out)
trim_ade.to_csv(file_ade_out)

#make a list of positions for each list
ade_positions = trim_ade["position"].to_numpy()
ah_positions = trim_ah["position"].to_numpy()

#check each line of ah data if ade position exists
for x in range(len(trim_ah)):
    cur_pos = trim_ah["position"].iloc[x]
    
    #if position is in each list then create list of the position, methylation in ade, and methylation in ah
    if cur_pos in ade_positions:
        ade_index = ade_positions.tolist().index(cur_pos)
        ade_meth = trim_ade["Methylated"].iloc[ade_index] / trim_ade["CpGs"].iloc[ade_index]
        ah_meth = trim_ah["Methylated"].iloc[x] / trim_ah["CpGs"].iloc[x]
        
        #attempt to append the list, if doesn't exist yet then initilize the list
        try:
            pos_list = np.append(pos_list,[[cur_pos,ade_meth,ah_meth]],axis=0)
        except:
            pos_list = np.array([[cur_pos,ade_meth,ah_meth]])
      
#convert the numpy array to a dataframe, and convert the position column back to an int
condense_data = pd.DataFrame(pos_list)
condense_data.convert_dtypes()

#write out the condensed dataset so 'gene_anme.xlsx' in the root directory
excel_out = directory+gene+"_chr"+chrom_num+".xlsx"

try:
    condense_data.to_excel(excel_out, index=False, header=["position","Meth_Ade","Meth_Ah"])
except:
    print("output file is open, please close the file")
    quit()
print("processing done")