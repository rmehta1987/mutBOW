import pandas as pd
from operator import itemgetter
import numpy as np
import argparse

def createOneDataFrame(path):
    ''' Creates a dataframe from a .maf file contaning mutation data
    Args:
        path: path to the maf file
    Returns:
        mutation_df a dataframe consiting of all the mutations
    '''
    used_col = ['Hugo_Symbol', 'Entrez_Gene_Id', 'Chromosome', 'Start_position',
                'End_position', 'Strand', 'Variant_Classification', 'Variant_Type',
                'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2',
                'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode', 'Transcript_Strand']
    # Note some maf files contain comments at the beginning noting version history, can ignore 
    # comments argument in pandas and passing '#'
    mutation_df=pd.read_csv(path,sep='\t',comment='#',encoding='utf-8',usecols=used_col)
    
    return mutation_df

def createAllDataFrame(path):
    
    used_col = ['Hugo_Symbol', 'Entrez_Gene_Id', 'Chromosome', 'Start_position',
                'End_position', 'Strand', 'Variant_Classification', 'Variant_Type',
                'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2',
                'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode', 'Transcript_Strand']
    mutation_df = pd.concat([pd.read_csv(os.path.join(result_dir, maf_file),
                                    sep='\t', 
                                    na_values='[Not Available]') 
                    for maf_file in maf_files])

def sortPatients(the_df):
    ''' Sorts the dataframe by patient respective of their mutations
    See: https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/ for column names

    Args:
        the_df - dataframe containing mutation data of all patients
    Returns:
        dataframe sorted by patients with mutations
    '''
    
    return the_df.sort_values(by=['Tumor_Sample_Barcode'], axis=0, ascending=True)

def createFiles(the_df):
    ''' Creates a csv file, where row is patient, and columns are the mutations
    then writes 5 files: 
    a list of all the mutations, an array containing all the unique mutations (the vocab),
    a list of numpy arrays, where columns represent the unique mutation index and each numpy array is the patient 
    a list of numpy arrays, where the columns represent the count of the associated mutation for each patient
    a CSV file with a matrix representing Patients x Mutations
    a csv file of the dataframe
    [Note this is a bag of words format for Patient Mutation Data]
    [Note we ignore variant of mutation, Misense, deletion, insertion as treated as the same]
    Args:
        the_df: a sorted dataframe by patient name, containing their mutation data
    '''
    temp_muts = the_df['Hugo_Symbol'].unique()
    temp_ids = the_df['Tumor_Sample_Barcode'].unique()
    mut_array = np.array(temp_muts) # Vocab

    # Create a dictionary such that the key is the mutation and the value is the index
    dict_mut = {x: i for i, x in enumerate(mut_array)}
    newpd = pd.DataFrame(columns=temp_muts, index=temp_ids)
    newpd.fillna(0,inplace=True)
    pat_idx = []  # list of numpy arrays, where columns represent the unique mutation index and each numpy array is the patient 
    pat_cnt = []  # list of numpy arrays, where the columns represent the count of the associated mutation for each patient

    for ind_id, id in enumerate(temp_ids):

        # Get the Count of Mutations for a specific patient
        filtered_counts = the_df[the_df['Tumor_Sample_Barcode']==id]['Hugo_Symbol'].value_counts()

        # See: https://stackoverflow.com/questions/18453566/python-dictionary-get-list-of-values-for-list-of-keys
        # Create a list of the indicies returned from filtered_counts
        pat_idx.append(np.array(itemgetter(*filtered_counts.index.values)(dict_mut)))

        # Create a list of the counts of the associate mutations for each patient
        pat_cnt.append(np.array(filtered_counts.values))

        # Create a dataframe of patients x mutations
        newpd.iloc[ind_id][filtered_counts.index.values] = filtered_counts.values

    # Is a stack like save, so the 2nd argument is the 1st argument when loading
    np.savez('patient_bow.npz', pat_idx, pat_cnt)
    np.savez('patient_vocab.npz', mut_array)

    the_df.to_csv('reduced_pat_mut.csv', sep=',',index=False)
    newpd.to_csv('pat_mut_mat.csv',sep=',')

    return newpd


"""Parse arguments."""
parser = argparse.ArgumentParser(
    description='Population random measure embedding model.')

# data configurations
parser.add_argument(
'--dataset',
default='broad.mit.edu_LUAD_IlluminaGA_DNASeq.maf',
type=str,
help='dataset folder or maf file name')

parser.add_argument(
'--multiple',
default=0,
type=int,
help='If there are multiple MAF files, 0 for No')

args = parser.parse_args()

if args.multiple:
    the_df = createAllDataFrame(args.dataset)
else:
    the_df = createOneDataFrame(args.dataset)

the_df = sortPatients(the_df)
createFiles(the_df)











        















    

    
    






