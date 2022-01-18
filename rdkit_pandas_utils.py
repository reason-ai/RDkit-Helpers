# -*- coding: utf-8 -*-
"""
Created on Mon Jan 17 13:00:28 2022

@author: Reilly
"""
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, PandasTools

class PandasUtils():
    '''A utility class for RDkit to make it quicker and easier to perform
    common tasks on a pandas dataframe. A pandas dataframe must be given as
    input when initializing the class
    
    Attributes
    ----------
    descriptors_factory: DescriptorsFactory
        a factory to get common 2d molecular descriptors and add their
        respective columns to the dataframe
    '''
    def __init__(self,df):
        self.df = df.copy()
        self.descriptors_factory = DescriptorsFactory()
    
    def add_descriptors_to_frame(self,descriptors=None,default=True):
        return self.descriptors_factory.get_descriptors(self.df,
                                                        descriptors,default)
    

class DescriptorsFactory():
    '''
    Options:
        'emw': Descriptors.ExactMolWt,
        'mw': Descriptors.MolWt,
        'n_ha': Descriptors.NumHAcceptors,
        'n_hd': Descriptors.NumHDonors,
        'mlp': Descriptors.MolLogP,
        'n_rb': Descriptors.NumRotatableBonds,
        'hac': Descriptors.HeavyAtomCount,
        'n_ve': Descriptors.NumValenceElectrons}
    '''
    funcs = {'emw': Descriptors.ExactMolWt,
             'amw': Descriptors.MolWt,
             'n_ha': Descriptors.NumHAcceptors,
             'n_hd': Descriptors.NumHDonors,
             'mlp': Descriptors.MolLogP,
             'n_rb': Descriptors.NumRotatableBonds,
             'hac': Descriptors.HeavyAtomCount,
             'n_ve': Descriptors.NumValenceElectrons}
    
    def __init__(self):
        pass
        
    def get_descriptors(self,df,args,default):
        '''
        Get desired descriptors from list and add them to frame
        ExactMolWt, NumHAcceptors, NumHDonors, MolLogp are most commonly
        obtained descriptors so they are set as default arguments if no list
        is passed in.

        Parameters
        ----------
        df : Pandas DataFrame
            DataFrame that columns are being added to.
        args : list
            List of descriptors to add to dataframe.
        default : bool
            Switch to determine whether default descriptors are returned.

        Returns
        -------
        df : Pandas DataFrame
            DataFrame with added columns.

        '''
        if args: default = False
        if default: args = ['emw','n_ha','n_hd','mlp'] 
        else: args = args
        
        if not self.check_args(args):
            print('No changes made')
            return df
        
        for arg in args:
            df[arg] = df['ROMol'].apply(
                DescriptorsFactory.funcs[arg])
        return df
    
    def check_args(self,args):
        for arg in args:
            if arg not in DescriptorsFactory.funcs:
                print(f'WARNING! Invalid keyword argument passed [{arg}]')
                return False
        return True


# Example usage          
datafile = 'some csv file with a smiles col'
df = pd.read_csv(datafile,index_col=0)
PandasTools.AddMoleculeColumnToFrame(df,'smiles')

desc_list = ['amw','n_ha','n_hd','mlp']

put = PandasUtils(df)
new_df = put.add_descriptors_to_frame(desc_list)



