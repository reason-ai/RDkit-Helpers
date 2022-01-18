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
    common tasks on a pandas dataframe. 
    
    Attributes
    ----------
    descriptors_factory: DescriptorsFactory
        a factory to get common 2d molecular descriptors and add their
        respective columns to the dataframe
        
    Returns
    ----------
    df: A modified dataframe with added columns. Dataframe must be assigned
    to a variable when calling any method
    '''
    def __init__(self):
        self.descriptors_factory = DescriptorsFactory()
    
    def add_descriptors_to_frame(self,df,descriptors=None,default=True):
        return self.descriptors_factory.get_descriptors(df,descriptors,
                                                      default)
        

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
            List of descriptors strings to add to dataframe.
        default : bool
            Switch to determine whether default descriptors are returned.

        Returns
        -------
        df : Pandas DataFrame
            DataFrame with added columns.

        '''
        df = df.copy()
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
datafile = 'some csv file with a smiles columns'
df = pd.read_csv(datafile)
PandasTools.AddMoleculeColumnToFrame(df,'smiles')

desc_list = ['amw','n_ve']

put = PandasUtils()
new_df = put.add_descriptors_to_frame(df,desc_list)
other_df = put.add_descriptors_to_frame(df) #no arguments so default values returned



