'''
so far just a one-time use script two compare two version of the vcf mapping table.
pandas df comparison bit could be useful later, and could be developed further
'''

from os import path

import pandas as pd

import base
DATA_DIR = path.join(base.ROOT_DIR, "data")

def compare_df(df1, df2):
    '''
    given two pandas dataframes, print out difference in row indexes, and values.
    '''

    intersection_rows = [field for field in df1.index if field in df2.index]
    for field in df1.index:
        if not field in intersection_rows:
            print("Only in df1: ", field)
    for field in df2.index:
        if not field in intersection_rows:
            print("Only in df2: ", field)

    intersection_cols = [field for field in df1.columns if field in df2.columns]
    for field in df1.columns:
        if not field in intersection_cols:
            print("Only in df1: ", field)
    for field in df2.columns:
        if not field in intersection_cols:
            print("Only in df2: ", field)

    df1 = df1[intersection_cols]
    df2 = df2[intersection_cols]

    for field in intersection_rows:
        comp = (df1.loc[field] == df2.loc[field])
        if not all(comp):
            mismatch_cols = list(comp[comp == 0].index)
            print("%s:\t%s" % (field, ", ".join(mismatch_cols)))


            
if __name__ == '__main__':
    fname_old = path.join(DATA_DIR, 'VCF Mapping Table - v0.4.8 variant table.tsv.bak')
    fname_new = path.join(DATA_DIR, 'VCF Mapping Table - v0.4.8 variant table.tsv')
#    fname_old = path.join(DATA_DIR, 'VCF Mapping Table - GeneTable v0.4.6.tsv.bak')
#    fname_new = path.join(DATA_DIR, 'VCF Mapping Table - GeneTable v0.4.6.tsv')

    KEY_FIELD = "field_name"

    df_old = pd.read_csv(fname_old, sep='\t', header=5, na_filter=False)
    df_old.set_index(KEY_FIELD, inplace=True)
    df_old.pop("no")
    df_new = pd.read_csv(fname_new, sep='\t', header=5, na_filter=False)
    df_new.set_index(KEY_FIELD, inplace=True)
    df_new.pop("no")

    compare_df(df_old, df_new)
