'''
so far just a one-time use script two compare two version of the vcf mapping table.
pandas df comparison bit could be useful later, and could be developed further
'''

from os import path

import pandas as pd

DATA_DIR = path.join(base.ROOT_DIR, "data")

def compare_df(df1, df2):
    '''
    given two pandas dataframes, print out difference in row indexes, and values.
    '''
    intersection = [field for field in df1.index if field in df2.index]
    for field in df1.index:
        if not field in intersection:
            print("Only in df1: ", field)
    for field in df2.index:
        if not field in intersection:
            print("Only in df2: ", field)

    for field in intersection:
        comp = (df1.loc[field] == df2.loc[field])
        if not all(comp):
            mismatch_cols = list(comp[comp == 0].index)
            print("%s:\t%s" % (field, ", ".join(mismatch_cols)))

if __name__ == '__main__':
    fname_old = path.join(DATA_DIR, 'VCF Mapping Table - v0.4 variant table - stable.tsv')
    fname_new = path.join(DATA_DIR, 'VCF Mapping Table - v0.4.7 variant table - stable.tsv')
    KEY_FIELD = "field_name"

    df_old = pd.read_csv(fname_old, sep='\t', header=5, na_filter=False)
    df_old.set_index(KEY_FIELD, inplace=True)
    df_old.pop("no")
    df_new = pd.read_csv(fname_new, sep='\t', header=5, na_filter=False)
    df_new.set_index(KEY_FIELD, inplace=True)
    df_new.pop("no")

    compare_df(df_old, df_new)
