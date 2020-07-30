'''
any utils I would want to add on pandas can go here.
'''

import pandas as pd

def print_full(df):
    '''
    prints a pandas data frame completely.
    '''
    if not isinstance(df, pd.core.frame.DataFrame):
        raise ValueError('df should be a pandas.core.frame.DataFrame')

    options_print = {
        'display.max_rows': len(df),
        'display.max_columns': None,
        'display.width': 2000,
        'display.float_format':'{:20,.2f}'.format,
        'display.max_colwidth': None
    }
    options_ori = {opt_name:pd.get_option(opt_name) for opt_name in options_print.keys()}

    for opt_name, opt_value in options_print.items():
        pd.set_option(opt_name, opt_value)

    print(df)

    for opt_name, opt_value in options_ori.items():
        pd.set_option(opt_name, opt_value)
