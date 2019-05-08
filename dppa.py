### Deep Protein Polarity Analyser (DPPA) v0.1
### Copyright (c) 2019 Jan Marans Agnella Justi & Mariana de França Costa
import argparse
import pandas as pd
import numpy as np
from tqdm import tqdm
from pathlib import Path
from src.helpers.auxfunc import AuxFuncPack

def main(args):

    print('Starting analysis for', args['target'], '@', args['input'])
    # create handler of custom functions
    auxf_handler = AuxFuncPack()
    # create list from .fasta target file
    target_list = auxf_handler.fasta_to_list(args['input'], args['target'])
    # create dfs
    df_pol_results = pd.DataFrame(columns=['ColNum','PossibleAminos','PossiblePols'])
    df_alert_results = pd.DataFrame(columns=['SeqName','ColNum','AlertType'])
    # execute deep searcher 
    number_columns = len(target_list[0][1]) # number of chars on sequence
    # on each column from .fasta target file
    for current_col in tqdm(range(0, number_columns)):
        root, df_alert_results = auxf_handler.deep_searcher(
            args['input'], target_list, current_col, df_alert_results
        )
        # get root list of elements
        amino_leaves = []
        for leaf in root.leaves:
            amino_leaves.append(leaf.amino)
        # remove unlisted chars
        blacklisted_str = ['-', 'X', 'B', 'Z', '?']
        amino_leaves = [x for x in amino_leaves if x not in blacklisted_str]
        # unify list
        unified_list = list(set(amino_leaves))
        # check length
        if len(unified_list) > 1:
            # get list of polarities
            filepath = Path.cwd() / 'src' / 'conv' / 'pols.csv'
            df_pols = pd.read_csv(filepath)
            curr_pol_list = []
            for curr_pol in range(0, df_pols.shape[0]):
                pol_amino_list = list(df_pols.Amino.at[curr_pol])
                for curr_amino in unified_list:
                    if curr_amino in pol_amino_list:
                        curr_pol_list.append(df_pols.Type.at[curr_pol])
            # unify polarity list
            curr_pol_list = list(set(curr_pol_list))
            # increment on df_pol_results
            df_pol_results = df_pol_results.append(
                {
                    'ColNum': current_col+1,
                    'PossibleAminos': unified_list,
                    'PossiblePols': curr_pol_list
                }
            , ignore_index=True)

    # export dfs
    filepath = Path.cwd() / 'batch' / args['input'] / args['target']
    alerts_filepath = str(filepath).replace('.fasta', '') + '-alerts.csv'
    pols_filepath = str(filepath).replace('.fasta', '') + '-pols.csv'
    print('Exporting alerts...')
    df_alert_results.to_csv(alerts_filepath, index=False)
    print('Exporting polarity results...')
    df_pol_results.to_csv(pols_filepath, index=False)     

    print('Done.')
    # end


if __name__ == '__main__':
    # get options from user
    parser = argparse.ArgumentParser(description='Analyse all protein alignment .fasta files.')
    parser.add_argument('-i', '--input', help='Input folder with .fasta files.', required=True)
    parser.add_argument('-t', '--target', help='Target .fasta file to be analysed.', required=True)
    user_args = vars(parser.parse_args())
    main(user_args)