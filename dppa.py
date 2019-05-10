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
            curr_pol_list = auxf_handler.get_polarities_list(unified_list)
            # get aminos percentages dict
            aminos_dict = auxf_handler.get_aminos_perc_dict(root, unified_list)
            # increment on df_pol_results
            df_pol_results = df_pol_results.append(
                {
                    'ColNum': current_col+1,
                    'PossibleAminos': aminos_dict,
                    'PossiblePols': curr_pol_list
                }
            , ignore_index=True)

    # export dfs to csv
    auxf_handler.export_dfs(
        args['input'], args['target'], df_alert_results, df_pol_results
    )        

    print('Done.')
    # end


if __name__ == '__main__':
    # get options from user
    parser = argparse.ArgumentParser(description='Analyse all protein alignment .fasta files.')
    parser.add_argument('-i', '--input', help='Input folder with .fasta files.', required=True)
    parser.add_argument('-t', '--target', help='Target .fasta file to be analysed.', required=True)
    user_args = vars(parser.parse_args())
    main(user_args)