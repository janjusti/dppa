### Deep Protein Polarity Analyser (DPPA) v0.1
### Copyright (c) 2019 Jan Marans Agnella Justi & Mariana de França Costa
import logging
import argparse
import pkg_resources
import pandas as pd
from tqdm import tqdm
from .src.helpers import auxfunc
from .src.conv import export_dfs

def main():
    # get options from user
    parser = argparse.ArgumentParser(description='Analyse all protein alignment .fasta files from a target.')
    parser.add_argument(
        '-t', '--target', help='Target .fasta file to be analysed.', required=True
    )
    parser.add_argument(
        '-r', '--report', help='Output report file type.', required=True, 
        choices=['csv', 'xls', 'all']
    )
    parser.add_argument(
        '--debug', help='Turn debug messages on.', action='store_true'
    )
    user_args = vars(parser.parse_args())
    # config logging
    set_debug_mode(user_args['debug'])
    # run and export results
    run_and_export(user_args)

def run(target_fn):
    # run analyser
    logging.debug(f'Starting analysis for {target_fn}')
    # create handler of custom functions
    auxf_handler = auxfunc.AuxFuncPack()
    # create list from .fasta target file
    target_list = auxf_handler.fasta_to_list(target_fn)
    # create dfs
    df_pol_results = pd.DataFrame(columns=['ColNum', 'PossibleAminos', 'PossiblePols', 'PolScore'])
    df_alert_results = pd.DataFrame(columns=['SeqName', 'ColNum', 'AlertType'])
    df_pols = pd.read_csv(pkg_resources.resource_filename(__name__, 'src/conv/pols.csv'))
    # get list of unknown and known aminos
    [unknown_aminos, known_aminos] = auxf_handler.get_lists_aminos(df_pols)
    # execute deep searcher
    number_columns = len(target_list[0][1]) # number of chars on sequence
    # on each column from .fasta target file
    isDebugModeNotActive = not logging.getLogger().isEnabledFor(logging.DEBUG)
    for current_col in tqdm(range(0, number_columns), disable=isDebugModeNotActive):
        root, df_alert_results = auxf_handler.deep_searcher(
            target_list, current_col, df_alert_results, unknown_aminos
        )
        # get root list of elements
        amino_leaves = []
        for leaf in root.leaves:
            if leaf.amino in known_aminos or leaf.amino in unknown_aminos:
                amino_leaves.append(leaf.amino)
            else:
                # invalid char, increment on df_alert_results
                df_alert_results = auxf_handler.add_alert_on_df(
                    df_alert_results,
                    leaf.name,
                    current_col,
                    'Unlisted Char ' + leaf.amino
                )
        # remove still unknown aminos
        amino_leaves = [x for x in amino_leaves if x not in unknown_aminos]
        # unify list
        unified_list = list(set(amino_leaves))
        # check length
        if len(unified_list) > 1:
            # get aminos percentages dict
            aminos_dict = auxf_handler.get_aminos_perc_dict(root, unified_list, unknown_aminos)
            # get polarities percentages dict
            pols_dict = auxf_handler.get_polarities_perc_dict(aminos_dict, df_pols)
            # get polarity score
            curr_pol_score = auxf_handler.get_pol_score(pols_dict, df_pols)
            # increment on df_pol_results
            df_pol_results = df_pol_results.append(
                {
                    'ColNum': current_col+1,
                    'PossibleAminos': aminos_dict,
                    'PossiblePols': pols_dict,
                    'PolScore': curr_pol_score
                }, ignore_index=True)
    return df_pol_results, df_alert_results

def export(target_name, report_type, df_pol_results, df_alert_results):
    # export dfs to csv
    dfexp_handler = export_dfs.DfExporter()
    dfexp_handler.export_dfs(
        target_name, report_type, df_pol_results, df_alert_results
    )

def run_and_export(args):
    # run
    [df_pol_results, df_alert_results] = run(args['target'])
    # export dfs to csv
    export(
        args['target'], args['report'], df_pol_results, df_alert_results
    )
    logging.debug('Done.')
    # end

def set_debug_mode(isActive):
    curr_level = logging.DEBUG if isActive else logging.INFO
    logging.basicConfig(level=curr_level, format='%(message)s')