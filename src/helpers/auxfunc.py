from pathlib import Path
from Bio.SeqIO import FastaIO
from anytree import AnyNode, LevelOrderGroupIter
from StyleFrame import StyleFrame, Styler, utils
import pandas as pd
import string

class AuxFuncPack:
    # package of custom functions

    def fasta_to_list(self, fastas_folder, fasta_fn):
        ### convert .fasta file into a list
        # list format: [(name, seq)]
        filepath = Path.cwd() / 'batch' / fastas_folder / fasta_fn
        fasta_list = list(FastaIO.SimpleFastaParser(open(filepath)))

        return fasta_list

    def add_alert_on_df(self, df_alert, seq_name, col_num, alert_type_str):
        df_alert = df_alert.append(
            {
                'SeqName': seq_name.replace('_', ' '),
                'ColNum': col_num+1,
                'AlertType': alert_type_str
            }
        , ignore_index=True)

        return df_alert
    
    def get_lists_aminos(self, df_pols):
        ### get list of unknown aminos, based on df_pols
        # get alphabet list
        unknown_aminos = list(string.ascii_uppercase)
        # add '?' and '-' char
        unknown_aminos.extend(['?', '-'])
        # start known aminos list
        known_aminos = []
        # remove all pol-defined letters
        for curr_pol in range(0, df_pols.shape[0]):
            list_curr_pol = list(df_pols.Amino.at[curr_pol])
            known_aminos.extend(list_curr_pol)
            unknown_aminos = [x for x in unknown_aminos if x not in list_curr_pol]

        return unknown_aminos, known_aminos

    def deep_searcher(self, fastas_folder, fasta_list, col_num, df_alert, unknown_aminos):
        # create local root
        loc_root = AnyNode(name='LocRoot', amino='LocRootAmino')
        # create nodes from column number
        for row in range(0, len(fasta_list)):
            AnyNode(
                name=fasta_list[row][0], 
                amino=fasta_list[row][1][col_num], 
                parent=loc_root
            )
        # run through nodes
        for node in loc_root.children:
            # check if amino is possibly deepable
            if node.amino in unknown_aminos:
                if 'consensus_sequence' in node.name and node.amino not in ['-', '?']: 
                    # non-gap deepable char
                    fasta_dict = dict(fasta_list)
                    # check previous number of gaps from sequence
                    num_gaps = fasta_dict[node.name][:col_num].count('-')
                    # shift column number to the correct one
                    shifted_col_num = col_num-num_gaps
                    # get fasta_list of this node's original sequence
                    deeper_fn = node.name.replace('_', ' ') 
                    deeper_fn = deeper_fn.replace(' consensus sequence','') + '.fasta'
                    deeper_list = self.fasta_to_list(fastas_folder, deeper_fn)
                    # go deeper
                    deeper_result = self.deep_searcher(fastas_folder, deeper_list, shifted_col_num, df_alert, unknown_aminos)
                    deeper_result[0].parent = node
                    df_alert = deeper_result[1]
                else:
                    # non-deepable char, increment on df_ale
                    df_alert = self.add_alert_on_df(
                        df_alert, 
                        node.name,
                        col_num, 
                        'Pure ' + node.amino
                    )

        return loc_root, df_alert

    def get_aminos_perc_dict(self, root, aminos_list, unknown_aminos):
        ### convert list of aminos to dict (amino -> its perc)
        # get list of aminos levels
        aminos_levels = [[node.amino for node in children] for children in LevelOrderGroupIter(
            root,
            filter_=lambda n: n.amino != 'LocRootAmino'
        )]
        # clear any empty lists (parent roots)
        aminos_levels = [x for x in aminos_levels if x]
        # run through levels
        deepables = [*unknown_aminos]
        deepables.remove('?') # '?' is not deepable
        deepables.remove('-') # '-' is not deepable
        aminos_dict = {}
        for target_amino in aminos_list:
            total_perc = 0
            deepness_perc = 1
            for level in aminos_levels:
                amino_perc_on_lvl = (level.count(target_amino)/len(level))
                total_perc = total_perc + (amino_perc_on_lvl)*deepness_perc
                num_of_deepables = sum(level.count(x) for x in deepables)
                deepables_perc_on_lvl = (num_of_deepables/len(level))
                deepness_perc = deepness_perc * deepables_perc_on_lvl
            aminos_dict[target_amino] = round(100*total_perc,3)

        return aminos_dict
    
    def get_polarities_perc_dict(self, aminos_dict, df_pols):
        ### get dict of polarities within aminos list
        # get list of tuples (pol_type -> perc)
        list_percs = []
        for curr_pol in range(0, df_pols.shape[0]):
            pol_amino_list = list(df_pols.Amino.at[curr_pol])
            for curr_amino in list(aminos_dict.items()):
                if curr_amino[0] in pol_amino_list:
                    list_percs.append((df_pols.Type.at[curr_pol], curr_amino[1]))
        # convert into unified dict
        pols_dict = {}
        for key, value in list_percs:
            pols_dict[key] = round(pols_dict.get(key, 0) + value, 2)

        return pols_dict

    def get_pol_score(self, pols_dict, df_pols):
        ### calculate polarity score
        # convert dict into list of tuples
        pols_percs = list(pols_dict.items())
        # get sum of current scores
        sum_scores = 0
        for curr_pol in pols_percs:
            sum_scores += df_pols.loc[df_pols['Type'] == curr_pol[0], 'Score'].values[0]
        # calculate
        pol_list_size = len(pols_percs)
        if pol_list_size is 1:
            pol_score = 0
        else:
            # get minimal value of perc
            min_pols = min(pols_percs, key = lambda t: t[1])[1]
            pol_score = round(sum_scores*3 + pol_list_size*0.6 + min_pols*0.01, 5)

        return pol_score

    def export_dfs(self, folder_name, target_name, report_type, df_alert_results, df_pol_results):
        ### export alerts and pols df to csv
        folder_path = Path.cwd() / 'batch' / folder_name / target_name
        # sorting by polscores
        df_pol_results = df_pol_results.sort_values(
            by=['PolScore'], ascending=False
        ).reset_index(drop=True)
        print('MaxScore:', df_pol_results.PolScore.at[0])
        print('Exporting', df_alert_results.shape[0], 'alerts...')
        print('Exporting polarity results...')
        # check user option
        if report_type == 'csv' or report_type == 'all':
            # option 1: export to csv format
            self.df_to_csv(df_alert_results, folder_path, 'alerts')
            self.df_to_csv(df_pol_results, folder_path, 'pols')
        if report_type == 'xls' or report_type == 'all':
            # option 2: export to xls format
            df_list = []
            sheetname_list = []
            # check if there is any pol result
            if df_pol_results.shape[0] > 0:
                df_list.append(df_pol_results)
                sheetname_list.append('Polarity')
            # chech if there is any alert
            if df_alert_results.shape[0] > 0:
                df_list.append(df_alert_results)
                sheetname_list.append('Alerts')
            self.df_to_xls(df_list, sheetname_list, folder_path)

    def df_to_csv(self, df, folder_path, fn_suffix):
        file_path = str(folder_path).replace('.fasta', '') + '-' + fn_suffix + '.csv'
        df.to_csv(file_path, index=False)

    def df_to_xls(self, list_dfs, list_sheet_names, folder_path):
        # create list of styleframes based on dfs
        list_sfs = []
        for df in list_dfs:
            list_sfs.append(StyleFrame(df))
        # apply cell patterns
        self.apply_highlights_pol_rows(list_sfs[0], list_dfs[0])
        # export
        file_path = str(folder_path).replace('.fasta', '') + '-report.xls'
        # best_fit factors
        StyleFrame.A_FACTOR = 6
        StyleFrame.P_FACTOR = 1.3
        with StyleFrame.ExcelWriter(file_path) as writer:
            for df_idx in range(0, len(list_sfs)):
                list_sfs[df_idx].to_excel(
                    writer, 
                    sheet_name=list_sheet_names[df_idx], 
                    index=False,
                    best_fit=list(list_dfs[df_idx].columns.values)
                )

    def apply_highlights_pol_rows(self, sf, df):
        # create masks based on row polarities
        pol_sets = df.PossiblePols.apply(lambda x: set(x.keys()))
        pp_mask = pol_sets.apply(lambda x: True if 'Pp' in x else False)
        pn_mask = pol_sets.apply(lambda x: True if 'Pn' in x else False)
        nc_mask = pol_sets.apply(lambda x: True if 'Nc' in x else False)  
        np_mask = pol_sets.apply(lambda x: True if 'Np' in x else False)      
        # list all cases and its colour code
        list_cases = []
        list_cases.append(tuple([sf[(pp_mask & pn_mask & nc_mask & np_mask)], 'a80000']))
        list_cases.append(tuple([sf[(pp_mask & pn_mask & (nc_mask ^ np_mask))], 'ff0000']))
        list_cases.append(tuple([sf[(pp_mask & pn_mask & ~(nc_mask | np_mask))], 'ff6400']))
        list_cases.append(tuple([sf[((pp_mask ^ pn_mask) & nc_mask & np_mask)], 'ff9200']))
        list_cases.append(tuple([sf[(pp_mask & ~pn_mask & (nc_mask ^ np_mask))], 'ffc800']))
        list_cases.append(tuple([sf[(~pp_mask & pn_mask & (nc_mask ^ np_mask))], 'e3b200']))
        list_cases.append(tuple([sf[~pp_mask & ~pn_mask & nc_mask & np_mask], 'fffe00']))
        list_cases.append(tuple([sf[sf.PolScore == 0], 'b3ff00']))
        # apply
        for case in list_cases:
            sf.apply_style_by_indexes(
                indexes_to_style=case[0],
                styler_obj=Styler(bg_color=case[1])
            )
        # end