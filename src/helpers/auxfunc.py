from pathlib import Path
from Bio.SeqIO import FastaIO
from anytree import AnyNode
import pandas as pd

class AuxFuncPack:
    # package of custom functions

    def fasta_to_list(self, fastas_folder, fasta_fn):
        # convert .fasta file into a list
        # list format: [(name, seq)]
        filepath = Path.cwd() / 'batch' / fastas_folder / fasta_fn
        fasta_list = list(FastaIO.SimpleFastaParser(open(filepath, 'rU')))
        return fasta_list

    def deep_searcher(self, fastas_folder, fasta_list, col_num, df_alert):
        # create local root
        loc_root = AnyNode(name='LocRoot')
        # create nodes from column number
        for row in range(0, len(fasta_list)):
            AnyNode(
                name=fasta_list[row][0], 
                amino=fasta_list[row][1][col_num], 
                parent=loc_root
            )
        # run through nodes
        for node in loc_root.children:
            # check amino cases
            if node.amino is '-':
                # increment on df_alert
                df_alert = df_alert.append(
                    {
                        'SeqName': node.name.replace('_', ' '),
                        'ColNum': col_num+1,
                        'AlertType': 'Gap'
                    }
                , ignore_index=True)
            if node.amino in ['X', 'B', 'Z', '?']:
                if 'consensus_sequence' in node.name: 
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
                    deeper_result = self.deep_searcher(fastas_folder, deeper_list, shifted_col_num, df_alert)
                    deeper_result[0].parent = node
                    df_alert = deeper_result[1]
                else:
                    # pure sequence, increment on df_alert
                    pure_str = 'Pure ' + node.amino
                    df_alert = df_alert.append(
                        {
                            'SeqName': node.name.replace('_', ' '),
                            'ColNum': col_num+1,
                            'AlertType': pure_str
                        }
                    , ignore_index=True)
                    pass

        return loc_root, df_alert

    def get_polarities_list(self, amino_list):
        # get list of polarities within aminos list
        filepath = Path.cwd() / 'src' / 'conv' / 'pols.csv'
        df_pols = pd.read_csv(filepath)
        curr_pol_list = []
        for curr_pol in range(0, df_pols.shape[0]):
            pol_amino_list = list(df_pols.Amino.at[curr_pol])
            for curr_amino in amino_list:
                if curr_amino in pol_amino_list:
                    curr_pol_list.append(df_pols.Type.at[curr_pol])
        # unify polarity list
        curr_pol_list = list(set(curr_pol_list))
        return curr_pol_list

    def export_dfs(self, folder_name, target_name, df_alert_results, df_pol_results):
        # export alerts and pols df to csv
        filepath = Path.cwd() / 'batch' / folder_name / target_name
        alerts_filepath = str(filepath).replace('.fasta', '') + '-alerts.csv'
        pols_filepath = str(filepath).replace('.fasta', '') + '-pols.csv'
        print('Exporting alerts...')
        df_alert_results.to_csv(alerts_filepath, index=False)
        print('Exporting polarity results...')
        df_pol_results.to_csv(pols_filepath, index=False) 