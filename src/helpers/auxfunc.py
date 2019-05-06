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
            if node.amino in ['X', 'B', 'Z']:
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