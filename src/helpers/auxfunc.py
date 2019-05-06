from pathlib import Path
from Bio.SeqIO import FastaIO
from anytree import AnyNode, RenderTree

class AuxFuncPack:
    # package of custom functions

    def fasta_to_list(self, fastas_folder, fasta_fn):
        # convert .fasta file into a list
        # list format: [(name, seq)]
        filepath = Path.cwd() / 'arq' / 'input' / fastas_folder / fasta_fn
        fasta_list = list(FastaIO.SimpleFastaParser(open(filepath, 'rU')))
        return fasta_list

    def deep_searcher(self, fastas_folder, fasta_list, col_num):
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
                pass
            if node.amino is 'X' or \
                node.amino is 'B' or \
                node.amino is 'Z':
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
                    self.deep_searcher(fastas_folder, deeper_list, shifted_col_num).parent = node
                else:
                    # pure sequence, alert
                    pass

        return loc_root