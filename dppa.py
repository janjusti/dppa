### Deep Protein Polarity Analyser (DPPA) v0.1
### Copyright (c) 2019 Jan Marans Agnella Justi & Mariana de França Costa
import argparse
import numpy as np
from src.helpers.timer import CustomTimer
from src.helpers.auxfunc import AuxFuncPack

def main(args):

    # initialize timer
    c_timer = CustomTimer()

    # create handler of custom functions
    auxf_handler = AuxFuncPack()
    # create list from .fasta target file
    target_list = auxf_handler.fasta_to_list(args['input'], args['target'])
    # execute deep searcher 
    number_columns = len(target_list[0][1]) # number of chars on sequence
    # on each column from .fasta target file
    for current_col in range(0, number_columns):
        root = auxf_handler.deep_searcher(args['input'], target_list, current_col)
        # get root list of elements
        amino_leaves = []
        for leaf in root.leaves:
            amino_leaves.append(leaf.amino)
        # remove unlisted chars
        blacklisted_str = ['-', 'X', 'B', 'Z']
        amino_leaves = [x for x in amino_leaves if x not in blacklisted_str]
        # unify list
        unified_list = list(set(amino_leaves))
        print(current_col, unified_list)
        # check length
        # if > 1, ++pol_result

    # test
    # from anytree import RenderTree
    # print(RenderTree(root))

    # export reports
    # end
    c_timer.end()


if __name__ == '__main__':
    # get options from user
    parser = argparse.ArgumentParser(description='Analyse all protein alignment .fasta files.')
    parser.add_argument('-i', '--input', help='Input folder with .fasta files.', required=True)
    parser.add_argument('-t', '--target', help='Target .fasta file to be analysed.', required=True)
    user_args = vars(parser.parse_args())
    main(user_args)