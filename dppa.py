### Deep Protein Polarity Analyser (DPPA) v0.1
### Copyright (c) 2019 Jan Marans Agnella Justi & Mariana de França Costa
import argparse
from src.helpers.timer import CustomTimer
from src.helpers.auxfunc import AuxFuncPack

def main(args):

    # initialize timer
    c_timer = CustomTimer()

    # create handler of custom functions
    auxf_handler = AuxFuncPack()
    # create iterator from .fasta target file
    target_iter = auxf_handler.fasta_to_iter(args['input'], args['target'])
    # execute deep searcher 
    # (?) better way to calculate this?
    number_columns = len(str(list(target_iter)[0].seq))
    # on each column from .fasta target file
    for current_col in range(0, number_columns):
        root = auxf_handler.deep_searcher(args['input'], target_iter, current_col)
        # get root list of elements
        # remove unlisted chars (-, X)
        # unify list
        # check length
        # if > 1, ++pol_result


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