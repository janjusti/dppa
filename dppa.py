### Deep Protein Polarity Analyser (DPPA) v0.1
### Copyright (c) 2019 Jan Marans Agnella Justi & Mariana de França Costa
import argparse
from src.obj.mps import MyProtSeq
from src.helpers.timer import CustomTimer

def main(args):

    # initialize timer
    c_timer = CustomTimer()

    # create array of MyProtSeqs from .fasta target file
    # execute protocol on each column from .fasta target file
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