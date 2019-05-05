### Deep Protein Polarity Analyser (DPPA) v0.1
import argparse
from src.helpers.timer import CustomTimer

def main(args):
    # initialize timer
    c_timer = CustomTimer()


if __name__ == '__main__':
    # get options from user
    parser = argparse.ArgumentParser(description='Analyse all protein alignment .fasta files.')
    parser.add_argument('-i', '--input', help='Input folder with .fasta files.', required=True)
    parser.add_argument('-t', '--target', help='Target .fasta file to be analysed.', required=True)
    user_args = vars(parser.parse_args())
    main(user_args)