from pathlib import Path
from Bio import SeqIO
import numpy as np

class AuxFuncPack:
    # package of custom functions

    def fasta_to_iter(self, fastas_folder, fasta_fn):
        # convert .fasta file into an interator
        filepath = Path.cwd() / 'arq' / 'input' / fastas_folder / fasta_fn
        fasta_iterator = SeqIO.parse(open(filepath, 'rU'), 'fasta')
        return fasta_iterator

    def deep_searcher(self, fastas_folder, fastas_iter, col_num):
        return 1