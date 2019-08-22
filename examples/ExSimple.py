from Bio import Entrez
from Bio.Align.Applications import MuscleCommandline
import dppa.core as solver
import time

def generate_fasta_from_ids(id_list, fasta_output):
    Entrez.email = 'example@example.com'
    f = open(fasta_output, 'w+')
    for curr_id in id_list:
        curr_req = Entrez.efetch(db='protein', id=curr_id, rettype='fasta')
        curr_seq = curr_req.read()
        f.write(curr_seq)
        print(curr_id, 'successfully fetched.')
        time.sleep(1)
    f.close()
    print(fasta_output, 'sucessfully created.')

def align_via_muscle(input_file, output_file):
    comm_muscle = MuscleCommandline(
        input=input_file, out=output_file
    )
    comm_muscle()
    print('Alignment', input_file, '>', output_file, 'done.')

def main_simple():
    id_list = ['AGM34409.1', 'AGM34408.1', 'AIQ82784.1', 'AIQ82783.1']
    generate_fasta_from_ids(id_list, 'simple-unaligned.fasta')

    align_via_muscle('simple-unaligned.fasta', 'simple.fasta')

    results = solver.run('simple.fasta')
    solver.export(results, 'csv', 'simple')