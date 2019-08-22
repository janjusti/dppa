from Bio import Entrez, AlignIO, SeqIO
from Bio.Align import AlignInfo
from Bio.SeqRecord import SeqRecord
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

def get_consensus_seq(fasta_file, q_id):
    align = AlignIO.read(fasta_file, 'fasta')
    seq = SeqRecord(
        AlignInfo.SummaryInfo(align).gap_consensus(), 
        id=q_id,
        description=''
    )
    print("'" + q_id + "' consensus sequence generated (original: " + fasta_file + ')')
    return seq

def write_fasta(fasta_file, seq_list):
    with open(fasta_file, 'w') as handle:
        SeqIO.write(seq_list, handle, 'fasta')
    print(fasta_file, 'successfully written.')

def main_group():
    groupA_ids = ['AGM34409.1', 'AGM34408.1']
    groupB_ids = ['AIQ82784.1', 'AIQ82783.1']
    generate_fasta_from_ids(groupA_ids, 'groupA-unaligned.fasta')
    generate_fasta_from_ids(groupB_ids, 'groupB-unaligned.fasta')

    align_via_muscle('groupA-unaligned.fasta', 'groupA.fasta')
    align_via_muscle('groupB-unaligned.fasta', 'groupB.fasta')

    seqA = get_consensus_seq('groupA.fasta', 'groupA_any_sentence')
    seqB = get_consensus_seq('groupB.fasta', 'groupB_any_sentence')
    write_fasta('groupAB-cons.fasta', [seqA, seqB])

    align_via_muscle('groupAB-cons.fasta', 'groups-target.fasta')

    results = solver.run('groups-target.fasta', searchable_keyphrase='any sentence')
    solver.export(results, 'csv', 'groups')