import os, sys
import time
import numpy as np
from Bio import SeqIO
from commands import getstatusoutput
import csv

import urllib

def fetch_pdb(id):
    url = 'http://www.rcsb.org/pdb/files/%s.pdb' % id
    return urllib.urlopen(url).read()

if __name__ == '__main__':

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("--family_csv", type=str, \
        default='pfam_structures.csv', help="Pfam structures csv file")
    parser.add_argument("--pfamid", type=str, help="PFam id run")
    parser.add_argument("--gap_open_cost", type=int, default=30,\
         help="Gap open cost for clustalw2")
    parser.add_argument("--exposure_cutoff", type=float, default=2.5,\
         help="Cutoff for surface exposure")
    args = parser.parse_args()


    with open(args.family_csv, 'rb') as fin:
        reader = csv.reader(fin)
        header = next(reader)
        print header
        rows = [row for row in reader]
    
    col2list = dict(zip(header, zip(*rows)))
   
    # Fetch the pdb for the pfamid
    if args.pfamid not in col2list['pfamid']:
        raise ValueError("Unknown pdamid {}".format(args.pfamid))

    #Fetch the pdb 
    # Find the index of the pfamid 
    pfamidx = next(idx for (idx,p) in enumerate(col2list['pfamid']) if p == args.pfamid)

    pdb_file =os.path.join(args.pfamid, col2list['pdbid'][pfamidx] + ".pdb") 
    if not os.path.exists(pdb_file):
        pdb_data = fetch_pdb(col2list['pdbid'][pfamidx]) 
        with open(pdb_file, "w") as fout:
            fout.write(pdb_data)
    else :
        print "Pdb {}  already exists".format(pdb_file);
    
    # Fetch the consensus sequence
    # Check of consensus already exists
    consensus_file = os.path.join(args.pfamid, 'baker.fasta')
    if not os.path.exists(consensus_file):
        from get_consensus import fetch_consensus
        fetch_consensus(args.pfamid) 
    else: 
        print "Consensus file : {} already exists".format(consensus_file)

    # Get the pdb sequence for the chain
    chain = col2list['chain'][pfamidx]
    pdb_fasta_file = os.path.join(args.pfamid, "pdb.fasta")
    if not os.path.exists(pdb_fasta_file):
        cmd = "python pdbtools/pdb_seq.py {0} -c {1}".format(pdb_file, chain)
        out = getstatusoutput(cmd)
        if out[0] > 0:
            raise EnvironmentError("pdb_seq.py cmd failed")
        
        with open(pdb_fasta_file, 'w') as fout:
            fout.write(out[1].strip()) 
    else:
        print "PDB fasta file : {} already exists".format(pdb_fasta_file)
        
    # concatenate the two files
    cat_file = os.path.join(args.pfamid, "joint.fasta")
    os.system("cat {} {} > {}".format(consensus_file, pdb_fasta_file, cat_file))
     
    # Runt the alignment code
    aligned_file = os.path.join(args.pfamid, "aligned.fasta")
    cmd = "./align_prot_baker.sh {} {} {}".format(cat_file, aligned_file, \
        args.gap_open_cost)
    print "Running {}".format(cmd)
    out = getstatusoutput(cmd)
    if out[0] > 0:
        raise EnvironmentError("pdb_seq.py cmd failed")

    # Fix the pdb numbering 
    cmd = "python pdbtools/pdb_residue_renumber.py {} -c {} -s 1 -g".format(pdb_file, chain)
    print "Running {}".format(cmd)
    out = getstatusoutput(cmd)
    if out[0] > 0:
        raise EnvironmentError("pdb_residue_renumber.py cmd failed")

    from_f, to_f = col2list['pdbid'][pfamidx] + "_res-renum.pdb", pdb_file
    cmd = "mv {} {}".format(from_f, to_f)
    print "Running {}".format(cmd)
    out = getstatusoutput(cmd)
    if out[0] > 0:
        raise EnvironmentError("mv failed")
    
    # Extract the mapping
    
    seq_gen = SeqIO.parse(aligned_file,"fasta")
    baker_seq = str(next(seq_gen).seq)
    pdb_seq = str(next(seq_gen).seq)
    baker_length = len(baker_seq.replace('-', ''))
    pdb_length = len(pdb_seq.replace('-', ''))
    
    assert len(baker_seq) == len(pdb_seq)
    
    overlaps = [idx for idx in \
        range(len(baker_seq)) if baker_seq[idx]!='-' and pdb_seq[idx]!='-']
    assert len(overlaps) <= min(baker_length, pdb_length)   

    baker_map = [idx for idx,e in enumerate(baker_seq) if e != '-'] 
    pdb_map = [idx for idx,e in enumerate(pdb_seq) if e != '-'] 
    
    baker_idxs = [i for (i,j) in enumerate(baker_map) if j in overlaps]
    pdb_idxs = [i for (i,j) in enumerate(pdb_map) if j in overlaps]

    # Launch pymol and make a list of calculations
    import pymol
    from pymol import cmd
    pymol.finish_launching()
    time.sleep(1)
    cmd.load(pdb_file)
    from surface_residues import findAtomExposure
    atoms = findAtomExposure('2ID5 and chain A')
    cmd.quit()

    time.sleep(1);
    from collections import defaultdict
    atomdict = defaultdict(float)
    for record in atoms:
        name, chain, residx, exposure = record
        atomdict[residx] = max(atomdict[residx], exposure)

    atomdict = dict(atomdict)
    core_res = [r for (r,exp) in atomdict.items() if exp < args.exposure_cutoff]
    surface_res = [r for (r,exp) in atomdict.items() if exp > args.exposure_cutoff]
    print core_res
    print surface_res
