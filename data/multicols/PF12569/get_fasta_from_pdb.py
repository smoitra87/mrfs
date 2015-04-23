import pymol
from pymol import cmd;
import sys
pymol.finish_launching()
import time ; time.sleep(1);
cmd.load('2XPI.pdb')
cmd.save('2XPI.fasta','chain A')
cmd.quit()
