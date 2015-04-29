from pymol import cmd

def ligandNeighbors(prot, chain, ligand):
   
    print ligand 
	# set some string names for temporary objects/selections
    tmpObj = cmd.get_unused_name("_tmp")
    tmpObj2 = cmd.get_unused_name("_bla")
 
	# operate on a new object & turn off the original
    sel = "{} and chain {}".format(prot, chain)
    cmd.create(tmpObj, sel)
    cmd.disable(sel)

    cmd.select(tmpObj2, "( resn {} around 10 ) and chain {}".format(ligand, chain))
 
    atoms = set()
    cmd.iterate(tmpObj2, "atoms.add((chain,resv))", space=locals())
    cmd.delete(tmpObj)

    return atoms
