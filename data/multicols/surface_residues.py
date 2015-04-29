'''
http://pymolwiki.org/index.php/FindSurfaceResidues
'''

from pymol import cmd


def findAtomExposure(selection="all"):
    """
DESCRIPTION

    Finds those atoms on the surface of a protein
    that have at least 'cutoff' exposed A**2 surface area.

USAGE

    findSurfaceAtoms [ selection, [ cutoff ]]

SEE ALSO

    findSurfaceResidues
    """
    tmpObj = cmd.get_unused_name("_tmp")
    cmd.create(tmpObj, "(" + selection + ") and polymer", zoom=0)

    cmd.set("dot_solvent", 1, tmpObj)
    cmd.get_area(selection=tmpObj, load_b=1)

    atoms = list()
    cmd.iterate(tmpObj, "atoms.append((name,chain,resv,b))", space=locals())
    cmd.delete(tmpObj)
    return atoms
