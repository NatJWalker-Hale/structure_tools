#! /usr/bin/python3

import sys
import argparse
from Bio import PDB


def parse_p2rank(inpath: str) -> list:
    outl = []
    with open(inpath, "r") as inf:
        for line in inf:
            line = line.strip().split(",")
            line = [x.strip() for x in line]  # clear bad whitespace
            if line[-1] == "1":
                outl.append((int(line[1]), line[2]))  # list of pos, res tuples
    return outl


def parse_p2rank_cent(inpath: str) -> tuple:
    with open(inpath, "r") as inf:
        next(inf)
        line = inf.readline() # only use top scoring pocket
        line = [x.strip() for x in line.split(",")]
        x, y, z = [float(v) for v in line[6:9]]

    return (x, y, z)


def parse_substitutions(inpath: str) -> list[list]:
    outl = []
    with open(inpath, "r") as inf:
        for line in inf:
            # next(inf)  # skip header
            line = line.strip().split()
            outl.append(line)
    return outl


def structure_distance_pocket(pdbpath: str, pocket: list,
                              substitutions: list[list], pocketcent=None):
    parser = PDB.PDBParser()
    struc = parser.get_structure(id="target", file=pdbpath)
    residues = [r for r in struc.get_residues()]  # all
    if pocketcent is None:
        pocketres = [residues[p] for p in [x[0]-1 for x in pocket]]  # just pocket
        # need to 1-index here for consistency with p2rank
        pocketCA = [x["CA"] for x in pocketres]  # alpha carbons of pocket
        x_coords = [x.coord[0] for x in pocketCA]
        y_coords = [x.coord[1] for x in pocketCA]
        z_coords = [x.coord[2] for x in pocketCA]
        # calculate centroid as mean
        centroid = (sum(x_coords) / len(x_coords),
                    sum(y_coords) / len(y_coords),
                    sum(z_coords) / len(z_coords))
    else:
        centroid = pocketcent
    # now let's make a fake atom to be our centroid
    myatom = PDB.Atom.Atom("CA", centroid, 0, 0, "null",
                               "null", 0, element="C")
    for s in substitutions:
        ca1 = residues[int(s[1]) - 1]["CA"]  # change to 0 idx
        ca2 = myatom
        s += [ca2 - ca1]
        

def structure_distance_pocket_all(pdbpath: str, pocket: list, pocketcent=None):
    parser = PDB.PDBParser()
    struc = parser.get_structure(id="target", file=pdbpath)
    residues = [r for r in struc.get_residues()]  # all
    if pocketcent is None:
        pocketres = [residues[p] for p in [x[0]-1 for x in pocket]]  # just pocket
        # need to 1-index here for consistency with p2rank
        pocketCA = [x["CA"] for x in pocketres]  # alpha carbons of pocket
        x_coords = [x.coord[0] for x in pocketCA]
        y_coords = [x.coord[1] for x in pocketCA]
        z_coords = [x.coord[2] for x in pocketCA]
        # calculate centroid as mean
        centroid = (sum(x_coords) / len(x_coords),
                    sum(y_coords) / len(y_coords),
                    sum(z_coords) / len(z_coords))
    else:
        centroid = pocketcent
    # now let's make a fake atom to be our centroid
    myatom = PDB.Atom.Atom("CA", centroid, 0, 0, "null",
                               "null", 0, element="C")
    dist = []
    for r in residues:
        ca1 = r["CA"]
        ca2 = myatom
        dist.append(ca2 - ca1)
    return dist


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("perrespocket", help="p2rank CSV output containing \
                        of each residue to binding pocket")
    parser.add_argument("-ps", "--pocketsumm", help="p2rank CSV output \
                        containing pocket summary")
    parser.add_argument("structure", help="PDB-formatted structure of protein")
    parser.add_argument("substitutions", help="TSV-formatted results from \
                        get_conv_subs_from_asr.py, with parent reference mode")
    parser.add_argument("-a", "--all", help="'All' mode: print distances to \
                        pocket centroid for all residues, not just subs",
                        default = False, action="store_true")
    args = parser.parse_args()

   
    pocket_res = parse_p2rank(args.perrespocket)
    subs = parse_substitutions(args.substitutions)

    if not args.all:
        if args.pocketsumm is not None:
            pocketcent = parse_p2rank_cent(args.pocketsumm)
            structure_distance_pocket(args.structure, pocket_res, subs,
                                      pocketcent)
            # mod subs in place
        else:
            structure_distance_pocket(args.structure, pocket_res, subs)
        print("pos_aln\tpos_par\tpar\tdesc\tsub\ttype\tdist")
        for s in subs:
            print("\t".join([str(x) for x in s]))
    else:
        if args.pocketsumm is not None:
            pocketcent = parse_p2rank_cent(args.pocketsumm)
            d = structure_distance_pocket_all(args.structure, pocket_res,
                                              pocketcent)
        else:
            d = structure_distance_pocket_all(args.structure, pocket_res)
        i = 1
        struc_name = args.structure.split("_")[0]
        print("\t".join(["ref", "pos", "dist"]))
        for distance in d:
            # 1-indexed
            print(f"{struc_name}\t{i}\t{distance}")
            i += 1

