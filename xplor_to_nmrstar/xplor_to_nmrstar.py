import pynmrstar
import sys

class xplor_to_nmrstar(object):
    """
    Generates xplor script and input files from NMR-STAR
    """

    def __init__(self,fname):
        seq,noe,angle,rdc = self.read_input_file(fname)
        self.write_seq_file(seq)
        self.write_noe_tbl(noe)
        self.write_dihedral_tbl(angle)

    @staticmethod
    def read_input_file(in_file):
        """
        Reads input NEF/NMR-STAR file
        :param in_file: input file name with proper path
        :return: (is file readable (True/False), Content type Entry/Saveframe/Loop, data object (data) )
        """
        is_ok = False
        try:
            in_data = pynmrstar.Entry.from_file(in_file)
            is_ok = True
            msg = "Entry"
        except ValueError:
            try:
                in_data = pynmrstar.Saveframe.from_file(in_file)
                is_ok = True
                msg = "Saveframe"
            except ValueError:
                try:
                    in_data = pynmrstar.Loop.from_file(in_file)
                    is_ok = True
                    msg = "Loop"
                except ValueError as e:
                    in_data = None
                    msg = "File contains no valid saveframe or loop. Invalid file PyNMRSTAR Error:{}".format(e)
        except IOError:
            in_data = None
            msg = "File not found"
        seq = in_data.get_loops_by_category('_Chem_comp_assembly')
        noe = in_data.get_loops_by_category('_Gen_dist_constraint')
        angle = in_data.get_loops_by_category('_Torsion_angle_constraint')
        rdc = in_data.get_loops_by_category('_RDC_constraint')
        return seq,noe,angle,rdc

    @staticmethod
    def write_seq_file(seq):
        with open('protein.seq','w') as fout:
            for dat in seq:
                s = dat.get_data_by_tag('Comp_ID')[0]
                for i in range(0,len(s),18):
                    if i!=0:
                        ss=" ".join(s[j:i])
                        fout.write('{}\n'.format(ss))
                        j=i
                    else:
                        j=i
                ss = " ".join(s[i:])
                fout.write('{}'.format(ss))

    def write_noe_tbl(self,noe):
        with open('noe.tbl','w') as fout:
            for dat in noe:
                cols = dat.get_tag_names()
                id = cols.index("_Gen_dist_constraint.ID")
                seq1 = cols.index("_Gen_dist_constraint.Comp_index_ID_1")
                atm1 = cols.index("_Gen_dist_constraint.Atom_ID_1")
                seq2 = cols.index("_Gen_dist_constraint.Comp_index_ID_2")
                atm2 = cols.index("_Gen_dist_constraint.Atom_ID_2")
                t_val = cols.index("_Gen_dist_constraint.Target_val")
                lb = cols.index("_Gen_dist_constraint.Distance_lower_bound_val")
                ub = cols.index("_Gen_dist_constraint.Distance_upper_bound_val")
                logic = cols.index("_Gen_dist_constraint.Member_logic_code")
                rid = ''
                atom1=[]
                atom2=[]
                d=0.0
                dm=0.0
                dp=0.0
                for row in dat.data:
                    if rid != row[id]:
                        #print (row[id], atom1, atom2, d, dm, dp)
                        if rid != '':
                            pseudo1 = self.get_pseudo_atom(atom1)
                            pseudo2 = self.get_pseudo_atom(atom2)
                            if pseudo1 == "H":
                                pseudo1 = "HN"
                            if pseudo2 == "H":
                                pseudo2 = "HN"
                            fout.write('assign (resid {} and name {}) '
                                       'resid {} and name {}) {} {} {}\n'.format(sq1,pseudo1,sq2,pseudo2,d,dm,dp))
                           # print (row[id], self.get_psudo_atom(atom1),self.get_psudo_atom(atom2),d,dm,dp)
                        atom1 = [row[atm1]]
                        atom2 = [row[atm2]]
                        sq1 = row[seq1]
                        sq2 = row[seq2]
                        d = row[t_val]
                        dp = round(float(row[ub])-float(row[t_val]),2)
                        dm = round(float(row[t_val]) - float(row[lb]),2)

                    else:
                        atom1.append(row[atm1])
                        atom2.append(row[atm2])
                    rid = row[id]

    @staticmethod
    def write_dihedral_tbl(angle):
        with open('dihedral.tbl', 'w') as fout:
            for dat in angle:
                cols = dat.get_tag_names()
                seq1 = cols.index("_Torsion_angle_constraint.Comp_index_ID_1")
                atm1 = cols.index("_Torsion_angle_constraint.Atom_ID_1")
                seq2 = cols.index("_Torsion_angle_constraint.Comp_index_ID_2")
                atm2 = cols.index("_Torsion_angle_constraint.Atom_ID_2")
                seq3 = cols.index("_Torsion_angle_constraint.Comp_index_ID_3")
                atm3 = cols.index("_Torsion_angle_constraint.Atom_ID_3")
                seq4 = cols.index("_Torsion_angle_constraint.Comp_index_ID_4")
                atm4 = cols.index("_Torsion_angle_constraint.Atom_ID_4")
                lb = cols.index("_Torsion_angle_constraint.Angle_lower_bound_val")
                ub = cols.index("_Torsion_angle_constraint.Angle_upper_bound_val")
                tg = cols.index("_Torsion_angle_constraint.Angle_target_val")
                for row in dat:
                    try:
                        t_value = float(row[tg])
                    except ValueError:
                        ang = abs(float(row[ub])-float(row[lb]))
                        t_value = round(float(row[lb])+ang/2.0,2)
                        tol = round(ang/2.0,2)
                    fout.write('assign (resid {} and name {}) (resid {} '
                               'and name {}\n'.format(row[seq1],row[atm1],row[seq2],row[atm2]))
                    fout.write('\t(resid {} and name {}) (resid {} '
                               'and name {}) 1.0 {} {} 2\n'.format(row[seq3],row[atm3],row[seq4],row[atm4],t_value,tol))

    @staticmethod
    def get_pseudo_atom(atom_list):
        if len(atom_list):
            if atom_list.count(atom_list[0])==len(atom_list):
                pseudo = atom_list[0]
            else:
                for i in range (1,len(atom_list)):
                    atm_list = [j[:-i] for j in atom_list]
                    if atm_list.count(atm_list[0])== len(atm_list):
                        pseudo = '{}#'.format(atm_list[0])
                        break
        else:
            pseudo = ''
        return pseudo

if __name__ =="__main__":
    p = xplor_to_nmrstar('nef_examples/2lci.str')

