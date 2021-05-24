##
# File: BMRBChemShiftStat.py
# Date: 20-Sep-2019
#

##
""" Wrapper class for retrieving BMRB chemical shift statistics.
    @author: Masashi Yokochi
"""
import sys
import os
import os.path
import csv
import re
import copy
import pickle
import collections

# from wwpdb.utils.config.ConfigInfo import ConfigInfo, getSiteId
# from wwpdb.utils.nmr.io.ChemCompIo import ChemCompReader

class BMRBChemShiftStat:
    """ Wrapper class for retrieving BMRB chemical shift statistics.
    """

    def __init__(self):
        # lazy import of others (non-standard residues)
        self.lazy_others = True

        # directory
        self.stat_dir = os.path.dirname(__file__) + '/bmrb_cs_stat/'

        # statistics objects
        self.aa_filt = []
        self.aa_full = []

        self.dna_filt = []
        self.dna_full = []

        self.rna_filt = []
        self.rna_full = []

        self.others = []
        self.extras = []

        self.__aa_comp_ids = set()
        self.__dna_comp_ids = set()
        self.__rna_comp_ids = set()
        self.__oth_comp_ids = set()
        self.__ext_comp_ids = set()

        self.__std_comp_ids = set()
        self.__all_comp_ids = set()

        self.aa_threshold = 0.1
        self.na_threshold = 0.3

        self.max_count_th = 10

        # CCD accessing utility
        self.__verbose = False
        self.__lfh = sys.stderr

        #self.__cI = ConfigInfo(getSiteId())
        self.__ccCvsPath = 'lib/'

        # self.__ccR = ChemCompReader(self.__verbose, self.__lfh)
        # self.__ccR.setCachePath(self.__ccCvsPath)

        self.__last_comp_id = None
        self.__last_comp_id_test = False
        self.__last_chem_comp_dict = None
        self.__last_chem_comp_atoms = None
        self.__last_chem_comp_bonds = None

        # taken from wwpdb.apps.ccmodule.io.ChemCompIo
        self.__chem_comp_atom_dict = [
                ('_chem_comp_atom.comp_id','%s','str',''),
                ('_chem_comp_atom.atom_id','%s','str',''),
                ('_chem_comp_atom.alt_atom_id','%s','str',''),
                ('_chem_comp_atom.type_symbol','%s','str',''),
                ('_chem_comp_atom.charge','%s','str',''),
                ('_chem_comp_atom.pdbx_align','%s','str',''),
                ('_chem_comp_atom.pdbx_aromatic_flag','%s','str',''),
                ('_chem_comp_atom.pdbx_leaving_atom_flag','%s','str',''),
                ('_chem_comp_atom.pdbx_stereo_config','%s','str',''),
                ('_chem_comp_atom.model_Cartn_x','%s','str',''),
                ('_chem_comp_atom.model_Cartn_y','%s','str',''),
                ('_chem_comp_atom.model_Cartn_z','%s','str',''),
                ('_chem_comp_atom.pdbx_model_Cartn_x_ideal','%s','str',''),
                ('_chem_comp_atom.pdbx_model_Cartn_y_ideal','%s','str',''),
                ('_chem_comp_atom.pdbx_model_Cartn_z_ideal','%s','str',''),
                ('_chem_comp_atom.pdbx_component_atom_id','%s','str',''),
                ('_chem_comp_atom.pdbx_component_comp_id','%s','str',''),
                ('_chem_comp_atom.pdbx_ordinal','%s','str','')
                ]

        atom_id = next(d for d in self.__chem_comp_atom_dict if d[0] == '_chem_comp_atom.atom_id')
        self.__cca_atom_id = self.__chem_comp_atom_dict.index(atom_id)

        aromatic_flag = next(d for d in self.__chem_comp_atom_dict if d[0] == '_chem_comp_atom.pdbx_aromatic_flag')
        self.__cca_aromatic_flag = self.__chem_comp_atom_dict.index(aromatic_flag)

        leaving_atom_flag = next(d for d in self.__chem_comp_atom_dict if d[0] == '_chem_comp_atom.pdbx_leaving_atom_flag')
        self.__cca_leaving_atom_flag = self.__chem_comp_atom_dict.index(leaving_atom_flag)

        type_symbol = next(d for d in self.__chem_comp_atom_dict if d[0] == '_chem_comp_atom.type_symbol')
        self.__cca_type_symbol = self.__chem_comp_atom_dict.index(type_symbol)

        # taken from wwpdb.apps.ccmodule.io.ChemCompIo
        self.__chem_comp_bond_dict = [
                ('_chem_comp_bond.comp_id','%s','str',''),
                ('_chem_comp_bond.atom_id_1','%s','str',''),
                ('_chem_comp_bond.atom_id_2','%s','str',''),
                ('_chem_comp_bond.value_order','%s','str',''),
                ('_chem_comp_bond.pdbx_aromatic_flag','%s','str',''),
                ('_chem_comp_bond.pdbx_stereo_config','%s','str',''),
                ('_chem_comp_bond.pdbx_ordinal','%s','str','')
                ]

        atom_id_1 = next(d for d in self.__chem_comp_bond_dict if d[0] == '_chem_comp_bond.atom_id_1')
        self.__ccb_atom_id_1 = self.__chem_comp_bond_dict.index(atom_id_1)

        atom_id_2 = next(d for d in self.__chem_comp_bond_dict if d[0] == '_chem_comp_bond.atom_id_2')
        self.__ccb_atom_id_2 = self.__chem_comp_bond_dict.index(atom_id_2)

        aromatic_flag = next(d for d in self.__chem_comp_bond_dict if d[0] == '_chem_comp_bond.pdbx_aromatic_flag')
        self.__ccb_aromatic_flag = self.__chem_comp_bond_dict.index(aromatic_flag)

        if not self.loadStatFromPickleFiles():
            self.loadStatFromCsvFiles()

    def isOk(self):
        """ Return whether all BMRB chemical shift statistics are available.
        """

        return len(self.aa_filt) > 0 and len(self.aa_full) > 0 and len(self.dna_filt) > 0 and len(self.dna_full) > 0 and len(self.rna_filt) > 0 and len(self.rna_full) and (len(self.others) > 0 or self.lazy_others)

    def hasCompId(self, comp_id):
        """ Return whether a given comp_id has BMRB chemical shift statistics.
        """

        if comp_id in self.__std_comp_ids:
            return True

        self.loadOtherStatFromCsvFiles()

        return comp_id in self.__all_comp_ids

    def getTypeOfCompId(self, comp_id):
        """ Return type of a given comp_id.
            @return: array of bool: peptide, nucleotide, carbohydrate
        """

        if comp_id in self.__aa_comp_ids:
            return True, False, False

        elif comp_id in self.__dna_comp_ids or comp_id in self.__rna_comp_ids:
            return False, True, False

        if self.__updateChemCompDict(comp_id):
            type = self.__last_chem_comp_dict['_chem_comp.type']

            if 'PEPTIDE' in type:
                return True, False, False
            elif 'DNA' in type or 'RNA' in type:
                return False, True, False
            elif 'SACCHARIDE' in type:
                return False, False, True

        peptide_like = len(self.getBackBoneAtoms(comp_id, True, True, False, False))
        nucleotide_like = len(self.getBackBoneAtoms(comp_id, True, False, True, False))
        carbohydrate_like = len(self.getBackBoneAtoms(comp_id, True, False, False, True))

        return peptide_like > nucleotide_like and peptide_like > carbohydrate_like,\
                nucleotide_like > peptide_like and nucleotide_like > carbohydrate_like,\
                carbohydrate_like > peptide_like and carbohydrate_like > nucleotide_like

    def hasEnoughStat(self, comp_id, primary=True):
        """ Return whether a given comp_id has enough chemical shift statistics.
        """

        if comp_id in self.__std_comp_ids:
            return True

        self.loadOtherStatFromCsvFiles()

        if not comp_id in self.__all_comp_ids:
            return False

        if primary:
            if any(i for i in self.others if i['comp_id'] == comp_id and i['primary']):
                return True
        else:
            if any(i for i in self.others if i['comp_id'] == comp_id and 'secondary' in i and i['secondary']):
                return True

        return False

    def get(self, comp_id, diamagnetic=True):
        """ Return BMRB chemical shift statistics for a given comp_id.
        """

        if comp_id in self.__aa_comp_ids:

            if diamagnetic:
                return [i for i in self.aa_filt if i['comp_id'] == comp_id]
            else:
                return [i for i in self.aa_full if i['comp_id'] == comp_id]

        elif comp_id in self.__dna_comp_ids:

            if diamagnetic:
                return [i for i in self.dna_filt if i['comp_id'] == comp_id]
            else:
                return [i for i in self.dna_full if i['comp_id'] == comp_id]

        elif comp_id in self.__rna_comp_ids:

            if diamagnetic:
                return [i for i in self.rna_filt if i['comp_id'] == comp_id]
            else:
                return [i for i in self.rna_full if i['comp_id'] == comp_id]

        self.loadOtherStatFromCsvFiles()

        if not comp_id in self.__all_comp_ids:
            return []

        return [i for i in self.others if i['comp_id'] == comp_id]

    def __get(self, comp_id, diamagnetic=True):
        """ Return atom list for a given comp_id.
        """

        if comp_id in self.__aa_comp_ids:

            if diamagnetic:
                return [i for i in self.aa_filt if i['comp_id'] == comp_id]
            else:
                return [i for i in self.aa_full if i['comp_id'] == comp_id]

        elif comp_id in self.__dna_comp_ids:

            if diamagnetic:
                return [i for i in self.dna_filt if i['comp_id'] == comp_id]
            else:
                return [i for i in self.dna_full if i['comp_id'] == comp_id]

        elif comp_id in self.__rna_comp_ids:

            if diamagnetic:
                return [i for i in self.rna_filt if i['comp_id'] == comp_id]
            else:
                return [i for i in self.rna_full if i['comp_id'] == comp_id]

        self.loadOtherStatFromCsvFiles()

        if not comp_id in self.__all_comp_ids:
            self.__appendExtraFromCcd(comp_id)

        if comp_id in self.__oth_comp_ids:
            return [i for i in self.others if i['comp_id'] == comp_id]

        elif comp_id in self.__ext_comp_ids:
            return [i for i in self.extras if i['comp_id'] == comp_id]

        return []

    def getMaxAmbigCodeWoSetId(self, comp_id, atom_id):
        """ Return maximum ambiguity code of a given atom that does not require declaration of ambiguity set ID.
            @return: one of (1, 2, 3), 0 for not found
        """

        if not comp_id in self.__std_comp_ids:
            self.loadOtherStatFromCsvFiles()

        if not comp_id in self.__all_comp_ids:
            self.__appendExtraFromCcd(comp_id)

        try:

            d = next(i['desc'] for i in self.__get(comp_id) if i['atom_id'] == atom_id)

            if 'geminal' in d:
                return 2
            elif d == 'aroma-opposite':
                return 3
            else:
                return 1

        except StopIteration:
            return 0

    def getGeminalAtom(self, comp_id, atom_id):
        """ Return geminal or aromatic opposite atom of a given atom.
        """

        if not comp_id in self.__std_comp_ids:
            self.loadOtherStatFromCsvFiles()

        if not comp_id in self.__all_comp_ids:
            self.__appendExtraFromCcd(comp_id)

        cs_stat = self.__get(comp_id)

        try:

            d = next(i['desc'] for i in cs_stat if i['atom_id'] == atom_id)

            if d == 'methyl-geminal' and atom_id[0] == 'H':
                return next(i['atom_id'] for i in cs_stat if i['desc'] == d and i['atom_id'] != atom_id and i['atom_id'][:-2] == atom_id[:-2] and i['atom_id'][-1] == atom_id[-1])
            elif 'geminal' in d or d == 'aroma-opposite':
                if not atom_id.endswith("'"):
                    return next(i['atom_id'] for i in cs_stat if i['desc'] == d and i['atom_id'] != atom_id and i['atom_id'][:-1] == atom_id[:-1])
                elif atom_id.endswith("''"):
                    return next(i['atom_id'] for i in cs_stat if i['desc'] == d and i['atom_id'] != atom_id and i['atom_id'] == atom_id[:-1])
                else:
                    return next(i['atom_id'] for i in cs_stat if i['desc'] == d and i['atom_id'] != atom_id and i['atom_id'] == atom_id + "'")
            else:
                return None

        except StopIteration:
            return None

    def getAllAtoms(self, comp_id, excl_minor_atom=False, primary=False):
        """ Return all atoms of a given comp_id.
        """

        if not comp_id in self.__std_comp_ids:
            self.loadOtherStatFromCsvFiles()

        if not comp_id in self.__all_comp_ids:
            self.__appendExtraFromCcd(comp_id)

        cs_stat = self.__get(comp_id)

        if comp_id in self.__std_comp_ids or primary:
            return [i['atom_id'] for i in cs_stat if\
                    (not excl_minor_atom or (excl_minor_atom and i['primary']))]

        return [i['atom_id'] for i in cs_stat if\
                (not excl_minor_atom or not 'secondary' in i or (excl_minor_atom and i['secondary']))]

    def getBackBoneAtoms(self, comp_id, excl_minor_atom=False, polypeptide_like=False, polynucleotide_like=False, carbohydrates_like=False):
        """ Return backbone atoms of a given comp_id.
        """

        if not comp_id in self.__std_comp_ids:
            self.loadOtherStatFromCsvFiles()

        if not comp_id in self.__all_comp_ids:
            self.__appendExtraFromCcd(comp_id)

        if polypeptide_like == False and polynucleotide_like == False and carbohydrates_like == False:
            polypeptide_like, polynucleotide_like, carbohydrates_like = self.getTypeOfCompId(comp_id)

        cs_stat = self.__get(comp_id)

        if comp_id in self.__aa_comp_ids:
            return [i['atom_id'] for i in cs_stat if i['atom_id'] in ['C', 'CA', 'CB', 'H', 'HA', 'HA2', 'HA3', 'N'] and
                    (not excl_minor_atom or (excl_minor_atom and i['primary']))]

        elif comp_id in self.__dna_comp_ids:
            return [i['atom_id'] for i in cs_stat if i['atom_id'] in ["C1'", "C2'", "C3'", "C4'", "C5'", "H1'", "H2'", "H2''", "H3'", "H4'", "H5'", "H5''", 'P'] and
                    (not excl_minor_atom or (excl_minor_atom and i['primary']))]

        elif comp_id in self.__rna_comp_ids:
            return [i['atom_id'] for i in cs_stat if i['atom_id'] in ["C1'", "C2'", "C3'", "C4'", "C5'", "H1'", "H2'", "H3'", "H4'", "H5'", "H5''", "HO2'", 'P'] and
                    (not excl_minor_atom or (excl_minor_atom and i['primary']))]

        elif polypeptide_like:
            return [i['atom_id'] for i in cs_stat if i['atom_id'] in ['C', 'CA', 'CB', 'H', 'HA', 'HA2', 'HA3', 'N'] and
                    (not excl_minor_atom or (excl_minor_atom and i['primary']))]

        elif polynucleotide_like:
            return [i['atom_id'] for i in cs_stat if i['atom_id'] in ["C1'", "C2'", "C3'", "C4'", "C5'", "H1'", "H2'", "H2''", "H3'", "H4'", "H5'", "H5''", 'P'] and
                    (not excl_minor_atom or not 'secondary' in i or (excl_minor_atom and i['secondary']))]

        elif carbohydrates_like:
            return [i['atom_id'] for i in cs_stat if i['atom_id'] in ["C1", "C2", "C3", "C4", "C5", "C6", "H61", "H62"] and
                    (not excl_minor_atom or not 'secondary' in i or (excl_minor_atom and i['secondary']))]

        return []

    def getAromaticAtoms(self, comp_id, excl_minor_atom=False, primary=False):
        """ Return aromatic atoms of a given comp_id.
        """

        if not comp_id in self.__std_comp_ids:
            self.loadOtherStatFromCsvFiles()

        if not comp_id in self.__all_comp_ids:
            self.__appendExtraFromCcd(comp_id)

        cs_stat = self.__get(comp_id)

        if comp_id in self.__std_comp_ids or primary:
            return [i['atom_id'] for i in cs_stat if 'aroma' in i['desc'] and
                    (not excl_minor_atom or (excl_minor_atom and i['primary']))]

        return [i['atom_id'] for i in cs_stat if 'aroma' in i['desc'] and
                (not excl_minor_atom or not 'secondary' in i or (excl_minor_atom and i['secondary']))]

    def getMethylAtoms(self, comp_id, excl_minor_atom=False, primary=False):
        """ Return atoms in methyl group of a given comp_id.
        """

        if not comp_id in self.__std_comp_ids:
            self.loadOtherStatFromCsvFiles()

        if not comp_id in self.__all_comp_ids:
            self.__appendExtraFromCcd(comp_id)

        cs_stat = self.__get(comp_id)

        if comp_id in self.__std_comp_ids or primary:
            return [i['atom_id'] for i in cs_stat if 'methyl' in i['desc'] and
                    (not excl_minor_atom or (excl_minor_atom and i['primary']))]

        return [i['atom_id'] for i in cs_stat if 'methyl' in i['desc'] and
                (not excl_minor_atom or not 'secondary' in i or (excl_minor_atom and i['secondary']))]

    def getRepresentativeMethylProtons(self, comp_id, excl_minor_atom=False, primary=False):
        """ Return representative protons in methyl group of a given comp_id.
        """

        ends_w_num = [a for a in self.getMethylAtoms(comp_id, excl_minor_atom, primary) if a.startswith('H') and a[-1].isdigit()]
        ends_w_alp = [a for a in self.getMethylAtoms(comp_id, excl_minor_atom, primary) if a.startswith('H') and not a[-1].isdigit()]

        atm_list = []

        if len(ends_w_num) > 0:
            atm_list.extend([a for a in ends_w_num if a.endswith('1')])

        if len(ends_w_alp) > 0:
            min_len = min([len(a) for a in ends_w_alp])
            atm_list.extend([a for a in ends_w_alp if len(a) == min_len])

        return atm_list

    def getNonRepresentativeMethylProtons(self, comp_id, excl_minor_atom=False, primary=False):
        """ Return non-representative protons in methyl group of a given comp_id.
        """

        rep_list = self.getRepresentativeMethylProtons(comp_id, excl_minor_atom, primary)

        return [a for a in self.getMethylAtoms(comp_id, excl_minor_atom, primary) if a.startswith('H') and not a in rep_list]

    def getSideChainAtoms(self, comp_id, excl_minor_atom=False, polypeptide_like=False, polynucleotide_like=False, carbohydrates_like=False):
        """ Return sidechain atoms of a given comp_id.
        """

        if not comp_id in self.__std_comp_ids:
            self.loadOtherStatFromCsvFiles()

        if not comp_id in self.__all_comp_ids:
            self.__appendExtraFromCcd(comp_id)

        if polypeptide_like == False and polynucleotide_like == False and carbohydrates_like == False:
            polypeptide_like, polynucleotide_like, carbohydrates_like = self.getTypeOfCompId(comp_id)

        bb_atoms = self.getBackBoneAtoms(comp_id, excl_minor_atom, polypeptide_like, polynucleotide_like, carbohydrates_like)

        try:
            if polypeptide_like:
                bb_atoms.remove('CB')
        except ValueError:
            pass

        cs_stat = self.__get(comp_id)

        if comp_id in self.__std_comp_ids or polypeptide_like:
            return [i['atom_id'] for i in cs_stat if not i['atom_id'] in bb_atoms and
                    (not excl_minor_atom or (excl_minor_atom and i['primary']))]

        return [i['atom_id'] for i in cs_stat if not i['atom_id'] in bb_atoms and
                (not excl_minor_atom or not 'secondary' in i or (excl_minor_atom and i['secondary']))]

    def printStat(self, atm_list):
        """ Print out BMRB chemical shift statistics.
        """

        for i in atm_list:
            print (i)

    def loadStatFromCsvFiles(self):
        """ Load all BMRB chemical shift statistics from CSV files.
        """

        file_name_list = [self.stat_dir + 'aa_filt.csv', self.stat_dir + 'aa_full.csv',
                          self.stat_dir + 'dna_filt.csv', self.stat_dir + 'dna_full.csv',
                          self.stat_dir + 'rna_filt.csv', self.stat_dir + 'rna_full.csv',
                          self.stat_dir + 'others.csv']

        if any(not os.path.exists(file_name) for file_name in file_name_list):
            return False

        self.aa_filt = self.loadStatFromCsvFile(self.stat_dir + 'aa_filt.csv', self.aa_threshold)
        self.aa_full = self.loadStatFromCsvFile(self.stat_dir + 'aa_full.csv', self.aa_threshold)

        self.dna_filt = self.loadStatFromCsvFile(self.stat_dir + 'dna_filt.csv', self.na_threshold)
        self.dna_full = self.loadStatFromCsvFile(self.stat_dir + 'dna_full.csv', self.na_threshold)

        self.rna_filt = self.loadStatFromCsvFile(self.stat_dir + 'rna_filt.csv', self.na_threshold)
        self.rna_full = self.loadStatFromCsvFile(self.stat_dir + 'rna_full.csv', self.na_threshold)

        if not self.lazy_others:
            self.others = self.loadStatFromCsvFile(self.stat_dir + 'others.csv', self.aa_threshold, self.na_threshold)

        self.__updateCompIdSet()

        return True

    def loadOtherStatFromCsvFiles(self):
        """ Load all BMRB chemical shift statistics from CSV files.
        """

        if len(self.others) > 0 or not self.lazy_others:
            return

        self.others = self.loadStatFromCsvFile(self.stat_dir + 'others.csv', self.aa_threshold, self.na_threshold)

        self.__updateCompIdSet()

        return True

    def loadStatFromCsvFile(self, file_name, primary_th, secondary_th=None):
        """ Load BMRB chemical shift statistics from a given CSV file.
        """

        atm_list = []

        with open(file_name, 'r') as f:
            reader = csv.DictReader(f)

            for row in reader:

                comp_id = row['comp_id']

                if not self.__updateChemCompDict(comp_id):
                    continue

                _atom_id = row['atom_id']

                # methyl proton group
                if _atom_id.startswith('M'):
                    _atom_id = re.sub(r'^M', 'H', _atom_id)

                    for i in range(1, 4):
                        _row = {}
                        _row['comp_id'] = comp_id
                        _row['atom_id'] = _atom_id + str(i)

                        if not self.__checkAtomNomenclature(_row['atom_id']):
                            continue

                        _row['count'] = int(row['count'])
                        _row['avg'] = float(row['avg'])
                        try:
                            _row['std'] = float(row['std'])
                        except ValueError:
                            _row['std'] = None
                        _row['min'] = float(row['min'])
                        _row['max'] = float(row['max'])
                        _row['desc'] = 'methyl'
                        _row['primary'] = False
                        _row['norm_freq'] = None

                        atm_list.append(_row)

                elif comp_id == 'HEM' and not re.match(r'^HM[A-D]$', _atom_id) is None: # others.csv dependent code

                    for i in ['', 'A', 'B']:
                        _row = {}
                        _row['comp_id'] = comp_id
                        _row['atom_id'] = _atom_id + i

                        if not self.__checkAtomNomenclature(_row['atom_id']):
                            continue

                        _row['count'] = int(row['count'])
                        _row['avg'] = float(row['avg'])
                        try:
                            _row['std'] = float(row['std'])
                        except ValueError:
                            _row['std'] = None
                        _row['min'] = float(row['min'])
                        _row['max'] = float(row['max'])
                        _row['desc'] = 'methyl'
                        _row['primary'] = False
                        _row['norm_freq'] = None

                        atm_list.append(_row)

                elif comp_id == 'HEB' and (not re.match(r'^HM[A-D]1$', _atom_id) is None or _atom_id == 'HBB1'): # others.csv dependent code

                    for i in range(1, 4):
                        _row = {}
                        _row['comp_id'] = comp_id
                        _row['atom_id'] = _atom_id[:-1] + str(i)

                        if not self.__checkAtomNomenclature(_row['atom_id']):
                            continue

                        _row['count'] = int(row['count'])
                        _row['avg'] = float(row['avg'])
                        try:
                            _row['std'] = float(row['std'])
                        except ValueError:
                            _row['std'] = None
                        _row['min'] = float(row['min'])
                        _row['max'] = float(row['max'])
                        _row['desc'] = 'methyl'
                        _row['primary'] = False
                        _row['norm_freq'] = None

                        atm_list.append(_row)

                elif comp_id == 'HEC' and (not re.match(r'^HM[A-D]$', _atom_id) is None or not re.match(r'^HB[BC]$', _atom_id) is None): # others.csv dependent code

                    for i in range(1, 4):
                        _row = {}
                        _row['comp_id'] = comp_id
                        _row['atom_id'] = _atom_id + str(i)

                        if not self.__checkAtomNomenclature(_row['atom_id']):
                            continue

                        _row['count'] = int(row['count'])
                        _row['avg'] = float(row['avg'])
                        try:
                            _row['std'] = float(row['std'])
                        except ValueError:
                            _row['std'] = None
                        _row['min'] = float(row['min'])
                        _row['max'] = float(row['max'])
                        _row['desc'] = 'methyl'
                        _row['primary'] = False
                        _row['norm_freq'] = None

                        atm_list.append(_row)

                # geminal proton group
                elif _atom_id.startswith('Q'):
                    _atom_id = re.sub(r'^Q', 'H', _atom_id)

                    for i in range(1, 4):
                        _row = {}
                        _row['comp_id'] = comp_id
                        _row['atom_id'] = _atom_id + str(i)

                        if not self.__checkAtomNomenclature(_row['atom_id']):
                            continue

                        _row['count'] = int(row['count'])
                        _row['avg'] = float(row['avg'])
                        try:
                            _row['std'] = float(row['std'])
                        except ValueError:
                            _row['std'] = None
                        _row['min'] = float(row['min'])
                        _row['max'] = float(row['max'])
                        _row['desc'] = 'methyl'
                        _row['primary'] = False
                        _row['norm_freq'] = None

                        atm_list.append(_row)

                elif not((comp_id == 'HEM' and not re.match(r'^HM[A-D][AB]$', _atom_id) is None) or\
                         (comp_id == 'HEB' and (not re.match(r'^HM[A-D][23]$', _atom_id) is None or not re.match(r'^HBB[23]', _atom_id) is None)) or\
                         (comp_id == 'HEC' and (not re.match(r'^HM[A-D][123]$', _atom_id) is None or not re.match(r'^HB[BC][123]$', _atom_id) is None))):
                    _row = {}
                    _row['comp_id'] = comp_id
                    _row['atom_id'] = _atom_id

                    if not self.__checkAtomNomenclature(_row['atom_id']):
                        continue

                    _row['count'] = int(row['count'])
                    _row['avg'] = float(row['avg'])
                    try:
                        _row['std'] = float(row['std'])
                    except ValueError:
                        _row['std'] = None
                    _row['min'] = float(row['min'])
                    _row['max'] = float(row['max'])
                    _row['desc'] = 'isolated'
                    _row['primary'] = False
                    _row['norm_freq'] = None

                    atm_list.append(_row)

        comp_ids = set([i['comp_id'] for i in atm_list])

        if not secondary_th is None: # extract rest of atoms for non-standard residues

            for comp_id in comp_ids:

                if self.__updateChemCompDict(comp_id):

                    for a in self.__last_chem_comp_atoms:

                        if a[self.__cca_leaving_atom_flag] == 'Y' or not a[self.__cca_type_symbol] in ['H', 'C', 'N', 'P']:
                            continue

                        if not any(i for i in atm_list if i['comp_id'] == comp_id and i['atom_id'] == a[self.__cca_atom_id]):

                            _row = {}
                            _row['comp_id'] = comp_id
                            _row['atom_id'] = a[self.__cca_atom_id]
                            _row['desc'] = 'isolated'
                            _row['primary'] = False
                            _row['norm_freq'] = None
                            _row['count'] = 0

                            atm_list.append(_row)

        self.__detectMethylProtonFromAtomNomenclature(comp_ids, atm_list)
        self.__detectGeminalProtonFromAtomNomenclature(comp_ids, atm_list)

        self.__detectGeminalCarbon(comp_ids, atm_list)
        self.__detectGeminalNitrogen(comp_ids, atm_list)

        self.__detectMajorResonance(comp_ids, atm_list, primary_th, secondary_th)

        return atm_list

    def __appendExtraFromCcd(self, comp_id):
        """ Append atom list as extra residue for a given comp_id.
        """

        if comp_id in self.__all_comp_ids or comp_id in self.__ext_comp_ids or not self.__updateChemCompDict(comp_id):
            return

        atm_list = []

        for a in self.__last_chem_comp_atoms:

            if a[self.__cca_leaving_atom_flag] == 'Y' or not a[self.__cca_type_symbol] in ['H', 'C', 'N', 'P']:
                continue

            _row = {}
            _row['comp_id'] = comp_id
            _row['atom_id'] = a[self.__cca_atom_id]
            _row['desc'] = 'isolated'
            _row['primary'] = False
            _row['norm_freq'] = None

            atm_list.append(_row)

        self.__ext_comp_ids.add(comp_id)

        if len(atm_list) == 0:
            return

        comp_ids = [comp_id]

        self.__detectMethylProtonFromAtomNomenclature(comp_ids, atm_list)
        self.__detectGeminalProtonFromAtomNomenclature(comp_ids, atm_list)

        self.__detectGeminalCarbon(comp_ids, atm_list)
        self.__detectGeminalNitrogen(comp_ids, atm_list)

        self.extras.extend(atm_list)

    def __updateChemCompDict(self, comp_id):
        """ Update CCD information for a given comp_id.
            @return: True for successfully update CCD information or False for the case a given comp_id does not exist in CCD
        """

        comp_id = comp_id.upper()

        if comp_id != self.__last_comp_id:
            self.__last_comp_id_test = False if '_' in comp_id else comp_id
            self.__last_comp_id = comp_id

            if self.__last_comp_id_test:
                self.__last_chem_comp_dict = self.__ccR.getChemCompDict()
                self.__last_chem_comp_atoms = self.__ccR.getAtomList()
                self.__last_chem_comp_bonds = self.__ccR.getBonds()

        return self.__last_comp_id_test

    def __checkAtomNomenclature(self, atom_id):
        """ Check atom nomenclature.
        """

        if any(a[self.__cca_atom_id] for a in self.__last_chem_comp_atoms if a[self.__cca_atom_id] == atom_id and a[self.__cca_leaving_atom_flag] != 'Y'):
            return True

        else:
            if self.__verbose:
                self.__lfh.write("+BMRBChemShiftStat.__checkAtomNomenclature() ++ Error  - Invalid atom nomenclature %s, comp_id %s\n" % (atom_id, self.__last_comp_id))
            return False

    def __detectMethylProtonFromAtomNomenclature(self, comp_ids, atm_list):
        """ Detect methyl proton from atom nomenclature.
        """

        for comp_id in comp_ids:
            _list = [i for i in atm_list if i['comp_id'] == comp_id]

            h_list = [i for i in _list if i['atom_id'].startswith('H') and i['desc'] == 'isolated']

            if self.__updateChemCompDict(comp_id):
                c_h_bonds = collections.Counter([b[self.__ccb_atom_id_1] for b in self.__last_chem_comp_bonds if b[self.__ccb_atom_id_1].startswith('C') and b[self.__ccb_atom_id_2].startswith('H')])

                for k, v in c_h_bonds.items():
                    if v == 3:

                        for i in _list:
                            if i['atom_id'] == k:
                                i['desc'] = 'methyl'

                        for i in h_list:
                            atom_id = i['atom_id']
                            if any(b for b in self.__last_chem_comp_bonds if b[self.__ccb_atom_id_1] == k and b[self.__ccb_atom_id_2] == atom_id):
                                i['desc'] = 'methyl'

            else:
                h_1 = [i['atom_id'][:-1] for i in h_list if i['atom_id'].endswith('1')]
                h_2 = [i['atom_id'][:-1] for i in h_list if i['atom_id'].endswith('2')]
                h_3 = [i['atom_id'][:-1] for i in h_list if i['atom_id'].endswith('3')]
                h_4 = [i['atom_id'][:-1] for i in h_list if i['atom_id'].endswith('4')]

                h_common = set(h_1) & set(h_2) & set(h_3) - set(h_4)

                for h in h_common:
                    for i in h_list:
                        atom_id = i['atom_id']
                        if atom_id == h + '1' or atom_id == h + '2' or atom_id == h + '3':
                            i['desc'] = 'methyl'

    def __detectGeminalProtonFromAtomNomenclature(self, comp_ids, atm_list):
        """ Detect geminal proton from atom nomenclature.
        """

        for comp_id in comp_ids:
            _list = [i for i in atm_list if i['comp_id'] == comp_id]

            h_list = [i for i in _list if i['atom_id'].startswith('H') and i['desc'] == 'isolated']

            if self.__updateChemCompDict(comp_id):
                aro_list = [a[self.__cca_atom_id] for a in self.__last_chem_comp_atoms if a[self.__cca_aromatic_flag] == 'Y']

                for i in _list:
                    if i['atom_id'] in aro_list:
                        i['desc'] = 'aroma'

                for aro in aro_list:
                    for i in h_list:
                        atom_id = i['atom_id']
                        if any(b for b in self.__last_chem_comp_bonds if b[self.__ccb_atom_id_1] == aro and b[self.__ccb_atom_id_2] == atom_id):
                            i['desc'] = 'aroma'

                leaving_atom_list = [a[self.__cca_atom_id] for a in self.__last_chem_comp_atoms if a[self.__cca_leaving_atom_flag] == 'Y']

                cn_h_bonds = collections.Counter([b[self.__ccb_atom_id_1] for b in self.__last_chem_comp_bonds if b[self.__ccb_atom_id_2].startswith('H') and not b[self.__ccb_atom_id_2] in leaving_atom_list])

                h_list = [i for i in _list if i['atom_id'].startswith('H') and i['desc'] == 'isolated']

                for k, v in cn_h_bonds.items():
                    if v == 2:
                        for i in h_list:
                            atom_id = i['atom_id']
                            if any(b for b in self.__last_chem_comp_bonds if b[self.__ccb_atom_id_1] == k and b[self.__ccb_atom_id_2] == atom_id):
                                i['desc'] = 'geminal'

                h_list = [i for i in _list if i['atom_id'].startswith('H') and i['desc'] == 'aroma']

                hvy_c_list = []

                pair = 0
                for h_1 in h_list:
                    if h_1['atom_id'][-1] in ['1', '2', '3']:
                        hvy_1 = next(b[self.__ccb_atom_id_1] for b in self.__last_chem_comp_bonds if b[self.__ccb_atom_id_2] == h_1['atom_id'])
                        for h_2 in h_list:
                            if h_2['atom_id'][-1] in ['1', '2', '3'] and h_list.index(h_1) < h_list.index(h_2):
                                hvy_2 = next(b[self.__ccb_atom_id_1] for b in self.__last_chem_comp_bonds if b[self.__ccb_atom_id_2] == h_2['atom_id'])
                                if hvy_1[:-1] == hvy_2[:-1]:
                                    hvy_1_c = set([b[self.__ccb_atom_id_1] for b in self.__last_chem_comp_bonds if b[self.__ccb_aromatic_flag] == 'Y' and b[self.__ccb_atom_id_2] == hvy_1]) |\
                                              set([b[self.__ccb_atom_id_2] for b in self.__last_chem_comp_bonds if b[self.__ccb_aromatic_flag] == 'Y' and b[self.__ccb_atom_id_1] == hvy_1])
                                    hvy_2_c = set([b[self.__ccb_atom_id_1] for b in self.__last_chem_comp_bonds if b[self.__ccb_aromatic_flag] == 'Y' and b[self.__ccb_atom_id_2] == hvy_2]) |\
                                              set([b[self.__ccb_atom_id_2] for b in self.__last_chem_comp_bonds if b[self.__ccb_aromatic_flag] == 'Y' and b[self.__ccb_atom_id_1] == hvy_2])
                                    set_hvy_c = hvy_1_c & hvy_2_c
                                    if len(set_hvy_c) == 1:
                                        for hvy_c in set_hvy_c:
                                            hvy_c_list.append(hvy_c)
                                        pair += 1

                if pair > 0 and pair % 2 == 0:
                    hvy_c_set_in_ring = set()
                    for hvy_c_1 in hvy_c_list:
                        for hvy_c_2 in hvy_c_list:
                            if hvy_c_1 < hvy_c_2:
                                hvy_set_1 = set([b[self.__ccb_atom_id_1] for b in self.__last_chem_comp_bonds if b[self.__ccb_aromatic_flag] == 'Y' and b[self.__ccb_atom_id_2] == hvy_c_1]) |\
                                            set([b[self.__ccb_atom_id_2] for b in self.__last_chem_comp_bonds if b[self.__ccb_aromatic_flag] == 'Y' and b[self.__ccb_atom_id_1] == hvy_c_1])
                                hvy_set_2 = set([b[self.__ccb_atom_id_1] for b in self.__last_chem_comp_bonds if b[self.__ccb_aromatic_flag] == 'Y' and b[self.__ccb_atom_id_2] == hvy_c_2]) |\
                                            set([b[self.__ccb_atom_id_2] for b in self.__last_chem_comp_bonds if b[self.__ccb_aromatic_flag] == 'Y' and b[self.__ccb_atom_id_1] == hvy_c_2])
                                in_ring = False
                                for hvy_1 in hvy_set_1:
                                    if in_ring:
                                        break
                                    for hvy_2 in hvy_set_2:
                                        if in_ring:
                                            break
                                        if any(b for b in self.__last_chem_comp_bonds if b[self.__ccb_aromatic_flag] == 'Y' and b[self.__ccb_atom_id_1] == hvy_1 and b[self.__ccb_atom_id_2] == hvy_2) or\
                                           any(b for b in self.__last_chem_comp_bonds if b[self.__ccb_aromatic_flag] == 'Y' and b[self.__ccb_atom_id_1] == hvy_2 and b[self.__ccb_atom_id_2] == hvy_1):
                                            in_ring = True
                                if in_ring:
                                    hvy_c_set_in_ring.add(hvy_c_1)
                                    hvy_c_set_in_ring.add(hvy_c_2)

                    for h_1 in h_list:
                        if h_1['atom_id'][-1] in ['1', '2', '3']:
                            hvy_1 = next(b[self.__ccb_atom_id_1] for b in self.__last_chem_comp_bonds if b[self.__ccb_atom_id_2] == h_1['atom_id'])
                            for h_2 in h_list:
                                if h_2['atom_id'][-1] in ['1', '2', '3'] and h_list.index(h_1) < h_list.index(h_2):
                                    hvy_2 = next(b[self.__ccb_atom_id_1] for b in self.__last_chem_comp_bonds if b[self.__ccb_atom_id_2] == h_2['atom_id'])
                                    if hvy_1[:-1] == hvy_2[:-1]:
                                        hvy_1_c = set([b[self.__ccb_atom_id_1] for b in self.__last_chem_comp_bonds if b[self.__ccb_aromatic_flag] == 'Y' and b[self.__ccb_atom_id_2] == hvy_1]) |\
                                                  set([b[self.__ccb_atom_id_2] for b in self.__last_chem_comp_bonds if b[self.__ccb_aromatic_flag] == 'Y' and b[self.__ccb_atom_id_1] == hvy_1])
                                        hvy_2_c = set([b[self.__ccb_atom_id_1] for b in self.__last_chem_comp_bonds if b[self.__ccb_aromatic_flag] == 'Y' and b[self.__ccb_atom_id_2] == hvy_2]) |\
                                                  set([b[self.__ccb_atom_id_2] for b in self.__last_chem_comp_bonds if b[self.__ccb_aromatic_flag] == 'Y' and b[self.__ccb_atom_id_1] == hvy_2])
                                        if len(hvy_1_c & hvy_2_c & hvy_c_set_in_ring) > 0:
                                            h_1['desc'] = 'aroma-opposite'
                                            h_2['desc'] = 'aroma-opposite'
                                            for i in _list:
                                                if i['atom_id'] == hvy_1 or i['atom_id'] == hvy_2:
                                                    i['desc'] = 'aroma-opposite'

            else:
                h_1 = [i['atom_id'][:-1] for i in h_list if i['atom_id'].endswith('1')]
                h_2 = [i['atom_id'][:-1] for i in h_list if i['atom_id'].endswith('2')]
                h_3 = [i['atom_id'][:-1] for i in h_list if i['atom_id'].endswith('3')]

                c_list = [i for i in _list if i['atom_id'].startswith('C')]

                c_1 = ['H' + i['atom_id'][1:-1] for i in c_list if i['atom_id'].endswith('1')]
                c_2 = ['H' + i['atom_id'][1:-1] for i in c_list if i['atom_id'].endswith('2')]
                c_3 = ['H' + i['atom_id'][1:-1] for i in c_list if i['atom_id'].endswith('3')]

                n_list = [i for i in _list if i['atom_id'].startswith('N')]

                n_1 = ['H' + i['atom_id'][1:-1] for i in n_list if i['atom_id'].endswith('1')]
                n_2 = ['H' + i['atom_id'][1:-1] for i in n_list if i['atom_id'].endswith('2')]
                n_3 = ['H' + i['atom_id'][1:-1] for i in n_list if i['atom_id'].endswith('3')]

                h_common = set(h_1) & set(h_2) - set(h_3)
                cn_common = set(c_1) & set(c_2) | set(c_1) & set(n_2) | set(n_1) & set(c_2)

                for h in h_common:
                    for i in h_list:
                        atom_id = i['atom_id']
                        if atom_id == h + '1' or atom_id == h + '2':
                            atom_id = 'N' + i['atom_id'][1:]
                            if not any(n for n in n_list if n['atom_id'] == atom_id):
                                i['desc'] = 'aroma' if h in cn_common and i['avg'] > 5.0 else 'geminal'

                h_common = set(h_2) & set(h_3) - set(h_1)
                cn_common = set(c_2) & set(c_3) | set(c_2) & set(n_3) | set(n_2) & set(c_3)

                for h in h_common:
                    for i in h_list:
                        atom_id = i['atom_id']
                        if atom_id == h + '2' or atom_id == h + '3':
                            atom_id = 'N' + i['atom_id'][1:]
                            if not any(n for n in n_list if n['atom_id'] == atom_id):
                                i['desc'] = 'aroma' if h in cn_common and i['avg'] > 5.0 else 'geminal'

                h_list = [i for i in _list if i['atom_id'].startswith('H') and i['desc'] == 'isolated']

                for h in h_list:
                    if h['avg'] > 5.0:
                        atom_id = 'C' + h['atom_id'][1:]
                        if any(c for c in c_list if c['atom_id'] == atom_id and c['avg'] > 95.0 and c['avg'] < 170.0):
                            h['desc'] = 'aroma'

                h_list = [i for i in _list if i['atom_id'].startswith('H') and i['desc'] == 'isolated']

                h_c = [i['atom_id'][:-1] for i in h_list if i['atom_id'].endswith("'") and not i['atom_id'].endswith("''")]
                h_cc = [i['atom_id'][:-2] for i in h_list if i['atom_id'].endswith("''")]

                c_c = ['H' + i['atom_id'][1:-1] for i in c_list if i['atom_id'].endswith("'") and not i['atom_id'].endswith("''")]
                c_cc = ['H' + i['atom_id'][1:-2] for i in c_list if i['atom_id'].endswith("''")]

                h_common = set(h_c) & set(h_cc) & set(c_c) - set(c_cc)

                for h in h_common:
                    for i in h_list:
                        atom_id = i['atom_id']
                        if atom_id == h + "'" or atom_id == h + "''":
                            atom_id = 'N' + i['atom_id'][1:]
                            if not any(n for n in n_list if n['atom_id'] == atom_id):
                                i['desc'] = 'geminal'

                h_list = [i for i in _list if i['atom_id'].startswith('H') and i['desc'] == 'aroma']

                h_1 = [i['atom_id'][:-1] for i in h_list if i['atom_id'].endswith('1')]
                h_2 = [i['atom_id'][:-1] for i in h_list if i['atom_id'].endswith('2')]
                h_3 = [i['atom_id'][:-1] for i in h_list if i['atom_id'].endswith('3')]

                h_common = set(h_1) & set(h_2)

                if len(h_common) > 0 and len(h_common) % 2 == 0:
                    for h in h_common:
                        for i in h_list:
                            atom_id = i['atom_id']
                            if atom_id == h + '1' or atom_id == h + '2':
                                i['desc'] = 'aroma-opposite'

                h_common = set(h_2) & set(h_3)

                if len(h_common) > 0 and len(h_common) % 2 == 0:
                    for h in h_common:
                        for i in h_list:
                            atom_id = i['atom_id']
                            if atom_id == h + '2' or atom_id == h + '3':
                                i['desc'] = 'aroma-opposite'

    def __detectGeminalCarbon(self, comp_ids, atm_list):
        """ Detect geminal carbon from atom nomenclature.
        """

        for comp_id in comp_ids:
            _list = [i for i in atm_list if i['comp_id'] == comp_id]

            if self.__updateChemCompDict(comp_id):
                methyl_c_list = [i['atom_id'] for i in _list if i['atom_id'].startswith('C') and i['desc'] == 'methyl']

                for methyl_c_1 in methyl_c_list:
                    for methyl_c_2 in methyl_c_list:
                        if methyl_c_list.index(methyl_c_1) < methyl_c_list.index(methyl_c_2):
                            hvy_1_c = set([b[self.__ccb_atom_id_2] for b in self.__last_chem_comp_bonds if b[self.__ccb_atom_id_1] == methyl_c_1 and not b[self.__ccb_atom_id_2].startswith('H')]) |\
                                      set([b[self.__ccb_atom_id_1] for b in self.__last_chem_comp_bonds if b[self.__ccb_atom_id_2] == methyl_c_1 and not b[self.__ccb_atom_id_1].startswith('H')])
                            hvy_2_c = set([b[self.__ccb_atom_id_2] for b in self.__last_chem_comp_bonds if b[self.__ccb_atom_id_1] == methyl_c_2 and not b[self.__ccb_atom_id_2].startswith('H')]) |\
                                      set([b[self.__ccb_atom_id_1] for b in self.__last_chem_comp_bonds if b[self.__ccb_atom_id_2] == methyl_c_2 and not b[self.__ccb_atom_id_1].startswith('H')])
                            hvy_common = hvy_1_c & hvy_2_c
                            if len(hvy_common) > 0:
                                for hvy_c in hvy_common:
                                    v = len([b for b in self.__last_chem_comp_bonds if b[self.__ccb_atom_id_1] == hvy_c and b[self.__ccb_atom_id_2] in methyl_c_list]) +\
                                        len([b for b in self.__last_chem_comp_bonds if b[self.__ccb_atom_id_2] == hvy_c and b[self.__ccb_atom_id_1] in methyl_c_list])

                                    if v == 2:
                                        for i in _list:
                                            if i['atom_id'] == methyl_c_1 or i['atom_id'] == methyl_c_2:
                                                i['desc'] = 'methyl-geminal'

                                                for methyl_h in [b[self.__ccb_atom_id_2] for b in self.__last_chem_comp_bonds if b[self.__ccb_atom_id_1] == methyl_c_1 and b[self.__ccb_atom_id_2].startswith('H')]:
                                                    for j in _list:
                                                        if j['atom_id'] == methyl_h:
                                                            j['desc'] = 'methyl-geminal'

                                                for methyl_h in [b[self.__ccb_atom_id_2] for b in self.__last_chem_comp_bonds if b[self.__ccb_atom_id_1] == methyl_c_2 and b[self.__ccb_atom_id_2].startswith('H')]:
                                                    for j in _list:
                                                        if j['atom_id'] == methyl_h:
                                                            j['desc'] = 'methyl-geminal'

            else:
                methyl_list = ['C' + i['atom_id'][1:-1] for i in _list if i['atom_id'].startswith('H') and i['desc'] == 'methyl' and i['atom_id'].endswith('1')]

                methyl_1 = [i[:-1] for i in methyl_list if i.endswith('1')]
                methyl_2 = [i[:-1] for i in methyl_list if i.endswith('2')]
                methyl_3 = [i[:-1] for i in methyl_list if i.endswith('3')]

                c_list = [i for i in _list if i['atom_id'].startswith('C')]

                for c in c_list:
                    if c['atom_id'] in methyl_list:
                        c['desc'] = 'methyl'
                    elif c['avg'] > 95.0 and c['avg'] < 170.0:
                        c['desc'] = 'aroma'

                methyl_common = set(methyl_1) & set(methyl_2) - set(methyl_3)

                for m in methyl_common:
                    for c in c_list:
                        atom_id = c['atom_id']
                        if atom_id == m + '1' or atom_id == m + '2':
                            c['desc'] = 'methyl-geminal'
                            for h in [i for i in _list if i['atom_id'].startswith('H') and i['desc'] == 'methyl' and i['atom_id'].startswith('H' + atom_id[1:])]:
                                h['desc'] = 'methyl-geminal'

                aroma_list = ['C' + i['atom_id'][1:] for i in _list if i['atom_id'].startswith('H') and i['desc'] == 'aroma']

                for c in c_list:
                    if c['atom_id'] in aroma_list:
                        c['desc'] = 'aroma'

                aroma_opposite_list = ['C' + i['atom_id'][1:] for i in _list if i['atom_id'].startswith('H') and i['desc'] == 'aroma-opposite']

                for c in c_list:
                    if c['atom_id'] in aroma_opposite_list:
                        c['desc'] = 'aroma-opposite'

    def __detectGeminalNitrogen(self, comp_ids, atm_list):
        """ Detect geminal nitrogen from atom nomenclature.
        """

        for comp_id in comp_ids:
            _list = [i for i in atm_list if i['comp_id'] == comp_id]

            geminal_n_list = ['N' + i['atom_id'][1:-1] for i in _list if i['atom_id'].startswith('H') and i['desc'] == 'geminal' and i['atom_id'].endswith('1')]

            if self.__updateChemCompDict(comp_id):
                for geminal_n_1 in geminal_n_list:
                    for geminal_n_2 in geminal_n_list:
                        if geminal_n_list.index(geminal_n_1) < geminal_n_list.index(geminal_n_2):
                            hvy_1_c = set([b[self.__ccb_atom_id_2] for b in self.__last_chem_comp_bonds if b[self.__ccb_atom_id_1] == geminal_n_1 and not b[self.__ccb_atom_id_2].startswith('H')]) |\
                                      set([b[self.__ccb_atom_id_1] for b in self.__last_chem_comp_bonds if b[self.__ccb_atom_id_2] == geminal_n_1 and not b[self.__ccb_atom_id_1].startswith('H')])
                            hvy_2_c = set([b[self.__ccb_atom_id_2] for b in self.__last_chem_comp_bonds if b[self.__ccb_atom_id_1] == geminal_n_2 and not b[self.__ccb_atom_id_2].startswith('H')]) |\
                                      set([b[self.__ccb_atom_id_1] for b in self.__last_chem_comp_bonds if b[self.__ccb_atom_id_2] == geminal_n_2 and not b[self.__ccb_atom_id_1].startswith('H')])
                            hvy_common = hvy_1_c & hvy_2_c
                            if len(hvy_common) > 0:
                                for hvy_c in hvy_common:
                                    v = len([b for b in self.__last_chem_comp_bonds if b[self.__ccb_atom_id_1] == hvy_c and b[self.__ccb_atom_id_2] in geminal_n_list]) +\
                                        len([b for b in self.__last_chem_comp_bonds if b[self.__ccb_atom_id_2] == hvy_c and b[self.__ccb_atom_id_1] in geminal_n_list])

                                    if v == 2:
                                        for i in _list:
                                            if i['atom_id'] == geminal_n_1 or i['atom_id'] == geminal_n_2:
                                                i['desc'] = 'geminal'

            else:
                geminal_n_1 = [i[:-1] for i in geminal_n_list if i.endswith('1')]
                geminal_n_2 = [i[:-1] for i in geminal_n_list if i.endswith('2')]
                geminal_n_3 = [i[:-1] for i in geminal_n_list if i.endswith('3')]

                n_list = [i for i in _list if i['atom_id'].startswith('N')]

                for n in n_list:
                    if n['avg'] > 125.0:
                        n['desc'] = 'aroma'
                        atom_id = 'H' + n['atom_id'][1:]
                        try:
                            h = next(i for i in _list if i['atom_id'] == atom_id and i['desc'] == 'isolated')
                            h['desc'] = 'aroma'
                        except StopIteration:
                            pass

                geminal_common = set(geminal_n_1) & set(geminal_n_2) - set(geminal_n_3)

                for g in geminal_common:
                    for n in n_list:
                        atom_id = n['atom_id']
                        if atom_id == g + '1' or atom_id == g + '2':
                            n['desc'] = 'geminal'

    def __detectMajorResonance(self, comp_ids, atm_list, primary_th, secondary_th=None):
        """ Detect major resonance based on count of assignments.
        """

        for comp_id in comp_ids:
            _list = [i for i in atm_list if i['comp_id'] == comp_id]

            max_count = max([i['count'] for i in _list])

            for i in _list:
                i['norm_freq'] = float("%.3f" % (float(i['count']) / max_count))
                if max_count >= self.max_count_th:
                    if i['count'] > max_count * primary_th:
                        i['primary'] = True
                    if (not secondary_th is None) and i['count'] > max_count * secondary_th:
                        i['secondary'] = True

    def writeStatAsPickleFiles(self):
        """ Write all BMRB chemical shift statistics as pickle files.
        """

        self.writeStatAsPickleFile(self.aa_filt, self.stat_dir + 'aa_filt.pkl')
        self.writeStatAsPickleFile(self.aa_full, self.stat_dir + 'aa_full.pkl')

        self.writeStatAsPickleFile(self.dna_filt, self.stat_dir + 'dna_filt.pkl')
        self.writeStatAsPickleFile(self.dna_full, self.stat_dir + 'dna_full.pkl')

        self.writeStatAsPickleFile(self.rna_filt, self.stat_dir + 'rna_filt.pkl')
        self.writeStatAsPickleFile(self.rna_full, self.stat_dir + 'rna_full.pkl')

        self.loadOtherStatFromCsvFiles()

        self.writeStatAsPickleFile(self.others, self.stat_dir + 'others.pkl')

    def writeStatAsPickleFile(self, atm_list, file_name):
        """ Write BMRB chemical shift statistics as pickle file.
        """

        with open(file_name, 'wb') as f:
            pickle.dump(atm_list, f)

    def loadStatFromPickleFiles(self):
        """ Load all BMRB chemical shift statistics from pickle files if possible.
        """

        file_name_list = [self.stat_dir + 'aa_filt.pkl', self.stat_dir + 'aa_full.pkl',
                          self.stat_dir + 'dna_filt.pkl', self.stat_dir + 'dna_full.pkl',
                          self.stat_dir + 'rna_filt.pkl', self.stat_dir + 'rna_full.pkl',
                          self.stat_dir + 'others.pkl']

        for file_name in file_name_list:
            if not os.path.exists(file_name):
                return False

        self.aa_filt = self.loadStatFromPickleFile(self.stat_dir + 'aa_filt.pkl')
        self.aa_full = self.loadStatFromPickleFile(self.stat_dir + 'aa_full.pkl')

        self.dna_filt = self.loadStatFromPickleFile(self.stat_dir + 'dna_filt.pkl')
        self.dna_full = self.loadStatFromPickleFile(self.stat_dir + 'dna_full.pkl')

        self.rna_filt = self.loadStatFromPickleFile(self.stat_dir + 'rna_filt.pkl')
        self.rna_full = self.loadStatFromPickleFile(self.stat_dir + 'rna_full.pkl')

        self.others = self.loadStatFromPickleFile(self.stat_dir + 'others.pkl')

        self.__updateCompIdSet()

        return True

    def loadStatFromPickleFile(self, file_name):
        """ Load BMRB chemical shift statistics from pickle file if possible.
        """

        if os.path.exists(file_name):

            with open(file_name, 'rb') as f:
                return pickle.load(f)

        return []

    def __updateCompIdSet(self):
        """ Update set of comp_id having BMRB chemical shift statistics
        """

        self.__aa_comp_ids = set([i['comp_id'] for i in self.aa_filt])
        self.__dna_comp_ids = set([i['comp_id'] for i in self.dna_filt])
        self.__rna_comp_ids = set([i['comp_id'] for i in self.rna_filt])

        self.__all_comp_ids |= self.__aa_comp_ids
        self.__all_comp_ids |= self.__dna_comp_ids
        self.__all_comp_ids |= self.__rna_comp_ids

        self.__std_comp_ids = copy.copy(self.__all_comp_ids)

        self.__oth_comp_ids = set([i['comp_id'] for i in self.others])

        self.__all_comp_ids |= self.__oth_comp_ids
