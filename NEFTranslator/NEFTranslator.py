##
# File: NEFTranslator.py
# Date:  02-May-2019  K. Baskaran
#

##
import sys
import os
import ntpath
import json
import logging
import re
import csv
import datetime
import pynmrstar
import itertools
import copy

#from pytz import utc
from packaging import version

from BMRBChemShiftStat import BMRBChemShiftStat


__version__ = '2.9.7'

__pynmrstar_v3__ = version.parse(pynmrstar.__version__) >= version.parse("3.0.0")

logging.getLogger().setLevel(logging.ERROR)
"""
def load_json_data(json_file):
    "" Load JSON data to dictionary.
        @param json_file: input JSON file path
        @return: dictionary object
    ""

    data_dict = []

    try:
        with open(json_file, 'r') as file:
            data_dict = json.loads(file.read())

    except Exception as e:
        logging.error(str(e))

    return data_dict
"""
def load_csv_data(csv_file, transpose=False):
    """ Load CSV data to list.
        @param cvs_file: input CSV file path
        @param transpose: transpose CSV data
        @return: list object
    """

    data_map = []

    try:
        data = []

        with open(csv_file, 'r') as file:
            csv_reader = csv.reader(file, delimiter=',')

            for r in csv_reader:
                if r[0][0] != '#':
                    data.append(r)

        data_map = list(map(list, zip(*data))) if transpose else data

    except Exception as e:
        logging.error(str(e))

    return data_map

def get_lp_tag(lp_data, tags):
    """ Return the selected loop tags by row as a list of lists.
    """

    return lp_data.get_tag(tags) if __pynmrstar_v3__ else lp_data.get_data_by_tag(tags)

def get_idx_msg(idx_tag_ids, tags, ent):
    """ Return description about current index.
        @author: Masashi Yokochi
        @return: description
    """

    try:

        idx_msg = ''

        if len(idx_tag_ids) > 0:
            for _j in idx_tag_ids:
                idx_msg += tags[_j] + " " + str(ent[tags[_j]]) + ", "

            idx_msg = '[Check row of ' + idx_msg[:-2] + '] '

        return idx_msg

    except KeyError:
        return ''

class NEFTranslator(object):
    """ Bi-directional translator between NEF and NMR-STAR
    """

    def __init__(self):
        self.lib_dir = os.path.dirname(__file__) + '/lib/'

        self.tagMap = load_csv_data(self.lib_dir + 'NEF_NMRSTAR_equivalence.csv', transpose=True)

        self.nef_mandatory_tag = load_csv_data(self.lib_dir + 'NEF_mandatory.csv')

        self.star_mandatory_tag = load_csv_data(self.lib_dir + 'NMR-STAR_mandatory.csv')
        """
        self.atomDict = load_json_data(self.lib_dir + 'atomDict.json')

        self.codeDict = load_json_data(self.lib_di + 'codeDict.json')
        """
        # whether to replace zero by empty if 'void-zero' is set
        self.replace_zero_by_null_in_case = False
        # whether to insert _Atom_chem_shift.Original_PDB_* items
        self.insert_original_pdb_cs_items = True

        # supported version
        self.nef_version = '1.1'
        self.star_version = '3.2.1.18'

        # format name
        self.nef_format_name = 'nmr_exchange_format'

        # empty value
        self.empty_value = (None, '', '.', '?')

        # true value
        self.true_value = ('true', 't', 'yes', 'y', '1')

        # NEF boolean values
        self.nef_boolean = ('true', 'false')

        # NMR-STAR boolean values
        self.star_boolean = ('yes', 'no')

        # paramagnetic elements, except for Oxygen
        self.paramag_elems = ('LI', 'NA', 'MG', 'AL', 'K', 'CA', 'SC', 'TI', 'V', 'MN', 'RB', 'SR', 'Y', 'ZR', 'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'SN', 'CS', 'BA', 'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU', 'HF', 'TA', 'W', 'RE', 'OS', 'IR', 'PT', 'FR', 'RA', 'AC')

        # ferromagnetic elements
        self.ferromag_elems = ('CR', 'FE', 'CO', 'NI')

        # non-metal elements
        self.non_metal_elems = ('H', 'C', 'N', 'O', 'P', 'S', 'SE')

        # ambiguity codes
        self.bmrb_ambiguity_codes = (1, 2, 3, 4, 5, 6, 9)

        # isotope numbers of NMR observable atoms
        self.atom_isotopes = {'H': [1, 2, 3],
                              'C': [13],
                              'N': [15, 14],
                              'O': [17],
                              'P': [31],
                              'S': [33],
                              'F': [19],
                              'CD': [113, 111],
                              'CA': [43]
                              }

        # list of atom_id of known residues
        self.atomDict = {"ALA": ["N", "CA", "C", "O", "CB", "OXT", "H", "H2", "HA", "HB1", "HB2", "HB3", "HXT"],
                         "ARG": ["N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2", "OXT", "H", "H2", "HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "HE", "HH11", "HH12", "HH21", "HH22", "HXT"],
                         "ASN": ["N", "CA", "C", "O", "CB", "CG", "OD1", "ND2", "OXT", "H", "H2", "HA", "HB2", "HB3", "HD21", "HD22", "HXT"],
                         "ASP": ["N", "CA", "C", "O", "CB", "CG", "OD1", "OD2", "OXT", "H", "H2", "HA", "HB2", "HB3", "HD2", "HXT"],
                         "ASX": ["N", "CA", "C", "O", "CB", "CG", "XD1", "XD2", "OXT", "H", "H2", "HA", "HB1", "HB2", "HXT"],
                         "CYS": ["N", "CA", "C", "O", "CB", "SG", "OXT", "H", "H2", "HA", "HB2", "HB3", "HG", "HXT"],
                         "GLN": ["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2", "OXT", "H", "H2", "HA", "HB2", "HB3", "HG2", "HG3", "HE21", "HE22", "HXT"],
                         "GLU": ["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2", "OXT", "H", "H2", "HA", "HB2", "HB3", "HG2", "HG3", "HE2", "HXT"],
                         "GLX": ["N", "CA", "C", "O", "CB", "CG", "CD", "XE1", "XE2", "HA", "OXT", "HXT", "HB1", "HB2", "HG1", "HG2", "H", "H2"],
                         "GLY": ["N", "CA", "C", "O", "OXT", "H", "H2", "HA2", "HA3", "HXT"],
                         "HIS": ["N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2", "OXT", "H", "H2", "HA", "HB2", "HB3", "HD1", "HD2", "HE1", "HE2", "HXT"],
                         "ILE": ["N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1", "OXT", "H", "H2", "HA", "HB", "HG12", "HG13", "HG21", "HG22", "HG23", "HD11", "HD12", "HD13", "HXT"],
                         "LEU": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "OXT", "H", "H2", "HA", "HB2", "HB3", "HG", "HD11", "HD12", "HD13", "HD21", "HD22", "HD23", "HXT"],
                         "LYS": ["N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ", "OXT", "H", "H2", "HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "HE2", "HE3", "HZ1", "HZ2", "HZ3", "HXT"],
                         "MET": ["N", "CA", "C", "O", "CB", "CG", "SD", "CE", "OXT", "H", "H2", "HA", "HB2", "HB3", "HG2", "HG3", "HE1", "HE2", "HE3", "HXT"],
                         "PHE": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OXT", "H", "H2", "HA", "HB2", "HB3", "HD1", "HD2", "HE1", "HE2", "HZ", "HXT"],
                         "PRO": ["N", "CA", "C", "O", "CB", "CG", "CD", "OXT", "H", "HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "HXT"],
                         "SER": ["N", "CA", "C", "O", "CB", "OG", "OXT", "H", "H2", "HA", "HB2", "HB3", "HG", "HXT"],
                         "THR": ["N", "CA", "C", "O", "CB", "OG1", "CG2", "OXT", "H", "H2", "HA", "HB", "HG1", "HG21", "HG22", "HG23", "HXT"],
                         "TRP": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2", "OXT", "H", "H2", "HA", "HB2", "HB3", "HD1", "HE1", "HE3", "HZ2", "HZ3", "HH2", "HXT"],
                         "TYR": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH", "OXT", "H", "H2", "HA", "HB2", "HB3", "HD1", "HD2", "HE1", "HE2", "HH", "HXT"],
                         "VAL": ["N", "CA", "C", "O", "CB", "CG1", "CG2", "OXT", "H", "H2", "HA", "HB", "HG11", "HG12", "HG13", "HG21", "HG22", "HG23", "HXT"],
                         "DA": ["OP3", "P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'", "N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4", "HOP3", "HOP2", "H5'", "H5''", "H4'", "H3'", "HO3'", "H2'", "H2''", "H1'", "H8", "H61", "H62", "H2"],
                         "DC": ["OP3", "P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'", "N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6", "HOP3", "HOP2", "H5'", "H5''", "H4'", "H3'", "HO3'", "H2'", "H2''", "H1'", "H41", "H42", "H5", "H6"],
                         "DG": ["OP3", "P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'", "N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4", "HOP3", "HOP2", "H5'", "H5''", "H4'", "H3'", "HO3'", "H2'", "H2''", "H1'", "H8", "H1", "H21", "H22"],
                         "DT": ["OP3", "P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'", "N1", "C2", "O2", "N3", "C4", "O4", "C5", "C7", "C6", "HOP3", "HOP2", "H5'", "H5''", "H4'", "H3'", "HO3'", "H2'", "H2''", "H1'", "H3", "H71", "H72", "H73", "H6"],
                         "DU": ["OP3", "P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'", "N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6", "HOP3", "HOP2", "H5'", "H5''", "H4'", "H3'", "HO3'", "H2'", "H2''", "H1'", "H3", "H5", "H6"],
                         "A": ["OP3", "P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", "N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4", "HOP3", "HOP2", "H5'", "H5''", "H4'", "H3'", "HO3'", "H2'", "HO2'", "H1'", "H8", "H61", "H62", "H2"],
                         "C": ["OP3", "P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", "N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6", "HOP3", "HOP2", "H5'", "H5''", "H4'", "H3'", "HO3'", "H2'", "HO2'", "H1'", "H41", "H42", "H5", "H6"],
                         "G": ["OP3", "P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", "N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4", "HOP3", "HOP2", "H5'", "H5''", "H4'", "H3'", "HO3'", "H2'", "HO2'", "H1'", "H8", "H1", "H21", "H22"],
                         "T": ["OP3", "P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'", "N1", "C2", "O2", "N3", "C4", "O4", "C5", "C7", "C6", "HOP3", "HOP2", "H5'", "H5''", "H4'", "H3'", "HO3'", "H2'", "H2''", "H1'", "H3", "H71", "H72", "H73", "H6"],
                         "U": ["OP3", "P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", "N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6", "HOP3", "HOP2", "H5'", "H5''", "H4'", "H3'", "HO3'", "H2'", "HO2'", "H1'", "H3", "H5", "H6"],
                         #
                         "PYL": ["CB2", "CG2", "CD2", "CE2", "N2", "CA2", "C2", "O2", "NZ", "CE", "CD", "CG", "CB", "CA", "C", "OXT", "O", "N", "HXT", "HA", "H", "H2", "HB3", "HB2", "HG3", "HG2", "HD3", "HD2", "HE3", "HE2", "HZ", "HA2", "HE22", "HD32", "HD22", "HG22", "HB12", "HB22", "HB32"],
                         "SEC": ["N", "CA", "CB", "SE", "C", "O", "OXT", "H", "H2", "HA", "HB2", "HB3", "HE", "HXT"],
                         "UNK": ["N", "CA", "C", "O", "CB", "CG", "OXT", "H", "H2", "HA", "HB1", "HB2", "HG1", "HG2", "HG3", "HXT"],
                         #
                         "N": ["OP3", "P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", "HOP3", "HOP2", "H5'", "H5''", "H4'", "H3'", "HO3'", "H2'", "HO2'", "H1'", "H1'2"],
                         }

        # taken from wwpdb.utils.align.SequenceReferenceData.py
        self.monDict3 = {'ALA': 'A',
                         'ARG': 'R',
                         'ASN': 'N',
                         'ASP': 'D',
                         'ASX': 'B',
                         'CYS': 'C',
                         'GLN': 'Q',
                         'GLU': 'E',
                         'GLX': 'Z',
                         'GLY': 'G',
                         'HIS': 'H',
                         'ILE': 'I',
                         'LEU': 'L',
                         'LYS': 'K',
                         'MET': 'M',
                         'PHE': 'F',
                         'PRO': 'P',
                         'SER': 'S',
                         'THR': 'T',
                         'TRP': 'W',
                         'TYR': 'Y',
                         'VAL': 'V',
                         'DA': 'A',
                         'DC': 'C',
                         'DG': 'G',
                         'DT': 'T',
                         'DU': 'U',
                         'DI': 'I',
                         'A': 'A',
                         'C': 'C',
                         'G': 'G',
                         'I': 'I',
                         'T': 'T',
                         'U': 'U'
                         }

        # limit number of dimensions
        self.lim_num_dim = 16

        # BMRB chemical shift statistics
        self.__csStat = BMRBChemShiftStat()

        # CCD accessing utility
        #self.__cI = ConfigInfo(getSiteId())
        #self.__ccCvsPath = self.__cI.get('SITE_CC_CVS_PATH')

        #self.__ccR = ChemCompReader(False, sys.stderr)
        #self.__ccR.setCachePath(self.__ccCvsPath)

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

        # readable item type
        self.readable_item_type = {'str': 'a string',
                                   'bool': 'a boolean value',
                                   'int': 'a integer',
                                   'index-int': 'a unique positive integer',
                                   'positive-int': 'a positive integer',
                                   'positive-int-as-str': 'a positive integer',
                                   'pointer-index': 'a integer acting as a pointer to the parent item',
                                   'float': 'a floating point number',
                                   'positive-float': 'a positive floating point number',
                                   'range-float': 'a floating point number in a specific range',
                                   'enum': 'a enumeration value',
                                   'enum-int': 'a enumeration value restricted to integers'}

        # alternative dictionary of constraint type
        self.dist_alt_constraint_type = {'nef': {'NOE': 'noe',
                                                 'NOE build-up': 'noe_build_up',
                                                 'noe build-up': 'noe_build_up',
                                                 'NOE buildup': 'noe_build_up',
                                                 'noe buildup': 'noe_build_up',
                                                 'NOE build up': 'noe_build_up',
                                                 'noe build up': 'noe_build_up',
                                                 'noe not seen': 'noe_not_seen',
                                                 'ROE': 'roe',
                                                 'roe build-up': 'roe_build_up',
                                                 'ROE buildup': 'roe_build_up',
                                                 'roe buildup': 'roe_build_up',
                                                 'ROE build up': 'roe_build_up',
                                                 'roe build up': 'roe_build_up',
                                                 'roe build-up': 'roe_build_up',
                                                 'hydrogen bond': 'hbond',
                                                 'Hbond': 'hbond',
                                                 'HBond': 'hbond',
                                                 'H-bond': 'hbond',
                                                 'h-bond': 'hbond',
                                                 'H-Bond': 'hbond',
                                                 'Hydrogen bond': 'hbond',
                                                 'disulfide bond': 'disulfide_bond',
                                                 'Disulfide bond': 'disulfide_bond',
                                                 'S-S bond': 'disulfide_bond',
                                                 'SS bond': 'disulfide_bond',
                                                 'SS-bond': 'disulfide_bond',
                                                 'disulfide bridge': 'disulfide_bond',
                                                 'Disulfide bridge': 'disulfide_bond',
                                                 'paramagnetic relaxation': 'pre',
                                                 'PRE': 'pre',
                                                 'Paramagnetic relaxation': 'pre',
                                                 'paramagnetic relaxation enhancement': 'pre',
                                                 'Paramagnetic relaxation enhancement': 'pre',
                                                 'general distance': 'undefined',
                                                 'distance': 'undefined',
                                                 'Mutation': 'mutation',
                                                 'chemical shift perturbation': 'shift_perturbation',
                                                 'shift perturbation': 'shift_perturbation',
                                                 'chem shift perturbation': 'shift_perturbation',
                                                 'CS perturbation': 'shift_perturbation',
                                                 'csp': 'shift_perturbation',
                                                 'CSP': 'shift_perturbation'
                                                 },
                                         'nmr-star': {'noe': 'NOE',
                                                      'noe_build_up': 'NOE build-up',
                                                      'noe build-up': 'NOE build-up',
                                                      'NOE buildup': 'NOE build-up',
                                                      'noe buildup': 'NOE build-up',
                                                      'NOE build up': 'NOE build-up',
                                                      'noe build up': 'NOE build-up',
                                                      'noe_not_seen': 'NOE not seen',
                                                      'noe not seen': 'NOE not seen',
                                                      'roe': 'ROE',
                                                      'roe_build_up': 'ROE build-up',
                                                      'roe build-up': 'ROE build-up',
                                                      'ROE buildup': 'ROE build-up',
                                                      'roe buildup': 'ROE build-up',
                                                      'ROE build up': 'ROE build-up',
                                                      'roe build up': 'ROE build-up',
                                                      'roe build-up': 'ROE build-up',
                                                      'hbond': 'hydrogen bond',
                                                      'Hbond': 'hydrogen bond',
                                                      'HBond': 'hydrogen bond',
                                                      'H-bond': 'hydrogen bond',
                                                      'h-bond': 'hydrogen bond',
                                                      'H-Bond': 'hydrogen bond',
                                                      'Hydrogen bond': 'hydrogen bond',
                                                      'disulfide_bond': 'disulfide bond',
                                                      'Disulfide bond': 'disulfide bond',
                                                      'S-S bond': 'disulfide bond',
                                                      'SS bond': 'disulfide bond',
                                                      'SS-bond': 'disulfide bond',
                                                      'disulfide bridge': 'disulfide bond',
                                                      'Disulfide bridge': 'disulfide bond',
                                                      'PRE': 'paramagnetic relaxation',
                                                      'pre': 'paramagnetic relaxation',
                                                      'Paramagnetic relaxation': 'paramagnetic relaxation',
                                                      'paramagnetic relaxation enhancement': 'paramagnetic relaxation',
                                                      'Paramagnetic relaxation enhancement': 'paramagnetic relaxation',
                                                      'Mutation': 'mutation',
                                                      'unknown': 'general distance',
                                                      'undefined': 'general distance',
                                                      'shift_perturbation': 'chemical shift perturbation',
                                                      'shift perturbation': 'chemical shift perturbation',
                                                      'chem shift perturbation': 'chemical shift perturbation',
                                                      'CS perturbation': 'chemical shift perturbation',
                                                      'csp': 'chemical shift perturbation',
                                                      'CSP': 'chemical shift perturbation'
                                                     }
                                         }

        self.dihed_alt_constraint_type = {'nef': {'J-couplings': 'jcoupling',
                                                  'j-couplings': 'jcoupling',
                                                  'J couplings': 'jcoupling',
                                                  'j couplings': 'jcoupling',
                                                  'Jcouplings': 'jcoupling',
                                                  'jcouplings': 'jcoupling',
                                                  'J-coupling': 'jcoupling',
                                                  'j-coupling': 'jcoupling',
                                                  'J coupling': 'jcoupling',
                                                  'j coupling': 'jcoupling',
                                                  'Jcoupling': 'jcoupling',
                                                  'chemical shift': 'chemical_shift',
                                                  'Chemical shift': 'chemical_shift',
                                                  'Chemical_shift': 'chemical_shift',
                                                  'chemical shifts': 'chemical_shift',
                                                  'Chemical shifts': 'chemical_shift',
                                                  'Chemical_shifts': 'chemical_shift',
                                                  'backbone chemical shifts': 'chemical_shift',
                                                  'Backbone chemical shifts': 'chemical_shift',
                                                  'Mainchain chemical shifts': 'chemical_shift',
                                                  'mainchain chemical shifts': 'chemical_shift',
                                                  'Main chain chemical shifts': 'chemical_shift',
                                                  'main chain chemical shifts': 'chemical_shift',
                                                  'bb chemical shifts': 'chemical_shift',
                                                  'BB chemical shifts': 'chemical_shift',
                                                  'backbone chemical shift': 'chemical_shift',
                                                  'Backbone chemical shift': 'chemical_shift',
                                                  'Mainchain chemical shift': 'chemical_shift',
                                                  'mainchain chemical shift': 'chemical_shift',
                                                  'Main chain chemical shift': 'chemical_shift',
                                                  'main chain chemical shift': 'chemical_shift',
                                                  'bb chemical shift': 'chemical_shift',
                                                  'BB chemical shift': 'chemical_shift',
                                                  'backbone chem shifts': 'chemical_shift',
                                                  'Backbone chem shifts': 'chemical_shift',
                                                  'Mainchain chem shifts': 'chemical_shift',
                                                  'mainchain chem shifts': 'chemical_shift',
                                                  'Main chain chem shifts': 'chemical_shift',
                                                  'main chain chem shifts': 'chemical_shift',
                                                  'bb chem shifts': 'chemical_shift',
                                                  'BB chem shifts': 'chemical_shift',
                                                  'backbone chem shift': 'chemical_shift',
                                                  'Backbone chem shift': 'chemical_shift',
                                                  'Mainchain chem shift': 'chemical_shift',
                                                  'mainchain chem shift': 'chemical_shift',
                                                  'Main chain chem shift': 'chemical_shift',
                                                  'main chain chem shift': 'chemical_shift',
                                                  'bb chem shift': 'chemical_shift',
                                                  'BB chem shift': 'chemical_shift',
                                                  'backbone cs': 'chemical_shift',
                                                  'Backbone cs': 'chemical_shift',
                                                  'Mainchain cs': 'chemical_shift',
                                                  'mainchain cs': 'chemical_shift',
                                                  'Main chain cs': 'chemical_shift',
                                                  'main chain cs': 'chemical_shift',
                                                  'bb cs': 'chemical_shift',
                                                  'BB cs': 'chemical_shift',
                                                  'backbone CS': 'chemical_shift',
                                                  'Backbone CS': 'chemical_shift',
                                                  'Mainchain CS': 'chemical_shift',
                                                  'mainchain CS': 'chemical_shift',
                                                  'Main chain CS': 'chemical_shift',
                                                  'main chain CS': 'chemical_shift',
                                                  'bb CS': 'chemical_shift',
                                                  'BB CS': 'chemical_shift',
                                                  'TALOS': 'chemical_shift',
                                                  'talos': 'chemical_shift',
                                                  'TALOS+': 'chemical_shift',
                                                  'talos+': 'chemical_shift',
                                                  'TALOS-N': 'chemical_shift',
                                                  'talos-n': 'chemical_shift'
                                                  },
                                          'nmr-star': {'jcoupling': 'J-couplings',
                                                       'Jcoupling': 'J-couplings',
                                                       'jcouplings': 'J-couplings',
                                                       'Jcouplings':'J-couplings',
                                                       'j-couplings': 'J-couplings',
                                                       'J couplings': 'J-couplings',
                                                       'j couplings': 'J-couplings',
                                                       'J-coupling': 'J-couplings',
                                                       'j-coupling': 'J-couplings',
                                                       'J coupling': 'J-couplings',
                                                       'j coupling': 'J-couplings',
                                                       'chemical_shift': 'backbone chemical shifts',
                                                       'Chemical_shift': 'backbone chemical shifts',
                                                       'chemical_shifts': 'backbone chemical shifts',
                                                       'Chemical_shifts': 'backbone chemical shifts',
                                                       'chemical shift': 'backbone chemical shifts',
                                                       'Chemical shift': 'backbone chemical shifts',
                                                       'chemical shifts': 'backbone chemical shifts',
                                                       'Chemical shifts': 'backbone chemical shifts',
                                                       'Backbone chemical shifts': 'backbone chemical shifts',
                                                       'Mainchain chemical shifts': 'backbone chemical shifts',
                                                       'mainchain chemical shifts': 'backbone chemical shifts',
                                                       'Main chain chemical shifts': 'backbone chemical shifts',
                                                       'main chain chemical shifts': 'backbone chemical shifts',
                                                       'bb chemical shifts': 'backbone chemical shifts',
                                                       'BB chemical shifts': 'backbone chemical shifts',
                                                       'backbone chemical shift': 'backbone chemical shifts',
                                                       'Backbone chemical shift': 'backbone chemical shifts',
                                                       'Mainchain chemical shift': 'backbone chemical shifts',
                                                       'mainchain chemical shift': 'backbone chemical shifts',
                                                       'Main chain chemical shift': 'backbone chemical shifts',
                                                       'main chain chemical shift': 'backbone chemical shifts',
                                                       'bb chemical shift': 'backbone chemical shifts',
                                                       'BB chemical shift': 'backbone chemical shifts',
                                                       'backbone chem shifts': 'backbone chemical shifts',
                                                       'Backbone chem shifts': 'backbone chemical shifts',
                                                       'Mainchain chem shifts': 'backbone chemical shifts',
                                                       'mainchain chem shifts': 'backbone chemical shifts',
                                                       'Main chain chem shifts': 'backbone chemical shifts',
                                                       'main chain chem shifts': 'backbone chemical shifts',
                                                       'bb chem shifts': 'backbone chemical shifts',
                                                       'BB chem shifts': 'backbone chemical shifts',
                                                       'backbone chem shift': 'backbone chemical shifts',
                                                       'Backbone chem shift': 'backbone chemical shifts',
                                                       'Mainchain chem shift': 'backbone chemical shifts',
                                                       'mainchain chem shift': 'backbone chemical shifts',
                                                       'Main chain chem shift': 'backbone chemical shifts',
                                                       'main chain chem shift': 'backbone chemical shifts',
                                                       'bb chem shift': 'backbone chemical shifts',
                                                       'BB chem shift': 'backbone chemical shifts',
                                                       'backbone cs': 'backbone chemical shifts',
                                                       'Backbone cs': 'backbone chemical shifts',
                                                       'Mainchain cs': 'backbone chemical shifts',
                                                       'mainchain cs': 'backbone chemical shifts',
                                                       'Main chain cs': 'backbone chemical shifts',
                                                       'main chain cs': 'backbone chemical shifts',
                                                       'bb cs': 'backbone chemical shifts',
                                                       'BB cs': 'backbone chemical shifts',
                                                       'backbone CS': 'backbone chemical shifts',
                                                       'Backbone CS': 'backbone chemical shifts',
                                                       'Mainchain CS': 'backbone chemical shifts',
                                                       'mainchain CS': 'backbone chemical shifts',
                                                       'Main chain CS': 'backbone chemical shifts',
                                                       'main chain CS': 'backbone chemical shifts',
                                                       'bb CS': 'backbone chemical shifts',
                                                       'BB CS': 'backbone chemical shifts',
                                                       'TALOS': 'backbone chemical shifts',
                                                       'talos': 'backbone chemical shifts',
                                                       'TALOS+': 'backbone chemical shifts',
                                                       'talos+': 'backbone chemical shifts',
                                                       'TALOS-N': 'backbone chemical shifts',
                                                       'talos-n': 'backbone chemical shifts'
                                                       }
                                          }

        self.rdc_alt_constraint_type = {'nef': {'RDC': 'measured',
                                                'rdc': 'measured'
                                                },
                                        'nmr-star': {'rdc': 'RDC',
                                                     'measured': 'RDC'
                                                     }
                                        }

        self.entity_del_atom_row = ['ID', 'Entity_assembly_ID', 'Comp_index_ID', 'Comp_ID', 'Atom_ID', 'Auth_entity_assembly_ID', 'Auth_seq_ID', 'Auth_comp_ID', 'Auth_atom_ID', 'Assembly_ID']

        self.int_pattern = re.compile(r'^([+-]?[1-9]\d*|0)$')

    def read_input_file(self, in_file):
        """ Read input NEF/NMR-STAR file.
            @param in_file: input NEF/NMR-STAR file path
            @return: status, Entry/Saveframe/Loop data type or message, data object
        """

        is_ok = True
        star_data = None

        try:
            star_data = pynmrstar.Entry.from_file(in_file)
            msg = 'Entry'

        except Exception as e1:

            try:
                star_data = pynmrstar.Saveframe.from_file(in_file)
                msg = 'Saveframe'

            except Exception as e2:

                try:
                    star_data = pynmrstar.Loop.from_file(in_file)
                    msg = 'Loop'

                except Exception as e3:

                    is_ok = False

                    if version.parse(pynmrstar.__version__) >= version.parse("2.6.5.1"):
                        if not "Invalid loop. Loops must start with the 'loop_' keyword." in str(e3):
                            msg = str(e3)
                        else:
                            msg = str(e1)
                    else:
                        if not 'internaluseyoushouldntseethis_frame' in str(e3):
                            msg = str(e3)
                        else:
                            msg = str(e1) # '%s contains no valid saveframe or loop. PyNMRSTAR ++ Error  - %s' % (os.path.basename(in_file), str(e))
        """
        except Exception as e:
            is_ok = False
            msg = str(e)
        """
        return is_ok, msg, star_data

    def get_one_letter_code(self, comp_id):
        """ Convert comp ID to 1-letter code.
            @change: support empty value by Masashi Yokochi
            @param comp_id: chemical compoent id
            @return: 1-letter code for a given comp_id
        """

        comp_id = comp_id.upper()

        if comp_id in self.monDict3:
            return self.monDict3[comp_id]
        elif comp_id in self.empty_value:
            return '.'
        else:
            return 'X'
        """
        if comp_id in self.codeDict:
            return self.codeDict[comp_id]
        elif comp_id in self.empty_value:
            return '.'
        else:
            return 'X'
        """
    '''
    def get_readable_time_stamp(self, time):
        """ Return time stamp in human readable format for logging.
            @param time: current system time via time.time()
            @return: readable time stamp
        """

        return datetime.datetime.fromtimestamp(time, tz=utc).strftime('%Y-%m-%d %H:%M:%S')
    '''
    def check_mandatory_tags(self, in_file=None, file_type=None):
        """ Returns list of missing mandatory saveframe/loop tags of the input file.
            @change: detect missing mandatory saveframe tags in Saveframe/Loop data as well as Entry data by Masashi Yokochi
            @param in_file: input NEF/NMR-STAR file path
            @param file_type: input file type either 'nef' or 'nmr-star'
            @return: list of missing mandatory saveframe tags, list of missing mandatory loop tags
        """

        mandatory_tag = self.nef_mandatory_tag if file_type == 'nef' else self.star_mandatory_tag

        missing_sf_tags = []
        missing_lp_tags = []

        try:
            star_data = pynmrstar.Entry.from_file(in_file)

            sf_list = [sf.category for sf in star_data.frame_list]

            for _tag in mandatory_tag:

                if _tag[0][0] == '_' and _tag[1] == 'yes':

                    try:
                        tags = star_data.get_tags([_tag[0]])
                        if len(tags[_tag[0]]) == 0 and _tag[0][1:].split('.')[0] in sf_list:
                            missing_sf_tags.append(_tag[0])
                    except: # ValueError:
                        missing_lp_tags.append(_tag[0])

        except: # ValueError:

            try:
                star_data = pynmrstar.Saveframe.from_file(in_file)

                for _tag in mandatory_tag:

                    if _tag[0][0] == '_' and _tag[1] == 'yes':

                        try:
                            tag = star_data.get_tag(_tag[0])
                            if len(tag) == 0 and _tag[0][1:].split('.')[0] == star_data.category:
                                missing_sf_tags.append(_tag[0])
                        except: # ValueError:
                            missing_lp_tags.append(_tag[0])

            except: # ValueError:

                try:
                    star_data = pynmrstar.Loop.from_file(in_file)

                    for _tag in mandatory_tag:

                        if _tag[0][0] == '_' and _tag[0][1:].split('.')[0] == star_data.category and _tag[1] == 'yes':

                            try:
                                get_lp_tag(star_data, _tag[0])
                            except: # ValueError:
                                missing_lp_tags.append(_tag[0])

                except: # ValueError:
                    pass

        return missing_sf_tags, missing_lp_tags

    def is_mandatory_tag(self, item, file_type):
        """ Check if a given tag is mandatory.
            @author: Masashi Yokochi
            @param item: item name
            @param file_type: input file type either 'nef' or 'nmr-star'
            @return: True for mandatory tag, False otherwise
        """

        mandatory_tag = self.nef_mandatory_tag if file_type == 'nef' else self.star_mandatory_tag

        return any(t[0] == item and t[1] == 'yes' for t in mandatory_tag)

    def validate_file(self, in_file, file_subtype='A'):
        """ Validate input NEF/NMR-STAR file.
            @param infile: input NEF/NMR-STAR file path
            @param file_subtype: should be 'A' or 'S' or 'R' where A for All in one file, S for chemical Shifts file, R for Restraints file
            @return: status, JSON message
        """

        is_valid = True
        info = []
        warning = []
        error = []

        file_type = 'unknown'

        err_template_for_missing_mandatory_loop = "The mandatory loop '%s' is missing. Deposition of %s is mandatory. Please re-upload the %s file."
        err_template_for_empty_mandatory_loop = "The mandatory loop '%s' is empty. Deposition of %s is mandatory. Please re-upload the %s file."

        try:

            is_done, data_type, star_data = self.read_input_file(in_file)

            if is_done:

                minimal_info_nef_a = ['_nef_chemical_shift', '_nef_distance_restraint']
                minimal_info_nef_s = ['_nef_chemical_shift']
                minimal_info_nef_r = ['_nef_distance_restraint']

                minimal_info_star_a = ['_Atom_chem_shift', '_Gen_dist_constraint']
                minimal_info_star_s = ['_Atom_chem_shift']
                minimal_info_star_r = ['_Gen_dist_constraint']

                sf_list, lp_list = self.get_data_content(star_data, data_type)

                info.append('{} saveframes and {} loops found'.format(len(sf_list), len(lp_list)))

                nef_sf_list = [i for i in sf_list if not i is None and 'nef' in i]
                nef_lp_list = [i for i in lp_list if not i is None and 'nef' in i]

                info.append('{} saveframes and {} loops found with NEF prefix'.format(len(nef_sf_list), len(nef_lp_list)))

                if len(nef_sf_list) > 0 or len(nef_lp_list) > 0:

                    is_nef_file = True
                    info.append('{} is a NEF file'.format(in_file))
                    file_type = 'nef'

                else:

                    is_nef_file = False
                    info.append('{} is an NMR-STAR file'.format(in_file))
                    file_type = 'nmr-star'

                if is_nef_file:
                    if file_subtype == 'A':

                        for lp_info in minimal_info_nef_a:
                            content_subtype = 'assigned chemical shifts' if 'shift' in lp_info else 'distance restraints'
                            if lp_info not in lp_list:
                                is_valid = False
                                error.append(err_template_for_missing_mandatory_loop % (lp_info, content_subtype, file_type.upper()))
                            else:
                                if self.is_empty_loop(star_data, lp_info, data_type):
                                    is_valid = False
                                    error.append(err_template_for_empty_mandatory_loop % (lp_info, content_subtype, file_type.upper()))

                    elif file_subtype == 'S':
                        content_subtype = 'assigned chemical shifts'

                        for lp_info in minimal_info_nef_s:
                            if lp_info not in lp_list:
                                is_valid = False
                                error.append(err_template_for_missing_mandatory_loop % (lp_info, content_subtype, file_type.upper()))
                            else:
                                if self.is_empty_loop(star_data, lp_info, data_type):
                                    is_valid = False
                                    error.append(err_template_for_empty_mandatory_loop % (lp_info, content_subtype, file_type.upper()))

                    elif file_subtype == 'R':
                        content_subtype = 'distance restraints'

                        for lp_info in minimal_info_nef_r:
                            if lp_info not in lp_list:
                                is_valid = False
                                error.append(err_template_for_missing_mandatory_loop % (lp_info, content_subtype, file_type.upper()))
                            else:
                                if self.is_empty_loop(star_data, lp_info, data_type):
                                    is_valid = False
                                    error.append(err_template_for_empty_mandatory_loop % (lp_info, content_subtype, file_type.upper()))

                    else:
                        is_valid = False
                        error.append('file_subtype flag should be A/S/R')

                else:
                    if file_subtype == 'A':

                        for lp_info in minimal_info_star_a:
                            content_subtype = 'assigned chemical shifts' if 'shift' in lp_info else 'distance restraints'
                            if lp_info not in lp_list:
                                is_valid = False
                                error.append(err_template_for_missing_mandatory_loop % (lp_info, content_subtype, file_type.upper()))
                            else:
                                if self.is_empty_loop(star_data, lp_info, data_type):
                                    is_valid = False
                                    error.append(err_template_for_empty_mandatory_loop % (lp_info, content_subtype, file_type.upper()))

                    elif file_subtype == 'S':
                        content_subtype = 'assigned chemical shifts'

                        for lp_info in minimal_info_star_s:
                            if lp_info not in lp_list:
                                is_valid = False
                                error.append(err_template_for_missing_mandatory_loop % (lp_info, content_subtype, file_type.upper()))
                            else:
                                if self.is_empty_loop(star_data, lp_info, data_type):
                                    is_valid = False
                                    error.append(err_template_for_empty_mandatory_loop % (lp_info, content_subtype, file_type.upper()))

                    elif file_subtype == 'R':
                        content_subtype = 'distance restraints'

                        for lp_info in minimal_info_star_r:
                            if lp_info not in lp_list:
                                is_valid = False
                                error.append(err_template_for_missing_mandatory_loop % (lp_info, content_subtype, file_type.upper()))
                            else:
                                if self.is_empty_loop(star_data, lp_info, data_type):
                                    is_valid = False
                                    error.append(err_template_for_empty_mandatory_loop % (lp_info, content_subtype, file_type.upper()))

                    else:
                        is_valid = False
                        error.append('file_subtype flag should be A/S/R')

            else:
                is_valid = False
                error.append(data_type)

        except Exception as e:
            is_valid = False
            error.append(str(e))

        return is_valid, json.dumps({'info': info, 'warning': warning, 'error': error, 'file_type': file_type})

    def is_empty_loop(self, star_data, lp_category, data_type):
        """ Check if a given loop is empty.
            @return: True for empty loop, False otherwise
        """

        if data_type == 'Entry':
            loops = star_data.get_loops_by_category(lp_category)

            return any(len(loop.data) == 0 for loop in loops)

        elif data_type == 'Saveframe':
            loop = star_data.get_loop_by_category(lp_category)

            return len(loop.data) == 0

        else:
            return len(star_data.data) == 0

    def is_empty(self, data):
        """ Check if given data has empty code.
            @author: Masashi Yokochi
            @return: True for empty data, False otherwise
        """

        return any(d in self.empty_value for d in data)

    def is_data(self, data):
        """ Check if given data has no empty code.
            @author: Masashi Yokochi
            @return: True for non-empty data, False for empty data
        """

        return not any(d in self.empty_value for d in data)

    def get_data_content(self, star_data, data_type):
        """ Extract saveframe categories and loop categories from star data object.
            @return: list of saveframe categories, list of loop categories
        """

        sf_list = []
        lp_list = []

        if data_type == 'Entry':

            for sf in star_data.frame_list:
                sf_list.append(sf.category)

                for lp in sf:
                    lp_list.append(lp.category)

        elif data_type == 'Saveframe':

            for lp in star_data:
                lp_list.append(lp.category)

        elif not star_data is None:
            lp_list.append(star_data.category)

        return sf_list, lp_list

    def get_seq_from_cs_loop(self, in_file):
        """ Extract sequence from chemical shift loop.
            @param in_file: NEF/NMR-STAR file
            @return: status, JSON message
        """

        is_valid, json_dumps = self.validate_file(in_file, 'S')

        data = json.loads(json_dumps)

        info = data['info']
        warning = data['warning']
        error = data['error']

        is_ok = False
        seq = []

        if is_valid:

            info.append('File successfully read ')
            in_data = self.read_input_file(in_file)[-1]

            if data['file_type'] == 'nmr-star':

                info.append('NMR-STAR')
                seq = self.get_star_seq(in_data)

                if len(seq[0]) > 0:
                    is_ok = True

                else:
                    error.append("Can't extract sequence from chemical shift loop")

            elif data['file_type'] == 'nef':

                info.append('NEF')
                seq = self.get_nef_seq(in_data)

                if len(seq[0]) > 0:
                    is_ok = True

                else:
                    error.append("Can't extract sequence from chemical shift loop")

            else:
                error.append("Can't identify file type, it is neither NEF nor NMR-STAR")

        else:
            error.append('File validation failed (or) File contains no chemical shift information')

        return is_ok, json.dumps({'info': info, 'warning': warning, 'error': error, 'file_type': data['file_type'], 'data': seq})

    def get_nef_seq(self, star_data, lp_category='nef_chemical_shift', seq_id='sequence_code', comp_id='residue_name',
                    chain_id='chain_code', allow_empty=False, allow_gap=False):
        """ Extract sequence from any given loops in an NEF file.
            @change: re-written by Masashi Yokochi
            @return: list of sequence information for each loop
        """

        user_warn_msg = ''

        try:
            loops = star_data.get_loops_by_category(lp_category)
        except: # AttributeError:
            try:
                loops = [star_data.get_loop_by_category(lp_category)]
            except: # AttributeError:
                loops = [star_data]

        data = [] # data of all loops

        tags = [seq_id, comp_id, chain_id]

        for loop in loops:
            cmp_dict = {}
            seq_dict = {}

            seq_data = []

            if set(tags) & set(loop.tags) == set(tags):
                seq_data = get_lp_tag(loop, tags)
                for i in seq_data:
                    if i[2] in self.empty_value:
                        i[2] = 1
            else:
                _tags_exist = False
                for i in range(1, self.lim_num_dim):
                    _tags = [seq_id + '_' + str(i), comp_id + '_' + str(i), chain_id + '_' + str(i)]
                    if set(_tags) & set(loop.tags) == set(_tags):
                        _tags_exist = True
                        seq_data += get_lp_tag(loop, _tags)

                if not _tags_exist:
                    missing_tags = list(set(tags) - set(loop.tags))
                    raise LookupError("Missing mandatory %s loop tag%s." % (missing_tags, 's' if len(missing_tags) > 1 else ''))

            if allow_empty:
                seq_data = list(filter(self.is_data, seq_data))
                if len(seq_data) == 0:
                    continue
            else:
                for l, i in enumerate(seq_data):
                    if self.is_empty(i) and l < len(loop.data):
                        r = {}
                        for j in range(len(loop.tags)):
                            r[loop.tags[j]] = loop.data[l][j]
                        user_warn_msg += '[Invalid data] Sequence must not be empty. #_of_row %s, data_of_row %s.\n' % (l + 1, r)
                        #raise ValueError("Sequence must not be empty. #_of_row %s, data_of_row %s." % (l + 1, r))

            for l, i in enumerate(seq_data):
                try:
                    int(i[0])
                except ValueError:
                    if l < len(loop.data):
                        r = {}
                        for j in range(len(loop.tags)):
                            r[loop.tags[j]] = loop.data[l][j]
                        user_warn_msg += '[Invalid data] %s must be int. #_of_row %s, data_of_row %s.\n' % (seq_id, l + 1, r)
                        #raise ValueError("%s must be int. #_of_row %s, data_of_row %s." % (seq_id, l + 1, r))

            if len(user_warn_msg) > 0:
                raise UserWarning(user_warn_msg)

            try:

                chains = sorted(set([i[2] for i in seq_data]))
                offset_seq_ids = {i[2]: 0 for i in seq_data}
                for c in chains:
                    min_seq_id = min([int(i[0]) for i in seq_data if i[2] == c])
                    if min_seq_id < 0:
                        offset_seq_ids[c] = min_seq_id * -1
                sorted_seq = sorted(set(['{} {:04d} {}'.format(i[2], int(i[0]) + offset_seq_ids[i[2]], i[1]) for i in seq_data]))

                chk_dict = {'{} {:04d}'.format(i[2], int(i[0])):i[1] for i in seq_data}

                for i in seq_data:
                    chk_key = '{} {:04d}'.format(i[2], int(i[0]))
                    if chk_dict[chk_key] != i[1]:
                        raise KeyError("Sequence must be unique. %s %s, %s %s, %s %s vs %s." % (chain_id, i[2], seq_id, i[0], comp_id, i[1], chk_dict[chk_key]))

                if len(sorted_seq[0].split(' ')[-1]) > 1:
                    if len(chains) > 1:
                        for c in chains:
                            cmp_dict[c] = [i.split(' ')[-1].upper() for i in sorted_seq if i.split(' ')[0] == c]
                            seq_dict[c] = [int(i.split(' ')[1]) - offset_seq_ids[c] for i in sorted_seq if i.split(' ')[0] == c]
                    else:
                        c = list(chains)[0]
                        cmp_dict[c] = [i.split(' ')[-1].upper() for i in sorted_seq]
                        seq_dict[c] = [int(i.split(' ')[1]) - offset_seq_ids[c] for i in sorted_seq]
                else:
                    if len(chains) > 1:
                        for c in chains:
                            cmp_dict[c] = [i.split(' ')[-1].upper() for i in sorted_seq if i.split(' ')[0] == c]
                            seq_dict[c] = [int(i.split(' ')[1]) - offset_seq_ids[c] for i in sorted_seq if i.split(' ')[0] == c]
                    else:
                        c = list(chains)[0]
                        cmp_dict[c] = [i.split(' ')[-1].upper() for i in sorted_seq]
                        seq_dict[c] = [int(i.split(' ')[1]) - offset_seq_ids[c] for i in sorted_seq]

                asm = [] # assembly of a loop

                for c in chains:
                    ent = {} # entity

                    ent['chain_id'] = str(c)

                    if allow_gap:
                        ent['seq_id'] = []
                        ent['comp_id'] = []

                        _seq_id_ = None

                        for _seq_id, _comp_id in zip(seq_dict[c], cmp_dict[c]):

                            if not _seq_id_ is None and _seq_id_ + 1 != _seq_id and _seq_id_ + 20 > _seq_id:
                                for s in range(_seq_id_ + 1, _seq_id):
                                    ent['seq_id'].append(s)
                                    ent['comp_id'].append('.')
                            ent['seq_id'].append(_seq_id)
                            ent['comp_id'].append(_comp_id)

                            _seq_id_ = _seq_id

                    else:
                        ent['seq_id'] = seq_dict[c]
                        ent['comp_id'] = cmp_dict[c]

                    if len(chains) > 1:
                        identity = []
                        for _c in chains:
                            if _c == c:
                                continue
                            if cmp_dict[_c] == cmp_dict[c]:
                                identity.append(_c)
                        if len(identity) > 0:
                            ent['identical_chain_id'] = identity

                    asm.append(ent)

                data.append(asm)

            except ValueError:
                #raise ValueError("%s must be int." % seq_id)
                pass

        if len(data) == 0:
            data.append([])

        return data

    def get_star_seq(self, star_data, lp_category='Atom_chem_shift', seq_id='Comp_index_ID', comp_id='Comp_ID',
                     chain_id='Entity_assembly_ID', allow_empty=False, allow_gap=False):
        """ Extract sequence from any given loops in an NMR-STAR file.
            @change: re-written by Masashi Yokochi
            @return: list of sequence information for each loop
        """

        user_warn_msg = ''

        try:
            loops = star_data.get_loops_by_category(lp_category)
        except: # AttributeError:
            try:
                loops = [star_data.get_loop_by_category(lp_category)]
            except: # AttributeError:
                loops = [star_data]

        data = [] # data of all loops

        tags = [seq_id, comp_id, chain_id]
        tags_ = [seq_id, comp_id]

        for loop in loops:
            cmp_dict = {}
            seq_dict = {}

            seq_data = []

            if set(tags) & set(loop.tags) == set(tags):
                seq_data = get_lp_tag(loop, tags)
                for i in seq_data:
                    if i[2] in self.empty_value:
                        i[2] = '1'
            elif set(tags_) & set(loop.tags) == set(tags_): # No Entity_assembly_ID tag case
                seq_data = get_lp_tag(loop, tags_)
                for i in seq_data:
                    i.append('1')
            else:
                _tags_exist = False
                for i in range(1, self.lim_num_dim):
                    _tags = [seq_id + '_' + str(i), comp_id + '_' + str(i), chain_id + '_' + str(i)]
                    _tags_ = [seq_id + '_' + str(i), comp_id + '_' + str(i)]
                    if set(_tags) & set(loop.tags) == set(_tags):
                        _tags_exist = True
                        seq_data += get_lp_tag(loop, _tags)
                    elif set(_tags_) & set(loop.tags) == set(_tags_):
                        _tags_exist = True
                        seq_data_ = get_lp_tag(loop, _tags_)
                        for i in seq_data_:
                            i.append('1')
                        seq_data += seq_data_

                if not _tags_exist:
                    missing_tags = list(set(tags) - set(loop.tags))
                    raise LookupError("Missing mandatory %s loop tag%s." % (missing_tags, 's' if len(missing_tags) > 1 else ''))

            if allow_empty:
                seq_data = list(filter(self.is_data, seq_data))
                if len(seq_data) == 0:
                    continue
            else:
                for l, i in enumerate(seq_data):
                    if self.is_empty(i) and l < len(loop.data):
                        r = {}
                        for j in range(len(loop.tags)):
                            r[loop.tags[j]] = loop.data[l][j]
                        user_warn_msg += '[Invalid data] Sequence must not be empty. #_of_row %s, data_of_row %s.\n' % (l + 1, r)
                        #raise ValueError("Sequence must not be empty. #_of_row %s, data_of_row %s." % (l + 1, r))

            for l, i in enumerate(seq_data):
                try:
                    int(i[0])
                except ValueError:
                    if l < len(loop.data):
                        r = {}
                        for j in range(len(loop.tags)):
                            r[loop.tags[j]] = loop.data[l][j]
                        user_warn_msg += '[Invalid data] %s must be int. #_of_row %s, data_of_row %s.\n' % (seq_id, l + 1, r)
                        #raise ValueError("%s must be int. #_of_row %s, data_of_row %s." % (seq_id, l + 1, r))

            if len(user_warn_msg) > 0:
                raise UserWarning(user_warn_msg)

            try:

                chains = sorted(set([i[2] for i in seq_data]))
                offset_seq_ids = {i[2]: 0 for i in seq_data}
                for c in chains:
                    min_seq_id = min([int(i[0]) for i in seq_data if i[2] == c])
                    if min_seq_id < 0:
                        offset_seq_ids[c] = min_seq_id * -1
                sorted_seq = sorted(set(['{} {:04d} {}'.format(i[2], int(i[0]) + offset_seq_ids[i[2]], i[1]) for i in seq_data]))

                chk_dict = {'{} {:04d}'.format(i[2], int(i[0])):i[1] for i in seq_data}

                for i in seq_data:
                    chk_key = '{} {:04d}'.format(i[2], int(i[0]))
                    if chk_dict[chk_key] != i[1]:
                        raise KeyError("Sequence must be unique. %s %s, %s %s, %s %s vs %s." % (chain_id, i[2], seq_id, i[0], comp_id, i[1], chk_dict[chk_key]))

                if len(sorted_seq[0].split(' ')[-1]) > 1:
                    if len(chains) > 1:
                        for c in chains:
                            cmp_dict[c] = [i.split(' ')[-1].upper() for i in sorted_seq if i.split(' ')[0] == c]
                            seq_dict[c] = [int(i.split(' ')[1]) - offset_seq_ids[c] for i in sorted_seq if i.split(' ')[0] == c]
                    else:
                        c = list(chains)[0]
                        cmp_dict[c] = [i.split(' ')[-1].upper() for i in sorted_seq]
                        seq_dict[c] = [int(i.split(' ')[1]) - offset_seq_ids[c] for i in sorted_seq]
                else:
                    if len(chains) > 1:
                        for c in chains:
                            cmp_dict[c] = [i.split(' ')[-1].upper() for i in sorted_seq if i.split(' ')[0] == c]
                            seq_dict[c] = [int(i.split(' ')[1]) - offset_seq_ids[c] for i in sorted_seq if i.split(' ')[0] == c]
                    else:
                        c = list(chains)[0]
                        cmp_dict[c] = [i.split(' ')[-1].upper() for i in sorted_seq]
                        seq_dict[c] = [int(i.split(' ')[1]) - offset_seq_ids[c] for i in sorted_seq]

                asm = [] # assembly of a loop

                for c in chains:
                    ent = {} # entity

                    ent['chain_id'] = str(c)

                    if allow_gap:
                        ent['seq_id'] = []
                        ent['comp_id'] = []

                        _seq_id_ = None

                        for _seq_id, _comp_id in zip(seq_dict[c], cmp_dict[c]):

                            if not _seq_id_ is None and _seq_id_ + 1 != _seq_id and _seq_id_ + 20 > _seq_id:
                                for s in range(_seq_id_ + 1, _seq_id):
                                    ent['seq_id'].append(s)
                                    ent['comp_id'].append('.')
                            ent['seq_id'].append(_seq_id)
                            ent['comp_id'].append(_comp_id)

                            _seq_id_ = _seq_id

                    else:
                        ent['seq_id'] = seq_dict[c]
                        ent['comp_id'] = cmp_dict[c]

                    if len(chains) > 1:
                        identity = []
                        for _c in chains:
                            if _c == c:
                                continue
                            if cmp_dict[_c] == cmp_dict[c]:
                                identity.append(_c)
                        if len(identity) > 0:
                            ent['identical_chain_id'] = identity

                    asm.append(ent)

                data.append(asm)

            except ValueError:
                #raise ValueError("%s must be int." % seq_id)
                pass

        if len(data) == 0:
            data.append([])

        return data

    def get_star_auth_seq(self, star_data, lp_category='Atom_chem_shift', aseq_id='Auth_seq_ID', acomp_id='Auth_comp_ID',
                          asym_id='Auth_asym_ID', seq_id='Comp_index_ID', chain_id='Entity_assembly_ID', allow_empty=True):
        """ Extract author sequence from any given loops in an NMR-STAR file.
            @author: Masashi Yokochi
            @return: list of author sequence information for each loop
        """

        user_warn_msg = ''

        try:
            loops = star_data.get_loops_by_category(lp_category)
        except: # AttributeError:
            try:
                loops = [star_data.get_loop_by_category(lp_category)]
            except: # AttributeError:
                loops = [star_data]

        data = [] # data of all loops

        tags = [aseq_id, acomp_id, asym_id, seq_id, chain_id]
        tags_ = [aseq_id, acomp_id, seq_id, asym_id]

        for loop in loops:
            seq_dict = {}
            acmp_dict = {}
            aseq_dict = {}
            asym_dict = {}

            seq_data = []

            if set(tags) & set(loop.tags) == set(tags):
                seq_data = get_lp_tag(loop, tags)
            elif set(tags_) & set(loop.tags) == set(tags_): # No Entity_assembly_ID tag case
                seq_data = get_lp_tag(loop, tags_)
                for i in seq_data:
                    i.append('1')
            else:
                _tags_exist = False
                for i in range(1, self.lim_num_dim):
                    _tags = [aseq_id + '_' + str(i), acomp_id + '_' + str(i), asym_id + '_' + str(i), seq_id + '_' + str(i), chain_id + '_' + str(i)]
                    _tags_ = [aseq_id + '_' + str(i), acomp_id + '_' + str(i), asym_id + '_' + str(i), seq_id + '_' + str(i)]
                    if set(_tags) & set(loop.tags) == set(_tags):
                        _tags_exist = True
                        seq_data += get_lp_tag(loop, _tags)
                    elif set(_tags_) & set(loop.tags) == set(_tags_):
                        _tags_exist = True
                        seq_data_ = get_lp_tag(loop, _tags_)
                        for i in seq_data_:
                            i.append('1')
                        seq_data += seq_data_

                if not _tags_exist:
                    missing_tags = list(set(tags) - set(loop.tags))
                    raise LookupError("Missing mandatory %s loop tag%s." % (missing_tags, 's' if len(missing_tags) > 1 else ''))

            if allow_empty:
                seq_data = list(filter(self.is_data, seq_data))
                if len(seq_data) == 0:
                    continue
            else:
                for l, i in enumerate(seq_data):
                    if self.is_empty(i) and l < len(loop.data):
                        r = {}
                        for j in range(len(loop.tags)):
                            r[loop.tags[j]] = loop.data[l][j]
                        user_warn_msg += '[Invalid data] Author sequence must not be empty. #_of_row %s, data_of_row %s.\n' % (l + 1, r)
                        #raise ValueError("Author sequence must not be empty. #_of_row %s, data_of_row %s." % (l + 1, r))

            for l, i in enumerate(seq_data):
                try:
                    int(i[3])
                except ValueError:
                    if l < len(loop.data):
                        r = {}
                        for j in range(len(loop.tags)):
                            r[loop.tags[j]] = loop.data[l][j]
                        user_warn_msg += '[Invalid data] %s must be int. #_of_row %s, data_of_row %s.\n' % (seq_id, l + 1, r)
                        #raise ValueError("%s must be int. #_of_row %s, data_of_row %s." % (seq_id, l + 1, r))

            if len(user_warn_msg) > 0:
                raise UserWarning(user_warn_msg)

            try:

                chains = sorted(set([i[4] for i in seq_data]))
                offset_seq_ids = {i[4]: 0 for i in seq_data}
                for c in chains:
                    min_seq_id = min([int(i[3]) for i in seq_data if i[4] == c])
                    if min_seq_id < 0:
                        offset_seq_ids[c] = min_seq_id * -1
                sorted_seq = sorted(set(['{}:{:04d}:{}:{: >4}:{}'.format(i[4], int(i[3]) + offset_seq_ids[i[4]], i[2], i[0], i[1]) for i in seq_data]))

                chk_dict = {'{}:{:04d}:{}:{: >4}'.format(i[4], int(i[3]), i[2], i[0]):i[1] for i in seq_data}

                for i in seq_data:
                    chk_key = '{}:{:04d}:{}:{: >4}'.format(i[4], int(i[3]), i[2], i[0])
                    if chk_dict[chk_key] != i[1]:
                        raise KeyError("Author sequence must be unique. %s %s, %s %s, %s %s, %s %s, %s %s vs %s." %\
                                       (chain_id, i[4], seq_id, i[3], asym_id, i[2], aseq_id, i[0], acomp_id, i[1], chk_dict[chk_key]))

                if len(sorted_seq[0].split(':')[-1]) > 1:
                    if len(chains) > 1:
                        for c in chains:
                            acmp_dict[c] = [i.split(':')[-1] for i in sorted_seq if i.split(':')[0] == c]
                            aseq_dict[c] = [i.split(':')[3].strip() for i in sorted_seq if i.split(':')[0] == c]
                            asym_dict[c] = [i.split(':')[2] for i in sorted_seq if i.split(':')[0] == c]
                            seq_dict[c] = [int(i.split(':')[1]) - offset_seq_ids[c] for i in sorted_seq if i.split(':')[0] == c]
                    else:
                        c = list(chains)[0]
                        acmp_dict[c] = [i.split(':')[-1] for i in sorted_seq]
                        aseq_dict[c] = [i.split(':')[3].strip() for i in sorted_seq]
                        asym_dict[c] = [i.split(':')[2] for i in sorted_seq]
                        seq_dict[c] = [int(i.split(':')[1]) - offset_seq_ids[c] for i in sorted_seq]
                else:
                    if len(chains) > 1:
                        for c in chains:
                            acmp_dict[c] = [i.split(':')[-1] for i in sorted_seq if i.split(':')[0] == c]
                            aseq_dict[c] = [i.split(':')[3].strip() for i in sorted_seq if i.split(':')[0] == c]
                            asym_dict[c] = [i.split(':')[2] for i in sorted_seq if i.split(':')[0] == c]
                            seq_dict[c] = [int(i.split(':')[1]) - offset_seq_ids[c] for i in sorted_seq if i.split(':')[0] == c]
                    else:
                        c = list(chains)[0]
                        acmp_dict[c] = [i.split(':')[-1] for i in sorted_seq]
                        aseq_dict[c] = [i.split(':')[3].strip() for i in sorted_seq]
                        asym_dict[c] = [i.split(':')[2] for i in sorted_seq]
                        seq_dict[c] = [int(i.split(':')[1]) - offset_seq_ids[c] for i in sorted_seq]

                asm = [] # assembly of a loop

                for c in chains:
                    ent = {} # entity

                    ent['chain_id'] = str(c)
                    ent['seq_id'] = seq_dict[c]
                    ent['auth_asym_id'] = asym_dict[c]
                    ent['auth_seq_id'] = aseq_dict[c]
                    ent['auth_comp_id'] = acmp_dict[c]

                    asm.append(ent)

                data.append(asm)

            except ValueError:
                #raise ValueError("%s must be int." % seq_id)
                pass

        if len(data) == 0:
            data.append([])

        return data

    def get_nef_comp_atom_pair(self, star_data, lp_category='nef_chemical_shift', comp_id='residue_name', atom_id='atom_name',
                               allow_empty=False):
        """ Wrapper function of get_comp_atom_pair() for an NEF file.
            @author: Masashi Yokochi
        """

        return self.get_comp_atom_pair(star_data, lp_category, comp_id, atom_id, allow_empty)

    def get_star_comp_atom_pair(self, star_data, lp_category='Atom_chem_shift', comp_id='Comp_ID', atom_id='Atom_ID',
                                allow_empty=False):
        """ Wrapper function of get_comp_atom_pair() for an NMR-STAR file.
            @author: Masashi Yokochi
        """

        return self.get_comp_atom_pair(star_data, lp_category, comp_id, atom_id, allow_empty)

    def get_star_auth_comp_atom_pair(self, star_data, lp_category='Atom_chem_shift', comp_id='Auth_comp_ID', atom_id='Auth_atom_ID',
                                     allow_empty=True):
        """ Wrapper function of get_comp_atom_pair() for pairs of author comp_id and author atom_id in an NMR-STAR file.
            @author: Masashi Yokochi
        """

        return self.get_comp_atom_pair(star_data, lp_category, comp_id, atom_id, allow_empty)

    def get_comp_atom_pair(self, star_data, lp_category, comp_id, atom_id, allow_empty):
        """ Extract unique pairs of comp_id and atom_id from any given loops in an NEF/NMR-STAR file.
            @author: Masashi Yokochi
            @return: list of unique pairs of comp_id and atom_id for each loop
        """

        user_warn_msg = ''

        try:
            loops = star_data.get_loops_by_category(lp_category)
        except: # AttributeError:
            try:
                loops = [star_data.get_loop_by_category(lp_category)]
            except: # AttributeError:
                loops = [star_data]

        data = [] # data of all loops

        tags = [comp_id, atom_id]

        for loop in loops:
            atm_dict = {}

            pair_data = []

            if set(tags) & set(loop.tags) == set(tags):
                pair_data = get_lp_tag(loop, tags)
            else:
                _tags_exist = False
                for i in range(1, self.lim_num_dim):
                    _tags = [comp_id + '_' + str(i), atom_id + '_' + str(i)]
                    if set(_tags) & set(loop.tags) == set(_tags):
                        _tags_exist = True
                        pair_data += get_lp_tag(loop, _tags)

                if not _tags_exist:
                    missing_tags = list(set(tags) - set(loop.tags))
                    raise LookupError("Missing mandatory %s loop tag%s." % (missing_tags, 's' if len(missing_tags) > 1 else ''))

            if allow_empty:
                pair_data = list(filter(self.is_data, pair_data))
                if len(pair_data) == 0:
                    continue
            else:
                for l, i in enumerate(pair_data):
                    if self.is_empty(i) and l < len(loop.data):
                        r = {}
                        for j in range(len(loop.tags)):
                            r[loop.tags[j]] = loop.data[l][j]
                        user_warn_msg += '[Invalid data] %s and %s must not be empty. #_of_row %s, data_of_row %s.\n' %\
                                         (comp_id, atom_id, l + 1, r)
                        #raise ValueError("%s and %s must not be empty. #_of_row %s, data_of_row %s." %\
                        #                 (comp_id, atom_id, l + 1, r))

            if len(user_warn_msg) > 0:
                raise UserWarning(user_warn_msg)

            comps = sorted(set([i[0].upper() for i in pair_data if not i[0] in self.empty_value]))
            sorted_comp_atom = sorted(set(['{} {}'.format(i[0].upper(), i[1]) for i in pair_data]))

            for c in comps:
                atm_dict[c] = [i.split(' ')[1] for i in sorted_comp_atom if i.split(' ')[0] == c]

            asm = [] # assembly of a loop

            for c in comps:
                ent = {} # entity

                ent['comp_id'] = c
                ent['atom_id'] = atm_dict[c]

                asm.append(ent)

            data.append(asm)

        if len(data) == 0:
            data.append([])

        return data

    def get_nef_atom_type_from_cs_loop(self, star_data, lp_category='nef_chemical_shift', atom_type='element', isotope_number='isotope_number', atom_id='atom_name',
                                       allow_empty=False):
        """ Wrapper function of get_atom_type_from_cs_loop() for an NEF file.
            @author: Masashi Yokochi
        """

        return self.get_atom_type_from_cs_loop(star_data, lp_category, atom_type, isotope_number, atom_id, allow_empty)

    def get_star_atom_type_from_cs_loop(self, star_data, lp_category='Atom_chem_shift', atom_type='Atom_type', isotope_number='Atom_isotope_number', atom_id='Atom_ID',
                                        allow_empty=False):
        """ Wrapper function of get_atom_type_from_cs_loop() for an NMR-SAR file.
            @author: Masashi Yokochi
        """

        return self.get_atom_type_from_cs_loop(star_data, lp_category, atom_type, isotope_number, atom_id, allow_empty)

    def get_atom_type_from_cs_loop(self, star_data, lp_category, atom_type, isotope_number, atom_id, allow_empty):
        """ Extract unique pairs of atom_type, isotope number, and atom_id from assigned chemical shifts in n NEF/NMR-SAR file.
            @author: Masashi Yokochi
            @return: list of unique pairs of atom_type, isotope number, and atom_id for each CS loop
        """

        user_warn_msg = ''

        try:
            loops = star_data.get_loops_by_category(lp_category)
        except: # AttributeError:
            try:
                loops = [star_data.get_loop_by_category(lp_category)]
            except: # AttributeError:
                loops = [star_data]

        data = [] # data of all loops

        tags = [atom_type, isotope_number, atom_id]

        for loop in loops:
            ist_dict = {}
            atm_dict = {}

            a_type_data = []

            if set(tags) & set(loop.tags) != set(tags):
                missing_tags = list(set(tags) - set(loop.tags))
                raise LookupError("Missing mandatory %s loop tag%s." % (missing_tags, 's' if len(missing_tags) > 1 else ''))

            a_type_data = get_lp_tag(loop, tags)

            if allow_empty:
                a_type_data = list(filter(self.is_data, a_type_data))
                if len(a_type_data) == 0:
                    continue
            else:
                for l, i in enumerate(a_type_data):
                    if self.is_empty(i) and l < len(loop.data):
                        r = {}
                        for j in range(len(loop.tags)):
                            r[loop.tags[j]] = loop.data[l][j]
                        user_warn_msg += '[Invalid data] %s, %s, and %s must not be empty. #_of_row %s, data_of_row %s.\n' %\
                                         (atom_type, isotope_number, atom_id, l + 1, r)
                        #raise ValueError("%s, %s, and %s must not be empty. #_of_row %s, data_of_row %s." %\
                        #                 (atom_type, isotope_number, atom_id, l + 1, r))

            for l, i in enumerate(a_type_data):
                try:
                    int(i[1])
                except ValueError:
                    if l < len(loop.data):
                        r = {}
                        for j in range(len(loop.tags)):
                            r[loop.tags[j]] = loop.data[l][j]
                        user_warn_msg += '[Invalid data] %s must be int. #_of_row %s, data_of_row %s.\n' % (isotope_number, l + 1, r)
                        #raise ValueError("%s must be int. #_of_row %s, data_of_row %s." % (isotope_number, l + 1, r))

            if len(user_warn_msg) > 0:
                raise UserWarning(user_warn_msg)

            try:

                a_types = sorted(set([i[0] for i in a_type_data]))
                sorted_ist = sorted(set(['{} {}'.format(i[0], i[1]) for i in a_type_data]))
                sorted_atm = sorted(set(['{} {}'.format(i[0], i[2]) for i in a_type_data]))

                for t in a_types:
                    ist_dict[t] = [int(i.split(' ')[1]) for i in sorted_ist if i.split(' ')[0] == t]
                    atm_dict[t] = [i.split(' ')[1] for i in sorted_atm if i.split(' ')[0] == t]

                asm = [] # assembly of a loop

                for t in a_types:
                    ent = {} # entity

                    ent['atom_type'] = t
                    ent['isotope_number'] = ist_dict[t]
                    ent['atom_id'] = atm_dict[t]

                    asm.append(ent)

                data.append(asm)

            except ValueError:
                #raise ValueError("%s must be int." % isotope_number)
                pass

        if len(data) == 0:
            data.append([])

        return data

    def get_star_ambig_code_from_cs_loop(self, star_data, lp_category='Atom_chem_shift', comp_id='Comp_ID', atom_id='Atom_ID', ambig_code='Ambiguity_code', ambig_set_id='Ambiguity_set_ID'):
        """ Extract unique pairs of comp_id, atom_id, and ambiguity code from assigned chemical shifts in an NMR-SAR file.
            @author: Masashi Yokochi
            @return: list of unique pairs of comp_id, atom_id, and ambiguity code for each CS loop
        """

        user_warn_msg = ''

        try:
            loops = star_data.get_loops_by_category(lp_category)
        except: # AttributeError:
            try:
                loops = [star_data.get_loop_by_category(lp_category)]
            except: # AttributeError:
                loops = [star_data]

        data = [] # data of all loops

        tags = [comp_id, atom_id, ambig_code, ambig_set_id]

        for loop in loops:
            atm_dict = {}

            ambig_data = []

            if set(tags) & set(loop.tags) != set(tags):
                missing_tags = list(set(tags) - set(loop.tags))
                raise LookupError("Missing mandatory %s loop tag%s." % (missing_tags, 's' if len(missing_tags) > 1 else ''))

            ambig_data = get_lp_tag(loop, tags)

            if len(ambig_data) == 0:
                data.append(None)
                continue

            for l, i in enumerate(ambig_data):
                # already checked elsewhere
                #if i[0] in self.empty_value:
                #   raise ValueError("%s must not be empty." % comp_id)
                #if i[1] in self.empty_value:
                #    raise ValueError("%s must not be empty." % atom_id)
                if not i[2] in self.empty_value:

                    try:
                       code = int(i[2])
                    except ValueError:
                        if l < len(loop.data):
                            r = {}
                            for j in range(len(loop.tags)):
                                r[loop.tags[j]] = loop.data[l][j]
                            user_warn_msg += '[Invalid data] %s must be one of %s. #_of_row %s, data_of_row %s.\n' % (ambig_code, list(self.bmrb_ambiguity_codes), l + 1, r)
                            #raise ValueError("%s must be one of %s. #_of_row %s, data_of_row %s." % (ambig_code, list(self.bmrb_ambiguity_codes), l + 1, r))

                    if not code in self.bmrb_ambiguity_codes:
                        if l < len(loop.data):
                            r = {}
                            for j in range(len(loop.tags)):
                                r[loop.tags[j]] = loop.data[l][j]
                            user_warn_msg += '[Invalid data] %s must be one of %s. #_of_row %s, data_of_row %s.\n' % (ambig_code, list(self.bmrb_ambiguity_codes), l + 1, r)
                            #raise ValueError("%s must be one of %s. #_of_row %s, data_of_row %s." % (ambig_code, list(self.bmrb_ambiguity_codes), l + 1, r))

                    if code >= 4:
                        if i[3] in self.empty_value and l < len(loop.data):
                            r = {}
                            for j in range(len(loop.tags)):
                                r[loop.tags[j]] = loop.data[l][j]
                            user_warn_msg += '[Invalid data] %s must not be empty for %s %s. #_of_row %s, data_of_row %s.\n' % (ambig_set_id, ambig_code, code, l + 1, r)
                            #raise ValueError("%s must not be empty for %s %s. #_of_row %s, data_of_row %s." % (ambig_set_id, ambig_code, code, l + 1, r))
                        else:
                            try:
                                int(i[3])
                            except ValueError:
                                if l < len(loop.data):
                                    r = {}
                                    for j in range(len(loop.tags)):
                                        r[loop.tags[j]] = loop.data[l][j]
                                    user_warn_msg += '[Invalid data] %s must be int. #_of_row %s, data_of_row %s.\n' % (ambig_set_id, l + 1, r)
                                    #raise ValueError("%s must be int. #_of_row %s, data_of_row %s." % (ambig_set_id, l + 1, r))

                if not i[3] in self.empty_value:

                    if i[2] in self.empty_value or not i[2] in ['4', '5', '6', '9']:
                        if l < len(loop.data):
                            r = {}
                            for j in range(len(loop.tags)):
                                r[loop.tags[j]] = loop.data[l][j]
                            user_warn_msg += '[Invalid data] %s must be empty for %s %s. #_of_row %s, data_of_row %s.\n' % (ambig_set_id, ambig_code, i[2], l + 1, r)
                            #raise ValueError("%s must be empty for %s %s. #_of_row %s, data_of_row %s." % (ambig_set_id, ambig_code, i[2], l + 1, r))

            if len(user_warn_msg) > 0:
                raise UserWarning(user_warn_msg)

            ambigs = sorted(set(['{}:{}'.format(i[0].upper(), i[2]) for i in ambig_data]))
            sorted_atm = sorted(set(['{}:{} {}'.format(i[0].upper(), i[2], i[1]) for i in ambig_data]))

            for a in ambigs:
                atm_dict[a] = [i.split(' ')[1] for i in sorted_atm if i.split(' ')[0] == a]

            asm = [] # assembly of a loop

            for a in ambigs:
                ent = {} # entity

                split_a = a.split(':')

                ent['comp_id'] = split_a[0]
                ent['ambig_code'] = None if split_a[1] in self.empty_value else int(split_a[1])
                ent['atom_id'] = atm_dict[a]

                asm.append(ent)

            data.append(asm)

        if len(data) == 0:
            data.append([])

        return data

    def get_nef_index(self, star_data, lp_category='nef_sequence', index_id='index'):
        """ Wrapper function of get_index() for an NEF file.
            @author: Masashi Yokochi
        """

        return self.get_index(star_data, lp_category, index_id)

    def get_star_index(self, star_data, lp_category='Chem_comp_assembly', index_id='NEF_index'):
        """ Wrapper function of get_index() for an NMR-STAR file.
            @author: Masashi Yokochi
        """

        return self.get_index(star_data, lp_category, index_id)

    def get_index(self, star_data, lp_category, index_id):
        """ Extract index_id from any given loops in an NEF/NMR-STAR file.
            @author: Masashi Yokochi
            @return: list of index for each loop
        """

        user_warn_msg = ''

        try:
            loops = star_data.get_loops_by_category(lp_category)
        except: # AttributeError:
            try:
                loops = [star_data.get_loop_by_category(lp_category)]
            except: # AttributeError:
                loops = [star_data]

        data = [] # data of all loops

        tags = [index_id]

        for loop in loops:
            index_data = []

            if set(tags) & set(loop.tags) == set(tags):
                index_data = get_lp_tag(loop, tags)
            else:
                raise LookupError("Missing mandatory %s loop tag." % index_id)

            for l, i in enumerate(index_data):
                if self.is_empty(i) and l < len(loop.data):
                    r = {}
                    for j in range(len(loop.tags)):
                        r[loop.tags[j]] = loop.data[l][j]
                    user_warn_msg += '[Invalid data] %s must not be empty. #_of_row %s, data_of_row %s.\n' % (index_id, l + 1, r)
                    #raise ValueError("%s must not be empty. #_of_row %s, data_of_row %s." % (index_id, l + 1, r))
                else:
                    try:
                        int(i[0])
                    except ValueError:
                        if l < len(loop.data):
                            r = {}
                            for j in range(len(loop.tags)):
                                r[loop.tags[j]] = loop.data[l][j]
                            user_warn_msg += '[Invalid data] %s must be int. #_of_row %s, data_of_row %s.\n' % (index_id, l + 1, r)
                            #raise ValueError("%s must be int. #_of_row %s, data_of_row %s." % (index_id, l + 1, r))

            if len(user_warn_msg) > 0:
                raise UserWarning(user_warn_msg)

            try:

                if __pynmrstar_v3__:
                    idxs = [int(i) for i in index_data]
                else:
                    idxs = [int(i) for i in index_data[0]]

                dup_idxs = [i for i in set(idxs) if idxs.count(i) > 1]

                if len(dup_idxs) > 0:
                    raise KeyError("%s must be unique in loop. %s are duplicated." % (index_id, dup_idxs))

                data.append(idxs)

            except ValueError:
                #raise ValueError("%s must be int." % index_id)
                pass

        if len(data) == 0:
            data.append([])

        return data

    def check_nef_data(self, star_data, lp_category='nef_chemical_shift',
                       key_items=[{'name': 'chain_code', 'type': 'str'},
                                  {'name': 'sequence_code', 'type': 'int'},
                                  {'name': 'residue_name', 'type': 'str'},
                                  {'name': 'atom_name', 'type': 'str'}],
                       data_items=[{'name': 'value', 'type': 'float', 'mandatory': True},
                                   {'name': 'value_uncertainty', 'type': 'positive-float', 'mandatory': False}]):
        """ Wrapper function of check_data() for an NEF file.
            @author: Masashi Yokochi
        """

        return self.check_data(star_data, lp_category, key_items, data_items)

    def check_star_data(self, star_data, lp_category='Atom_chem_shift',
                        key_items=[{'name': 'Entity_assembly_ID', 'type': 'int'},
                                   {'name': 'Comp_index_ID', 'type': 'int'},
                                   {'name': 'Comp_ID', 'type': 'str'},
                                   {'name': 'Atom_ID', 'type': 'str'}],
                        data_items=[{'name': 'Val', 'type': 'float', 'mandatory': True},
                                    {'name': 'Val_err', 'type': 'positive-float', 'mandatory': False}]):
        """ Wrapper function of check_data() for an NMR-STAR file.
            @author: Masashi Yokochi
        """

        return self.check_data(star_data, lp_category, key_items, data_items)

    def check_data(self, star_data, lp_category, key_items, data_items, allowed_tags=None, disallowed_tags=None,
                   test_on_index=False, enforce_non_zero=False, enforce_sign=False, enforce_range=False, enforce_enum=False,
                   excl_missing_data=False):
        """ Extract data with sanity check from any given loops in an NEF/NMR-STAR file.
            @author: Masashi Yokochi
            @return: list of extracted data for each loop
        """

        try:
            loops = star_data.get_loops_by_category(lp_category)
        except: # AttributeError:
            try:
                loops = [star_data.get_loop_by_category(lp_category)]
            except: # AttributeError:
                loops = [star_data]

        user_warn_msg = ''

        data = [] # data of all loops

        item_types = ('str', 'bool', 'int', 'index-int', 'positive-int', 'positive-int-as-str', 'pointer-index', 'float', 'positive-float', 'range-float', 'enum', 'enum-int')

        key_names = [k['name'] for k in key_items]
        data_names = [d['name'] for d in data_items]
        mand_data_names = [d['name'] for d in data_items if d['mandatory']]
        _mand_data_names = [d['name'] for d in data_items if d['mandatory'] and 'default-from' in d]

        key_len = len(key_items)

        for k in key_items:
            if not k['type'] in item_types:
                raise TypeError("Type %s of data item %s must be one of %s." % (k['type'], k['name'], item_types))

        for d in data_items:
            if not d['type'] in item_types:
                raise TypeError("Type %s of data item %s must be one of %s." % (d['type'], d['name'], item_types))

        if not allowed_tags is None:

            if len(key_names) > 0 and (set(key_names) | set(allowed_tags)) != set(allowed_tags):
                raise LookupError("Key items %s must not exists." % ((set(key_names) | set(allowed_tags)) - set(allowed_tags)))

            if len(data_names) > 0 and (set(data_names) | set(allowed_tags)) != set(allowed_tags):
                raise LookupError("Data items %s must not exists." % ((set(data_names) | set(allowed_tags)) - set(allowed_tags)))

            for d in data_items:
                if 'group-mandatory' in d and d['group-mandatory']:
                    group = d['group']
                    if not group['member-with'] is None:
                        for m in group['member-with']:
                            if not m in allowed_tags:
                                raise ValueError("Member data item %s of %s must exists in allowed tags." % (m, d['name']))
                    if not group['coexist-with'] is None:
                        for c in group['coexist-with']:
                            if not c in allowed_tags:
                                raise ValueError("Coexisting data item %s of %s must exists in allowed tags." % (c, d['name']))
                    if 'smaller-than' in group and group['smaller-than']:
                        for s in group['smaller-than']:
                            if not s in allowed_tags:
                                raise ValueError("Smaller data item %s of %s must exists in allowed tags." % (s, d['name']))
                    if 'larger-than' in group and group['larger-than']:
                        for l in group['larger-than']:
                            if not l in allowed_tags:
                                raise ValueError("Larger data item %s of %s must exists in allowed tags." % (l, d['name']))
                    if 'not-equal-to' in group and group['not-equal-to']:
                        for l in group['not-equal-to']:
                            if not l in allowed_tags:
                                raise ValueError("Nonequal data item %s of %s must exists in allowed tags." % (l, d['name']))

        for loop in loops:
            tag_data = []

            if len(key_names) > 0 and set(key_names) & set(loop.tags) != set(key_names):
                missing_tags = list(set(key_names) - set(loop.tags))
                for k in key_items:
                    if k['name'] in missing_tags:
                        if 'default' in k:
                            for row in loop.data:
                                row.append(k['default'])
                            loop.add_tag(k['name'])
                        else:
                            raise LookupError("Missing mandatory %s loop tag%s." % (missing_tags, 's' if len(missing_tags) > 1 else ''))

            if len(mand_data_names) > 0 and set(mand_data_names) & set(loop.tags) != set(mand_data_names):
                missing_tags = list(set(mand_data_names) - set(loop.tags))
                for k in key_items:
                    if k['name'] in missing_tags:
                        if 'default' in k:
                            for row in loop.data:
                                row.append(k['default'])
                            loop.add_tag(k['name'])
                        else:
                            raise LookupError("Missing mandatory %s loop tag%s." % (missing_tags, 's' if len(missing_tags) > 1 else ''))

            if len(_mand_data_names) > 0 and set(_mand_data_names) & set(loop.tags) != set(_mand_data_names):
                missing_tags = list(set(_mand_data_names) - set(loop.tags))
                for d in data_items:
                    if d['name'] in missing_tags:
                        if 'default-from' in d:
                            if d['type'] == 'index-int':
                                for l, row in enumerate(loop.data, start=1):
                                    row.append(l)
                                loop.add_tag(d['name'])
                            elif d['name'] == 'element' or d['name'] == 'Atom_type':
                                from_col = loop.tags.index(d['default-from'])
                                for row in loop.data:
                                    ref = row[from_col]
                                    if ref.startswith('H') or ref.startswith('Q') or ref.startswith('M'):
                                        row.append('H')
                                    else:
                                        row.append(ref[0])
                                loop.add_tag(d['name'])
                            elif d['name'] == 'isotope_number' or d['name'] == 'Atom_isotope_number':
                                from_col = loop.tags.index(d['default-from'])
                                for row in loop.data:
                                    ref = row[from_col]
                                    if ref.startswith('H') or ref.startswith('Q') or ref.startswith('M'):
                                        row.append(1)
                                    else:
                                        row.append(self.atom_isotopes[ref[0]][0])
                                loop.add_tag(d['name'])

            if not disallowed_tags is None:
                if len(set(loop.tags) & set(disallowed_tags)) > 0:
                    disallow_tags = list(set(loop.tags) & set(disallowed_tags))
                    raise LookupError("Disallowed %s loop tag%s exist%s." % (disallow_tags, 's' if len(disallow_tags) > 1 else '', '' if len(disallow_tags) > 1 else 's'))

            for d in data_items:
                if 'group-mandatory' in d and d['group-mandatory']:
                    name = d['name']
                    group = d['group']
                    if name in loop.tags:
                        if not group['coexist-with'] is None:
                            for c in group['coexist-with']:
                                if not c in loop.tags:
                                    missing_tags = list(set(group['coexist-with']).add(name) - set(loop.tags))
                                    raise LookupError("Missing mandatory %s loop tag%s." % (missing_tags, 's' if len(missing_tags) > 1 else ''))

                    elif not group['member-with'] is None:
                        has_member = False
                        for m in group['member-with']:
                            if m in loop.tags:
                                has_member = True
                                break
                        if not has_member:
                            missing_tags = list(set(group['member-with']).add(name) - set(loop.tags))
                            raise LookupError("Missing mandatory %s loop tag%s." % (missing_tags, 's' if len(missing_tags) > 1 else ''))

            tags = [k['name'] for k in key_items]
            for data_name in set(data_names) & set(loop.tags):
                tags.append(data_name)

            tag_len = len(tags)

            static_val = {}
            for name in tags:
                static_val[name] = None

            idx_tag_ids = set()
            for j in range(tag_len):
                name = tags[j]
                if name in key_names:
                    for k in key_items:
                        if k['name'] == name and k['type'] == 'index-int':
                            idx_tag_ids.add(j)
                else:
                    for d in data_items:
                        if d['name'] == name and d['type']== 'index-int':
                            idx_tag_ids.add(j)

            relax_key_ids = set()
            for j in range(tag_len):
                name = tags[j]
                for d in data_items:
                    if d['name'] == name and 'relax-key-if-exist' in d and d['relax-key-if-exist']:
                        relax_key_ids.add(j)

            tag_data = get_lp_tag(loop, tags)

            if test_on_index and len(idx_tag_ids) > 0:

                for l, _j in enumerate(idx_tag_ids):

                    try:
                        idxs = [int(i[_j]) for i in tag_data]

                        dup_idxs = [i for i in set(idxs) if idxs.count(i) > 1]

                        if len(dup_idxs) > 0:
                            raise KeyError("%s must be unique in loop. %s are duplicated." % (tags[_j], dup_idxs))

                    except ValueError:
                        r = {}
                        for j in range(len(loop.tags)):
                            r[loop.tags[j]] = loop.data[l][j]
                        raise ValueError("%s must be int. #_of_row %s, data_of_row %s." % (tags[_j], l + 1, r))

            if not excl_missing_data:
                for l, i in enumerate(tag_data):
                    for j in range(tag_len):
                        if i[j] in self.empty_value:
                            name = tags[j]
                            if name in key_names:
                                r = {}
                                for _j in range(len(loop.tags)):
                                    r[loop.tags[_j]] = loop.data[l][_j]
                                raise ValueError("%s must not be empty. #_of_row %s, data_of_row %s." % (name, l + 1, r))
                            else:
                                for d in data_items:
                                    if d['name'] == name and d['mandatory'] and not 'default' in d:
                                        r = {}
                                        for _j in range(len(loop.tags)):
                                            r[loop.tags[_j]] = loop.data[l][_j]
                                        raise ValueError("%s must not be empty. #_of_row %s, data_of_row %s." % (name, l + 1, r))

            if test_on_index and key_len > 0:
                keys = set()

                rechk = False

                for l, i in enumerate(tag_data):

                    key = ''
                    for j in range(key_len):
                        key += ' ' + i[j]
                    key.rstrip()

                    if key in keys:

                        relax_key = False

                        if len(relax_key_ids) > 0:
                            for j in relax_key_ids:
                                if not i[j] is self.empty_value:
                                    relax_key = True
                                    break

                        if relax_key:
                            rechk = True

                        else:
                            msg = ''
                            for j in range(key_len):
                                msg += key_names[j] + ' %s, ' % i[j]

                            idx_msg = ''

                            if len(idx_tag_ids) > 0:
                                for _j in idx_tag_ids:
                                    idx_msg += tags[_j] + ' '

                                    for _i in tag_data:
                                        _key = ''
                                        for j in range(key_len):
                                            _key += " " + _i[j]
                                            _key.rstrip()

                                        if key == _key:
                                            idx_msg += _i[_j] + ' vs '

                                    idx_msg = idx_msg[:-4] + ', '

                                idx_msg = '[Check rows of ' + idx_msg[:-2] + '] '

                            user_warn_msg += '[Multiple data] %sDuplicated rows having the following values %s exist in a loop.\n' % (idx_msg, msg.rstrip().rstrip(','))

                    else:
                        keys.add(key)

                if rechk:
                    keys = set()

                    for l, i in enumerate(tag_data):

                        key = ''
                        for j in range(key_len):
                            key += ' ' + i[j]
                        for j in relax_key_ids:
                            key += ' ' + i[j]
                        key.rstrip()

                        if key in keys:

                            msg = ''
                            for j in range(key_len):
                                msg += key_names[j] + ' %s, ' % i[j]
                            for j in relax_key_ids:
                                if not i[j] in self.empty_value:
                                    msg += tags[j] + ' %s, ' % i[j]

                            idx_msg = ''

                            if len(idx_tag_ids) > 0:
                                for _j in idx_tag_ids:
                                    idx_msg += tags[_j] + ' '

                                    for _i in tag_data:
                                        _key = ''
                                        for j in range(key_len):
                                            _key += " " + _i[j]
                                        for j in relax_key_ids:
                                            _key += " " + _i[j]
                                            _key.rstrip()

                                        if key == _key:
                                            idx_msg += _i[_j] + ' vs '

                                    idx_msg = idx_msg[:-4] + ', '

                                idx_msg = '[Check rows of ' + idx_msg[:-2] + '] '

                            user_warn_msg += '[Multiple data] %sDuplicated rows having the following values %s exist in a loop.\n' % (idx_msg, msg.rstrip().rstrip(','))

                        else:
                            keys.add(key)

            asm = [] # assembly of a loop

            for l, i in enumerate(tag_data):
                ent = {} # entity

                missing_mandatory_data = False

                for j in range(tag_len):
                    name = tags[j]
                    val = i[j]
                    if j < key_len:
                        k = key_items[j]
                        type = k['type']
                        if type == 'bool':
                            try:
                                ent[name] = val.lower() in self.true_value
                            except:
                                if excl_missing_data:
                                    missing_mandatory_data = True
                                    continue
                                else:
                                    raise ValueError("%s%s '%s' must be %s." % (get_idx_msg(idx_tag_ids, tags, ent), name, val, self.readable_item_type[type]))
                        elif type == 'int':
                            try:
                                ent[name] = int(val)
                            except:
                                if 'default-from' in k and k['default-from'] == 'self':
                                    ent[name] = self.letter_to_int(val)
                                elif 'default' in k:
                                    ent[name] = int(k['default'])
                                elif excl_missing_data:
                                    missing_mandatory_data = True
                                    continue
                                else:
                                    raise ValueError("%s%s '%s' must be %s." % (get_idx_msg(idx_tag_ids, tags, ent), name, val, self.readable_item_type[type]))
                        elif type == 'index-int' or type == 'positive-int' or type == 'positive-int-as-str':
                            try:
                                ent[name] = int(val)
                            except:
                                if 'default-from' in k and k['default-from'] == 'self':
                                    ent[name] = self.letter_to_int(val, 1)
                                elif 'default' in k:
                                    ent[name] = int(k['default'])
                                elif excl_missing_data:
                                    missing_mandatory_data = True
                                    continue
                                else:
                                    raise ValueError("%s%s '%s' must be %s." % (get_idx_msg(idx_tag_ids, tags, ent), name, val, self.readable_item_type[type]))
                            if (type == 'index-int' and ent[name] <= 0) or (type == 'positive-int' and (ent[name] < 0 or (ent[name] == 0 and 'enforce-non-zero' in k and k['enforce-non-zero']))):
                                raise ValueError("%s%s '%s' must be %s." % (get_idx_msg(idx_tag_ids, tags, ent), name, val, self.readable_item_type[type]))
                            elif ent[name] == 0 and enforce_non_zero:
                                if 'void-zero' in k:
                                    if self.replace_zero_by_null_in_case:
                                        loop.data[l][loop.tags.index(name)] = None
                                    ent[name] = None
                                else:
                                    user_warn_msg += "[Zero value error] %s%s '%s' should not be zero, as defined by %s.\n" % (get_idx_msg(idx_tag_ids, tags, ent), name, val, self.readable_item_type[type])
                            if type == 'positive-int-as-str':
                                ent[name] = str(ent[name])
                        elif type == 'pointer-index':
                            try:
                                ent[name] = int(val)
                            except:
                                if 'default-from' in k and k['default-from'] == 'self':
                                    ent[name] = self.letter_to_int(val, 1)
                                elif 'default' in k:
                                    ent[name] = int(k['default'])
                                elif excl_missing_data:
                                    missing_mandatory_data = True
                                    continue
                                else:
                                    raise ValueError("%s%s '%s' must be %s." % (get_idx_msg(idx_tag_ids, tags, ent), name, val, self.readable_item_type[type]))
                            if ent[name] <= 0:
                                raise ValueError("%s%s '%s' must be %s." % (get_idx_msg(idx_tag_ids, tags, ent), name, val, self.readable_item_type[type]))
                            elif static_val[name] is None:
                                static_val[name] = val
                            elif val != static_val[name] and test_on_index:
                                raise ValueError("%s%s %s vs %s must be %s." % (get_idx_msg(idx_tag_ids, tags, ent), name, val, static_val[name], type))
                        elif type == 'float':
                            try:
                                ent[name] = float(val)
                            except:
                                if excl_missing_data:
                                    missing_mandatory_data = True
                                    continue
                                else:
                                    raise ValueError("%s%s '%s' must be %s." % (get_idx_msg(idx_tag_ids, tags, ent), name, val, self.readable_item_type[type]))
                        elif type == 'positive-float':
                            try:
                                ent[name] = float(val)
                            except:
                                if excl_missing_data:
                                    missing_mandatory_data = True
                                    continue
                                else:
                                    raise ValueError("%s%s '%s' must be %s." % (get_idx_msg(idx_tag_ids, tags, ent), name, val, self.readable_item_type[type]))
                            if ent[name] < 0.0 or (ent[name] == 0.0 and 'enforce-non-zero' in k and k['enforce-non-zero']):
                                raise ValueError("%s%s '%s' must be %s." % (get_idx_msg(idx_tag_ids, tags, ent), name, val, self.readable_item_type[type]))
                            elif ent[name] == 0.0 and enforce_non_zero:
                                if 'void-zero' in k:
                                    if self.replace_zero_by_null_in_case:
                                        loop.data[l][loop.tags.index(name)] = None
                                    ent[name] = None
                                else:
                                    user_warn_msg += "[Zero value error] %s%s '%s' should not be zero, as defined by %s.\n" % (get_idx_msg(idx_tag_ids, tags, ent), name, val, self.readable_item_type[type])
                        elif type == 'range-float':
                            try:
                                _range = k['range']
                                ent[name] = float(val)
                            except KeyError:
                                raise ValueError('Range of key item %s is not defined' % name)
                            except:
                                if excl_missing_data:
                                    missing_mandatory_data = True
                                    continue
                                elif not enforce_range:
                                    ent[name] = None
                                    continue
                                else:
                                    user_warn_msg += "[Range value error] %s%s '%s' must be %s.\n" % (get_idx_msg(idx_tag_ids, tags, ent), name, val, self.readable_item_type[type])
                            if ('min_exclusive' in _range and _range['min_exclusive'] == 0.0 and ent[name] <= 0.0) or ('min_inclusive' in _range and _range['min_inclusive'] == 0.0 and ent[name] < 0):
                                if ent[name] < 0.0:
                                    if ('max_inclusive' in _range and abs(ent[name]) > _range['max_inclusive']) or ('max_exclusive' in _range and abs(ent[name]) >= _range['max_exclusive']) or ('enforce-sign' in k and k['enforce-sign']):
                                        if not enforce_range:
                                            ent[name] = None
                                        else:
                                            user_warn_msg += "[Range value error] %s%s '%s' must be within range %s.\n" % (get_idx_msg(idx_tag_ids, tags, ent), name, val, _range)
                                    elif enforce_sign:
                                        user_warn_msg += "[Negative value error] %s%s '%s' should not have negative value for %s, %s.\n" % (get_idx_msg(idx_tag_ids, tags, ent), name, val, self.readable_item_type[type], _range)
                                elif ent[name] == 0.0 and 'enforce-non-zero' in k and k['enforce-non-zero']:
                                    if not enforce_range:
                                        ent[name] = None
                                    else:
                                        user_warn_msg += "[Range value error] %s%s '%s' must be within range %s.\n" % (get_idx_msg(idx_tag_ids, tags, ent), name, val, _range)
                                elif ent[name] == 0.0 and enforce_non_zero:
                                    if 'void-zero' in k:
                                        if self.replace_zero_by_null_in_case:
                                            loop.data[l][loop.tags.index(name)] = None
                                        ent[name] = None
                                    else:
                                        user_warn_msg += "[Zero value error] %s%s '%s' should not be zero, as defined by %s, %s.\n" % (get_idx_msg(idx_tag_ids, tags, ent), name, val, self.readable_item_type[type], _range)
                            elif ('min_exclusive' in _range and ent[name] <= _range['min_exclusive']) or ('min_inclusive' in _range and ent[name] < _range['min_inclusive']) or ('max_inclusive' in _range and ent[name] > _range['max_inclusive']) or ('max_exclusive' in _range and ent[name] >= _range['max_exclusive']):
                                if 'void-zero' in k and ent[name] == 0.0:
                                    if self.replace_zero_by_null_in_case:
                                        loop.data[l][loop.tags.index(name)] = None
                                    ent[name] = None
                                elif not enforce_range:
                                    ent[name] = None
                                else:
                                    user_warn_msg += "[Range value error] %s%s '%s' must be within range %s.\n" % (get_idx_msg(idx_tag_ids, tags, ent), name, val, _range)
                        elif type == 'enum':
                            try:
                                enum = k['enum']
                                if not val in enum:
                                    if 'enum-alt' in k and val in k['enum-alt']:
                                        val = k['enum-alt'][val]
                                        i[j] = val
                                    elif 'enforce-enum' in k and k['enforce-enum']:
                                        if excl_missing_data:
                                            missing_mandatory_data = True
                                            continue
                                        else:
                                            raise ValueError("%s%s '%s' must be one of %s." % (get_idx_msg(idx_tag_ids, tags, ent), name, val, enum))
                                    elif enforce_enum:
                                        user_warn_msg += "[Enumeration error] %s%s '%s' should be one of %s.\n" % (get_idx_msg(idx_tag_ids, tags, ent), name, val, enum)
                                ent[name] = val
                            except KeyError:
                                raise ValueError('Enumeration of key item %s is not defined' % name)
                        elif type == 'enum-int':
                            try:
                                enum = k['enum']
                                if not int(val) in enum:
                                    if 'enforce-enum' in k and k['enforce-enum']:
                                        raise ValueError("%s%s '%s' must be one of %s." % (get_idx_msg(idx_tag_ids, tags, ent), name, val, enum))
                                    elif enforce_enum:
                                        user_warn_msg += "[Enumeration error] %s%s '%s' should be one of %s.\n" % (get_idx_msg(idx_tag_ids, tags, ent), name, val, enum)
                                ent[name] = int(val)
                            except KeyError:
                                raise ValueError('Enumeration of key item %s is not defined' % name)
                            except:
                                if excl_missing_data:
                                    missing_mandatory_data = True
                                    continue
                                else:
                                    raise ValueError("%s%s '%s' must be %s." % (get_idx_msg(idx_tag_ids, tags, ent), name, val, self.readable_item_type[type]))
                        else:
                            if val in self.empty_value:
                                missing_mandatory_data = True
                                continue
                            if 'uppercase' in k and k['uppercase']:
                                ent[name] = val.upper()
                            else:
                                ent[name] = val

                    else:
                        for d in data_items:
                            if d['name'] == name:
                                type = d['type']
                                if val in self.empty_value and ('enum' in type or (not 'default-from' in d and not 'default' in d)):
                                    if d['mandatory'] and excl_missing_data:
                                        missing_mandatory_data = True
                                        continue
                                    ent[name] = None
                                elif type == 'bool':
                                    try:
                                        ent[name] = val in self.true_value
                                    except:
                                        if excl_missing_data:
                                            ent[name] = None
                                            continue
                                        else:
                                            raise ValueError("%s%s '%s' must be %s." % (get_idx_msg(idx_tag_ids, tags, ent), name, val, self.readable_item_type[type]))
                                elif type == 'int':
                                    try:
                                        ent[name] = int(val)
                                    except:
                                        if 'default-from' in d and d['default-from'] == 'self':
                                            ent[name] = self.letter_to_int(val)
                                        elif 'default' in d:
                                            ent[name] = int(d['default'])
                                        elif excl_missing_data:
                                            ent[name] = None
                                            continue
                                        else:
                                            raise ValueError("%s%s '%s' must be %s." % (get_idx_msg(idx_tag_ids, tags, ent), name, val, self.readable_item_type[type]))
                                elif type == 'index-int' or type == 'positive-int' or type == 'positive-int-as-str':
                                    try:
                                        ent[name] = int(val)
                                    except:
                                        if 'default-from' in d and d['default-from'] == 'self':
                                            ent[name] = self.letter_to_int(val, 1)
                                        elif 'default' in d:
                                            ent[name] = int(d['default'])
                                        elif excl_missing_data:
                                            ent[name] = None
                                            continue
                                        else:
                                            raise ValueError("%s%s '%s' must be %s." % (get_idx_msg(idx_tag_ids, tags, ent), name, val, self.readable_item_type[type]))
                                    if (type == 'index-int' and ent[name] <= 0) or (type == 'positive-int' and (ent[name] < 0 or (ent[name] == 0 and 'enforce-non-zero' in d and d['enforce-non-zero']))):
                                        raise ValueError("%s%s '%s' must be %s." % (get_idx_msg(idx_tag_ids, tags, ent), name, val, self.readable_item_type[type]))
                                    elif ent[name] == 0 and enforce_non_zero:
                                        if 'void-zero' in d:
                                            if self.replace_zero_by_null_in_case:
                                                loop.data[l][loop.tags.index(name)] = None
                                            ent[name] = None
                                        else:
                                            user_warn_msg += "[Zero value error] %s%s '%s' should not be zero, as defined by %s.\n" % (get_idx_msg(idx_tag_ids, tags, ent), name, val, self.readable_item_type[type])
                                    if type == 'positive-int-as-str':
                                        ent[name] = str(ent[name])
                                elif type == 'pointer-index':
                                    try:
                                        ent[name] = int(val)
                                    except:
                                        if 'default-from' in d and d['default-from'] == 'self':
                                            ent[name] = self.letter_to_int(val, 1)
                                        elif 'default' in d:
                                            ent[name] = int(d['default'])
                                        elif excl_missing_data:
                                            ent[name] = None
                                            continue
                                        else:
                                            raise ValueError("%s%s '%s' must be %s." % (get_idx_msg(idx_tag_ids, tags, ent), name, val, self.readable_item_type[type]))
                                    if ent[name] <= 0:
                                        raise ValueError("%s%s '%s' must be %s." % (get_idx_msg(idx_tag_ids, tags, ent), name, val, self.readable_item_type[type]))
                                    elif static_val[name] is None:
                                        static_val[name] = val
                                    elif val != static_val[name] and test_on_index:
                                        raise ValueError("%s%s %s vs %s must be %s." % (get_idx_msg(idx_tag_ids, tags, ent), name, val, static_val[name], type))
                                elif type == 'float':
                                    try:
                                        ent[name] = float(val)
                                    except:
                                        if excl_missing_data:
                                            ent[name] = None
                                            continue
                                        else:
                                            raise ValueError("%s%s '%s' must be %s." % (get_idx_msg(idx_tag_ids, tags, ent), name, val, self.readable_item_type[type]))
                                elif type == 'positive-float':
                                    try:
                                        ent[name] = float(val)
                                    except:
                                        if excl_missing_data:
                                            ent[name] = None
                                            continue
                                        else:
                                            raise ValueError("%s%s '%s' must be %s." % (get_idx_msg(idx_tag_ids, tags, ent), name, val, self.readable_item_type[type]))
                                    if ent[name] < 0.0 or (ent[name] == 0.0 and 'enforce-non-zero' in d and d['enforce-non-zero']):
                                        raise ValueError("%s%s '%s' must be %s." % (get_idx_msg(idx_tag_ids, tags, ent), name, val, self.readable_item_type[type]))
                                    elif ent[name] == 0.0 and enforce_non_zero:
                                        if 'void-zero' in d:
                                            if self.replace_zero_by_null_in_case:
                                                loop.data[l][loop.tags.index(name)] = None
                                            ent[name] = None
                                        else:
                                            user_warn_msg += "[Zero value error] %s%s '%s' should not be zero, as defined by %s.\n" % (get_idx_msg(idx_tag_ids, tags, ent), name, val, self.readable_item_type[type])
                                elif type == 'range-float':
                                    try:
                                        _range = d['range']
                                        ent[name] = float(val)
                                    except KeyError:
                                        raise ValueError('Range of data item %s is not defined' % name)
                                    except:
                                        if excl_missing_data:
                                            ent[name] = None
                                            continue
                                        elif not enforce_range:
                                            ent[name] = None
                                            continue
                                        else:
                                            user_warn_msg += "[Range value error] %s%s '%s' must be %s.\n" % (get_idx_msg(idx_tag_ids, tags, ent), name, val, self.readable_item_type[type])
                                    if ('min_exclusive' in _range and _range['min_exclusive'] == 0.0 and ent[name] <= 0.0) or ('min_inclusive' in _range and _range['min_inclusive'] == 0.0 and ent[name] < 0):
                                        if ent[name] < 0.0:
                                            if ('max_inclusive' in _range and abs(ent[name]) > _range['max_inclusive']) or ('max_exclusive' in _range and abs(ent[name]) >= _range['max_exclusive']) or ('enforce-sign' in d and d['enforce-sign']):
                                                if not enforce_range:
                                                    ent[name] = None
                                                else:
                                                    user_warn_msg += "[Range value error] %s%s '%s' must be within range %s.\n" % (get_idx_msg(idx_tag_ids, tags, ent), name, val, _range)
                                            elif enforce_sign:
                                                user_warn_msg += "[Negative value error] %s%s '%s' should not have negative value for %s, %s.\n" % (get_idx_msg(idx_tag_ids, tags, ent), name, val, self.readable_item_type[type], _range)
                                        elif ent[name] == 0.0 and 'enforce-non-zero' in d and d['enforce-non-zero']:
                                            if not enforce_range:
                                                ent[name] = None
                                            else:
                                                user_warn_msg += "[Range value error] %s%s '%s' must be within range %s.\n" % (get_idx_msg(idx_tag_ids, tags, ent), name, val, _range)
                                        elif ent[name] == 0.0 and enforce_non_zero:
                                            if 'void-zero' in d:
                                                if self.replace_zero_by_null_in_case:
                                                    loop.data[l][loop.tags.index(name)] = None
                                                ent[name] = None
                                            else:
                                                user_warn_msg += "[Zero value error] %s%s '%s' should not be zero, as defined by %s, %s.\n" % (get_idx_msg(idx_tag_ids, tags, ent), name, val, self.readable_item_type[type], _range)
                                    elif ('min_exclusive' in _range and ent[name] <= _range['min_exclusive']) or ('min_inclusive' in _range and ent[name] < _range['min_inclusive']) or ('max_inclusive' in _range and ent[name] > _range['max_inclusive']) or ('max_exclusive' in _range and ent[name] >= _range['max_exclusive']):
                                        if 'void-zero' in d and ent[name] == 0.0:
                                            if self.replace_zero_by_null_in_case:
                                                loop.data[l][loop.tags.index(name)] = None
                                            ent[name] = None
                                        elif not enforce_range:
                                            ent[name] = None
                                        else:
                                            user_warn_msg += "[Range value error] %s%s '%s' must be within range %s.\n" % (get_idx_msg(idx_tag_ids, tags, ent), name, val, _range)
                                elif type == 'enum':
                                    try:
                                        enum = d['enum']
                                        if not val in enum:
                                            if 'enum-alt' in d and val in d['enum-alt']:
                                                val = d['enum-alt'][val]
                                                i[j] = val
                                            elif 'enforce-enum' in d and d['enforce-enum']:
                                                if excl_missing_data:
                                                    ent[name] = None
                                                    continue
                                                else:
                                                    raise ValueError("%s%s '%s' must be one of %s." % (get_idx_msg(idx_tag_ids, tags, ent), name, val, enum))
                                            elif enforce_enum:
                                                user_warn_msg += "[Enumeration error] %s%s '%s' should be one of %s.\n" % (get_idx_msg(idx_tag_ids, tags, ent), name, val, enum)
                                        ent[name] = val
                                    except KeyError:
                                        raise ValueError('Enumeration of data item %s is not defined' % name)
                                elif type == 'enum-int':
                                    try:
                                        enum = d['enum']
                                        if not int(val) in enum:
                                            if 'enforce-enum' in d and d['enforce-enum']:
                                                raise ValueError("%s%s '%s' must be one of %s." % (get_idx_msg(idx_tag_ids, tags, ent), name, val, enum))
                                            elif enforce_enum:
                                                user_warn_msg += "[Enumeration error] %s%s '%s' should be one of %s.\n" % (get_idx_msg(idx_tag_ids, tags, ent), name, val, enum)
                                        ent[name] = int(val)
                                    except KeyError:
                                        raise ValueError('Enumeration of data item %s is not defined' % name)
                                    except:
                                        if excl_missing_data:
                                            ent[name] = None
                                            continue
                                        else:
                                            raise ValueError("%s%s '%s' must be %s." % (get_idx_msg(idx_tag_ids, tags, ent), name, val, self.readable_item_type[type]))
                                else:
                                    if 'uppercase' in d and d['uppercase']:
                                        ent[name] = val.upper()
                                    else:
                                        ent[name] = val

                for d in data_items:
                    if 'group-mandatory' in d and d['group-mandatory']:
                        name = d['name']
                        group = d['group']
                        if name in ent and not ent[name] is None:
                            if not group['coexist-with'] is None:
                                has_coexist = True
                                for c in group['coexist-with']:
                                    if not c in ent or ent[c] is None:
                                        raise ValueError("%sOne of data item %s must not be empty for a row having %s '%s'." % (get_idx_msg(idx_tag_ids, tags, ent), c, name, ent[name]))

                            if 'smaller-than' in group and not group['smaller-than'] is None:
                                for s in group['smaller-than']:
                                    if s in ent and not ent[s] is None:
                                        if ent[name] < ent[s]:
                                            if 'circular-shift' in group:
                                                ent[s] -= abs(group['circular-shift'])
                                            if ent[name] < ent[s]:
                                                raise ValueError("%sData item %s '%s' must be larger than %s '%s'." % (get_idx_msg(idx_tag_ids, tags, ent), name, ent[name], s, ent[s]))

                            if 'larger-than' in group and not group['larger-than'] is None:
                                for l in group['larger-than']:
                                    if l in ent and not ent[l] is None:
                                        if ent[name] > ent[l]:
                                            if 'circular-shift' in group:
                                                ent[l] += abs(group['circular-shift'])
                                            if ent[name] > ent[l]:
                                                raise ValueError("%sData item %s '%s' must be smaller than %s '%s'." % (get_idx_msg(idx_tag_ids, tags, ent), name, ent[name], l, ent[l]))

                            if 'not-equal-to' in group and not group['not-equal-to'] is None:
                                for n in group['not-equal-to']:
                                    if n in ent and not ent[n] is None:
                                        if ent[name] == ent[n]:
                                            raise ValueError("%sData item %s '%s' must not be equal to %s '%s'." % (get_idx_msg(idx_tag_ids, tags, ent), name, ent[name], n, ent[n]))

                        elif not group['member-with'] is None:
                            has_member = False
                            for m in group['member-with']:
                                if m in ent and not ent[m] is None:
                                    has_member = True
                                    break
                            if not has_member:
                                member = set(group['member-with'])
                                member.add(name)
                                raise ValueError("%sOne of data items %s must not be empty." % (get_idx_msg(idx_tag_ids, tags, ent), member))

                if missing_mandatory_data:
                    continue

                asm.append(ent)

            data.append(asm)

        if len(user_warn_msg) > 0:
            raise UserWarning(user_warn_msg)

        if len(data) == 0:
            data.append([])

        return data

    def letter_to_int(self, code, min=0):
        """ Return digit from a given chain code.
        """

        alphabet = 'abcdefghijklmnopqrstuvwxyz'

        unit = 1
        ret = 0

        for c in ''.join(reversed(code.lower())):

            if c.isdigit():
                ret += unit * int(c)
            elif c.isalpha():
                ret += unit * (alphabet.index(c) + 1)
            else:
                continue

            unit *= 27

        return ret if ret >= min else min

    def index_to_letter(self, index):
        """ Return chain code from a given index (0 based).
        """

        if index < 0:
            return '.'
        if index > 19683:
            index = index % 19683

        if index < 27:
            return str(chr(65 + index))
        elif index < 729:
            return str(chr(64 + (index // 27))) + str(chr(65 + (index % 27)))
        else:
            return str(chr(64 + (index // 729))) + str(chr(64 + ((index % 729) // 27))) + str(chr(65 + (index % 27)))

    def get_conflict_id(self, star_data, lp_category, key_items):
        """ Return list of row ID of duplicated/conflicted rows except for rows of the first occurrence.
            @author: Masashi Yokochi
            @return: list of duplicated/conflicted row IDs in reverse order for each loop
        """

        try:
            loops = star_data.get_loops_by_category(lp_category)
        except: # AttributeError:
            try:
                loops = [star_data.get_loop_by_category(lp_category)]
            except: # AttributeError:
                loops = [star_data]

        data = [] # data of all loops

        key_names = [k['name'] for k in key_items]
        uppercases = [('uppercase' in k and k['uppercase']) for k in key_items]

        key_len = len(key_items)

        for loop in loops:

            if len(key_names) > 0 and set(key_names) & set(loop.tags) != set(key_names):
                missing_tags = list(set(key_names) - set(loop.tags))
                for k in key_items:
                    if k['name'] in missing_tags:
                        if 'default' in k:
                            for row in loop.data:
                                row.append(k['default'])
                            loop.add_tag(k['name'])
                        else:
                            raise LookupError("Missing mandatory %s loop tag%s." % (missing_tags, 's' if len(missing_tags) > 1 else ''))

            keys = set()
            dup_ids = set()

            for l, i in enumerate(get_lp_tag(loop, key_names)):

                key = ''
                for j in range(key_len):
                    key += ' ' + (i[j].upper() if uppercases[j] else i[j])
                key.rstrip()

                if key in keys:
                    dup_ids.add(l)

                else:
                    keys.add(key)

            data.append(sorted(list(dup_ids), reverse=True))

        return data

    def get_conflict_id_set(self, star_data, lp_category, key_items):
        """ Return list of row ID set of redundant/inconsistent rows.
            @author: Masashi Yokochi
            @return: list of redundant/inconsistent row ID set for each loop
        """

        try:
            loops = star_data.get_loops_by_category(lp_category)
        except: # AttributeError:
            try:
                loops = [star_data.get_loop_by_category(lp_category)]
            except: # AttributeError:
                loops = [star_data]

        data = [] # data of all loops

        key_names = [k['name'] for k in key_items]
        uppercases = [('uppercase' in k and k['uppercase']) for k in key_items]

        key_len = len(key_items)

        for loop in loops:

            if len(key_names) > 0 and set(key_names) & set(loop.tags) != set(key_names):
                missing_tags = list(set(key_names) - set(loop.tags))
                raise LookupError("Missing mandatory %s loop tag%s." % (missing_tags, 's' if len(missing_tags) > 1 else ''))

            tag_data = get_lp_tag(loop, key_names)

            keys = set()
            dup_ids = set()

            for l, i in enumerate(tag_data):

                key = ''
                for j in range(key_len):
                    key += ' ' + (i[j].upper() if uppercases[j] else i[j])
                key.rstrip()

                if key in keys:
                    dup_ids.add(l)

                else:
                    keys.add(key)

            conflict_id = sorted(list(dup_ids), reverse=True)

            if len(conflict_id) == 0:
                data.append(None)

            else:
                conflict_id_set = []

                for l in conflict_id:

                    key = ''
                    for j in range(key_len):
                        key += ' ' + tag_data[l][j]
                        key.rstrip()

                    id_set = [l]

                    for m in range(l):

                        _key = ''
                        for j in range(key_len):
                            _key += ' ' + tag_data[m][j]
                            _key.rstrip()

                        if key == _key:
                            id_set.append(m)

                            if m in conflict_id:
                                conflict_id.remove(m)

                    if len(id_set) > 1:
                        conflict_id_set.insert(0, sorted(id_set))

                data.append(conflict_id_set)

        return data

    def get_conflict_atom_id(self, star_data, file_type, lp_category, key_items):
        """ Return list of row ID of rows includes self atoms.
            @author: Masashi Yokochi
            @return: list of duplicated/conflicted row IDs in reverse order for each loop
        """

        try:
            loops = star_data.get_loops_by_category(lp_category)
        except: # AttributeError:
            try:
                loops = [star_data.get_loop_by_category(lp_category)]
            except: # AttributeError:
                loops = [star_data]

        data = [] # data of all loops

        key_names = [k['name'] for k in key_items]
        uppercases = [('uppercase' in k and k['uppercase']) for k in key_items]

        key_len = len(key_items)

        for loop in loops:

            if len(key_names) > 0 and set(key_names) & set(loop.tags) != set(key_names):
                missing_tags = list(set(key_names) - set(loop.tags))
                for k in key_items:
                    if k['name'] in missing_tags:
                        if 'default' in k:
                            for row in loop.data:
                                row.append(k['default'])
                            loop.add_tag(k['name'])
                        else:
                            raise LookupError("Missing mandatory %s loop tag%s." % (missing_tags, 's' if len(missing_tags) > 1 else ''))

            atom_keys = self.get_atom_keys(loop.get_tag_names(), file_type)

            len_atom_keys = len(atom_keys)

            dup_ids = set()

            for l, i in enumerate(get_lp_tag(loop, key_names)):

                for j in range(0, len_atom_keys - 2):
                    atom_key_j = atom_keys[j]
                    chain_id_j = i[key_names.index(atom_key_j['chain_tag'])]
                    seq_id_j = i[key_names.index(atom_key_j['seq_tag'])]
                    atom_id_j = i[key_names.index(atom_key_j['atom_tag'])]

                    for k in range(j + 1, len_atom_keys - 1):
                        atom_key_k = atom_keys[k]
                        chain_id_k = i[key_names.index(atom_key_k['chain_tag'])]
                        seq_id_k = i[key_names.index(atom_key_k['seq_tag'])]
                        atom_id_k = i[key_names.index(atom_key_k['atom_tag'])]

                        if chain_id_j == chain_id_k and seq_id_j == seq_id_k and atom_id_j == atom_id_k:
                            dup_ids.add(l)
                            break

            data.append(sorted(list(dup_ids), reverse=True))

        return data

    def check_sf_tag(self, star_data, file_type, category, tag_items, allowed_tags=None,
                     enforce_non_zero=False, enforce_sign=False, enforce_range=False, enforce_enum=False):
        """ Extract saveframe tags with sanity check.
            @author: Masashi Yokochi
            @return: list of extracted saveframe tags
        """

        user_warn_msg = ''

        item_types = ('str', 'bool', 'int', 'positive-int', 'positive-int-as-str', 'float', 'positive-float', 'range-float', 'enum', 'enum-int')

        tag_names = [t['name'] for t in tag_items]
        mand_tag_names = [t['name'] for t in tag_items if t['mandatory']]

        for t in tag_items:
            if not t['type'] in item_types:
                raise TypeError("Type %s of tag item %s must be one of %s." % (t['type'], t['name'], item_types))

        if not allowed_tags is None:

            if (set(tag_names) | set(allowed_tags)) != set(allowed_tags):
                raise LookupError("Tag items %s must not exists." % ((set(tag_names) | set(allowed_tags)) - set(allowed_tags)))

            for t in tag_items:
                if 'group-mandatory' in t and t['group-mandatory']:
                    group = t['group']
                    if not group['member-with'] is None:
                        for m in group['member-with']:
                            if not m in allowed_tags:
                                raise ValueError("Member tag item %s of %s must exists in allowed tags." % (m, t['name']))
                    if not group['coexist-with'] is None:
                        for c in group['coexist-with']:
                            if not c in allowed_tags:
                                raise ValueError("Coexisting tag item %s of %s must exists in allowed tags." % (c, t['name']))
                    if 'smaller-than' in group and group['smaller-than']:
                        for s in group['smaller-than']:
                            if not s in allowed_tags:
                                raise ValueError("Smaller tag item %s of %s must exists in allowed tags." % (s, t['name']))
                    if 'larger-than' in group and group['larger-than']:
                        for l in group['larger-than']:
                            if not l in allowed_tags:
                                raise ValueError("Larger tag item %s of %s must exists in allowed tags." % (l, t['name']))
                    if 'not-equal-to' in group and group['not-equal-to']:
                        for l in group['not-equal-to']:
                            if not l in allowed_tags:
                                raise ValueError("Nonequal tag item %s of %s must exists in allowed tags." % (l, t['name']))

        sf_tags = {i[0]:i[1] for i in star_data.tags}

        if len(mand_tag_names) > 0 and set(mand_tag_names) & set(sf_tags.keys()) != set(mand_tag_names):
            missing_tags = list(set(mand_tag_names) - set(sf_tags.keys()))
            raise LookupError("Missing mandatory %s saveframe tag%s." % (missing_tags, 's' if len(missing_tags) > 1 else ''))

        for t in tag_items:
            if 'group-mandatory' in t and t['group-mandatory']:
                name = t['name']
                group = t['group']
                if name in sf_tags.keys():
                    if not group['coexist-with'] is None:
                        for c in group['coexist-with']:
                            if not c in sf_tags.keys():
                                missing_tags = list(set(group['coexist-with']).add(name) - set(sf_tags.keys()))
                                raise LookupError("Missing mandatory %s saveframe tag%s." % (missing_tags, 's' if len(missing_tags) > 1 else ''))

                elif not group['member-with'] is None:
                    has_member = False
                    for m in group['member-with']:
                        if m in sf_tags.keys():
                            has_member = True
                            break
                    if not has_member:
                        missing_tags = list(set(group['member-with']).add(name) - set(sf_tags.keys()))
                        raise LookupError("Missing mandatory %s saveframe tag%s." % (missing_tags, 's' if len(missing_tags) > 1 else ''))

        for name, val in sf_tags.items():
            if val in self.empty_value:
                for t in tag_items:
                    if t['name'] == name and t['mandatory'] and not 'default' in t:
                        raise ValueError("%s must not be empty." % name)

        ent = {} # entity

        for name, val in sf_tags.items():
            for t in tag_items:
                if t['name'] == name:
                    type = t['type']
                    if val in self.empty_value and not 'enum' in type and not 'default-from' in t and not 'default' in t:
                       ent[name] = None
                    elif type == 'bool':
                        try:
                            ent[name] = val in self.true_value
                        except:
                            raise ValueError("%s '%s' must be %s." % (name, val, self.readable_item_type[type]))
                    elif type == 'int':
                        try:
                            ent[name] = int(val)
                        except:
                            if 'default-from' in t and t['default-from'] == 'self':
                                ent[name] = self.letter_to_int(val)
                            elif 'default' in t:
                                ent[name] = int(t['default'])
                            else:
                                raise ValueError("%s '%s' must be %s." % (name, val, self.readable_item_type[type]))
                    elif type == 'positive-int' or type == 'positive-int-as-str':
                        try:
                            ent[name] = int(val)
                        except:
                            if 'default-from' in t and t['default-from'] == 'self':
                                ent[name] = self.letter_to_int(val, 1)
                            elif 'default' in t:
                                ent[name] = int(t['default'])
                            else:
                                raise ValueError("%s '%s' must be %s." % (name, val, self.readable_item_type[type]))
                        if ent[name] < 0 or (ent[name] == 0 and 'enforce-non-zero' in t and t['enforce-non-zero']):
                            raise ValueError("%s '%s' must be %s." % (name, val, self.readable_item_type[type]))
                        elif ent[name] == 0 and enforce_non_zero:
                            if 'void-zero' in t:
                                if self.replace_zero_by_null_in_case:
                                    star_data.tags[sf_tags.keys().index(name)][1] = None
                                ent[name] = None
                            else:
                                user_warn_msg += "[Zero value error] %s '%s' should not be zero, as defined by %s.\n" % (name, val, self.readable_item_type[type])
                        if type == 'positive-int-as-str':
                            ent[name] = str(ent[name])
                    elif type == 'float':
                        try:
                            ent[name] = float(val)
                        except:
                            raise ValueError("%s '%s' must be %s." % (name, val, self.readable_item_type[type]))
                    elif type == 'positive-float':
                        try:
                            ent[name] = float(val)
                        except:
                            raise ValueError("%s '%s' must be %s." % (name, val, self.readable_item_type[type]))
                        if ent[name] < 0.0 or (ent[name] == 0.0 and 'enforce-non-zero' in t and t['enforce-non-zero']):
                            raise ValueError("%s '%s' must be %s." % (name, val, self.readable_item_type[type]))
                        elif ent[name] == 0.0 and enforce_non_zero:
                            if 'void-zero' in t:
                                if self.replace_zero_by_null_in_case:
                                    star_data.tags[sf_tags.keys().index(name)][1] = None
                                ent[name] = None
                            else:
                                user_warn_msg += "[Zero value error] %s '%s' should not be zero, as defined by %s.\n" % (name, val, self.readable_item_type[type])
                    elif type == 'range-float':
                        try:
                            _range = t['range']
                            ent[name] = float(val)
                        except KeyError:
                            raise ValueError('Range of tag item %s is not defined.' % name)
                        except:
                            if not enforce_range:
                                ent[name] = None
                                continue
                            else:
                                user_warn_msg += "[Range value error] %s '%s' must be %s.\n" % (name, val, self.readable_item_type[type])
                        if ('min_exclusive' in _range and _range['min_exclusive'] == 0.0 and ent[name] <= 0.0) or ('min_inclusive' in _range and _range['min_inclusive'] == 0.0 and ent[name] < 0):
                            if ent[name] < 0.0:
                                if ('max_inclusive' in _range and abs(ent[name]) > _range['max_inclusive']) or ('max_exclusive' in _range and abs(ent[name]) >= _range['max_exclusive']) or ('enforce-sign' in t and t['enforce-sign']):
                                    if not enforce_range:
                                        ent[name] = None
                                    else:
                                        user_warn_msg += "[Range value error] %s '%s' must be within range %s.\n" % (name, val, _range)
                                elif enforce_sign:
                                    user_warn_msg += "[Negative value error] %s '%s' should not have negative value for %s, %s.\n" % (name, val, self.readable_item_type[type], _range)
                            elif ent[name] == 0.0 and 'enforce-non-zero' in t and t['enforce-non-zero']:
                                if not enforce_range:
                                    ent[name] = None
                                else:
                                    user_warn_msg += "[Range value error] %s '%s' must be within range %s.\n" % (name, val, _range)
                            elif ent[name] == 0.0 and enforce_non_zero:
                                if 'void-zero' in t:
                                    if self.replace_zero_by_null_in_case:
                                        star_data.tags[sf_tags.keys().index(name)][1] = None
                                    ent[name] = None
                                else:
                                    user_warn_msg += "[Zero value error] %s '%s' should not be zero, as defined by %s, %s.\n" % (name, val, self.readable_item_type[type], _range)
                        elif ('min_exclusive' in _range and ent[name] <= _range['min_exclusive']) or ('min_inclusive' in _range and ent[name] < _range['min_inclusive']) or ('max_inclusive' in _range and ent[name] > _range['max_inclusive']) or ('max_exclusive' in _range and ent[name] >= _range['max_exclusive']):
                            if 'void-zero' in t and ent[name] == 0.0:
                                if self.replace_zero_by_null_in_case:
                                    star_data.tags[sf_tags.keys().index(name)][1] = None
                                ent[name] = None
                            elif not enforce_range:
                                ent[name] = None
                            else:
                                user_warn_msg += "[Range value error] %s '%s' must be within range %s.\n" % (name, val, _range)
                    elif type == 'enum':
                        if val in self.empty_value:
                            val = '?' # '.' raises internal error in NmrDpUtility
                        try:
                            enum = t['enum']
                            if not val in enum:
                                if 'enum-alt' in t and val in t['enum-alt']:
                                    tagNames = [_t[0] for _t in star_data.tags]
                                    itCol = tagNames.index(name)
                                    itName = '_' + category + '.' + t['name']
                                    if val == '?' and enforce_enum:
                                        if self.is_mandatory_tag(itName, file_type):
                                            user_warn_msg += "[Enumeration error] The mandatory type %s '%s' is missing and the type must be one of %s. '%s' will be given unless you would like to fix the type and re-upload the %s file.\n" % (itName, val, enum, t['enum-alt'][val], file_type.upper())
                                            val = t['enum-alt'][val]
                                            star_data.tags[itCol][1] = val
                                        else:
                                            user_warn_msg += "[Enumeration error] %s '%s' should be one of %s. The type may be filled with either 'undefined' or estimated value unless you would like to fix the type and re-upload the %s file.\n" % (name, val, enum, file_type.upper())
                                    else:
                                        val = t['enum-alt'][val]
                                        star_data.tags[itCol][1] = val
                                elif 'enforce-enum' in t and t['enforce-enum']:
                                    raise ValueError("%s '%s' must be one of %s." % (name, val, enum))
                                elif enforce_enum:
                                    user_warn_msg += "[Enumeration error] %s '%s' should be one of %s.\n" % (name, val, enum)
                            ent[name] = None if val in self.empty_value else val
                        except KeyError:
                            raise ValueError('Enumeration of tag item %s is not defined.' % name)
                    elif type == 'enum-int':
                        try:
                            enum = t['enum']
                            if not int(val) in enum:
                                if 'enforce-enum' in t and t['enforce-enum']:
                                    raise ValueError("%s '%s' must be one of %s." % (name, val, enum))
                                elif enforce_enum:
                                    user_warn_msg += "[Enumeration error] %s '%s' should be one of %s.\n" % (name, val, enum)
                            ent[name] = int(val)
                        except KeyError:
                            raise ValueError('Enumeration of tag item %s is not defined.' % name)
                        except:
                            raise ValueError("%s '%s' must be %s." % (name, val, self.readable_item_type[type]))
                    else:
                        ent[name] = val

            for t in tag_items:
                if 'group-mandatory' in t and t['group-mandatory']:
                    name = t['name']
                    group = t['group']
                    if name in ent and not ent[name] is None:
                        if not group['coexist-with'] is None:
                            has_coexist = True
                            for c in group['coexist-with']:
                                if not c in ent or ent[c] is None:
                                    raise ValueError("One of tag item %s must not be empty due to %s '%s'." % (c, name, ent[name]))

                        if 'smaller-than' in group and not group['smaller-than'] is None:
                            for s in group['smaller-than']:
                                if s in ent and not ent[s] is None:
                                    if ent[name] < ent[s]:
                                        if 'circular-shift' in group:
                                            ent[s] -= abs(group['circular-shift'])
                                        if ent[name] < ent[s]:
                                            raise ValueError("Tag item %s '%s' must be larger than %s '%s'." % (name, ent[name], s, ent[s]))

                        if 'larger-than' in group and not group['larger-than'] is None:
                            for l in group['larger-than']:
                                if l in ent and not ent[l] is None:
                                    if ent[name] > ent[l]:
                                        if 'circular-shift' in group:
                                            ent[l] += abs(group['circular-shift'])
                                        if ent[name] > ent[l]:
                                            raise ValueError("Tag item %s '%s' must be smaller than %s '%s'." % (name, ent[name], l, ent[l]))

                        if 'not-equal-to' in group and not group['not-equal-to'] is None:
                            for n in group['not-equal-to']:
                                if n in ent and not ent[n] is None:
                                    if ent[name] == ent[n]:
                                        raise ValueError("Tag item %s '%s' must not be equal to %s '%s'." % (name, ent[name], n, ent[n]))

                    elif not group['member-with'] is None:
                        has_member = False
                        for m in group['member-with']:
                            if m in ent and not ent[m] is None:
                                has_member = True
                                break
                        if not has_member:
                            member = set(group['member-with'])
                            member.add(name)
                            raise ValueError("One of tag items %s must not be empty." % member)

        if len(user_warn_msg) > 0:
            raise UserWarning(user_warn_msg)

        return ent

    def __updateChemCompDict(self, comp_id):
        """ Update CCD information for a given comp_id.
            @return: True for successfully update CCD information or False for the case a given comp_id does not exist in CCD
        """

        if comp_id != self.__last_comp_id:
            self.__last_comp_id_test = False if '_' in comp_id else self.__ccR.setCompId(comp_id)
            self.__last_comp_id = comp_id

            if self.__last_comp_id_test:
                self.__last_chem_comp_dict = self.__ccR.getChemCompDict()
                self.__last_chem_comp_atoms = self.__ccR.getAtomList()
                self.__last_chem_comp_bonds = self.__ccR.getBonds()

        return self.__last_comp_id_test

    def validate_comp_atom(self, comp_id, atom_id):
        """ Validate atom_id of comp_id.
            @change: support non-standard residue by Masashi Yokochi
            @return: True for valid atom_id of comp_id, False otherwise
        """

        comp_id = comp_id.upper()

        if comp_id in self.empty_value:
            return False

        atoms = []

        if comp_id in self.atomDict:
            atoms = self.atomDict[comp_id]

        else:
            self.__updateChemCompDict(comp_id)

            if self.__last_comp_id_test: # matches with comp_id in CCD
                atoms = [a[self.__cca_atom_id] for a in self.__last_chem_comp_atoms]
            else:
                return False

        return atom_id in atoms

    def validate_atom(self, star_data, lp_category='Atom_chem_shift', comp_id='Comp_ID', atom_id='Atom_ID'):
        """ Validate atom_id in a given loop against CCD.
            @change: support non-standard residue by Masashi Yokochi
            @return: list of valid row data
        """

        try:
            loop_data = star_data.get_loops_by_category(lp_category)
        except: # AttributeError:
            try:
                loop_data = [star_data.get_loop_by_category(lp_category)]
            except: # AttributeError:
                loop_data = [star_data]

        valid_row = []

        for i in loop_data:

            try:

                data = get_lp_tag(i, [comp_id, atom_id])

                for j in data:

                    _comp_id = j[0].upper()
                    _atom_id = j[1].upper()

                    if _comp_id in self.empty_value:
                        valid_row.append(j)

                    elif _comp_id in self.atomDict:
                        if not _atom_id in self.atomDict[_comp_id]:
                            valid_row.append(j)

                    else:
                        self.__updateChemCompDict(_comp_id)

                        if self.__last_comp_id_test: # matches with comp_id in CCD
                            if not _atom_id in [a[self.__cca_atom_id] for a in self.__last_chem_comp_atoms]:
                                valid_row.append(j)
                        else:
                            valid_row.append(j)

            except ValueError:
                logging.error('Missing one of data items %s' % (comp_id, atom_id))

        return valid_row

    def get_star_tag(self, nef_tag):
        """ Return NMR-STAR saveframe/loop tags corresponding to NEF tag.
            @change: rename the original get_nmrstar_tag() to get_star_tag() by Masashi Yokochi
            @return: NMR-STAR's author tag and data tag corresponding to a given NEF tag
        """

        try:

            n = self.tagMap[0].index(nef_tag)

            return self.tagMap[1][n], self.tagMap[2][n], n # author tag, data tag

        except ValueError:
            return None, None, None

    def get_nef_tag(self, star_tag):
        """ Return NEF saveframe/loop tags corresponding to NMR-STAR tag.
            @author: Masashi Yokochi
            @return: NEF data tag corresponding to a given NMR-STAR tag
        """

        try:

            n = self.tagMap[2].index(star_tag) # NEF has no auth tags

            return self.tagMap[0][n], n

        except ValueError: # None for NMR-STAR specific tag
            return None, None

    def get_star_loop_tags(self, nef_loop_tags):
        """ Return list of NMR-STAR loop tags corresponding to NEF loop tags.
            @change: rename the original get_nmrstar_loop_tags() to get_star_loop_tags() by Masashi Yokochi
            @return: list of NMR-STAR loop tags corresponding to given NEF loop tags
        """

        out_tag_w_ordinal = {}
        out_tag = []

        for j in nef_loop_tags:

            auth_tag, data_tag, ordinal = self.get_star_tag(j)

            if auth_tag is None:
                continue

            out_tag_w_ordinal[data_tag] = ordinal

            if auth_tag != data_tag:
                out_tag_w_ordinal[auth_tag] = ordinal + 100

        out_tag = [t[0] for t in sorted(out_tag_w_ordinal.items(), key=lambda x:x[1])]

        lp_category = nef_loop_tags[0].split('.')[0]

        if lp_category == '_nef_sequence':
            out_tag.append('_Chem_comp_assembly.Assembly_ID')

        elif lp_category == '_nef_covalent_links':
            out_tag.insert(0, '_Bond.ID')
            out_tag.insert(1, '_Bond.Type')
            out_tag.insert(2, '_Bond.Value_order')
            out_tag.append('_Bond.Assembly_ID')

        elif lp_category == '_nef_chemical_shift':
            out_tag.insert(0, '_Atom_chem_shift.ID')
            out_tag.insert(9, '_Atom_chem_shift.Ambiguity_code')
            out_tag.insert(10, '_Atom_chem_shift.Ambiguity_set_ID')

            if self.insert_original_pdb_cs_items:
                out_tag.append('_Atom_chem_shift.Original_PDB_strand_ID')
                out_tag.append('_Atom_chem_shift.Original_PDB_residue_no')
                out_tag.append('_Atom_chem_shift.Original_PDB_residue_name')
                out_tag.append('_Atom_chem_shift.Original_PDB_atom_name')

            out_tag.append('_Atom_chem_shift.Details')
            out_tag.append('_Atom_chem_shift.Assigned_chem_shift_list_ID')

        elif lp_category == '_nef_distance_restraint':
            out_tag.append('_Gen_dist_constraint.Member_logic_code')
            # out_tag.append('_Gen_dist_constraint.Details')
            out_tag.append('_Gen_dist_constraint.Gen_dist_constraint_list_ID')

        elif lp_category == '_nef_dihedral_restraint':
            # out_tag.append('_Torsion_angle_constraint.Details')
            out_tag.append('_Torsion_angle_constraint.Torsion_angle_constraint_list_ID')

        elif lp_category == '_nef_rdc_restraint':
            # out_tag.append('_RDC_constraint.Details'
            out_tag.append('_RDC_constraint.RDC_constraint_list_ID')

        elif lp_category == '_nef_peak':
            out_tag.append('_Peak_row_format.Details')
            out_tag.append('_Peak_row_format.Spectral_peak_list_ID')

        elif lp_category == '_nef_spectrum_dimension':
            out_tag.append('_Spectral_dim.Spectral_peak_list_ID')

        elif lp_category == '_nef_spectrum_dimension_transfer':
            out_tag.append('_Spectral_dim_transfer.Spectral_peak_list_ID')

        return out_tag

    def get_nef_loop_tags(self, star_loop_tags):
        """ Return list of NEF loop tags corresponding to NMR-STAR loop tags.
            @author: Masashi Yokochi
            @return: list of NEF loop tags corresponding to given NMR-STAR loop tags
        """

        out_tag_w_ordinal = {}
        out_tag = []

        for j in star_loop_tags:

            data_tag, ordinal = self.get_nef_tag(j)

            if not data_tag is None:
                out_tag_w_ordinal[data_tag] = ordinal

        out_tag = [t[0] for t in sorted(out_tag_w_ordinal.items(), key=lambda x:x[1])]

        if len(out_tag) > 0 and out_tag[0].startswith('_nef_sequence.') and not '_nef_sequence.index' in out_tag:
            out_tag.insert(0, '_nef_sequence.index')

        return out_tag

    def get_star_atom(self, comp_id, nef_atom, details=None, leave_unmatched=True):
        """ Return list of instanced atom_id of a given NEF atom (including wildcard codes) and its ambiguity code.
            @change: support non-standard residue by Masashi Yokochi
            @change: rename the original get_nmrstar_atom() to get_star_atom() by Masashi Yokochi
            @return: list of instanced atom_id of a given NEF atom, ambiguity_code, and description
        """

        comp_id = comp_id.upper()

        if comp_id in self.empty_value:
            return [], None, None

        comp_code = self.get_one_letter_code(comp_id)

        atom_list = []
        ambiguity_code = 1

        atoms = []

        if comp_id in self.atomDict:
            atoms = self.atomDict[comp_id]

        else:
            self.__updateChemCompDict(comp_id)

            if self.__last_comp_id_test: # matches with comp_id in CCD
                atoms = [a[self.__cca_atom_id] for a in self.__last_chem_comp_atoms]

            else:

                if leave_unmatched:
                    details = 'Unknown non-standard residue %s found.' % comp_id
                else:
                    logging.critical('Unknown non-standard residue %s found.' % comp_id)

        try:

            ref_atom = re.findall(r'(\S+)([xyXY])([%*])$|(\S+)([%*])$|(\S+)([xyXY]$)', nef_atom)[0]

            atm_set = [ref_atom.index(i) for i in ref_atom if i != '']

            pattern = None

            if atm_set == [0, 1, 2]: # endswith [xyXY][%*]

                atom_type = ref_atom[0]
                xy_code = ref_atom[1].lower()

                len_atom_type = len(atom_type)

                pattern = re.compile(r'%s\S\d+' % (atom_type))

                alist2 = [i for i in atoms if re.search(pattern, i)]

                xid = sorted(set([int(i[len_atom_type]) for i in alist2]))

                if xy_code == 'x':
                    atom_list = [i for i in alist2 if int(i[len_atom_type]) == xid[0]]
                else:
                    atom_list = [i for i in alist2 if int(i[len_atom_type]) == xid[1]]

                ambiguity_code = self.__csStat.getMaxAmbigCodeWoSetId(comp_id, atom_list[0])


            elif atm_set == [3, 4]: # endswith [%*] but neither [xyXY][%*]

                atom_type = ref_atom[3]
                wc_code = ref_atom[4]

                if wc_code == '%':
                    if comp_code == 'X':
                        pattern = re.compile(r'%s\S?$' % atom_type)
                    else:
                        pattern = re.compile(r'%s\d+' % atom_type)
                elif wc_code == '*':
                    pattern = re.compile(r'%s\S+' % atom_type)
                else:
                    logging.critical('Invalid NEF atom nomenclature %s found.' % nef_atom)

                atom_list = [i for i in atoms if re.search(pattern, i)]

                methyl_atoms = self.__csStat.getMethylAtoms(comp_id)

                ambiguity_code = 1 if atom_list[0] in methyl_atoms else self.__csStat.getMaxAmbigCodeWoSetId(comp_id, atom_list[0])

            elif atm_set == [5, 6]: # endswith [xyXY]

                atom_type = ref_atom[5]
                xy_code = ref_atom[6].lower()

                pattern = re.compile(r'%s[^\']+' % atom_type)

                atom_list = [i for i in atoms if re.search(pattern, i)]

                atom_list_len = len(atom_list)

                if atom_list_len != 2:
                    if atom_list_len > 2:
                        atom_list = []
                elif xy_code == 'y':
                    atom_list = atom_list[-1:]
                elif xy_code == 'x':
                    atom_list = atom_list[:1]
                else:
                    logging.critical('Invalid NEF atom nomenclature %s found.' % nef_atom)

                ambiguity_code = self.__csStat.getMaxAmbigCodeWoSetId(comp_id, atom_list[0])

            else:
                logging.critical('Invalid NEF atom nomenclature %s found.' % nef_atom)

        except IndexError:
            pass

        if len(atom_list) == 0:

            if nef_atom == 'HN' and self.__csStat.getTypeOfCompId(comp_id)[0]:
                return self.get_star_atom(comp_id, 'H', 'HN converted to H.' if leave_unmatched else None, leave_unmatched)

            if comp_id in self.atomDict:

                methyl_atoms = self.__csStat.getMethylAtoms(comp_id)

                if comp_code != 'X' and not nef_atom.endswith('%') and comp_code != 'X' and not nef_atom.endswith('*') and nef_atom + '1' in methyl_atoms:
                    return self.get_star_atom(comp_id, nef_atom + '%', ('%s converted to %s%%.' % (nef_atom, nef_atom)) if leave_unmatched else None, leave_unmatched)

                if nef_atom[-1].lower() == 'x' or nef_atom[-1].lower() == 'y' and nef_atom[:-1] + '1' in methyl_atoms:
                    return self.get_star_atom(comp_id, nef_atom[:-1] + '%', ('%s converted to %s%%.' % (nef_atom, nef_atom[:-1])) if leave_unmatched else None, leave_unmatched)

                if ((comp_code != 'X' and nef_atom[-1] == '%') or nef_atom[-1] == '*') and not (nef_atom[:-1] + '1' in methyl_atoms) and\
                    len(nef_atom) > 2 and (nef_atom[-2].lower() == 'x' or nef_atom[-2].lower() == 'y'):
                    return self.get_star_atom(comp_id, nef_atom[:-2] + ('1' if nef_atom[-2].lower() == 'x' else '2') + '%', ('%s converted to %s%%.' % (nef_atom, nef_atom[:-2] + ('1' if nef_atom[-2].lower() == 'x' else '2'))) if leave_unmatched else None, leave_unmatched)

            if nef_atom in atoms:
                atom_list.append(nef_atom)

            elif leave_unmatched:
                atom_list.append(nef_atom)
                ambiguity_code = None
                if details is None:
                    details = '%s is invalid atom_id in comp_id %s.' % (nef_atom, comp_id)

        return atom_list, ambiguity_code, details

    def get_nef_atom(self, comp_id, star_atom_list, details={}, leave_unmatched=True):
        """ Return list of all instanced atom_id of given NMR-STAR atoms with ambiguity code and CS value in a given comp_id.
            @author: Masashi Yokochi
            @return: list of instanced atom_id of given NMR-STAR atoms, descriptions, and atom conversion dictionary for conversion of other loops
        """

        comp_id = comp_id.upper()

        if comp_id in self.empty_value:
            return [], None, None

        atom_list = []
        atom_id_map = {}

        self.__updateChemCompDict(comp_id)

        atoms = []

        if comp_id in self.atomDict:
            atoms = self.atomDict[comp_id]

        elif self.__last_comp_id_test: # matches with comp_id in CCD
            atoms = [a[self.__cca_atom_id] for a in self.__last_chem_comp_atoms]

        methyl_atoms = self.__csStat.getMethylAtoms(comp_id)

        proc_atom_list = set()

        try:

            for a in star_atom_list:

                atom_id = a['atom_id']
                _ambig_code = a['ambig_code']
                _value = a['value']

                if atom_id in proc_atom_list:
                    continue

                if not atom_id in atoms:

                    _atom_id = None

                    if atom_id == 'HN' and self.__csStat.getTypeOfCompId(comp_id)[0]:
                        _atom_id = 'H'
                    elif atom_id.startswith('QQ'):
                        _atom_id = 'H' + atom_id[2:] + '%'
                    elif atom_id.startswith('QR'):
                        qr_atoms = sorted(set([atom_id[:-1] + '%' for atom_id in self.__csStat.getAromaticAtoms(comp_id) if atom_id[0] == 'H' and self.__csStat.getMaxAmbigCodeWoSetId(comp_id, atom_id) == 3]))
                        if len(qr_atoms) > 0:
                            for qr_atom in qr_atoms:
                                _a = copy.copy(a)
                                _a['atom_id'] = qr_atom
                                star_atom_list.append(_a)
                            star_atom_list.remove(a)
                            return self.get_nef_atom(comp_id, star_atom_list, details, leave_unmatched)
                    elif atom_id.startswith('Q') or atom_id.startswith('M'):
                        _atom_id = 'H' + atom_id[1:] + '%'
                    elif atom_id + '2' in self.__csStat.getAllAtoms(comp_id):
                        _atom_id = atom_id + '%'

                    if _atom_id is None:

                        if leave_unmatched:
                            atom_list.append(atom_id)
                            if not self.__last_comp_id_test:
                                details[atom_id] = 'Unknown non-standard residue %s found.' % comp_id
                            else:
                                details[atom_id] = '%s is invalid atom_id in comp_id %s.' % (atom_id, comp_id)
                            atom_id_map[atom_id] = atom_id
                        else:
                            if not self.__last_comp_id_test:
                                logging.critical('Unknown non-standard residue %s found.' % comp_id)
                            else:
                                logging.critical('Invalid atom nomenclature %s found.' % atom_id)

                        continue

                    atom_id = _atom_id

                if _ambig_code in self.empty_value:
                    ambig_code = self.__csStat.getMaxAmbigCodeWoSetId(comp_id, atom_id)
                    if atom_id.endswith("'"):
                        ambig_code = 1
                else:
                    ambig_code = int(_ambig_code)

                if not ambig_code in self.bmrb_ambiguity_codes:

                    if leave_unmatched:
                        atom_list.append(atom_id)
                        details[atom_id] = '%s has invalid ambiguity code %s.' % (atom_id, ambig_code)
                        atom_id_map[atom_id] = atom_id
                    else:
                        logging.critical('Invalid ambiguity code %s for atom_id %s found.' % (ambig_code, atom_id))

                    continue

                elif ambig_code == 1:

                    if atom_id[0] == 'H' and atom_id in methyl_atoms:

                        methyl_c, methyl_h_list = self.get_group(comp_id, atom_id)

                        if methyl_h_list is None:

                            if leave_unmatched:
                                atom_list.append(atom_id)
                                details[atom_id] = '%s is invalid atom_id in comp_id %s.' % (atom_id, comp_id)
                                atom_id_map[atom_id] = atom_id
                            else:
                                logging.critical('Invalid atom nomenclature %s found.' % atom_id)

                        else:

                            has_methyl_proton = len([_a['atom_id'] for _a in star_atom_list if _a['atom_id'] in methyl_h_list and (_a['value'] is None or _a['value'] == _value)]) == 3

                            if has_methyl_proton:

                                name_len = [len(n) for n in methyl_h_list]
                                max_len = max(name_len)
                                min_len = min(name_len)

                                if max_len == min_len or len(atom_id) == max_len:
                                    nef_atom = atom_id[:-1] + '%'
                                else: # For example, HEM HM[A-D]
                                    nef_atom = atom_id + '%'

                                atom_list.append(nef_atom)
                                details[nef_atom] = None
                                for b in methyl_h_list:
                                    atom_id_map[b] = nef_atom
                                    proc_atom_list.add(b)

                            else:
                                atom_list.append(atom_id)
                                details[atom_id] = None
                                atom_id_map[atom_id] = atom_id

                    else:
                        atom_list.append(atom_id)
                        details[atom_id] = None
                        atom_id_map[atom_id] = atom_id

                elif ambig_code == 2:

                    if self.__csStat.getMaxAmbigCodeWoSetId(comp_id, atom_id) != 2:

                        if leave_unmatched:
                            atom_list.append(atom_id)
                            details[atom_id] = '%s has invalid ambiguity code %s.' % (atom_id, ambig_code)
                            atom_id_map[atom_id] = atom_id
                        else:
                            logging.critical('Invalid ambiguity code %s for atom_id %s found.' % (ambig_code, atom_id))

                    elif atom_id[0] == 'H':

                        if atom_id in methyl_atoms:

                            methyl_c, methyl_h_list = self.get_group(comp_id, atom_id)
                            methyl_c_2, methyl_h_list_2 = self.get_geminal_group(comp_id, methyl_c)

                            if methyl_h_list is None:

                                if leave_unmatched:
                                    atom_list.append(atom_id)
                                    details[atom_id] = '%s is invalid atom_id in comp_id %s.' % (atom_id, comp_id)
                                    atom_id_map[atom_id] = atom_id
                                else:
                                    logging.critical('Invalid atom nomenclature %s found.' % atom_id)

                            else:

                                has_methyl_proton = len([_a['atom_id'] for _a in star_atom_list if _a['atom_id'] in methyl_h_list and (_a['value'] is None or _a['value'] == _value)]) == 3
                                has_methyl_proton_2 = not methyl_h_list_2 is None and len([_a['atom_id'] for _a in star_atom_list if _a['atom_id'] in methyl_h_list_2 and (_a['value'] is None or _a['value'] == _value)]) == 3

                                nef_atom_prefix = 'x'
                                nef_atom_prefix_2 = 'y'

                                if has_methyl_proton and has_methyl_proton_2:

                                    methyl_proton_value = next((_a['value'] for _a in star_atom_list if _a['atom_id'] in methyl_h_list), None)
                                    methyl_proton_value_2 = next((_a['value'] for _a in star_atom_list if _a['atom_id'] in methyl_h_list_2), None)

                                    if not methyl_proton_value is None and not methyl_proton_value_2 is None and float(methyl_proton_value_2) < float(methyl_proton_value):
                                        nef_atom_prefix = 'y'
                                        nef_atom_prefix_2 = 'x'

                                    nef_atom = methyl_h_list[0][:-2] + nef_atom_prefix + '%'

                                    atom_list.append(nef_atom)
                                    details[nef_atom] = None
                                    for b in methyl_h_list:
                                        atom_id_map[b] = nef_atom
                                        proc_atom_list.add(b)

                                    nef_atom_2 = methyl_h_list_2[0][:-2] + nef_atom_prefix_2 + '%'

                                    atom_list.append(nef_atom_2)
                                    details[nef_atom_2] = None
                                    for b in methyl_h_list_2:
                                        atom_id_map[b] = nef_atom_2
                                        proc_atom_list.add(b)

                                elif has_methyl_proton:

                                    nef_atom_prefix = 'x' if methyl_h_list[0] < methyl_h_list_2[0] else 'y'
                                    nef_atom = methyl_h_list[0][:-2] + nef_atom_prefix + '%'

                                    atom_list.append(nef_atom)
                                    details[nef_atom] = None
                                    for b in methyl_h_list:
                                        atom_id_map[b] = nef_atom
                                        proc_atom_list.add(b)

                                else:
                                    atom_list.append(atom_id)
                                    details[atom_id] = None
                                    atom_id_map[atom_id] = atom_id

                        else:

                            geminal_dummy, geminal_h_list = self.get_group(comp_id, atom_id)

                            if geminal_h_list is None:

                                if leave_unmatched:
                                    atom_list.append(atom_id)
                                    details[atom_id] = '%s is invalid atom_id in comp_id %s.' % (atom_id, comp_id)
                                    atom_id_map[atom_id] = atom_id
                                else:
                                    logging.critical('Invalid atom nomenclature %s found.' % atom_id)

                            else:

                                has_geminal_proton = len([_a['atom_id'] for _a in star_atom_list if _a['atom_id'] in geminal_h_list]) == 2

                                nef_atom_prefix = 'x'
                                nef_atom_prefix_2 = 'y'

                                if atom_id.endswith("'"):

                                    atom_list.append(atom_id)
                                    details[atom_id] = None
                                    atom_id_map[atom_id] = atom_id

                                elif has_geminal_proton:

                                    geminal_proton_value = next((_a['value'] for _a in star_atom_list if _a['atom_id'] == geminal_h_list[0]), None)
                                    geminal_proton_value_2 = next((_a['value'] for _a in star_atom_list if _a['atom_id'] == geminal_h_list[1]), None)

                                    if not geminal_proton_value is None and not geminal_proton_value_2 is None and float(geminal_proton_value_2) < float(geminal_proton_value):
                                        nef_atom_prefix = 'y'
                                        nef_atom_prefix_2 = 'x'

                                    nef_atom = geminal_h_list[0][:-1] + nef_atom_prefix

                                    atom_list.append(nef_atom)
                                    details[nef_atom] = None
                                    atom_id_map[geminal_h_list[0]] = nef_atom

                                    nef_atom_2 = geminal_h_list[1][:-1] + nef_atom_prefix_2

                                    atom_list.append(nef_atom_2)
                                    details[nef_atom_2] = None
                                    atom_id_map[geminal_h_list[1]] = nef_atom_2

                                    for b in geminal_h_list:
                                        proc_atom_list.add(b)

                                else:

                                    nef_atom_prefix = 'x' if (atom_id == geminal_h_list[0] and geminal_h_list[0] < geminal_h_list[1]) or (atom_id == geminal_h_list[1] and geminal_h_list[1] < geminal_h_list[0]) else 'y'
                                    nef_atom = atom_id[:-1] + nef_atom_prefix

                                    atom_list.append(nef_atom)
                                    details[nef_atom] = None
                                    atom_id_map[atom_id] = nef_atom

                    else:

                        atom_id_2, h_list_dummy = self.get_geminal_group(comp_id, atom_id)

                        if atom_id_2 is None:

                            if leave_unmatched:
                                atom_list.append(atom_id)
                                details[atom_id] = '%s is invalid atom_id in comp_id %s.' % (atom_id, comp_id)
                                atom_id_map[atom_id] = atom_id
                            else:
                                logging.critical('Invalid atom nomenclature %s found.' % atom_id)

                        else:

                            has_atom_id_2 = len([_a['atom_id'] for _a in star_atom_list if _a['atom_id'] == atom_id_2]) == 1

                            nef_atom_prefix = 'x'
                            nef_atom_prefix_2 = 'y'

                            if has_atom_id_2:

                                atom_id_value = next((_a['value'] for _a in star_atom_list if _a['atom_id'] == atom_id), None)
                                atom_id_value_2 = next((_a['value'] for _a in star_atom_list if _a['atom_id'] == atom_id_2), None)

                                if not atom_id_value is None and not atom_id_value_2 is None and float(atom_id_value_2) < float(atom_id_value):
                                    nef_atom_prefix = 'y'
                                    nef_atom_prefix_2 = 'x'

                                nef_atom = atom_id[:-1] + nef_atom_prefix

                                atom_list.append(nef_atom)
                                details[nef_atom] = None
                                atom_id_map[atom_id] = nef_atom
                                proc_atom_list.add(atom_id)

                                nef_atom_2 = atom_id_2[:-1] + nef_atom_prefix_2

                                atom_list.append(nef_atom_2)
                                details[nef_atom_2] = None
                                atom_id_map[atom_id_2] = nef_atom_2
                                proc_atom_list.add(atom_id_2)

                            else:

                                nef_atom_prefix = 'x' if atom_id < atom_id_2 else 'y'
                                nef_atom = atom_id[:-1] + nef_atom_prefix

                                atom_list.append(nef_atom)
                                details[nef_atom] = None
                                atom_id_map[atom_id] = nef_atom

                elif ambig_code == 3:

                    if self.__csStat.getMaxAmbigCodeWoSetId(comp_id, atom_id) != 3:

                        if leave_unmatched:
                            atom_list.append(atom_id)
                            details[atom_id] = '%s has invalid ambiguity code %s.' % (atom_id, ambig_code)
                            atom_id_map[atom_id] = atom_id
                        else:
                            logging.critical('Invalid ambiguity code %s for atom_id %s found.' % (ambig_code, atom_id))

                    else:

                        atom_id_2 = self.__csStat.getGeminalAtom(comp_id, atom_id)

                        if atom_id_2 is None:

                            if leave_unmatched:
                                atom_list.append(atom_id)
                                details[atom_id] = '%s is invalid atom_id in comp_id %s.' % (atom_id, comp_id)
                                atom_id_map[atom_id] = atom_id
                            else:
                                logging.critical('Invalid atom nomenclature %s found.' % atom_id)

                        else:

                            has_atom_id_2 = len([_a['atom_id'] for _a in star_atom_list if _a['atom_id'] == atom_id_2]) == 1

                            nef_atom_prefix = 'x'
                            nef_atom_prefix_2 = 'y'

                            if has_atom_id_2:

                                atom_id_value = next((_a['value'] for _a in star_atom_list if _a['atom_id'] == atom_id), None)
                                atom_id_value_2 = next((_a['value'] for _a in star_atom_list if _a['atom_id'] == atom_id_2), None)

                                if not atom_id_value is None and not atom_id_value_2 is None and float(atom_id_value_2) < float(atom_id_value):
                                    nef_atom_prefix = 'y'
                                    nef_atom_prefix_2 = 'x'

                                nef_atom = atom_id[:-1] + nef_atom_prefix

                                atom_list.append(nef_atom)
                                details[nef_atom] = None
                                atom_id_map[atom_id] = nef_atom
                                proc_atom_list.add(atom_id)

                                nef_atom_2 = atom_id_2[:-1] + nef_atom_prefix_2

                                atom_list.append(nef_atom_2)
                                details[nef_atom_2] = None
                                atom_id_map[atom_id_2] = nef_atom_2
                                proc_atom_list.add(atom_id_2)

                            else:

                                nef_atom_prefix = 'x' if atom_id < atom_id_2 else 'y'
                                nef_atom = atom_id[:-1] + nef_atom_prefix

                                atom_list.append(nef_atom)
                                details[nef_atom] = None
                                atom_id_map[atom_id] = nef_atom

                else:

                    atom_list.append(atom_id)
                    details[atom_id] = None
                    atom_id_map[atom_id] = atom_id

                proc_atom_list.add(atom_id)

        except KeyError:
            pass

        return atom_list, details, atom_id_map

    def get_group(self, comp_id, atom_id):
        """ Return heavy atom name and list of proton names bonded to the heavy atom.
            @author: Masashi Yokochi
            @return: heavy atom name and list of proton names
        """

        self.__updateChemCompDict(comp_id)

        if not self.__last_comp_id_test or atom_id is None or not atom_id[0] in ['H', 'C', 'N', 'O']:
            return None, None

        try:

            ccb = next(b for b in self.__last_chem_comp_bonds
                       if (b[self.__ccb_atom_id_1] == atom_id and (atom_id[0] == 'H' or b[self.__ccb_atom_id_2][0] == 'H')) or
                          (b[self.__ccb_atom_id_2] == atom_id and (atom_id[0] == 'H' or b[self.__ccb_atom_id_1][0] == 'H')))

            hvy_col = self.__ccb_atom_id_1 if ccb[self.__ccb_atom_id_2 if atom_id[0] == 'H' else self.__ccb_atom_id_1] == atom_id else self.__ccb_atom_id_2
            pro_col = self.__ccb_atom_id_2 if self.__ccb_atom_id_1 == hvy_col else self.__ccb_atom_id_1

            hvy = ccb[hvy_col]

            return hvy, [b[pro_col] for b in self.__last_chem_comp_bonds if b[hvy_col] == hvy and b[pro_col][0] == 'H']

        except StopIteration:
            return None, None

    def get_geminal_group(self, comp_id, atom_id):
        """ Return geminal heavy atom and list of proton names bonded to the geminal heavy atom.
            @author: Masashi Yokochi
            @return: geminal heavy atom name and list of geminal proton names
        """

        self.__updateChemCompDict(comp_id)

        if not self.__last_comp_id_test or atom_id is None or not atom_id[0] in ['H', 'C', 'N', 'O']:
            return None, None

        atom_id, h_list = self.get_group(comp_id, atom_id)

        if atom_id is None:
            return None, None

        h_list_len = len(h_list)

        try:

            ccb = next(b for b in self.__last_chem_comp_bonds if (b[self.__ccb_atom_id_2] == atom_id and b[self.__ccb_atom_id_1][0] != 'H') or
                                                                 (b[self.__ccb_atom_id_1] == atom_id and b[self.__ccb_atom_id_2][0] != 'H'))

            hvy_conn = ccb[self.__ccb_atom_id_1 if ccb[self.__ccb_atom_id_2] == atom_id else self.__ccb_atom_id_2]

            hvy_2 = next(c[self.__ccb_atom_id_1 if c[self.__ccb_atom_id_2] == hvy_conn else self.__ccb_atom_id_2]
                         for c in self.__last_chem_comp_bonds\
                         if (c[self.__ccb_atom_id_2] == hvy_conn and c[self.__ccb_atom_id_1] != atom_id and c[self.__ccb_atom_id_1][0] != 'H' and\
                             not self.get_group(comp_id, c[self.__ccb_atom_id_1])[1] is None and\
                             len(self.get_group(comp_id, c[self.__ccb_atom_id_1])[1]) == h_list_len) or
                            (c[self.__ccb_atom_id_1] == hvy_conn and c[self.__ccb_atom_id_2] != atom_id and c[self.__ccb_atom_id_2][0] != 'H' and\
                             not self.get_group(comp_id, c[self.__ccb_atom_id_2])[1] is None and\
                             len(self.get_group(comp_id, c[self.__ccb_atom_id_2])[1]) == h_list_len))

            return self.get_group(comp_id, hvy_2)

        except StopIteration:
            return None, None

    def nef2star_seq_row(self, nef_tags, star_tags, loop_data, report=None):
        """ Translate rows of data in sequence loop from NEF into NMR-STAR.
            @change: rename the original translate_seq_row() to nef2star_seq_row() by Masashi Yokochi
            @param nef_tags: list of NEF tags
            @param star_tags: list of NMR-STAR tags
            @param loop_data: loop data of NEF
            @param report: NMR data processing report
            @return: rows of NMR-STAR, rows of _Entity_deleted_atom loop (residue variants), otherwise None
        """

        out_row = []
        aux_row = []

        chain_index = nef_tags.index('_nef_sequence.chain_code')
        seq_index = nef_tags.index('_nef_sequence.sequence_code')
        comp_index = nef_tags.index('_nef_sequence.residue_name')

        for _star_chain, nef_chain in enumerate(self.authChainId, start=1):

            seq_list = sorted(set([int(i[seq_index]) for i in loop_data if i[chain_index] == nef_chain]))

            if len(seq_list) == 0:
                continue

            cif_chain = None
            if not report is None:
                seq_align = report.getSequenceAlignmentWithNmrChainId(nef_chain)
                if not seq_align is None:
                    cif_chain = seq_align['test_chain_id']

            if not cif_chain is None:
                _star_chain = str(self.letter_to_int(cif_chain))

            offset = None

            for _nef_seq in seq_list:

                _cif_seq = None
                if not cif_chain is None:
                    try:
                        _cif_seq = seq_align['test_seq_id'][seq_align['ref_seq_id'].index(_nef_seq)]
                    except:
                        pass

                nef_seq = str(_nef_seq)

                if offset is None:
                    offset = 1 - _nef_seq

                _star_seq = _nef_seq + offset

                i = next(i for i in loop_data if i[chain_index] == nef_chain and i[seq_index] == nef_seq)

                out = [None] * len(star_tags)

                variant = None

                for j in nef_tags:

                    auth_tag, data_tag, _ = self.get_star_tag(j)

                    if auth_tag is None:
                        continue

                    data = i[nef_tags.index(j)]

                    if j == '_nef_sequence.chain_code':
                        if cif_chain is None:
                            out[star_tags.index(auth_tag)] = data
                        else:
                            out[star_tags.index(auth_tag)] = cif_chain
                    elif j == '_nef_sequence.sequence_code':
                        if _cif_seq is None:
                            out[star_tags.index(auth_tag)] = data
                        else:
                            out[star_tags.index(auth_tag)] = str(_cif_seq)
                    else:
                        out[star_tags.index(auth_tag)] = data
                        if j == '_nef_sequence.residue_variant' and not data in self.empty_value and '-' in data:
                            variant = data

                    if auth_tag != data_tag:

                        data_index = star_tags.index(data_tag)

                        if j == '_nef_sequence.chain_code':
                            out[data_index] = _star_chain
                        elif j == '_nef_sequence.sequence_code':
                            out[data_index] = _star_seq
                        elif data in self.nef_boolean:
                            out[data_index] = 'yes' if data in self.true_value else 'no'
                        else:
                            out[data_index] = data

                out_row.append(out)

                self.authSeqMap[(nef_chain, _nef_seq)] = (_star_chain, _star_seq)
                self.selfSeqMap[(nef_chain, _nef_seq)] = (nef_chain if cif_chain is None else cif_chain,
                                                          _nef_seq if _cif_seq is None else _cif_seq)

                if not variant is None:
                    aux = [None] * len(self.entity_del_atom_row)
                    for l, aux_tag in enumerate(self.entity_del_atom_row):
                        if aux_tag == 'Entity_assembly_ID':
                            aux[l] = _star_chain
                        elif aux_tag == 'Comp_index_ID':
                            aux[l] = _star_seq
                        elif aux_tag == 'Auth_entity_assembly_ID':
                            aux[l] = nef_chain
                        elif aux_tag == 'Auth_seq_ID':
                            aux[l] = _nef_seq
                        elif aux_tag == 'Comp_ID' or aux_tag == 'Auth_comp_ID':
                            aux[l] = i[comp_index]

                    for _variant in variant.split(','):
                        _variant_ = _variant.strip(' ')
                        if _variant_.startswith('-'):
                            atom_list = self.get_star_atom(i[comp_index], _variant_[1:])[0]
                            if len(atom_list) > 0:
                                for atom in atom_list:
                                    _aux = copy.copy(aux)
                                    for l, aux_tag in enumerate(self.entity_del_atom_row):
                                        if aux_tag == 'ID':
                                            _aux[l] = len(aux_row) + 1
                                        elif aux_tag == 'Atom_ID':
                                            _aux[l] = atom
                                        elif aux_tag == 'Auth_atom_ID':
                                            _aux[l] = _variant_[1:]

                                    aux_row.append(_aux)

        return out_row, None if len(aux_row) == 0 else aux_row

    def star2nef_seq_row(self, star_tags, nef_tags, loop_data, report=None, entity_del_atom_loop=None):
        """ Translate rows of data in sequence loop from NMR-STAR into NEF.
            @author: Masashi Yokochi
            @param star_tags: list of NMR-STAR tags
            @param nef_tags: list of NEF tags
            @param loop_data: loop data of NMR-STAR
            @param report: NMR data processing report
            @return: rows of NEF
        """

        out_row = []

        has_nef_index = '_Chem_comp_assembly.NEF_index' in star_tags

        chain_index = star_tags.index('_Chem_comp_assembly.Entity_assembly_ID')
        seq_index = star_tags.index('_Chem_comp_assembly.Comp_index_ID')
        comp_index = star_tags.index('_Chem_comp_assembly.Comp_ID')

        if not entity_del_atom_loop is None:
            aux_tags = entity_del_atom_loop.get_tag_names()
            aux_chain_index = aux_tags.index('_Entity_deleted_atom.Entity_assembly_ID')
            aux_seq_index = aux_tags.index('_Entity_deleted_atom.Comp_index_ID')
            aux_atom_index = aux_tags.index('_Entity_deleted_atom.Atom_ID')

        self.star2nef_chain_mapping = {}
        self.star2cif_chain_mapping = {}

        seq_list = {}

        for cid, star_chain in enumerate(self.authChainId):

            seq_list[star_chain] = sorted(set([int(i[seq_index]) for i in loop_data if i[chain_index] == star_chain]))

            if len(seq_list[star_chain]) == 0:
                continue

            nef_chain = self.index_to_letter(cid)

            self.star2nef_chain_mapping[star_chain] = nef_chain

        if not report is None:

            for star_chain in self.authChainId:

                if len(seq_list[star_chain]) == 0:
                    continue

                cif_chain = None
                if not report is None:
                    seq_align = report.getSequenceAlignmentWithNmrChainId(star_chain)
                    if not seq_align is None:
                        cif_chain = seq_align['test_chain_id']

                self.star2cif_chain_mapping[star_chain] = cif_chain

            for star_chain in self.star2nef_chain_mapping.keys():
                if self.star2cif_chain_mapping[star_chain] is None:
                    for _star_chain in [k for k, v in self.star2nef_chain_mapping.items() if v != self.star2nef_chain_mapping[star_chain]]:
                        if not self.star2nef_chain_mapping[_star_chain] in self.star2cif_chain_mapping.values():
                            self.star2cif_chain_mapping[star_chain] = self.star2nef_chain_mapping[_star_chain]
                            break

            for star_chain in self.star2nef_chain_mapping.keys():
                if self.star2cif_chain_mapping[star_chain] is None:
                    self.star2cif_chain_mapping[star_chain] = self.star2nef_chain_mapping[star_chain]

        index = 1

        for star_chain in self.authChainId:

            _star_chain = int(star_chain)

            if len(seq_list[star_chain]) == 0:
                continue

            nef_chain = self.star2nef_chain_mapping[star_chain]

            cif_chain = None
            if not report is None:
                cif_chain = self.star2cif_chain_mapping[star_chain]
                seq_align = report.getSequenceAlignmentWithNmrChainId(star_chain)

            offset = None

            for _star_seq in seq_list[star_chain]:

                _cif_seq = None
                if not cif_chain is None:
                    try:
                        _cif_seq = seq_align['test_seq_id'][seq_align['ref_seq_id'].index(_star_seq)]
                    except:
                        pass

                star_seq = str(_star_seq)

                if offset is None:
                    offset = 1 - _star_seq

                _nef_seq = (_star_seq + offset) if _cif_seq is None else _cif_seq

                i = next(i for i in loop_data if i[chain_index] == star_chain and i[seq_index] == star_seq)

                out = [None] * len(nef_tags)

                if not has_nef_index:
                    out[0] = index

                for j in star_tags:

                    nef_tag, _ = self.get_nef_tag(j)

                    if nef_tag is None:
                        continue

                    data = i[star_tags.index(j)]

                    data_index = nef_tags.index(nef_tag)

                    if nef_tag == '_nef_sequence.chain_code':
                        if cif_chain is None:
                            out[data_index] = nef_chain
                        else:
                            out[data_index] = cif_chain
                    elif nef_tag == '_nef_sequence.sequence_code':
                        if _cif_seq is None:
                            out[data_index] = _nef_seq
                        else:
                            out[data_index] = _cif_seq
                    elif data in self.star_boolean:
                        out[data_index] = 'true' if data in self.true_value else 'false'
                    elif nef_tag == '_nef_sequence.residue_variant':
                        if not data in self.empty_value or entity_del_atom_loop is None:
                            out[data_index] = data
                        else:
                            star_atom_list = [{'atom_id': k[aux_atom_index], 'ambig_code': None, 'value': None} for k in entity_del_atom_loop if k[aux_chain_index] == star_chain and k[aux_seq_index] == star_seq]
                            if len(star_atom_list) == 0:
                                out[data_index] = data
                            else:
                                atom_list = self.get_nef_atom(i[comp_index], star_atom_list)[0]
                                if len(atom_list) == 0:
                                    out[data_index] = data
                                else:
                                    _atom_list = ['-' + a for a in atom_list]
                                    out[data_index] = ','.join(_atom_list)

                    else:
                        out[data_index] = data

                index += 1

                out_row.append(out)

                self.authSeqMap[(_star_chain, _star_seq)] = (nef_chain, _nef_seq)
                self.selfSeqMap[(_star_chain, _star_seq)] = (_star_chain if cif_chain is None else cif_chain,
                                                             _star_seq if _cif_seq is None else _cif_seq)

        return out_row

    def nef2star_bond_row(self, nef_tags, star_tags, loop_data):
        """ Translate rows of data in bond loop from NEF into NMR-STAR.
            @author: Masashi Yokochi
            @param nef_tags: list of NEF tags
            @param star_tags: list of NMR-STAR tags
            @param loop_data: loop data of NEF
            @return: rows of NMR-STAR
        """

        out_row = []

        nef_comp_index_1 = nef_tags.index('_nef_covalent_links.residue_name_1')
        nef_atom_index_1 = nef_tags.index('_nef_covalent_links.atom_name_1')
        nef_comp_index_2 = nef_tags.index('_nef_covalent_links.residue_name_2')
        nef_atom_index_2 = nef_tags.index('_nef_covalent_links.atom_name_2')

        if '_Bond.ID' in star_tags:
            star_id_index = star_tags.index('_Bond.ID')
        else:
            star_id_index = -1

        star_type_index = star_tags.index('_Bond.Type')
        star_value_order_index = star_tags.index('_Bond.Value_order')

        seq_ident_tags = self.get_seq_ident_tags(nef_tags, 'nef')

        for tag in seq_ident_tags:
            chain_tag = tag['chain_tag']
            seq_tag = tag['seq_tag']

            if chain_tag.endswith('_1'):
                chain_tag_1 = chain_tag
                seq_tag_1 = seq_tag
            else:
                chain_tag_2 = chain_tag
                seq_tag_2 = seq_tag

        key_indices = [star_tags.index(j) for j in ['_Bond.Entity_assembly_ID_1', '_Bond.Comp_index_ID_1', '_Bond.Atom_ID_1',
                                                    '_Bond.Entity_assembly_ID_2', '_Bond.Comp_index_ID_2', '_Bond.Atom_ID_2']]

        index = 1

        for i in loop_data:

            buf_row = []

            tag_map = {}
            self_tag_map = {}

            for tag in seq_ident_tags:
                chain_tag = tag['chain_tag']
                seq_tag = tag['seq_tag']

                nef_chain = i[nef_tags.index(chain_tag)]
                _nef_seq = i[nef_tags.index(seq_tag)]
                if type(_nef_seq) == str and not _nef_seq in self.empty_value:
                    _nef_seq = int(_nef_seq)

                seq_key = (nef_chain, _nef_seq)

                try:
                    tag_map[chain_tag], tag_map[seq_tag] = self.authSeqMap[seq_key]
                    self_tag_map[chain_tag], self_tag_map[seq_tag] = self.selfSeqMap[seq_key]
                except KeyError:
                    tag_map[chain_tag] = nef_chain
                    tag_map[seq_tag] = _nef_seq
                    self_tag_map[chain_tag] = nef_chain
                    self_tag_map[seq_tag] = _nef_seq

            intra_residue = i[nef_tags.index(chain_tag_1)] == i[nef_tags.index(chain_tag_2)] and i[nef_tags.index(seq_tag_1)] == i[nef_tags.index(seq_tag_2)]

            atom_list_1, ambiguity_code_1, details_1 = self.get_star_atom(i[nef_comp_index_1], i[nef_atom_index_1])
            atom_list_2, ambiguity_code_2, details_2 = self.get_star_atom(i[nef_comp_index_2], i[nef_atom_index_2])

            for k in atom_list_1:

                for l in atom_list_2:

                    if intra_residue and l == k:
                        continue

                    buf = [None] * len(star_tags)

                    for j in nef_tags:

                        auth_tag, data_tag, _ = self.get_star_tag(j)

                        if auth_tag is None:
                            continue

                        data = i[nef_tags.index(j)]

                        if 'chain_code' in j or 'sequence_code' in j:
                            buf[star_tags.index(auth_tag)] = self_tag_map[j]
                        else:
                            buf[star_tags.index(auth_tag)] = data

                        if auth_tag != data_tag:

                            data_index = star_tags.index(data_tag)

                            if 'chain_code' in j or 'sequence_code' in j:
                                buf[data_index] = tag_map[j]
                            elif j == '_nef_covalent_links.atom_name_1':
                                buf[data_index] = k
                            elif j == '_nef_covalent_links.atom_name_2':
                                buf[data_index] = l
                            else:
                                buf[data_index] = data

                    # 'amide': (C=O)-N
                    # 'directed': N-O
                    # 'disulfide': S-S
                    # 'ester': (C=O)-O
                    # 'ether': -O-
                    # 'hydrogen': (O/N/F)-H-(O/N/F)
                    # 'metal coordination': (N/O/S)-Metal
                    # 'peptide': (C=O)-N
                    # 'thioether': -S-
                    # 'oxime': >C=N-OH
                    # 'thioester': -(C=O)-S-
                    # 'phosphoester': ?
                    # 'phosphodiester': -(PO4)-
                    # 'diselenide': Se-Se
                    buf[star_type_index] = 'covalent'
                    if k == 'SG' and l == 'SG':
                        buf[star_type_index] = 'disulfide'
                    elif k == 'SE' and l == 'SE':
                        buf[star_type_index] = 'diselenide'
                    elif (k in self.non_metal_elems and (l in self.paramag_elems or l in self.ferromag_elems)) or\
                         (l in self.non_metal_elems and (k in self.paramag_elems or k in self.ferromag_elems)):
                        buf[star_type_index] = 'metal coordination'
                    elif ((k == 'C' and l == 'N') or (l == 'C' and k == 'N')) and i[nef_tags.index(chain_tag_1)] == i[nef_tags.index(chain_tag_2)] and i[nef_tags.index(seq_tag_1)] != i[nef_tags.index(seq_tag_2)]:
                        buf[star_type_index] = 'peptide'
                    buf[star_value_order_index] = 'sing'

                    buf_row.append(buf)

            keys = set()

            for i in buf_row:

                key = ''
                for j in key_indices:
                    key += ' ' + str(i[j])
                key.rstrip()

                if key in keys:
                    continue

                keys.add(key)

                if star_id_index >= 0:
                    i[star_id_index] = index

                index += 1

                out_row.append(i)

        return out_row

    def nef2star_cs_row(self, nef_tags, star_tags, loop_data):
        """ Translate data in chemical shift loop from NEF into NMR-STAR.
            @change: rename from original translate_cs_row() to nef2star_cs_row() by Masashi Yokochi
            @param nef_tags: list of NEF tags
            @param star_tags: list of NMR-STAR tags
            @param loop_data: loop data of NEF
            @return: rows of NMR-STAR
        """

        out_row = []

        chain_index = nef_tags.index('_nef_chemical_shift.chain_code')
        seq_index = nef_tags.index('_nef_chemical_shift.sequence_code')
        comp_index = nef_tags.index('_nef_chemical_shift.residue_name')
        atom_index = nef_tags.index('_nef_chemical_shift.atom_name')

        if '_Atom_chem_shift.ID' in star_tags:
            star_id_index = star_tags.index('_Atom_chem_shift.ID')
        else:
            star_id_index = -1

        star_ambig_code_index = star_tags.index('_Atom_chem_shift.Ambiguity_code')
        star_ambig_set_id_index = star_tags.index('_Atom_chem_shift.Ambiguity_set_ID')
        star_details_index = star_tags.index('_Atom_chem_shift.Details')

        if self.insert_original_pdb_cs_items:
            star_original_chani_id = star_tags.index('_Atom_chem_shift.Original_PDB_strand_ID')
            star_original_seq_id = star_tags.index('_Atom_chem_shift.Original_PDB_residue_no')
            star_original_comp_id = star_tags.index('_Atom_chem_shift.Original_PDB_residue_name')
            star_original_atom_id = star_tags.index('_Atom_chem_shift.Original_PDB_atom_name')

        index = 1

        for nef_chain in self.authChainId:

            mapped_seq_id = [s for c, s in self.authSeqMap.keys() if c == nef_chain]
            unmapped_seq_id = sorted(set([int(i[seq_index]) for i in loop_data if i[chain_index] == nef_chain and not i[seq_index] in self.empty_value and not self.int_pattern.match(i[seq_index]) is None and not int(i[seq_index]) in mapped_seq_id]))

            if len(unmapped_seq_id) > 0:
                mapped_seq_id.extend(unmapped_seq_id)

            for _nef_seq in sorted(mapped_seq_id):

                nef_seq = str(_nef_seq)

                try:
                    star_chain, _star_seq = self.authSeqMap[(nef_chain, _nef_seq)]
                    cif_chain, _cif_seq = self.selfSeqMap[(nef_chain, _nef_seq)]
                except KeyError:
                    star_chain  = nef_chain
                    _star_seq = _nef_seq
                    cif_chain = nef_chain
                    _cif_seq = _nef_seq

                in_row = [i for i in loop_data if i[chain_index] == nef_chain and i[seq_index] == nef_seq]

                if len(in_row) == 0:
                    continue

                for i in in_row:

                    atom_list, ambiguity_code, details = self.get_star_atom(i[comp_index], i[atom_index])

                    for atom in atom_list:

                        out = [None] * len(star_tags)

                        for j in nef_tags:

                            auth_tag, data_tag, _ = self.get_star_tag(j)

                            if auth_tag is None:
                                continue

                            data = i[nef_tags.index(j)]

                            if j == '_nef_chemical_shift.chain_code':
                                out[star_tags.index(auth_tag)] = cif_chain
                            elif j == '_nef_chemical_shift.sequence_code':
                                out[star_tags.index(auth_tag)] = _cif_seq
                            else:
                                out[star_tags.index(auth_tag)] = data

                            if auth_tag != data_tag:

                                if j == '_nef_chemical_shift.atom_name':
                                    out[star_tags.index(data_tag)] = atom
                                    if self.insert_original_pdb_cs_items:
                                        out[star_original_atom_id] = data
                                elif j == '_nef_chemical_shift.chain_code':
                                    out[star_tags.index(data_tag)] = star_chain
                                    if self.insert_original_pdb_cs_items:
                                        out[star_original_chani_id] = data
                                elif j == '_nef_chemical_shift.sequence_code':
                                    out[star_tags.index(data_tag)] = _star_seq
                                    if self.insert_original_pdb_cs_items:
                                        out[star_original_seq_id] = data
                                elif j == '_nef_chemical_shift.residue_name':
                                    out[star_tags.index(data_tag)] = data
                                    if self.insert_original_pdb_cs_items:
                                        out[star_original_comp_id] = data
                                else:
                                    out[star_tags.index(data_tag)] = data

                        if star_id_index >= 0:
                            out[star_id_index] = index

                        out[star_ambig_code_index] = ambiguity_code
                        out[star_ambig_set_id_index] = None
                        out[star_details_index] = details

                        index += 1

                        out_row.append(out)

        return out_row

    def star2nef_cs_row(self, star_tags, nef_tags, loop_data):
        """ Translate data in chemical shift loop from NMR-STAR into NEF.
            @author: Masashi Yokochi
            @param star_tags: list of NMR-STAR tags
            @param nef_tags: list of NEF tags
            @param loop_data: loop data of NMR-STAR
            @return: rows of NEF
        """

        out_row = []

        chain_index = star_tags.index('_Atom_chem_shift.Entity_assembly_ID')
        seq_index = star_tags.index('_Atom_chem_shift.Comp_index_ID')
        comp_index = star_tags.index('_Atom_chem_shift.Comp_ID')
        atom_index = star_tags.index('_Atom_chem_shift.Atom_ID')
        ambig_index = star_tags.index('_Atom_chem_shift.Ambiguity_code')
        value_index = star_tags.index('_Atom_chem_shift.Val')

        for star_chain in self.authChainId:

            _star_chain = int(star_chain)
            _cif_chain = None if not star_chain in self.star2cif_chain_mapping else self.star2cif_chain_mapping[star_chain]

            mapped_seq_id = [s for c, s in self.authSeqMap.keys() if c == _star_chain]
            unmapped_seq_id = set([int(i[seq_index]) for i in loop_data if i[chain_index] == star_chain and not i[seq_index] in self.empty_value and not self.int_pattern.match(i[seq_index]) is None and not int(i[seq_index]) in mapped_seq_id])

            if len(unmapped_seq_id) > 0:
                mapped_seq_id.extend(unmapped_seq_id)

            for _star_seq in sorted(mapped_seq_id):

                star_seq = str(_star_seq)

                in_row = [i for i in loop_data if i[chain_index] == star_chain and i[seq_index] == star_seq]

                if len(in_row) == 0:
                    continue

                star_atom_list = [{'atom_id': i[atom_index], 'ambig_code': i[ambig_index], 'value': i[value_index]} for i in in_row]

                atom_list, details, atom_id_map = self.get_nef_atom(in_row[0][comp_index], star_atom_list)

                if len(atom_list) == 0:
                    continue

                seq_key = (_star_chain, _star_seq)

                try:
                    cif_chain, _cif_seq = self.authSeqMap[seq_key]
                except KeyError:
                    try:
                        cif_chain = self.selfSeqMap[(_star_chain, 1)][0]
                        _cif_seq = _star_seq
                    except KeyError:
                        if _star_chain in self.empty_value or not _star_chain in self.authChainId:
                            nef_chain = _star_chain
                        else:
                            cid = self.authChainId.index(_star_chain)
                            nef_chain = self.index_to_letter(cid)
                        cif_chain = nef_chain
                        _cif_seq = _star_seq

                if not _cif_chain is None:
                    cif_chain = _cif_chain

                if not seq_key in self.atomIdMap:
                    self.atomIdMap[seq_key] = {}

                self.atomIdMap[seq_key].update(atom_id_map)

                for atom in atom_list:

                    out = [None] * len(nef_tags)

                    star_atom = next((k for k, v in atom_id_map.items() if v == atom), atom)
                    i = next((l for l in in_row if l[atom_index] == star_atom), None)

                    if i is None:

                        if star_atom.endswith('%'):
                            star_atom = star_atom.replace('%', '')
                            i = next((l for l in in_row if l[atom_index] == star_atom), None)

                        if i is None:
                            continue

                    for j in star_tags:

                        nef_tag, _ = self.get_nef_tag(j)

                        if nef_tag is None:
                            continue

                        data_index = nef_tags.index(nef_tag)

                        if nef_tag == '_nef_chemical_shift.atom_name':
                            out[data_index] = atom
                        elif nef_tag == '_nef_chemical_shift.chain_code':
                            out[data_index] = cif_chain
                        elif nef_tag == '_nef_chemical_shift.sequence_code':
                            out[data_index] = _cif_seq
                        else:
                            out[data_index] = i[star_tags.index(j)]

                    out_row.append(out)

        return out_row

    def get_seq_ident_tags(self, in_tags, file_type):
        """ Return list of tags utilized for sequence identification.
            @change: rename from original get_residue_identifier() to get_seq_ident_tags() by Masashi Yokochi
            @change: return list of dictionary
            @param in_tags: list of tags
            @param file_type: input file type either 'nef' or 'nmr-star'
            @return: list of tags utilized for sequence identification in given tags
        """

        out_tags = []

        for j in range(1, self.lim_num_dim):

            if file_type == 'nef':
                chain_tag_suffix = '.chain_code_%s' % j
            else:
                chain_tag_suffix = '.Entity_assembly_ID_%s' % j

            try:
                chain_tag = next(i for i in in_tags if i.endswith(chain_tag_suffix))
            except StopIteration:
                break

            if file_type == 'nef':
                seq_tag_suffix = '.sequence_code_%s' % j
            else:
                seq_tag_suffix = '.Comp_index_ID_%s' % j

            try:
                seq_tag = next(i for i in in_tags if i.endswith(seq_tag_suffix))
            except StopIteration:
                break

            out_tags.append({'chain_tag': chain_tag, 'seq_tag': seq_tag})

        return out_tags

    def get_atom_keys(self, in_tags, file_type):
        """ Return list of keys utilized for atom identification.
            @change: return list of dictionary
            @param in_tags: list of tags
            @param file_type: input file type either 'nef' or 'nmr-star'
            @return: list of keys utilized for atom identification in given tags
        """

        out_tags = []

        for j in range(1, self.lim_num_dim):

            if file_type == 'nef':
                chain_tag_suffix = '.chain_code_%s' % j
            else:
                chain_tag_suffix = '.Entity_assembly_ID_%s' % j

            try:
                chain_tag = next(i for i in in_tags if i.endswith(chain_tag_suffix))
            except StopIteration:
                break

            if file_type == 'nef':
                seq_tag_suffix = '.sequence_code_%s' % j
            else:
                seq_tag_suffix = '.Comp_index_ID_%s' % j

            try:
                seq_tag = next(i for i in in_tags if i.endswith(seq_tag_suffix))
            except StopIteration:
                break

            if file_type == 'nef':
                atom_tag_suffix = '.atom_name_%s' % j
            else:
                atom_tag_suffix = '.Atom_ID_%s' % j

            try:
                atom_tag = next(i for i in in_tags if i.endswith(atom_tag_suffix))
            except StopIteration:
                break

            out_tags.append({'chain_tag': chain_tag.split('.')[1], 'seq_tag': seq_tag.split('.')[1], 'atom_tag': atom_tag.split('.')[1]})

        return out_tags

    def nef2star_dist_row(self, nef_tags, star_tags, loop_data):
        """ Translate rows of data in distance restraint loop from NEF into NMR-STAR.
            @change: rename the original translate_restraint_row() to nef2star_dist_row() by Masashi Yokochi
            @param nef_tags: list of NEF tags
            @param star_tags: list of NMR-STAR tags
            @param loop_data: loop data of NEF
            @return: rows of NMR-STAR
        """

        out_row = []

        nef_comp_index_1 = nef_tags.index('_nef_distance_restraint.residue_name_1')
        nef_atom_index_1 = nef_tags.index('_nef_distance_restraint.atom_name_1')
        nef_comp_index_2 = nef_tags.index('_nef_distance_restraint.residue_name_2')
        nef_atom_index_2 = nef_tags.index('_nef_distance_restraint.atom_name_2')

        seq_ident_tags = self.get_seq_ident_tags(nef_tags, 'nef')

        for tag in seq_ident_tags:
            chain_tag = tag['chain_tag']
            seq_tag = tag['seq_tag']

            if chain_tag.endswith('_1'):
                chain_tag_1 = chain_tag
                seq_tag_1 = seq_tag
            else:
                chain_tag_2 = chain_tag
                seq_tag_2 = seq_tag

        try:
            index_index = star_tags.index('_Gen_dist_constraint.Index_ID')
        except ValueError:
            index_index = -1

        key_indices = [star_tags.index(j) for j in ['_Gen_dist_constraint.Entity_assembly_ID_1', '_Gen_dist_constraint.Comp_index_ID_1', '_Gen_dist_constraint.Atom_ID_1',
                                                    '_Gen_dist_constraint.Entity_assembly_ID_2', '_Gen_dist_constraint.Comp_index_ID_2', '_Gen_dist_constraint.Atom_ID_2']]

        member_code_index = star_tags.index('_Gen_dist_constraint.Member_logic_code')

        id_index = nef_tags.index('_nef_distance_restraint.restraint_id')

        id_list = sorted(set([int(i[id_index]) for i in loop_data]))

        index = 1

        for id in id_list:

            in_row = [i for i in loop_data if i[id_index] == str(id)]

            buf_row = []

            for i in in_row:

                tag_map = {}
                self_tag_map = {}

                for tag in seq_ident_tags:
                    chain_tag = tag['chain_tag']
                    seq_tag = tag['seq_tag']

                    nef_chain = i[nef_tags.index(chain_tag)]
                    _nef_seq = i[nef_tags.index(seq_tag)]
                    if type(_nef_seq) == str and not _nef_seq in self.empty_value:
                        _nef_seq = int(_nef_seq)

                    seq_key = (nef_chain, _nef_seq)

                    try:
                        tag_map[chain_tag], tag_map[seq_tag] = self.authSeqMap[seq_key]
                        self_tag_map[chain_tag], self_tag_map[seq_tag] = self.selfSeqMap[seq_key]
                    except KeyError:
                        tag_map[chain_tag] = nef_chain
                        tag_map[seq_tag] = _nef_seq
                        self_tag_map[chain_tag] = nef_chain
                        self_tag_map[seq_tag] = _nef_seq

                intra_residue = i[nef_tags.index(chain_tag_1)] == i[nef_tags.index(chain_tag_2)] and i[nef_tags.index(seq_tag_1)] == i[nef_tags.index(seq_tag_2)]

                atom_list_1, ambiguity_code_1, details_1 = self.get_star_atom(i[nef_comp_index_1], i[nef_atom_index_1])
                atom_list_2, ambiguity_code_2, details_2 = self.get_star_atom(i[nef_comp_index_2], i[nef_atom_index_2])

                or_code = len(atom_list_1) * len(atom_list_2) > 1

                for k in atom_list_1:

                    for l in atom_list_2:

                        if intra_residue and l == k:
                            continue

                        buf = [None] * len(star_tags)

                        for j in nef_tags:

                            auth_tag, data_tag, _ = self.get_star_tag(j)

                            if auth_tag is None:
                                continue

                            data = i[nef_tags.index(j)]

                            if 'chain_code' in j or 'sequence_code' in j:
                                buf[star_tags.index(auth_tag)] = self_tag_map[j]
                            else:
                                buf[star_tags.index(auth_tag)] = data

                            if auth_tag != data_tag:

                                data_index = star_tags.index(data_tag)

                                if 'chain_code' in j or 'sequence_code' in j:
                                    buf[data_index] = tag_map[j]
                                elif j == '_nef_distance_restraint.atom_name_1':
                                    buf[data_index] = k
                                elif j == '_nef_distance_restraint.atom_name_2':
                                    buf[data_index] = l
                                else:
                                    buf[data_index] = data
                            """
                            if details_1 is None and details_2 is None:
                                pass

                            else:

                                details_index = star_tags.index('_Gen_dist_constraint.Details')

                                buf[details_index] = ' '.join(filter(None, [details_1, detail_2]))
                            """
                        if or_code:
                            buf[member_code_index] = 'OR'

                        buf_row.append(buf)

            keys = set()

            for i in buf_row:

                key = ''
                for j in key_indices:
                    key += ' ' + str(i[j])
                key.rstrip()

                if key in keys:
                    continue

                keys.add(key)

                if index_index >= 0:
                    i[index_index] = index

                index += 1

                out_row.append(i)

        return out_row

    def star2nef_dist_row(self, star_tags, nef_tags, loop_data):
        """ Translate rows of data in distance restraint loop from NMR-STAR into NEF.
            @author: Masashi Yokochi
            @param star_tags: list of NMR-STAR tags
            @param nef_tags: list of NEF tags
            @param loop_data: loop data of NMR-STAR
            @return: rows of NEF
        """

        out_row = []

        comp_1_index = star_tags.index('_Gen_dist_constraint.Comp_ID_1')
        comp_2_index = star_tags.index('_Gen_dist_constraint.Comp_ID_2')

        seq_ident_tags = self.get_seq_ident_tags(star_tags, 'nmr-star')

        for tag in seq_ident_tags:
            chain_tag = tag['chain_tag']
            seq_tag = tag['seq_tag']

            if chain_tag.endswith('_1'):
                chain_tag_1 = chain_tag
                #seq_tag_1 = seq_tag
            else:
                pass
                #chain_tag_2 = chain_tag
                #seq_tag_2 = seq_tag

        try:
            index_index = nef_tags.index('_nef_distance_restraint.index')
        except ValueError:
            index_index = -1

        key_indices = [nef_tags.index(j) for j in ['_nef_distance_restraint.chain_code_1', '_nef_distance_restraint.sequence_code_1', '_nef_distance_restraint.atom_name_1',
                                                   '_nef_distance_restraint.chain_code_2', '_nef_distance_restraint.sequence_code_2', '_nef_distance_restraint.atom_name_2']]

        id_index = star_tags.index('_Gen_dist_constraint.ID')

        id_list = sorted(set([int(i[id_index]) for i in loop_data]))

        index = 1

        for id in id_list:

            in_row = [i for i in loop_data if i[id_index] == str(id)]

            buf_row = []

            for i in in_row:

                tag_map = {}

                for tag in seq_ident_tags:
                    chain_tag = tag['chain_tag']
                    seq_tag = tag['seq_tag']

                    star_chain = i[star_tags.index(chain_tag)]
                    _star_chain = star_chain
                    if type(star_chain) == str and not star_chain in self.empty_value:
                        _star_chain = int(star_chain)

                    _star_seq = i[star_tags.index(seq_tag)]
                    if type(_star_seq) == str and not _star_seq in self.empty_value:
                        _star_seq = int(_star_seq)

                    seq_key = (_star_chain, _star_seq)

                    try:
                        tag_map[chain_tag], tag_map[seq_tag] = self.authSeqMap[seq_key]
                    except KeyError:
                        try:
                            nef_chain = self.selfSeqMap[(_star_chain, 1)][0]
                        except KeyError:
                            if _star_chain in self.empty_value or not _star_chain in self.authChainId:
                                nef_chain = _star_chain
                            else:
                                cid = self.authChainId.index(_star_chain)
                                nef_chain = self.index_to_letter(cid)
                        tag_map[chain_tag] = nef_chain
                        tag_map[seq_tag] = _star_seq

                    if star_chain in self.star2cif_chain_mapping:
                        tag_map[chain_tag] = self.star2cif_chain_mapping[star_chain]

                    if chain_tag == chain_tag_1:
                        seq_key_1 = seq_key
                    else:
                        seq_key_2 = seq_key

                buf = [None] * len(nef_tags)

                for j in star_tags:

                    nef_tag, _ = self.get_nef_tag(j)

                    if nef_tag is None:
                        continue

                    data = i[star_tags.index(j)]

                    data_index = nef_tags.index(nef_tag)

                    if 'chain_code' in nef_tag or 'sequence_code' in nef_tag:
                        buf[data_index] = tag_map[j]
                    elif nef_tag == '_nef_distance_restraint.atom_name_1':
                        try:
                            if not self.atomIdMap is None:
                                buf[data_index] = self.atomIdMap[seq_key_1][data]
                        except KeyError:
                            atom_list = self.get_nef_atom(i[comp_1_index], [{'atom_id': data, 'ambig_code': None, 'value': None}])[0]
                            if len(atom_list) > 0:
                                buf[data_index] = atom_list[0]
                            else:
                                buf[data_index] = data
                    elif nef_tag == '_nef_distance_restraint.atom_name_2':
                        try:
                            if not self.atomIdMap is None:
                                buf[data_index] = self.atomIdMap[seq_key_2][data]
                        except KeyError:
                            atom_list = self.get_nef_atom(i[comp_2_index], [{'atom_id': data, 'ambig_code': None, 'value': None}])[0]
                            if len(atom_list) > 0:
                                buf[data_index] = atom_list[0]
                            else:
                                buf[data_index] = data
                    else:
                        buf[data_index] = data

                buf_row.append(buf)

            keys = set()

            for i in buf_row:

                key = ''
                for j in key_indices:
                    key += ' ' + str(i[j])
                key.rstrip()

                if key in keys:
                    continue

                keys.add(key)

                if index_index >= 0:
                    i[index_index] = index

                index += 1

                out_row.append(i)

        return out_row

    def nef2star_dihed_row(self, nef_tags, star_tags, loop_data):
        """ Translate rows of data in dihedral restraint loop from NEF into NMR-STAR.
            @author: Masashi Yokochi
            @param nef_tags: list of NEF tags
            @param star_tags: list of NMR-STAR tags
            @param loop_data: loop data of NEF
            @return: rows of NMR-STAR
        """

        out_row = []

        nef_comp_index_1 = nef_tags.index('_nef_dihedral_restraint.residue_name_1')
        nef_atom_index_1 = nef_tags.index('_nef_dihedral_restraint.atom_name_1')
        nef_comp_index_2 = nef_tags.index('_nef_dihedral_restraint.residue_name_2')
        nef_atom_index_2 = nef_tags.index('_nef_dihedral_restraint.atom_name_2')
        nef_comp_index_3 = nef_tags.index('_nef_dihedral_restraint.residue_name_3')
        nef_atom_index_3 = nef_tags.index('_nef_dihedral_restraint.atom_name_3')
        nef_comp_index_4 = nef_tags.index('_nef_dihedral_restraint.residue_name_4')
        nef_atom_index_4 = nef_tags.index('_nef_dihedral_restraint.atom_name_4')

        seq_ident_tags = self.get_seq_ident_tags(nef_tags, 'nef')

        for tag in seq_ident_tags:
            chain_tag = tag['chain_tag']
            seq_tag = tag['seq_tag']

            if chain_tag.endswith('_1'):
                chain_tag_1 = chain_tag
                seq_tag_1 = seq_tag
            elif chain_tag.endswith('_2'):
                chain_tag_2 = chain_tag
                seq_tag_2 = seq_tag
            elif chain_tag.endswith('_3'):
                chain_tag_3 = chain_tag
                seq_tag_3 = seq_tag
            else:
                chain_tag_4 = chain_tag
                seq_tag_4 = seq_tag

        try:
            index_index = star_tags.index('_Torsion_angle_constraint.Index_ID')
        except ValueError:
            index_index = -1

        key_indices = [star_tags.index(j) for j in ['_Torsion_angle_constraint.Entity_assembly_ID_1', '_Torsion_angle_constraint.Comp_index_ID_1', '_Torsion_angle_constraint.Atom_ID_1',
                                                    '_Torsion_angle_constraint.Entity_assembly_ID_2', '_Torsion_angle_constraint.Comp_index_ID_2', '_Torsion_angle_constraint.Atom_ID_2',
                                                    '_Torsion_angle_constraint.Entity_assembly_ID_3', '_Torsion_angle_constraint.Comp_index_ID_3', '_Torsion_angle_constraint.Atom_ID_3',
                                                    '_Torsion_angle_constraint.Entity_assembly_ID_4', '_Torsion_angle_constraint.Comp_index_ID_4', '_Torsion_angle_constraint.Atom_ID_4']]

        id_index = nef_tags.index('_nef_dihedral_restraint.restraint_id')

        id_list = sorted(set([int(i[id_index]) for i in loop_data]))

        index = 1

        for id in id_list:

            in_row = [i for i in loop_data if i[id_index] == str(id)]

            buf_row = []

            for i in in_row:

                tag_map = {}
                self_tag_map = {}

                for tag in seq_ident_tags:
                    chain_tag = tag['chain_tag']
                    seq_tag = tag['seq_tag']

                    nef_chain = i[nef_tags.index(chain_tag)]
                    _nef_seq = i[nef_tags.index(seq_tag)]
                    if type(_nef_seq) == str and not _nef_seq in self.empty_value:
                        _nef_seq = int(_nef_seq)

                    seq_key = (nef_chain, _nef_seq)

                    try:
                        tag_map[chain_tag], tag_map[seq_tag] = self.authSeqMap[seq_key]
                        self_tag_map[chain_tag], self_tag_map[seq_tag] = self.selfSeqMap[seq_key]
                    except KeyError:
                        tag_map[chain_tag] = nef_chain
                        tag_map[seq_tag] = _nef_seq
                        self_tag_map[chain_tag] = nef_chain
                        self_tag_map[seq_tag] = _nef_seq

                intra_residue_12 = i[nef_tags.index(chain_tag_1)] == i[nef_tags.index(chain_tag_2)] and i[nef_tags.index(seq_tag_1)] == i[nef_tags.index(seq_tag_2)]
                intra_residue_13 = i[nef_tags.index(chain_tag_1)] == i[nef_tags.index(chain_tag_3)] and i[nef_tags.index(seq_tag_1)] == i[nef_tags.index(seq_tag_3)]
                intra_residue_14 = i[nef_tags.index(chain_tag_1)] == i[nef_tags.index(chain_tag_4)] and i[nef_tags.index(seq_tag_1)] == i[nef_tags.index(seq_tag_4)]
                intra_residue_23 = i[nef_tags.index(chain_tag_2)] == i[nef_tags.index(chain_tag_3)] and i[nef_tags.index(seq_tag_2)] == i[nef_tags.index(seq_tag_3)]
                intra_residue_24 = i[nef_tags.index(chain_tag_2)] == i[nef_tags.index(chain_tag_4)] and i[nef_tags.index(seq_tag_2)] == i[nef_tags.index(seq_tag_4)]
                intra_residue_34 = i[nef_tags.index(chain_tag_3)] == i[nef_tags.index(chain_tag_4)] and i[nef_tags.index(seq_tag_3)] == i[nef_tags.index(seq_tag_4)]

                atom_list_1, ambiguity_code_1, details_1 = self.get_star_atom(i[nef_comp_index_1], i[nef_atom_index_1])
                atom_list_2, ambiguity_code_2, details_2 = self.get_star_atom(i[nef_comp_index_2], i[nef_atom_index_2])
                atom_list_3, ambiguity_code_3, details_3 = self.get_star_atom(i[nef_comp_index_3], i[nef_atom_index_3])
                atom_list_4, ambiguity_code_4, details_4 = self.get_star_atom(i[nef_comp_index_4], i[nef_atom_index_4])

                for k in atom_list_1:

                    for l in atom_list_2:

                        if intra_residue_12 and l == k:
                            continue

                        for m in atom_list_3:

                            if (intra_residue_13 and m == k) or (intra_residue_23 and m == l):
                                continue

                            for n in atom_list_4:

                                if (intra_residue_14 and n == k) or (intra_residue_24 and n == l) or (intra_residue_34 and n == m):
                                    continue

                                buf = [None] * len(star_tags)

                                for j in nef_tags:

                                    auth_tag, data_tag, _ = self.get_star_tag(j)

                                    if auth_tag is None:
                                        continue

                                    data = i[nef_tags.index(j)]

                                    if 'chain_code' in j or 'sequence_code' in j:
                                        buf[star_tags.index(auth_tag)] = self_tag_map[j]
                                    else:
                                        buf[star_tags.index(auth_tag)] = data

                                    if auth_tag != data_tag:

                                        data_index = star_tags.index(data_tag)

                                        if 'chain_code' in j or 'sequence_code' in j:
                                            buf[data_index] = tag_map[j]
                                        elif j == '_nef_dihedral_restraint.atom_name_1':
                                            buf[data_index] = k
                                        elif j == '_nef_dihedral_restraint.atom_name_2':
                                            buf[data_index] = l
                                        elif j == '_nef_dihedral_restraint.atom_name_3':
                                            buf[data_index] = m
                                        elif j == '_nef_dihedral_restraint.atom_name_4':
                                            buf[data_index] = n
                                        else:
                                            buf[data_index] = data
                                    """
                                    if details_1 is None and details_2 is None and details_3 is None and details_4 is None:
                                        pass

                                    else:

                                        details_index = star_tags.index('_Torsion_angle_constraint.Details')

                                        buf[details_index] = ' '.join(filter(None, [details_1, detail_2, details_3, details_4]))
                                    """
                                buf_row.append(buf)

            keys = set()

            for i in buf_row:

                key = ''
                for j in key_indices:
                    key += ' ' + str(i[j])
                key.rstrip()

                if key in keys:
                    continue

                keys.add(key)

                if index_index >= 0:
                    i[index_index] = index

                index += 1

                out_row.append(i)

        return out_row

    def nef2star_rdc_row(self, nef_tags, star_tags, loop_data):
        """ Translate rows of data in rdc restraint loop from NEF into NMR-STAR.
            @author: Masashi Yokochi
            @param nef_tags: list of NEF tags
            @param star_tags: list of NMR-STAR tags
            @param loop_data: loop data of NEF
            @return: rows of NMR-STAR
        """

        out_row = []

        nef_comp_index_1 = nef_tags.index('_nef_rdc_restraint.residue_name_1')
        nef_atom_index_1 = nef_tags.index('_nef_rdc_restraint.atom_name_1')
        nef_comp_index_2 = nef_tags.index('_nef_rdc_restraint.residue_name_2')
        nef_atom_index_2 = nef_tags.index('_nef_rdc_restraint.atom_name_2')

        seq_ident_tags = self.get_seq_ident_tags(nef_tags, 'nef')

        for tag in seq_ident_tags:
            chain_tag = tag['chain_tag']
            seq_tag = tag['seq_tag']

            if chain_tag.endswith('_1'):
                chain_tag_1 = chain_tag
                seq_tag_1 = seq_tag
            else:
                chain_tag_2 = chain_tag
                seq_tag_2 = seq_tag

        try:
            index_index = star_tags.index('_RDC_constraint.Index_ID')
        except ValueError:
            index_index = -1

        key_indices = [star_tags.index(j) for j in ['_RDC_constraint.Entity_assembly_ID_1', '_RDC_constraint.Comp_index_ID_1', '_RDC_constraint.Atom_ID_1',
                                                    '_RDC_constraint.Entity_assembly_ID_2', '_RDC_constraint.Comp_index_ID_2', '_RDC_constraint.Atom_ID_2']]

        id_index = nef_tags.index('_nef_rdc_restraint.restraint_id')

        id_list = sorted(set([int(i[id_index]) for i in loop_data]))

        index = 1

        for id in id_list:

            in_row = [i for i in loop_data if i[id_index] == str(id)]

            buf_row = []

            for i in in_row:

                tag_map = {}
                self_tag_map = {}

                for tag in seq_ident_tags:
                    chain_tag = tag['chain_tag']
                    seq_tag = tag['seq_tag']

                    nef_chain = i[nef_tags.index(chain_tag)]
                    _nef_seq = i[nef_tags.index(seq_tag)]
                    if type(_nef_seq) == str and not _nef_seq in self.empty_value:
                        _nef_seq = int(_nef_seq)

                    seq_key = (nef_chain, _nef_seq)

                    try:
                        tag_map[chain_tag], tag_map[seq_tag] = self.authSeqMap[seq_key]
                        self_tag_map[chain_tag], self_tag_map[seq_tag] = self.selfSeqMap[seq_key]
                    except KeyError:
                        tag_map[chain_tag] = nef_chain
                        tag_map[seq_tag] = _nef_seq
                        self_tag_map[chain_tag] = nef_chain
                        self_tag_map[seq_tag] = _nef_seq

                intra_residue = i[nef_tags.index(chain_tag_1)] == i[nef_tags.index(chain_tag_2)] and i[nef_tags.index(seq_tag_1)] == i[nef_tags.index(seq_tag_2)]

                atom_list_1, ambiguity_code_1, details_1 = self.get_star_atom(i[nef_comp_index_1], i[nef_atom_index_1])
                atom_list_2, ambiguity_code_2, details_2 = self.get_star_atom(i[nef_comp_index_2], i[nef_atom_index_2])

                for k in atom_list_1:

                    for l in atom_list_2:

                        if intra_residue and l == k:
                            continue

                        buf = [None] * len(star_tags)

                        for j in nef_tags:

                            auth_tag, data_tag, _ = self.get_star_tag(j)

                            if auth_tag is None:
                                continue

                            data = i[nef_tags.index(j)]

                            if 'chain_code' in j or 'sequence_code' in j:
                                buf[star_tags.index(auth_tag)] = self_tag_map[j]
                            else:
                                buf[star_tags.index(auth_tag)] = data

                            if auth_tag != data_tag:

                                data_index = star_tags.index(data_tag)

                                if 'chain_code' in j or 'sequence_code' in j:
                                    buf[data_index] = tag_map[j]
                                elif j == '_nef_rdc_restraint.atom_name_1':
                                    buf[data_index] = k
                                elif j == '_nef_rdc_restraint.atom_name_2':
                                    buf[data_index] = l
                                else:
                                    buf[data_index] = data
                            """
                            if details_1 is None and details_2 is None:
                                pass

                            else:

                                details_index = star_tags.index('_RDC_constraint.Details')

                                buf[details_index] = ' '.join(filter(None, [details_1, detail_2]))
                            """
                        buf_row.append(buf)

            keys = set()

            for i in buf_row:

                key = ''
                for j in key_indices:
                    key += ' ' + str(i[j])
                key.rstrip()

                if key in keys:
                    continue

                keys.add(key)

                if index_index >= 0:
                    i[index_index] = index

                index += 1

                out_row.append(i)

        return out_row

    def nef2star_peak_row(self, nef_tags, star_tags, loop_data):
        """ Translate rows of data in spectral peak loop from NEF into NMR-STAR.
            @param nef_tags: list of NEF tags
            @param star_tags: list of NMR-STAR tags
            @param loop_data: loop data of NEF
            @return: rows of NMR-STAR
        """

        out_row = []

        nef_comp_indices = []
        nef_atom_indices = []

        seq_ident_tags = self.get_seq_ident_tags(nef_tags, 'nef')

        num_dim = len(seq_ident_tags)

        for d in range(1, num_dim + 1):
            nef_comp_indices.append(nef_tags.index('_nef_peak.residue_name_%s' % d))
            nef_atom_indices.append(nef_tags.index('_nef_peak.atom_name_%s' % d))

        try:
            index_index = star_tags.index('_Peak_row_format.Index_ID')
        except ValueError:
            index_index = -1

        key_indices = [star_tags.index(j) for j in [k for k in star_tags if k.startswith('_Peak_row_format.Entity_assembly_ID') or\
                                                                            k.startswith('_Peak_row_format.Comp_index_ID') or\
                                                                            k.startswith('_Peak_row_format.Atom_ID')]]

        id_index = nef_tags.index('_nef_peak.peak_id')

        id_list = sorted(set([int(i[id_index]) for i in loop_data]))

        index = 1

        for id in id_list:

            in_row = [i for i in loop_data if i[id_index] == str(id)]

            buf_row = []

            for i in in_row:

                tag_map = {}
                self_tag_map = {}

                for tag in seq_ident_tags:
                    chain_tag = tag['chain_tag']
                    seq_tag = tag['seq_tag']

                    nef_chain = i[nef_tags.index(chain_tag)]
                    _nef_seq = i[nef_tags.index(seq_tag)]
                    if type(_nef_seq) == str and not _nef_seq in self.empty_value:
                        _nef_seq = int(_nef_seq)

                    seq_key = (nef_chain, _nef_seq)

                    try:
                        tag_map[chain_tag], tag_map[seq_tag] = self.authSeqMap[seq_key]
                        self_tag_map[chain_tag], self_tag_map[seq_tag] = self.selfSeqMap[seq_key]
                    except KeyError:
                        tag_map[chain_tag] = nef_chain
                        tag_map[seq_tag] = _nef_seq
                        self_tag_map[chain_tag] = nef_chain
                        self_tag_map[seq_tag] = _nef_seq

                details = ''
                a = []

                for nef_comp_index, nef_atom_index  in zip(nef_comp_indices, nef_atom_indices):
                    atom_list = self.get_star_atom(i[nef_comp_index], i[nef_atom_index])[0]
                    len_atom_list = len(atom_list)
                    if len_atom_list == 0:
                        atom_list.append('.')
                    elif len_atom_list > 1:
                        details += '%s -> %s, ' % (i[nef_atom_index], atom_list)
                    a.append(atom_list)

                if num_dim == 1:

                    for comb in itertools.product(a[0]):
                        buf_row.append(self.__nef2star_peak_row(nef_tags, star_tags, tag_map, self_tag_map, i, details, comb))

                elif num_dim == 2:

                    for comb in itertools.product(a[0], a[1]):
                        buf_row.append(self.__nef2star_peak_row(nef_tags, star_tags, tag_map, self_tag_map, i, details, comb))

                elif num_dim == 3:

                    for comb in itertools.product(a[0], a[1], a[2]):
                        buf_row.append(self.__nef2star_peak_row(nef_tags, star_tags, tag_map, self_tag_map, i, details, comb))

                elif num_dim == 4:

                    for comb in itertools.product(a[0], a[1], a[2], a[3]):
                        buf_row.append(self.__nef2star_peak_row(nef_tags, star_tags, tag_map, self_tag_map, i, details, comb))

                elif num_dim == 5:

                    for comb in itertools.product(a[0], a[1], a[2], a[3], a[4]):
                        buf_row.append(self.__nef2star_peak_row(nef_tags, star_tags, tag_map, self_tag_map, i, details, comb))

                elif num_dim == 6:

                    for comb in itertools.product(a[0], a[1], a[2], a[3], a[4], a[5]):
                        buf_row.append(self.__nef2star_peak_row(nef_tags, star_tags, tag_map, self_tag_map, i, details, comb))

                elif num_dim == 7:

                    for comb in itertools.product(a[0], a[1], a[2], a[3], a[4], a[5], a[6]):
                        buf_row.append(self.__nef2star_peak_row(nef_tags, star_tags, tag_map, self_tag_map, i, details, comb))

                elif num_dim == 8:

                    for comb in itertools.product(a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7]):
                        buf_row.append(self.__nef2star_peak_row(nef_tags, star_tags, tag_map, self_tag_map, i, details, comb))

                elif num_dim == 9:

                    for comb in itertools.product(a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8]):
                        buf_row.append(self.__nef2star_peak_row(nef_tags, star_tags, tag_map, self_tag_map, i, details, comb))

                elif num_dim == 10:

                    for comb in itertools.product(a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9]):
                        buf_row.append(self.__nef2star_peak_row(nef_tags, star_tags, tag_map, self_tag_map, i, details, comb))

                elif num_dim == 11:

                    for comb in itertools.product(a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10]):
                        buf_row.append(self.__nef2star_peak_row(nef_tags, star_tags, tag_map, self_tag_map, i, details, comb))

                elif num_dim == 12:

                    for comb in itertools.product(a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10], a[11]):
                        buf_row.append(self.__nef2star_peak_row(nef_tags, star_tags, tag_map, self_tag_map, i, details, comb))

                elif num_dim == 13:

                    for comb in itertools.product(a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10], a[11], a[12]):
                        buf_row.append(self.__nef2star_peak_row(nef_tags, star_tags, tag_map, self_tag_map, i, details, comb))

                elif num_dim == 14:

                    for comb in itertools.product(a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10], a[11], a[12], a[13]):
                        buf_row.append(self.__nef2star_peak_row(nef_tags, star_tags, tag_map, self_tag_map, i, details, comb))

                elif num_dim == 15:

                    for comb in itertools.product(a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10], a[11], a[12], a[13], a[14]):
                        buf_row.append(self.__nef2star_peak_row(nef_tags, star_tags, tag_map, self_tag_map, i, details, comb))

            keys = set()

            for i in buf_row:

                key = ''
                for j in key_indices:
                    key += ' ' + str(i[j])
                key.rstrip()

                if key in keys:
                    continue

                keys.add(key)

                if index_index >= 0:
                    i[index_index] = index

                index += 1

                out_row.append(i)

        return out_row

    def __nef2star_peak_row(self, nef_tags, star_tags, tag_map, self_tag_map, row, details, comb):
        """ Translate rows of data in spectral peak loop from NEF into NMR-STAR.
        """

        buf = [None] * len(star_tags)

        for j in nef_tags:

            auth_tag, data_tag, _ = self.get_star_tag(j)

            if auth_tag is None:
                continue

            data = row[nef_tags.index(j)]

            if 'chain_code' in j or 'sequence_code' in j:
                buf[star_tags.index(auth_tag)] = self_tag_map[j]
            else:
                buf[star_tags.index(auth_tag)] = data

            if auth_tag != data_tag:

                data_index = star_tags.index(data_tag)

                if 'chain_code' in j or 'sequence_code' in j:
                    buf[data_index] = tag_map[j]
                elif j.startswith('_nef_peak.atom_name'):
                    buf[data_index] = comb[int(j[20:]) - 1]
                else:
                    buf[data_index] = data

            if len(details) > 0:
                details_index = star_tags.index('_Peak_row_format.Details')

                buf[details_index] = details[:-2]

        return buf

    def star2nef_peak_row(self, star_tags, nef_tags, loop_data):
        """ Translate rows of data in spectral peak loop from NMR-STAR into NEF.
            @author: Masashi Yokochi
            @param star_tags: list of NMR-STAR tags
            @param nef_tags: list of NEF tags
            @param loop_data: loop data of NMR-STAR
            @return: rows of NEF
        """

        out_row = []

        comp_indices = []

        seq_ident_tags = self.get_seq_ident_tags(star_tags, 'nmr-star')

        num_dim = len(seq_ident_tags)

        for d in range(1, num_dim + 1):
            comp_indices.append(star_tags.index('_Peak_row_format.Comp_ID_%s' % d))

        try:
            index_index = nef_tags.index('_nef_peak.index')
        except ValueError:
            index_index = -1

        key_indices = [nef_tags.index(j) for j in [k for k in nef_tags if k.startswith('_nef_peak.chain_code') or\
                                                                          k.startswith('_nef_peak.sequence_code') or\
                                                                          k.startswith('_nef_peak.atom_name')]]

        id_index = star_tags.index('_Peak_row_format.ID')

        id_list = sorted(set([int(i[id_index]) for i in loop_data]))

        index = 1

        for id in id_list:

            in_row = [i for i in loop_data if i[id_index] == str(id)]

            buf_row = []

            for i in in_row:

                tag_map = {}

                s = []

                for tag in seq_ident_tags:
                    chain_tag = tag['chain_tag']
                    seq_tag = tag['seq_tag']

                    star_chain = i[star_tags.index(chain_tag)]
                    _star_chain = star_chain
                    if type(star_chain) == str and not star_chain in self.empty_value:
                        _star_chain = int(star_chain)

                    _star_seq = i[star_tags.index(seq_tag)]
                    if type(_star_seq) == str and not _star_seq in self.empty_value:
                        _star_seq = int(_star_seq)

                    seq_key = (_star_chain, _star_seq)

                    try:
                        tag_map[chain_tag], tag_map[seq_tag] = self.authSeqMap[seq_key]
                    except KeyError:
                        try:
                            nef_chain = self.selfSeqMap[(_star_chain, 1)][0]
                        except KeyError:
                            if _star_chain in self.empty_value or not _star_chain in self.authChainId:
                                nef_chain = _star_chain
                            else:
                                cid = self.authChainId.index(_star_chain)
                                nef_chain = self.index_to_letter(cid)
                        tag_map[chain_tag] = nef_chain
                        tag_map[seq_tag] = _star_seq

                    if star_chain in self.star2cif_chain_mapping:
                        tag_map[chain_tag] = self.star2cif_chain_mapping[star_chain]

                    s.append(seq_key)

                buf = [None] * len(nef_tags)

                for j in star_tags:

                    nef_tag, _ = self.get_nef_tag(j)

                    if nef_tag is None:
                        continue

                    data = i[star_tags.index(j)]

                    data_index = nef_tags.index(nef_tag)

                    if 'chain_code' in nef_tag or 'sequence_code' in nef_tag:
                        buf[data_index] = tag_map[j]
                    elif nef_tag.startswith('_nef_peak.atom_name'):
                        try:
                            if not self.atomIdMap is None:
                                buf[data_index] = self.atomIdMap[s[int(nef_tag[20:]) - 1]][data]
                        except KeyError:
                            atom_list = self.get_nef_atom(i[comp_indices[int(nef_tag[20:]) - 1]], [{'atom_id': data, 'ambig_code': None, 'value': None}])[0]
                            if len(atom_list) > 0:
                                buf[data_index] = atom_list[0]
                            else:
                                buf[data_index] = data
                    else:
                        buf[data_index] = data

                buf_row.append(buf)

            keys = set()

            for i in buf_row:

                key = ''
                for j in key_indices:
                    key += ' ' + str(i[j])
                key.rstrip()

                if key in keys:
                    continue

                keys.add(key)

                if index_index >= 0:
                    i[index_index] = index

                index += 1

                out_row.append(i)

        return out_row

    def nef2star_row(self, nef_tags, star_tags, in_row):
        """ Translate rows of data in a loop from NEF into NMR-STAR.
            @change: rename from original translate_row() to nef2star_row() by Masashi Yokochi
            @param nef_tags: list of NEF tags
            @param star_tags: list of NMR-STAR tags
            @param in_row: rows of NEF
            @return: rows of NMR-STAR
        """

        out_row = []

        tag_map = {}
        self_tag_map = {}

        for tag in self.get_seq_ident_tags(nef_tags, 'nef'):
            chain_tag = tag['chain_tag']
            seq_tag = tag['seq_tag']

            nef_chain = in_row[nef_tags.index(chain_tag)]
            _nef_seq = in_row[nef_tags.index(seq_tag)]
            if type(_nef_seq) == str and not _nef_seq in self.empty_value:
                _nef_seq = int(_nef_seq)

            seq_key = (nef_chain, _nef_seq)

            try:
                tag_map[chain_tag], tag_map[seq_tag] = self.authSeqMap[seq_key]
                self_tag_map[chain_tag], self_tag_map[seq_tag] = self.selfSeqMap[seq_key]
            except KeyError:
                tag_map[chain_tag] = nef_chain
                tag_map[seq_tag] = _nef_seq
                self_tag_map[chain_tag] = nef_chain
                self_tag_map[seq_tag] = _nef_seq

        if len(nef_tags) != len(star_tags):

            out = [None] * len(star_tags)

            for j in nef_tags:

                auth_tag, data_tag, _ = self.get_star_tag(j)

                if auth_tag is None:
                    continue

                data = in_row[nef_tags.index(j)]

                if 'chain_code' in j or 'sequence_code' in j:
                    out[star_tags.index(auth_tag)] = self_tag_map[j]
                else:
                    out[star_tags.index(auth_tag)] = data

                if auth_tag != data_tag:

                    data_index = star_tags.index(data_tag)

                    if 'chain_code' in j or 'sequence_code' in j:
                        out[data_index] = tag_map[j]
                    elif data in self.nef_boolean:
                        out[data_index] = 'yes' if data in self.true_value else 'no'
                    else:
                        out[data_index] = data

            out_row.append(out)

        else:
            out_row.append(in_row)

        return out_row

    def star2nef_row(self, star_tags, nef_tags, in_row):
        """ Translate rows of data in a loop from NMR-STAR into NEF.
            @author: Masashi Yokochi
            @param star_tags: list of NMR-STAR tags
            @param nef_tags: list of NEF tags
            @param in_row: rows of NMR-STAR
            @return: rows of NEF
        """

        out_row = []

        tag_map = {}

        for tag in self.get_seq_ident_tags(star_tags, 'nmr-star'):
            chain_tag = tag['chain_tag']
            seq_tag = tag['seq_tag']

            star_chain = in_row[star_tags.index(chain_tag)]
            _star_chain = star_chain
            if type(star_chain) == str and not star_chain in self.empty_value:
                _star_chain = int(star_chain)

            _star_seq = in_row[star_tags.index(seq_tag)]
            if type(_star_seq) == str and not _star_seq in self.empty_value:
                _star_seq = int(_star_seq)

            try:
                tag_map[chain_tag], tag_map[seq_tag] = self.authSeqMap[(_star_chain, _star_seq)]
            except KeyError:
                try:
                    nef_chain = self.selfSeqMap[(_star_chain, 1)][0]
                except KeyError:
                    if _star_chain in self.empty_value or not _star_chain in self.authChainId:
                        nef_chain = _star_chain
                    else:
                        cid = self.authChainId.index(_star_chain)
                        nef_chain = self.index_to_letter(cid)
                tag_map[chain_tag] = nef_chain
                tag_map[seq_tag] = _star_seq

            if star_chain in self.star2cif_chain_mapping:
                tag_map[chain_tag] = self.star2cif_chain_mapping[star_chain]

        out = [None] * len(nef_tags)

        for j in star_tags:

            nef_tag, _ = self.get_nef_tag(j)

            if nef_tag is None:
                continue

            data = in_row[star_tags.index(j)]

            data_index = nef_tags.index(nef_tag)

            if 'chain_code' in nef_tag or 'sequence_code' in nef_tag:
                out[data_index] = tag_map[j]
            elif data in self.star_boolean:
                out[data_index] = 'true' if data in self.true_value else 'false'
            else:
                out[data_index] = data

        out_row.append(out)

        return out_row

    def star2nef_peak_can(self, in_sf, out_sf):
        """ Translate NMR-STAR canonical spectral peak loops (_Peak, _Peak_general_char, _Peak.char, _Assigned_peak_chem_shift) into NEF.
            @author: Masashi Yokochi
            @param in_sf: input saveframe
            @param out_sf: output saveframe
            @return: assigned chemical shift list id
        """

        try:
            num_dim = int(in_sf.get_tag('Number_of_spectral_dimensions')[0])
        except ValueError:
            return None

        try:
            pk_loop = in_sf.get_loop_by_category('_Peak')
        except KeyError:
            return None

        try:
            pk_gen_char_loop = in_sf.get_loop_by_category('_Peak_general_char')
        except KeyError:
            pk_gen_char_loop = None

        try:
            pk_char_loop = in_sf.get_loop_by_category('_Peak_char')
        except KeyError:
            return None

        try:
            pk_assign_loop = in_sf.get_loop_by_category('_Assigned_peak_chem_shift')
        except KeyError:
            pk_assign_loop = None

        out_lp = pynmrstar.Loop.from_scratch()
        out_tags = ['_nef_peak.index', '_nef_peak.peak_id', '_nef_peak.volume', '_nef_peak.volume_uncertainty', '_nef_peak.height', '_nef_peak.height_uncertainty']
        for d in range(1, num_dim + 1):
            out_tags.append('_nef_peak.position_%s' % d)
            out_tags.append('_nef_peak.position_uncertainty_%s' % d)
        for d in range(1, num_dim + 1):
            out_tags.append('_nef_peak.chain_code_%s' % d)
            out_tags.append('_nef_peak.sequence_code_%s' % d)
            out_tags.append('_nef_peak.residue_name_%s' % d)
            out_tags.append('_nef_peak.atom_name_%s' % d)

        for tag in out_tags:
            out_lp.add_tag(tag)

        pk_tags = pk_loop.get_tag_names()

        pk_id_col = pk_tags.index('_Peak.ID')

        if not pk_gen_char_loop is None:
            pk_gen_char_tags = pk_gen_char_loop.get_tag_names()
            pk_gen_char_id_col = pk_gen_char_tags.index('_Peak_general_char.Peak_ID')
            pk_gen_char_val_col = pk_gen_char_tags.index('_Peak_general_char.Intensity_val')
            pk_gen_char_val_err_col = pk_gen_char_tags.index('_Peak_general_char.Intensity_val_err')
            pk_gen_char_type_col = pk_gen_char_tags.index('_Peak_general_char.Measurement_method')

        pk_char_tags = pk_char_loop.get_tag_names()

        pk_char_id_col = pk_char_tags.index('_Peak_char.Peak_ID')
        pk_char_dim_id_col = pk_char_tags.index('_Peak_char.Spectral_dim_ID')
        pk_char_pos_col = pk_char_tags.index('_Peak_char.Chem_shift_val')
        pk_char_pos_err_col = pk_char_tags.index('_Peak_char.Chem_shift_val_err')

        if not pk_assign_loop is None:
            pk_assign_tags = pk_assign_loop.get_tag_names()
            pk_assign_id_col = pk_assign_tags.index('_Assigned_peak_chem_shift.Peak_ID')
            pk_assign_dim_id_col = pk_assign_tags.index('_Assigned_peak_chem_shift.Spectral_dim_ID')
            pk_assign_chain_id_col = pk_assign_tags.index('_Assigned_peak_chem_shift.Entity_assembly_ID')
            pk_assign_seq_id_col = pk_assign_tags.index('_Assigned_peak_chem_shift.Comp_index_ID')
            pk_assign_comp_id_col = pk_assign_tags.index('_Assigned_peak_chem_shift.Comp_ID')
            pk_assign_atom_id_col = pk_assign_tags.index('_Assigned_peak_chem_shift.Atom_ID')
            pk_assign_cs_list_id_col = pk_assign_tags.index('_Assigned_peak_chem_shift.Assigned_chem_shift_list_ID')

        cs_list_id_set = set()

        index = 1

        for pk in pk_loop.data:

            out = [None] * len(out_tags)

            pk_id = pk[pk_id_col]

            out[0] = index
            out[1] = pk_id

            if not pk_gen_char_loop is None:
                out[2], out[3] = next(((pk_gen_char[pk_gen_char_val_col], pk_gen_char[pk_gen_char_val_err_col])\
                                       for pk_gen_char in pk_gen_char_loop.data\
                                       if pk_gen_char[pk_gen_char_id_col] == pk_id and pk_gen_char[pk_gen_char_type_col] == 'volume'), (None, None))

                out[4], out[5] = next(((pk_gen_char[pk_gen_char_val_col], pk_gen_char[pk_gen_char_val_err_col])\
                                       for pk_gen_char in pk_gen_char_loop.data\
                                       if pk_gen_char[pk_gen_char_id_col] == pk_id and pk_gen_char[pk_gen_char_type_col] == 'height'), (None, None))

            l = 6

            for d in range(1, num_dim + 1):
                out[l], out[l + 1] = next(((pk_char[pk_char_pos_col], pk_char[pk_char_pos_err_col]) for pk_char in pk_char_loop.data\
                                           if pk_char[pk_char_id_col] == pk_id and int(pk_char[pk_char_dim_id_col]) == d), (None, None))
                l += 2

            if not pk_assign_loop is None:
                for d in range(1, num_dim + 1):
                    pk_assign = next((pk_assign for pk_assign in pk_assign_loop.data\
                                      if pk_assign[pk_assign_id_col] == pk_id and int(pk_assign[pk_assign_dim_id_col]) == d), None)

                    if not pk_assign is None:

                        _star_chain = pk_assign[pk_assign_chain_id_col]
                        if type(_star_chain) == str and not _star_chain in self.empty_value:
                            _star_chain = int(_star_chain)

                        _star_seq = pk_assign[pk_assign_seq_id_col]
                        if type(_star_seq) == str and not _star_seq in self.empty_value:
                            _star_seq = int(_star_seq)

                        seq_key = (_star_chain, _star_seq)

                        try:
                            nef_chain, nef_seq = self.authSeqMap[seq_key]
                        except KeyError:
                            try:
                                nef_chain = self.selfSeqMap[(_star_chain, 1)][0]
                            except KeyError:
                                if _star_chain in self.empty_value or not _star_chain in self.authChainId:
                                    nef_chain = _star_chain
                                else:
                                    cid = self.authChainId.index(_star_chain)
                                    nef_chain = self.index_to_letter(cid)
                            nef_seq = _star_seq

                        if _star_chain in self.star2cif_chain_mapping:
                            nef_chain = self.star2cif_chain_mapping[_star_chain]

                        out[l] = nef_chain
                        out[l + 1] = nef_seq

                        comp_id = pk_assign[pk_assign_comp_id_col]
                        atom_id = pk_assign[pk_assign_atom_id_col]

                        try:
                            if not self.atomIdMap is None:
                                _atom_id = self.atomIdMap[seq_key][atom_id]
                        except KeyError:
                            atom_list = self.get_nef_atom(comp_id, [{'atom_id': atom_id, 'ambig_code': None, 'value': None}])[0]
                            if len(atom_list) > 0:
                                _atom_id = atom_list[0]
                            else:
                                _atom_id = atom_id

                        out[l + 2] = comp_id
                        out[l + 3] = _atom_id

                        cs_list_id = pk_assign[pk_assign_cs_list_id_col]

                        if not cs_list_id in self.empty_value:
                            cs_list_id_set.add(cs_list_id)

                    l = l + 4

            out_lp.add_data(out)

            index += 1

        out_sf.add_loop(out_lp)

        return None if len(cs_list_id_set) == 0 else list(cs_list_id_set)[0]

    def nef_to_nmrstar(self, nef_file, star_file=None, report=None):
        """ Convert NEF file to NMR-STAR file.
            @param nef_file: input NEF file path
            @param star_file: output NMR-STAR file path
            @param report: NMR data processing report object (optional)
        """

        (file_path, file_name) = ntpath.split(os.path.realpath(nef_file))

        is_done = True
        info = []
        warning = []
        error = []

        if star_file is None:
            star_file = file_path + '/' + file_name.split('.')[0] + '.str'

        is_readable, dat_content, nef_data = self.read_input_file(nef_file)

        try:
            star_data = pynmrstar.Entry.from_scratch(nef_data.entry_id)
        except: # AttributeError:
            star_data = pynmrstar.Entry.from_scratch(file_name.split('.')[0])
            warning.append('Not a complete Entry')

        if is_readable:

            if dat_content == 'Entry':
                if len(nef_data.get_loops_by_category('nef_sequence')) == 0: # DAOTHER-6694
                    error.append("Missing mandatory '_nef_sequence' category.")
                    return False, json.dumps({'info': info, 'warning': warning, 'error': error})
                self.authChainId = sorted(list(set(nef_data.get_loops_by_category('nef_sequence')[0].get_tag('chain_code'))))
            elif dat_content == 'Saveframe':
                self.authChainId = sorted(list(set(nef_data[0].get_tag('chain_code'))))
            elif dat_content == 'Loop':
                self.authChainId = sorted(list(set(nef_data.get_tag('chain_code'))))
            else:
                is_done = False
                error.append('File content unknown.')

            self.authSeqMap = None
            self.selfSeqMap = None

            asm_id = 0
            cs_list_id = 0
            dist_list_id = 0
            dihed_list_id = 0
            rdc_list_id = 0
            peak_list_id = 0

            if dat_content == 'Entry':

                for saveframe in nef_data:
                    sf = pynmrstar.Saveframe.from_scratch(saveframe.name)

                    if saveframe.category == 'nef_nmr_meta_data':
                        sf.set_tag_prefix('Entry')

                    elif saveframe.category == 'nef_molecular_system':
                        asm_id += 1
                        sf.set_tag_prefix('Assembly')

                    elif saveframe.category == 'nef_chemical_shift_list':
                        cs_list_id += 1
                        sf.set_tag_prefix('Assigned_chem_shift_list')

                    elif saveframe.category == 'nef_distance_restraint_list':
                        dist_list_id += 1
                        sf.set_tag_prefix('Gen_dist_constraint_list')

                    elif saveframe.category == 'nef_dihedral_restraint_list':
                        dihed_list_id += 1
                        sf.set_tag_prefix('Torsion_angle_constraint_list')

                    elif saveframe.category == 'nef_rdc_restraint_list':
                        rdc_list_id += 1
                        sf.set_tag_prefix('RDC_constraint_list')

                    elif saveframe.category == 'nef_nmr_spectrum':
                        peak_list_id += 1
                        sf.set_tag_prefix('Spectral_peak_list')

                    else:
                        continue

                    for tag in saveframe.tags:

                        if tag[0].lower() == 'sf_category':
                            auth_tag, data_tag, _ = self.get_star_tag(saveframe.category)
                            if not auth_tag is None:
                                sf.add_tag('Sf_category', auth_tag)
                        elif saveframe.category == 'nef_distance_restraint_list' and tag[0] == 'restraint_origin':
                            nef_tag = '{}.{}'.format(saveframe.tag_prefix, tag[0])
                            auth_tag, data_tag, _ = self.get_star_tag(nef_tag)
                            if not auth_tag is None:
                                sf.add_tag(auth_tag, tag[1] if not tag[1] in self.dist_alt_constraint_type['nmr-star'] else self.dist_alt_constraint_type['nmr-star'][tag[1]])
                        elif saveframe.category == 'nef_dihedral_restraint_list' and tag[0] == 'restraint_origin':
                            nef_tag = '{}.{}'.format(saveframe.tag_prefix, tag[0])
                            auth_tag, data_tag, _ = self.get_star_tag(nef_tag)
                            if not auth_tag is None:
                                sf.add_tag(auth_tag, tag[1] if not tag[1] in self.dihed_alt_constraint_type['nmr-star'] else self.dihed_alt_constraint_type['nmr-star'][tag[1]])
                        elif saveframe.category == 'nef_rdc_restraint_list' and tag[0] == 'restraint_origin':
                            nef_tag = '{}.{}'.format(saveframe.tag_prefix, tag[0])
                            auth_tag, data_tag, _ = self.get_star_tag(nef_tag)
                            if not auth_tag is None:
                                sf.add_tag(auth_tag, tag[1] if not tag[1] in self.rdc_alt_constraint_type['nmr-star'] else self.rdc_alt_constraint_type['nmr-star'][tag[1]])
                        else:
                            nef_tag = '{}.{}'.format(saveframe.tag_prefix, tag[0])
                            auth_tag, data_tag, _ = self.get_star_tag(nef_tag)
                            if not auth_tag is None:
                                sf.add_tag(auth_tag, tag[1])

                    has_covalent_links = any(loop for loop in saveframe if loop.category == '_nef_covalent_links')
                    aux_rows = None

                    for loop in saveframe:

                        lp = pynmrstar.Loop.from_scratch()
                        tags = self.get_star_loop_tags(loop.get_tag_names())

                        if len(tags) == 0:
                            continue

                        for tag in tags:
                            lp.add_tag(tag)

                        if loop.category == '_nef_sequence':
                            if self.authSeqMap is None:
                                self.authSeqMap = {}
                                self.selfSeqMap = {}
                            rows, aux_rows = self.nef2star_seq_row(loop.get_tag_names(), lp.get_tag_names(), loop.data, report)
                            for d in rows:
                                d[lp.get_tag_names().index('_Chem_comp_assembly.Assembly_ID')] = asm_id
                                lp.add_data(d)

                        elif loop.category == '_nef_covalent_links':
                            rows = self.nef2star_bond_row(loop.get_tag_names(), lp.get_tag_names(), loop.data)
                            for d in rows:
                                d[lp.get_tag_names().index('_Bond.Assembly_ID')] = asm_id
                                lp.add_data(d)

                        elif loop.category == '_nef_chemical_shift':
                            rows = self.nef2star_cs_row(loop.get_tag_names(), lp.get_tag_names(), loop.data)
                            for d in rows:
                                d[lp.get_tag_names().index('_Atom_chem_shift.Assigned_chem_shift_list_ID')] = cs_list_id
                                lp.add_data(d)

                        elif loop.category == '_nef_distance_restraint':
                            rows = self.nef2star_dist_row(loop.get_tag_names(), lp.get_tag_names(), loop.data)
                            for d in rows:
                                d[lp.get_tag_names().index('_Gen_dist_constraint.Gen_dist_constraint_list_ID')] = dist_list_id
                                lp.add_data(d)

                        elif loop.category == '_nef_dihedral_restraint':
                            rows = self.nef2star_dihed_row(loop.get_tag_names(), lp.get_tag_names(), loop.data)
                            for d in rows:
                                d[lp.get_tag_names().index('_Torsion_angle_constraint.Torsion_angle_constraint_list_ID')] = dihed_list_id
                                lp.add_data(d)

                        elif loop.category == '_nef_rdc_restraint':
                            rows = self.nef2star_rdc_row(loop.get_tag_names(), lp.get_tag_names(), loop.data)
                            for d in rows:
                                d[lp.get_tag_names().index('_RDC_constraint.RDC_constraint_list_ID')] = rdc_list_id
                                lp.add_data(d)

                        elif loop.category == '_nef_peak':
                            rows = self.nef2star_peak_row(loop.get_tag_names(), lp.get_tag_names(), loop.data)
                            for d in rows:
                                d[lp.get_tag_names().index('_Peak_row_format.Spectral_peak_list_ID')] = peak_list_id
                                lp.add_data(d)

                        else:

                            for data in loop.data:

                                if loop.category == '_nef_spectrum_dimension':
                                    rows = self.nef2star_row(loop.get_tag_names(), lp.get_tag_names(), data)
                                    for d in rows:
                                        d[lp.get_tag_names().index('_Spectral_dim.Spectral_peak_list_ID')] = peak_list_id
                                        lp.add_data(d)

                                elif loop.category == '_nef_spectrum_dimension_transfer':
                                    rows = self.nef2star_row(loop.get_tag_names(), lp.get_tag_names(), data)
                                    for d in rows:
                                        d[lp.get_tag_names().index('_Spectral_dim_transfer.Spectral_peak_list_ID', )] = peak_list_id
                                        lp.add_data(d)

                                else:
                                    rows = self.nef2star_row(loop.get_tag_names(), lp.get_tag_names(), data)
                                    for d in rows:
                                        lp.add_data(d)

                        sf.add_loop(lp)

                        if not aux_rows is None and ((loop.category == '_nef_sequence' and not has_covalent_links) or\
                                                     (loop.category == '_nef_covalent_links' and has_covalent_links)):
                            lp = pynmrstar.Loop.from_scratch()
                            for _tag in self.entity_del_atom_row:
                                lp.add_tag('_Entity_deleted_atom.%s' % _tag)
                            for d in aux_rows:
                                d[lp.get_tag_names().index('_Entity_deleted_atom.Assembly_ID')] = asm_id
                                lp.add_data(d)
                            sf.add_loop(lp)

                    if saveframe.category == 'nef_nmr_meta_data':
                        sf.add_tag('NMR_STAR_version', self.star_version)

                        try:
                            loop = sf.get_loop_by_category('_Software_applied_methods')
                            row = []
                            for t in loop.tags:
                                if t == 'Software_name':
                                    row.append('NEFTranslator')
                                elif t == 'Script_name':
                                    row.append('nef_to_nmrstar')
                                else:
                                    row.append('.')
                            loop.add_data(row)
                        except KeyError:
                            pass

                    elif saveframe.category == 'nef_molecular_system':
                        sf.add_tag('ID', asm_id)

                    elif saveframe.category == 'nef_chemical_shift_list':
                        sf.add_tag('ID', cs_list_id)

                    elif saveframe.category == 'nef_distance_restraint_list':
                        sf.add_tag('ID', dist_list_id)

                    elif saveframe.category == 'nef_dihedral_restraint_list':
                        sf.add_tag('ID', dihed_list_id)

                    elif saveframe.category == 'nef_rdc_restraint_list':
                        sf.add_tag('ID', rdc_list_id)

                    elif saveframe.category == 'nef_nmr_spectrum':
                        sf.add_tag('ID', peak_list_id)

                    else:
                        continue

                    star_data.add_saveframe(sf)

            elif dat_content == 'Saveframe' or dat_content == 'Loop':

                if dat_content == 'Saveframe':
                    saveframe = nef_data
                    sf = pynmrstar.Saveframe.from_scratch(saveframe.name)

                    if saveframe.category == 'nef_nmr_meta_data':
                        sf.set_tag_prefix('Entry')

                    elif saveframe.category == 'nef_molecular_system':
                        asm_id += 1
                        sf.set_tag_prefix('Assembly')

                    elif saveframe.category == 'nef_chemical_shift_list':
                        cs_list_id += 1
                        sf.set_tag_prefix('Assigned_chem_shift_list')

                    elif saveframe.category == 'nef_distance_restraint_list':
                        dist_list_id += 1
                        sf.set_tag_prefix('Gen_dist_constraint_list')

                    elif saveframe.category == 'nef_dihedral_restraint_list':
                        dihed_list_id += 1
                        sf.set_tag_prefix('Torsion_angle_constraint_list')

                    elif saveframe.category == 'nef_rdc_restraint_list':
                        rdc_list_id += 1
                        sf.set_tag_prefix('RDC_constraint_list')

                    elif saveframe.category == 'nef_nmr_spectrum':
                        peak_list_id += 1
                        sf.set_tag_prefix('Spectral_peak_list')

                    for tag in saveframe.tags:

                        if tag[0].lower() == 'sf_category':
                            auth_tag, data_tag, _ = self.get_star_tag(saveframe.category)
                            if not auth_tag is None:
                                sf.add_tag('Sf_category', auth_tag)
                        elif saveframe.category == 'nef_distance_restraint_list' and tag[0] == 'restraint_origin':
                            nef_tag = '{}.{}'.format(saveframe.tag_prefix, tag[0])
                            auth_tag, data_tag, _ = self.get_star_tag(nef_tag)
                            if not auth_tag is None:
                                sf.add_tag(auth_tag, tag[1] if not tag[1] in self.dist_alt_constraint_type['nmr-star'] else self.dist_alt_constraint_type['nmr-star'][tag[1]])
                        elif saveframe.category == 'nef_dihedral_restraint_list' and tag[0] == 'restraint_origin':
                            nef_tag = '{}.{}'.format(saveframe.tag_prefix, tag[0])
                            auth_tag, data_tag, _ = self.get_star_tag(nef_tag)
                            if not auth_tag is None:
                                sf.add_tag(auth_tag, tag[1] if not tag[1] in self.dihed_alt_constraint_type['nmr-star'] else self.dihed_alt_constraint_type['nmr-star'][tag[1]])
                        elif saveframe.category == 'nef_rdc_restraint_list' and tag[0] == 'restraint_origin':
                            nef_tag = '{}.{}'.format(saveframe.tag_prefix, tag[0])
                            auth_tag, data_tag, _ = self.get_star_tag(nef_tag)
                            if not auth_tag is None:
                                sf.add_tag(auth_tag, tag[1] if not tag[1] in self.rdc_alt_constraint_type['nmr-star'] else self.rdc_alt_constraint_type['nmr-star'][tag[1]])
                        else:
                            nef_tag = '{}.{}'.format(saveframe.tag_prefix, tag[0])
                            auth_tag, data_tag, _ = self.get_star_tag(nef_tag)
                            if not auth_tag is None:
                                sf.add_tag(auth_tag, tag[1])

                else:
                    sf = pynmrstar.Saveframe.from_scratch(nef_data.category)

                    if nef_data.category == '_nef_program_script':
                        sf.set_tag_prefix('Entry')
                        sf.add_tag('Sf_category', 'entry_information')

                    elif nef_data.category == '_audit': # DAOTHER-6327
                        sf.set_tag_prefix('Entry')
                        sf.add_tag('Sf_category', 'entry_information')

                    elif nef_data.category == '_nef_sequence':
                        asm_id += 1
                        sf.set_tag_prefix('Assembly')
                        sf.add_tag('Sf_category', 'assembly')

                    elif nef_data.category == '_nef_chemical_shift':
                        cs_list_id += 1
                        sf.set_tag_prefix('Assigned_chem_shift_list')
                        sf.add_tag('Sf_category', 'assigned_chemical_shifts')

                    elif nef_data.category == '_nef_distance_restraint':
                        dist_list_id += 1
                        sf.set_tag_prefix('Gen_dist_constraint_list')
                        sf.add_tag('Sf_category', 'general_distance_constraints')

                    elif nef_data.category == '_nef_dihedral_restraint':
                        dihed_list_id += 1
                        sf.set_tag_prefix('Torsion_angle_constraint_list')
                        sf.add_tag('Sf_category', 'torsion_angle_constraints')

                    elif nef_data.category == '_nef_rdc_restraint':
                        rdc_list_id += 1
                        sf.set_tag_prefix('RDC_constraint_list')
                        sf.add_tag('Sf_category', 'RDC_constraints')

                    elif nef_data.category == '_nef_peak':
                        peak_list_id += 1
                        sf.set_tag_prefix('Spectral_peak_list')
                        sf.add_tag('Sf_category', 'spectral_peak_list')

                    sf.add_tag('Sf_framecode', nef_data.category)

                    saveframe = [nef_data]

                has_covalent_links = any(loop for loop in saveframe if loop.category == '_nef_covalent_links')
                aux_rows = None

                for loop in saveframe:

                    lp = pynmrstar.Loop.from_scratch()
                    tags = self.get_star_loop_tags(loop.get_tag_names())

                    if len(tags) == 0:
                        continue

                    for tag in tags:
                        lp.add_tag(tag)

                    if loop.category == '_nef_sequence':
                        if self.authSeqMap is None:
                            self.authSeqMap = {}
                            self.selfSeqMap = {}
                        rows, aux_rows = self.nef2star_seq_row(loop.get_tag_names(), lp.get_tag_names(), loop.data, report)
                        for d in rows:
                            d[lp.get_tag_names().index('_Chem_comp_assembly.Assembly_ID')] = asm_id
                            lp.add_data(d)

                    elif loop.category == '_nef_covalent_links':
                        rows = self.nef2star_bond_row(loop.get_tag_names(), lp.get_tag_names(), loop.data)
                        for d in rows:
                            d[lp.get_tag_names().index('_Bond.Assembly_ID')] = asm_id
                            lp.add_data(d)

                    elif loop.category == '_nef_chemical_shift':
                        rows = self.nef2star_cs_row(loop.get_tag_names(), lp.get_tag_names(), loop.data)
                        for d in rows:
                            d[lp.get_tag_names().index('_Atom_chem_shift.Assigned_chem_shift_list_ID')] = cs_list_id
                            lp.add_data(d)

                    elif loop.category == '_nef_distance_restraint':
                        rows = self.nef2star_dist_row(loop.get_tag_names(), lp.get_tag_names(), loop.data)
                        for d in rows:
                            d[lp.get_tag_names().index('_Gen_dist_constraint.Gen_dist_constraint_list_ID')] = dist_list_id
                            lp.add_data(d)

                    elif loop.category == '_nef_dihedral_restraint':
                        rows = self.nef2star_dihed_row(loop.get_tag_names(), lp.get_tag_names(), loop.data)
                        for d in rows:
                            d[lp.get_tag_names().index('_Torsion_angle_constraint.Torsion_angle_constraint_list_ID')] = dihed_list_id
                            lp.add_data(d)

                    elif loop.category == '_nef_rdc_restraint':
                        rows = self.nef2star_rdc_row(loop.get_tag_names(), lp.get_tag_names(), loop.data)
                        for d in rows:
                            d[lp.get_tag_names().index('_RDC_constraint.RDC_constraint_list_ID')] = rdc_list_id
                            lp.add_data(d)

                    elif loop.category == '_nef_peak':
                        rows = self.nef2star_peak_row(loop.get_tag_names(), lp.get_tag_names(), loop.data)
                        for d in rows:
                            d[lp.get_tag_names().index('_Peak_row_format.Spectral_peak_list_ID')] = peak_list_id
                            lp.add_data(d)

                    else:

                        for data in loop.data:

                            if loop.category == '_nef_spectrum_dimension':
                                rows = self.nef2star_row(loop.get_tag_names(), lp.get_tag_names(), data)
                                for d in rows:
                                    d[lp.get_tag_names().index('_Spectral_dim.Spectral_peak_list_ID')] = peak_list_id
                                    lp.add_data(d)

                            elif loop.category == '_nef_spectrum_dimension_transfer':
                                rows = self.nef2star_row(loop.get_tag_names(), lp.get_tag_names(), data)
                                for d in rows:
                                    d[lp.get_tag_names().index('_Spectral_dim_transfer.Spectral_peak_list_ID', )] = peak_list_id
                                    lp.add_data(d)

                            else:
                                rows = self.nef2star_row(loop.get_tag_names(), lp.get_tag_names(), data)
                                for d in rows:
                                    lp.add_data(d)

                    sf.add_loop(lp)

                    if not aux_rows is None and ((loop.category == '_nef_sequence' and not has_covalent_links) or\
                                                 (loop.category == '_nef_covalent_links' and has_covalent_links)):
                        lp = pynmrstar.Loop.from_scratch()
                        for _tag in self.entity_del_atom_row:
                            lp.add_tag('_Entity_deleted_atom.%s' % _tag)
                        for d in aux_rows:
                            d[lp.get_tag_names().index('_Entity_deleted_atom.Assembly_ID')] = asm_id
                            lp.add_data(d)
                        sf.add_loop(lp)

                if saveframe.category == 'nef_nmr_meta_data':
                    sf.add_tag('NMR_STAR_version', self.star_version)

                    try:
                        loop = sf.get_loop_by_category('_Software_applied_methods')
                        row = []
                        for t in loop.tags:
                            if t == 'Software_name':
                                row.append('NEFTranslator')
                            elif t == 'Script_name':
                                row.append('nef_to_nmrstar')
                            else:
                                row.append('.')
                        loop.add_data(row)
                    except KeyError:
                        pass

                elif saveframe.category == 'nef_molecular_system':
                    sf.add_tag('ID', asm_id)

                elif saveframe.category == 'nef_chemical_shift_list':
                    sf.add_tag('ID', cs_list_id)

                elif saveframe.category == 'nef_distance_restraint_list':
                    sf.add_tag('ID', dist_list_id)

                elif saveframe.category == 'nef_dihedral_restraint_list':
                    sf.add_tag('ID', dihed_list_id)

                elif saveframe.category == 'nef_rdc_restraint_list':
                    sf.add_tag('ID', rdc_list_id)

                elif saveframe.category == 'nef_nmr_spectrum':
                    sf.add_tag('ID', peak_list_id)

                star_data.add_saveframe(sf)

            if is_done:
                #star_data.normalize() # do not invoke normalize() to preserve author provided Peak_row_format.ID using pynmrstar v3 library

                if __pynmrstar_v3__:
                    star_data.write_to_file(star_file, skip_empty_tags=False)
                else:
                    star_data.write_to_file(star_file)
                info.append('File {} successfully written.'.format(star_file))

        else:
            is_done = False
            error.append('Input file not readable.')

        return is_done, json.dumps({'info': info, 'warning': warning, 'error': error})

    def nmrstar_to_nef(self, star_file, nef_file=None, report=None):
        """ Convert NMR-STAR file to NEF file.
            @author: Masashi Yokochi
            @param star_file: input NMR-STAR file path
            @param nef_file: output NEF file path
            @param report: NMR data processing report object (optional)
        """

        (file_path, file_name) = ntpath.split(os.path.realpath(star_file))

        is_done = True
        info = []
        warning = []
        error = []

        if nef_file is None:
            nef_file = file_path + '/' + file_name.split('.')[0] + '.nef'

        is_readable, dat_content, star_data = self.read_input_file(star_file)

        try:
            nef_data = pynmrstar.Entry.from_scratch(star_data.entry_id)
        except: # AttributeError:
            nef_data = pynmrstar.Entry.from_scratch(file_name.split('.')[0])
            warning.append('Not a complete Entry.')

        if is_readable:

            if dat_content == 'Entry':
                if len(star_data.get_loops_by_category('Chem_comp_assembly')) == 0: # DAOTHER-6694
                    error.append("Missing mandatory '_Chem_comp_assembly' category.")
                    return False, json.dumps({'info': info, 'warning': warning, 'error': error})
                self.authChainId = sorted(list(set(star_data.get_loops_by_category('Chem_comp_assembly')[0].get_tag('Entity_assembly_ID'))), key=lambda x:float(re.sub('[^\d]+', '', x)))
            elif dat_content == 'Saveframe':
                self.authChainId = sorted(list(set(star_data[0].get_tag('Entity_assembly_ID'))), key=lambda x:float(re.sub('[^\d]+', '', x)))
            elif dat_content == 'Loop':
                self.authChainId = sorted(list(set(star_data.get_tag('Entity_assembly_ID'))), key=lambda x:float(re.sub('[^\d]+', '', x)))
            else:
                is_done = False
                error.append('File content unknown.')

            self.authSeqMap = None
            self.selfSeqMap = None
            self.atomIdMap = None

            if dat_content == 'Entry':

                for saveframe in star_data:
                    sf = pynmrstar.Saveframe.from_scratch(saveframe.name)

                    nef_tag, _ = self.get_nef_tag(saveframe.category)

                    if nef_tag is None:
                        continue

                    sf.set_tag_prefix(nef_tag)

                    sf.add_tag('sf_category', nef_tag)
                    sf.add_tag('sf_framecode', saveframe.name)

                    for tag in saveframe.tags:
                        tag_name = tag[0].lower()
                        if tag_name == 'sf_category' or tag_name == 'sf_framecode':
                            continue
                        elif saveframe.category == 'entry_information':
                            if tag_name == 'source_data_format':
                                sf.add_tag('format_name', self.nef_format_name)
                            elif tag_name == 'source_data_format_version':
                                sf.add_tag('format_version', self.nef_version)
                            else:
                                nef_tag, _ = self.get_nef_tag(saveframe.tag_prefix + '.' + tag[0])
                                if not nef_tag is None:
                                    sf.add_tag(nef_tag, tag[1])
                        elif saveframe.category == 'general_distance_constraints' and tag_name == 'constraint_type':
                            nef_tag, _ = self.get_nef_tag(saveframe.tag_prefix + '.' + tag[0])
                            if not nef_tag is None:
                                sf.add_tag(nef_tag, tag[1] if not tag[1] in self.dist_alt_constraint_type['nef'] else self.dist_alt_constraint_type['nef'][tag[1]])
                        elif saveframe.category == 'torsion_angle_constraints' and tag_name == 'constraint_type':
                            nef_tag, _ = self.get_nef_tag(saveframe.tag_prefix + '.' + tag[0])
                            if not nef_tag is None:
                                sf.add_tag(nef_tag, tag[1] if not tag[1] in self.dihed_alt_constraint_type['nef'] else self.dihed_alt_constraint_type['nef'][tag[1]])
                        elif saveframe.category == 'RDC_constraints' and tag_name == 'constraint_type':
                            nef_tag, _ = self.get_nef_tag(saveframe.tag_prefix + '.' + tag[0])
                            if not nef_tag is None:
                                sf.add_tag(nef_tag, tag[1] if not tag[1] in self.rdc_alt_constraint_type['nef'] else self.rdc_alt_constraint_type['nef'][tag[1]])
                        else:
                            nef_tag, _ = self.get_nef_tag(saveframe.tag_prefix + '.' + tag[0])
                            if not nef_tag is None:
                                sf.add_tag(nef_tag, tag[1])

                    entity_del_atom_loop = next((loop for loop in saveframe if loop.category == '_Entity_deleted_atom'), None)

                    has_pk_can_format = False
                    has_pk_row_format = False

                    for loop in saveframe:

                        try:

                            lp = pynmrstar.Loop.from_scratch()
                            tags = self.get_nef_loop_tags(loop.get_tag_names())
                            tag_set = set(tags)
                            if len(tags) > len(tag_set):
                                tags = list(tag_set)
                            for tag in tags:
                                lp.add_tag(tag)

                            if loop.category == '_Chem_comp_assembly':
                                if self.authSeqMap is None:
                                    self.authSeqMap = {}
                                    self.selfSeqMap = {}
                                rows = self.star2nef_seq_row(loop.get_tag_names(), lp.get_tag_names(), loop.data, report, None if entity_del_atom_loop is None else entity_del_atom_loop)
                                for d in rows:
                                    lp.add_data(d)

                            elif loop.category == '_Atom_chem_shift':
                                if self.atomIdMap is None:
                                    self.atomIdMap = {}

                                rows = self.star2nef_cs_row(loop.get_tag_names(), lp.get_tag_names(), loop.data)
                                for d in rows:
                                    lp.add_data(d)

                            elif loop.category == '_Gen_dist_constraint':
                                rows = self.star2nef_dist_row(loop.get_tag_names(), lp.get_tag_names(), loop.data)
                                for d in rows:
                                    lp.add_data(d)

                            elif loop.category == '_Peak_row_format':
                                rows = self.star2nef_peak_row(loop.get_tag_names(), lp.get_tag_names(), loop.data)
                                for d in rows:
                                    lp.add_data(d)
                                has_pk_row_format = True

                            elif len(tags) > 0:
                                for data in loop.data:
                                    rows = self.star2nef_row(loop.get_tag_names(), lp.get_tag_names(), data)
                                    for d in rows:
                                        lp.add_data(d)

                            elif loop.category == '_Peak':
                                has_pk_can_format = True

                            sf.add_loop(lp)

                        except ValueError:
                            pass

                    if saveframe.category == 'entry_information':
                        has_format_name = False
                        has_format_ver = False
                        for tags in sf.tags:
                            if tags[0] == 'format_name':
                                has_format_name = True
                            elif tags[0] == 'format_version':
                                has_format_ver = True

                        if not has_format_name:
                            sf.add_tag('format_name', self.nef_format_name)
                        if not has_format_ver:
                            sf.add_tag('format_version', self.nef_version)

                        try:
                            loop = sf.get_loop_by_category('_nef_program_script')
                            row = []
                            for t in loop.tags:
                                if t == 'program_name':
                                    row.append('NEFTranslator')
                                elif t == 'script_name':
                                    row.append('nmrstar_to_nef')
                                else:
                                    row.append('.')
                            loop.add_data(row)
                        except KeyError:
                            pass

                    if saveframe.category == 'spectral_peak_list' and has_pk_can_format and not has_pk_row_format:
                        cs_list_id = self.star2nef_peak_can(saveframe, sf)
                        if not cs_list_id is None and (len(sf.get_tag('chemical_shift_list')) == 0 or sf.get_tag('chemical_shift_list') in self.empty_value):
                            for cs_sf in star_data:
                                if cs_sf.get_tag('Sf_category')[0] == 'assigned_chemical_shifts' and cs_sf.get_tag('ID')[0] == cs_list_id and not cs_sf.name in self.empty_value:
                                    if len(sf.get_tag('chemical_shift_list')) == 0:
                                        sf.add_tag('chemical_shift_list', cs_sf.name)
                                    else:
                                        sf.tags[sf.tags.index('chemical_shift_list')][1] = cs_sf.name
                                    break

                    nef_data.add_saveframe(sf)

            elif dat_content == 'Saveframe' or dat_content == 'Loop':

                if dat_content == 'Saveframe':
                    saveframe = star_data
                    sf = pynmrstar.Saveframe.from_scratch(saveframe.name)

                    nef_tag, _ = self.get_nef_tag(saveframe.category)

                    if not nef_tag is None:

                        sf.set_tag_prefix(nef_tag)

                        sf.add_tag('sf_category', nef_tag)
                        sf.add_tag('sf_framecode', saveframe.name)

                        for tag in saveframe.tags:
                            tag_name = tag[0].lower()
                            if tag_name == 'sf_category' or tag_name == 'sf_framecode':
                                continue
                            elif saveframe.category == 'entry_information':
                                if tag_name == 'source_data_format':
                                    sf.add_tag('format_name', self.nef_format_name)
                                elif tag_name == 'source_data_format_version':
                                    sf.add_tag('format_version', self.nef_version)
                                else:
                                    nef_tag, _ = self.get_nef_tag(saveframe.tag_prefix + '.' + tag[0])
                                    if not nef_tag is None:
                                        sf.add_tag(nef_tag, tag[1])
                            elif saveframe.category == 'general_distance_constraints' and tag_name == 'constraint_type':
                                nef_tag, _ = self.get_nef_tag(saveframe.tag_prefix + '.' + tag[0])
                                if not nef_tag is None:
                                    sf.add_tag(nef_tag, tag[1] if not tag[1] in self.dist_alt_constraint_type['nef'] else self.dist_alt_constraint_type['nef'][tag[1]])
                            elif saveframe.category == 'torsion_angle_constraints' and tag_name == 'constraint_type':
                                nef_tag, _ = self.get_nef_tag(saveframe.tag_prefix + '.' + tag[0])
                                if not nef_tag is None:
                                    sf.add_tag(nef_tag, tag[1] if not tag[1] in self.dihed_alt_constraint_type['nef'] else self.dihed_alt_constraint_type['nef'][tag[1]])
                            elif saveframe.category == 'RDC_constraints' and tag_name == 'constraint_type':
                                nef_tag, _ = self.get_nef_tag(saveframe.tag_prefix + '.' + tag[0])
                                if not nef_tag is None:
                                    sf.add_tag(nef_tag, tag[1] if not tag[1] in self.rdc_alt_constraint_type['nef'] else self.rdc_alt_constraint_type['nef'][tag[1]])
                            else:
                                nef_tag, _ = self.get_nef_tag(saveframe.tag_prefix + '.' + tag[0])
                                if not nef_tag is None:
                                    sf.add_tag(nef_tag, tag[1])

                else:
                    sf = pynmrstar.Saveframe.from_scratch(star_data.category)

                    if star_data.category == '__Software_applied_methods':
                        sf.set_tag_prefix('nef_nmr_meta_data')
                        sf.add_tag('sf_category', 'nef_molecular_system')

                    elif star_data.category == '_Assembly':
                        sf.set_tag_prefix('nef_molecular_system')
                        sf.add_tag('sf_category', 'nef_molecular_system')

                    elif star_data.category == '_Atom_chem_shift':
                        sf.set_tag_prefix('nef_chemical_shift_list')
                        sf.add_tag('sf_category', 'nef_chemical_shift_list')

                    elif star_data.category == '_Gen_dist_constraint':
                        sf.set_tag_prefix('nef_distance_restraint_list')
                        sf.add_tag('sf_category', 'nef_distance_restraint_list')

                    elif star_data.category == '_RDC_constraint':
                        sf.set_tag_prefix('nef_rdc_restraint_list')
                        sf.add_tag('sf_category', 'nef_rdc_restraint_list')

                    elif star_data.category == '_Peak_row_format':
                        sf.set_tag_prefix('nef_nmr_spectrum')
                        sf.add_tag('sf_category', 'nef_nmr_spectrum')

                    sf.add_tag('sf_framecode', star_data.category)

                    saveframe = [star_data]

                entity_del_atom_loop = next((loop for loop in saveframe if loop.category == '_Entity_deleted_atom'), None)

                has_pk_can_format = False
                has_pk_row_format = False

                for loop in saveframe:

                    try:

                        lp = pynmrstar.Loop.from_scratch()
                        tags = self.get_nef_loop_tags(loop.get_tag_names())
                        tag_set = set(tags)
                        if len(tags) > len(tag_set):
                            tags = list(tag_set)
                        for tag in tags:
                            lp.add_tag(tag)

                        if loop.category == '_Chem_comp_assembly':
                            if self.authSeqMap is None:
                                self.authSeqMap = {}
                                self.selfSeqMap = {}
                            rows = self.star2nef_seq_row(loop.get_tag_names(), lp.get_tag_names(), loop.data, report, None if entity_del_atom_loop is None else entity_del_atom_loop)
                            for d in rows:
                                lp.add_data(d)

                        elif loop.category == '_Atom_chem_shift':
                            if self.atomIdMap is None:
                                self.atomIdMap = {}

                            rows = self.star2nef_cs_row(loop.get_tag_names(), lp.get_tag_names(), loop.data)
                            for d in rows:
                                lp.add_data(d)

                        elif loop.category == '_Gen_dist_constraint':
                            rows = self.star2nef_dist_row(loop.get_tag_names(), lp.get_tag_names(), loop.data)
                            for d in rows:
                                lp.add_data(d)

                        elif loop.category == '_Peak_row_format':
                            rows = self.star2nef_peak_row(loop.get_tag_names(), lp.get_tag_names(), loop.data)
                            for d in rows:
                                lp.add_data(d)
                            has_pk_row_format = True

                        elif len(tags) > 0:
                            for data in loop.data:
                                rows = self.star2nef_row(loop.get_tag_names(), lp.get_tag_names(), data)
                                for d in rows:
                                    lp.add_data(d)

                        elif loop.category == '_Peak':
                             has_pk_can_format = True

                        sf.add_loop(lp)

                    except ValueError:
                        pass

                if saveframe.category == 'entry_information':
                    has_format_name = False
                    has_format_ver = False
                    for tags in sf.tags:
                        if tags[0] == 'format_name':
                            has_format_name = True
                        elif tags[0] == 'format_version':
                            has_format_ver = True

                    if not has_format_name:
                        sf.add_tag('format_name', self.nef_format_name)
                    if not has_format_ver:
                        sf.add_tag('format_version', self.nef_version)

                    try:
                        loop = sf.get_loop_by_category('_nef_program_script')
                        row = []
                        for t in loop.tags:
                            if t == 'program_name':
                                row.append('NEFTranslator')
                            elif t == 'script_name':
                                row.append('nmrstar_to_nef')
                            else:
                                row.append('.')
                        loop.add_data(row)
                    except KeyError:
                        pass

                if saveframe.category == 'spectral_peak_list' and has_pk_can_format and not has_pk_row_format:
                    cs_list_id = self.star2nef_peak_can(saveframe, sf)
                    if not cs_list_id is None and (len(sf.get_tag('chemical_shift_list')) == 0 or sf.get_tag('chemical_shift_list') in self.empty_value):
                        for cs_sf in star_data:
                            if cs_sf.get_tag('Sf_category')[0] == 'assigned_chemical_shifts' and cs_sf.get_tag('ID')[0] == cs_list_id and not cs_sf.name in self.empty_value:
                                if len(sf.get_tag('chemical_shift_list')) == 0:
                                    sf.add_tag('chemical_shift_list', cs_sf.name)
                                else:
                                    sf.tags[sf.tags.index('chemical_shift_list')][1] = cs_sf.name
                                break

                nef_data.add_saveframe(sf)

            if is_done:
                if __pynmrstar_v3__:
                    nef_data.write_to_file(nef_file, skip_empty_tags=False)
                else:
                    nef_data.write_to_file(nef_file)
                info.append('File {} successfully written.'.format(nef_file))

        else:
            is_done = False
            error.append('Input file not readable.')

        return is_done, json.dumps({'info': info, 'warning': warning, 'error': error})

if __name__ == "__main__":
    _nefT = NEFTranslator()
    _nefT.nef_to_nmrstar('data/2l9r.nef')
    print(_nefT.validate_file('data/2l9r.str','A'))
