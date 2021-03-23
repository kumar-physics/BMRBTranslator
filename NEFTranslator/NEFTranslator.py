#!/usr/bin/env python

"""
This module does the following jobs
1. Validate the NEF and NMR-STAR files
2. Extract the sequence information
3. Format conversion

@author: Kumaran Baskaran

"""

# Make sure print function work in python2 and python3
from __future__ import print_function

import sys
import os
import ntpath
import json
import logging
import re
import csv
import datetime
import pynmrstar
from pytz import utc

PY3 = (sys.version_info[0] == 3)

(scriptPath, scriptName) = ntpath.split(os.path.realpath(__file__))

__version__ = "v1.2.0"


class NEFTranslator(object):
    """
    NEF to NMR-STAR translator object
    """
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s\t%(levelname)s\t%(message)s')

    mapFile = scriptPath + '/lib/NEF_NMRSTAR_equivalence.csv'
    NEFinfo = scriptPath + '/lib/NEF_mandatory.csv'
    NMRSTARinfo = scriptPath + '/lib/NMR-STAR_mandatory.csv'
    atmFile = scriptPath + '/lib/atomDict.json'
    codeFile = scriptPath + '/lib/codeDict.json'

    def __init__(self):
        ch = logging.StreamHandler()
        ch.setFormatter(self.formatter)
        self.logger.addHandler(ch)
        (isOk, msg, self.tagMap) = self.load_csv_data(self.mapFile, transpose=True)
        if not isOk:
            self.logger.error(msg)
        (isOk, msg, self.atomDict) = self.load_json_data(self.atmFile)
        if not isOk:
            self.logger.error(msg)
        (isOk, msg, self.codeDict) = self.load_json_data(self.codeFile)
        if not isOk:
            self.logger.error(msg)
        (isOk, msg, self.NEFinfo) = self.load_csv_data(self.NEFinfo)
        if not isOk:
            self.logger.error(msg)
        ch.flush()

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
        return is_ok, msg, in_data

    @staticmethod
    def load_json_data(json_file):
        """
        Loads json data files from lib folder
        :param json_file: json file
        :return: dictionay
        """
        try:
            with open(json_file, 'r') as jsonF:
                data_dict = json.loads(jsonF.read())
            is_ok = True
            msg = "{} file is read!".format(json_file)
        except IOError:
            msg = "{} file is missing!".format(json_file)
            is_ok = False
            data_dict = []
        return is_ok, msg, data_dict

    @staticmethod
    def load_csv_data(csv_file, transpose=False):
        """
        Loads csv data files from lib
        :param csv_file: csv file
        :param transpose: transpose multidimensional csv lists
        :return: list
        """
        try:
            with open(csv_file, 'r') as csv_file:
                spam_reader = csv.reader(csv_file, delimiter=',')
                csv_dat = []
                for r in spam_reader:
                    # print r
                    if r[0][0] != '#':
                        csv_dat.append(r)
            if transpose:
                csv_map = list(map(list, zip(*csv_dat)))
            else:
                csv_map = csv_dat
            msg = "{} file is read!".format(csv_file)
            is_ok = True
        except IOError:
            msg = "{} file is missing!".format(csv_file)
            csv_map = []
            is_ok = False
        return is_ok, msg, csv_map

    def get_one_letter_code(self, res):
        """
        Returns one leter amino acid code
        :param res: Three letter code
        :return: One leter code
        """
        try:
            ol = self.codeDict[res.upper()]
        except KeyError:
            ol = '?'
        return ol

    @staticmethod
    def time_stamp(ts):
        """
        Returns time stamp in human readable format for logging
        :param ts: current system time from time.time()
        :return: returns '%Y-%m-%d %H:%M:%S'
        """
        return datetime.datetime.fromtimestamp(ts, tz=utc).strftime('%Y-%m-%d %H:%M:%S')

    def check_mandatory_tags(self,in_file= None ,file_type = None):
        """
        Returns list of missing mandatory tags for the input file. List of mandaotry tags for NEF file is defined
        :param in_file: NEF / NMR-STAR file
        :return: list of missing mandaotry tags
        """
        if file_type is None or in_file is None:
            print ('Please specify the file and file type (NEF or NMR-STAR')
            exit(1)
        else:
            ent = pynmrstar.Entry.from_file(in_file)
            missing_tags = []
            if file_type == "NEF":
                mand_tag = self.NEFinfo
            else:
                mand_tag = self.NMRSTARinfo
            for tag in mand_tag:
                if tag[0][0]=="_" and tag[1] == 'yes':
                    try:
                        ent.get_tags([tag[0]])
                    except ValueError:
                        missing_tags.append(tag[0])
        return missing_tags



    def validate_file(self, in_file, file_type='A'):
        """
        Validates input NEF/NMR-STAR file.
        file_type flags can be 'A' or 'S' or 'R'.
        A for  All in one file,
        S for chemical Shifts file,
        R for 'Restraints file
        """
        info = []
        warning = []
        error = []
        file_format = 'UNKNOWN'
        try:
            file_info = self.read_input_file(in_file)
            if file_info[0]:
                in_data = file_info[2]
                minimal_info_nef_a = ['_nef_chemical_shift', '_nef_distance_restraint']
                minimal_info_nef_s = ['_nef_chemical_shift']
                minimal_info_nef_r = ['_nef_distance_restraint']
                minimal_info_nmrstar_a = ['_Atom_chem_shift', '_Gen_dist_constraint']
                minimal_info_nmrstar_s = ['_Atom_chem_shift']
                minimal_info_nmrstar_r = ['_Gen_dist_constraint']
                (sf_list, lp_list) = self.get_data_content(in_data, file_info[1])
                msg = "{} saveframes and {} loops found".format(len(sf_list), len(lp_list))
                info.append(msg)
                nef_sf_list = [i for i in sf_list if 'nef' in i]
                nef_lp_list = [i for i in lp_list if 'nef' in i]
                msg = "{} saveframes and {} loops found with nef prefix".format(len(nef_sf_list), len(nef_lp_list))
                info.append(msg)
                if len(nef_sf_list) > 0 or len(nef_lp_list) > 0:
                    is_nef_file = True
                    msg = "{} is a NEF file".format(in_file)
                    file_format = 'NEF'
                else:
                    is_nef_file = False
                    msg = "{} is a NMR-STAR file".format(in_file)
                    file_format = "NMR-STAR"
                info.append(msg)
                if is_nef_file:
                    is_valid = True
                    if file_type == "A":
                        for lp_info in minimal_info_nef_a:
                            if lp_info not in lp_list:
                                is_valid = False
                                error.append('{} loop not found'.format(lp_info))
                            else:
                                if self.is_empty_loop(in_data, lp_info, file_info[1]):
                                    is_valid = False
                                    error.append('{} loop is empty'.format(lp_info))
                    elif file_type == "S":
                        for lp_info in minimal_info_nef_s:
                            if lp_info not in lp_list:
                                is_valid = False
                                error.append('{} loop not found'.format(lp_info))
                            else:
                                if self.is_empty_loop(in_data, lp_info, file_info[1]):
                                    is_valid = False
                                    error.append('{} loop is empty'.format(lp_info))
                    elif file_type == "R":
                        for lp_info in minimal_info_nef_r:
                            if lp_info not in lp_list:
                                is_valid = False
                                error.append('{} loop not found'.format(lp_info))
                            else:
                                if self.is_empty_loop(in_data, lp_info, file_info[1]):
                                    is_valid = False
                                    error.append('{} loop is empty'.format(lp_info))
                    else:
                        msg = "file_type flag should be A/S/R"
                        error.append(msg)
                        is_valid = False
                else:
                    is_valid = True
                    if file_type == "A":
                        for lp_info in minimal_info_nmrstar_a:
                            if lp_info not in lp_list:
                                is_valid = False
                                error.append('{} loop not found'.format(lp_info))
                            else:
                                if self.is_empty_loop(in_data, lp_info, file_info[1]):
                                    is_valid = False
                                    error.append('{} loop is empty'.format(lp_info))
                    elif file_type == "S":
                        for lp_info in minimal_info_nmrstar_s:
                            if lp_info not in lp_list:
                                is_valid = False
                                error.append('{} loop not found'.format(lp_info))
                            else:
                                if self.is_empty_loop(in_data, lp_info, file_info[1]):
                                    is_valid = False
                                    error.append('{} loop is empty'.format(lp_info))
                    elif file_type == "R":
                        for lp_info in minimal_info_nmrstar_r:
                            if lp_info not in lp_list:
                                is_valid = False
                                error.append('{} loop not found'.format(lp_info))
                            else:
                                if self.is_empty_loop(in_data, lp_info, file_info[1]):
                                    is_valid = False
                                    error.append('{} loop is empty'.format(lp_info))
                    else:
                        msg = "file_type flag should be A/S/R"
                        error.append(msg)
                        is_valid = False

                # info.append(is_nef_file)
            else:
                error.append(file_info[1])
                is_valid = False
        except IOError:
            msg = "File not found {}".format(in_file)
            error.append(msg)
            is_valid = False

        return is_valid, json.dumps({'info': info, 'warning': warning, 'error': error, 'FILE': file_format})

    @staticmethod
    def is_empty_loop(star_data, lp_category, data_flag):
        """
        Check if a given loop is empty
        """
        if data_flag == "Entry":
            loops = star_data.get_loops_by_category(lp_category)
            for loop in loops:
                if len(loop.data) == 0:
                    return True
            return False
        elif data_flag == "Saveframe":
            loop = star_data.get_loop_by_category(lp_category)
            return len(loop.data) == 0
        else:
            return len(star_data.data) == 0

    @staticmethod
    def is_empty_data(data):
        """
        Check if given data has empty code
        """
        for d in data:
            if d in (None, '', '.', '?'):
                return True
        return False

    @staticmethod
    def is_data(data):
        """
        Check if given data has no empty code
        """
        for d in data:
            if d in (None, '', '.', '?'):
                return False
        return True

    @staticmethod
    def get_data_content(star_data, data_flag):
        sf_list = []
        lp_list = []
        if data_flag == "Entry":
            for sf in star_data.frame_list:
                sf_list.append(sf.category)
                for lp in sf:
                    lp_list.append(lp.category)
        elif data_flag == "Saveframe":
            for lp in star_data:
                lp_list.append(lp.category)
        else:
            lp_list.append(star_data.category)
        return sf_list, lp_list

    def get_seq_from_cs_loop(self, in_file):
        """
        Extracts sequence from checmial shift loop
        :param in_file: NEF/NMR-STAR file
        :return: status flag,json data
        """
        (flg, json_data) = self.validate_file(in_file, 'S')
        dat = json.loads(json_data)
        info = dat['info']
        warning = dat['warning']
        error = dat['error']
        is_ok = False
        seq = []
        if flg:
            info.append('File successfully read')
            in_dat = self.read_input_file(in_file)[-1]
            if dat['FILE'] == "NMR-STAR":
                info.append('NMR-STAR')
                seq = self.get_nmrstar_seq(in_dat)[0]
                if len(seq):
                    is_ok = True
                else:
                    error.append("Can't extract sequence from chemical shift loop")
            elif dat['FILE'] == "NEF":
                info.append('NEF')
                seq = self.get_nef_seq(in_dat)[0]
                if len(seq):
                    is_ok = True
                else:
                    error.append("Can't extract sequence from chemical shift loop")
            else:
                error.append("Can't identify file type, it is neither NEF nor NMR-STAR")
        else:
            error.append('File validation failed (or) File contains no chemical shift information')
        return is_ok, json.dumps({'info': info, 'warning': warning, 'error': error, 'FILE': dat['FILE'], 'DATA': seq})

    @staticmethod
    def get_nef_seq(str_data, lp_category='nef_chemical_shift', seq_id='sequence_code', res_id='residue_name',
                    chain_id='chain_code', allow_empty=False):
        """Extracts sequence from any given loop from a NEF file"""
        try:
            loops = str_data.get_loops_by_category(lp_category)
        except AttributeError:
            try:
                loops = [str_data.get_loop_by_category(lp_category)]
            except AttributeError:
                loops = [str_data]

        tags = [seq_id, res_id, chain_id]

        seq = []
        sid = []

        for loop in loops:
            seq_dict = {}
            sid_dict = {}

            seq_dat = []

            if set(tags) & set(loop.tags) == set(tags):
                seq_dat = loop.get_data_by_tag(tags)
            else:
                _tags_exist = False
                for i in range(1, 5): # expand up to 4 dimensions
                    _tags = [seq_id + '_' + str(i), res_id + '_' + str(i), chain_id + '_' + str(i)]
                    if set(_tags) & set(loop.tags) == set(_tags):
                        _tags_exist = True
                        seq_dat += loop.get_data_by_tag(_tags)
                    else:
                        break
                if not _tags_exist:
                    seq_dat = loop.get_data_by_tag(tags) # raise ValueError

            if allow_empty:
                seq_dat = list(filter(NEFTranslator.is_data, seq_dat))
            else:
                for i in seq_dat:
                    if NEFTranslator.is_empty_data(i):
                        raise ValueError('Sequence must not be empty. loop_category=%s, chain_id=%s, seq_id=%s, res_id=%s' % (lp_category, i[2], i[0], i[1]))

            chains = (set([i[2] for i in seq_dat]))
            seq_concat = (sorted(set(['{}-{:04d}-{}'.format(i[2], int(i[0]), i[1]) for i in seq_dat])))

            chk_dict = {'{}-{:04d}'.format(i[2], int(i[0])):i[1] for i in seq_dat}

            for i in seq_dat:
                chk_key = '{}-{:04d}'.format(i[2], int(i[0]))
                if chk_dict[chk_key] != i[1]:
                    raise KeyError('Sequence must be unique. loop_category=%s, chain_id=%s, seq_id=%s, res_id=%s vs %s' % (lp_category, i[2], i[0], i[1], chk_dict[chk_key]))

            if len(seq_concat[0].split("-")[-1]) > 1:
                if len(chains) > 1:
                    for c in chains:
                        # seq_array = "".join([self.getOneLetter(i.split("-")[-1]) for i in seq_concat if i.split("-")[0] == c])
                        seq_array = [i.split("-")[-1] for i in seq_concat if i.split("-")[0] == c]
                        sid_array = [int(i.split("-")[1]) for i in seq_concat if i.split("-")[0] == c]
                        seq_dict[c] = seq_array
                        sid_dict[c] = sid_array
                else:
                    # seq_array = "".join([self.getOneLetter(i.split("-")[-1]) for i in seq_concat])
                    seq_array = [i.split("-")[-1] for i in seq_concat]
                    sid_array = [int(i.split("-")[1]) for i in seq_concat]
                    seq_dict[list(chains)[0]] = seq_array
                    sid_dict[list(chains)[0]] = sid_array
            else:
                if len(chains) > 1:
                    for c in chains:
                        # seq_array = "".join([i.split("-")[-1] for i in seq_concat if i.split("-")[0] == c])
                        seq_array = [i.split("-")[-1] for i in seq_concat if i.split("-")[0] == c]
                        sid_array = [int(i.split("-")[1]) for i in seq_concat if i.split("-")[0] == c]
                        seq_dict[c] = seq_array
                        sid_dict[c] = sid_array
                else:
                    # seq_array = "".join([i.split("-")[-1] for i in seq_concat])
                    seq_array = [i.split("-")[-1] for i in seq_concat]
                    sid_array = [int(i.split("-")[1]) for i in seq_concat]
                    seq_dict[list(chains)[0]] = seq_array
                    sid_dict[list(chains)[0]] = sid_array
            seq.append(seq_dict)
            sid.append(sid_dict)
        return seq, sid

    @staticmethod
    def get_nmrstar_seq(str_data, lp_category='Atom_chem_shift', seq_id='Comp_index_ID', res_id='Comp_ID',
                        chain_id='Entity_assembly_ID', allow_empty=False):
        """Extracts sequence from any given NMR-STAR file"""
        try:
            loops = str_data.get_loops_by_category(lp_category)
        except AttributeError:
            try:
                loops = [str_data.get_loop_by_category(lp_category)]
            except AttributeError:
                loops = [str_data]

        seq = []
        sid = []

        tags = [seq_id, res_id, chain_id]
        tags_ = [seq_id, res_id]

        for loop in loops:
            seq_dict = {}
            sid_dict = {}

            seq_dat = []

            if set(tags) & set(loop.tags) == set(tags):
                seq_dat = loop.get_data_by_tag(tags)
            elif set(tags_) & set(loop.tags) == set(tags_): # No Entity_assembly_ID tag case
                seq_dat = loop.get_data_by_tag(tags_)
                for i in seq_dat:
                    i.append('1')
            else:
                _tags_exist = False
                for i in range(1, 5): # expand up to 4 dimensions
                    _tags = [seq_id + '_' + str(i), res_id + '_' + str(i), chain_id + '_' + str(i)]
                    _tags_ = [seq_id + '_' + str(i), res_id + '_' + str(i)]
                    if set(_tags) & set(loop.tags) == set(_tags):
                        _tags_exist = True
                        seq_dat += loop.get_data_by_tag(_tags)
                    elif set(_tags_) & set(loop.tags) == set(_tags_):
                        _tags_exist = True
                        seq_dat_ = loop.get_data_by_tag(_tags_)
                        for i in seq_dat_:
                            i.append('1')
                        seq_dat += seq_dat_
                    else:
                        break
                if not _tags_exist:
                    seq_dat = loop.get_data_by_tag(tags) # raise ValueError

            if allow_empty:
                seq_dat = list(filter(NEFTranslator.is_data, seq_dat))
            else:
                for i in seq_dat:
                    if NEFTranslator.is_empty_data(i):
                        raise ValueError('Sequence must not be empty. loop_category=%s, chain_id=%s, seq_id=%s, res_id=%s' % (lp_category, i[2], i[0], i[1]))

            chains = (set([i[2] for i in seq_dat]))
            seq_concat = (sorted(set(['{}-{:04d}-{}'.format(i[2], int(i[0]), i[1]) for i in seq_dat])))

            chk_dict = {'{}-{:04d}'.format(i[2], int(i[0])):i[1] for i in seq_dat}

            for i in seq_dat:
                chk_key = '{}-{:04d}'.format(i[2], int(i[0]))
                if chk_dict[chk_key] != i[1]:
                    raise KeyError('Sequence must be unique. loop_category=%s, chain_id=%s, seq_id=%s, res_id=%s vs %s' % (lp_category, i[2], i[0], i[1], chk_dict[chk_key]))

            if len(seq_concat[0].split("-")[-1]) > 1:
                if len(chains) > 1:
                    for c in chains:
                        # seq_array = "".join([self.getOneLetter(i.split("-")[-1]) for i in seq_concat if i.split("-")[0] == c])
                        seq_array = [i.split("-")[-1] for i in seq_concat if i.split("-")[0] == c]
                        sid_array = [int(i.split("-")[1]) for i in seq_concat if i.split("-")[0] == c]
                        seq_dict[c] = seq_array
                        sid_dict[c] = sid_array
                else:
                    # seq_array = "".join([self.getOneLetter(i.split("-")[-1]) for i in seq_concat])
                    seq_array = [i.split("-")[-1] for i in seq_concat]
                    sid_array = [int(i.split("-")[1]) for i in seq_concat]
                    seq_dict[list(chains)[0]] = seq_array
                    sid_dict[list(chains)[0]] = sid_array
            else:
                if len(chains) > 1:
                    for c in chains:
                        # seq_array = "".join([i.split("-")[-1] for i in seq_concat if i.split("-")[0] == c])
                        seq_array = [i.split("-")[-1] for i in seq_concat if i.split("-")[0] == c]
                        sid_array = [int(i.split("-")[1]) for i in seq_concat if i.split("-")[0] == c]
                        seq_dict[c] = seq_array
                        seeq_id_dict[c] = sid_array
                else:
                    # seq_array = "".join([i.split("-")[-1] for i in seq_concat])
                    seq_array = [i.split("-")[-1] for i in seq_concat]
                    sid_array = [int(i.split("-")[1]) for i in seq_concat]
                    seq_dict[list(chains)[0]] = seq_array
                    sid_dict[list(chains)[0]] = sid_array
            seq.append(seq_dict)
            sid.append(sid_dict)
        return seq, sid

    def validate_atom(self, star_data, lp_category='Atom_chem_shift', seq_id='Comp_index_ID', res_id='Comp_ID',
                      atom_id='Atom_ID'):
        """
        Validates the atoms in a given loop against IUPAC standard

        """
        try:
            loop_data = star_data.get_loops_by_category(lp_category)
        except AttributeError:
            try:
                loop_data = [star_data.get_loop_by_category(lp_category)]
            except AttributeError:
                loop_data = [star_data]

        ns = []
        for lp in loop_data:
            try:
                atm_data = lp.get_data_by_tag([seq_id, res_id, atom_id])
                for i in atm_data:
                    try:
                        if i[2] not in self.atomDict[i[1].upper()]:
                            ns.append(i)
                    except KeyError:
                        ns.append(i)
            except ValueError:
                print("One of the following tag is missing ", seq_id, res_id, atom_id)

            # nonStandard = [i for i in atm_data if i[2] not in self.atomDict[i[1].upper]]
            # ns.append(nonStandard)
        return ns

    def get_nmrstar_tag(self, tag):
        try:
            n = self.tagMap[0].index(tag)
        except ValueError:
            return []
        return [self.tagMap[1][n], self.tagMap[2][n]]
    def get_nef_tag(self, tag):
        try:
            n = self.tagMap[1].index(tag)
            s_tag = self.tagMap[0][n]
        except KeyError:
            try:
                n = self.tagMap[2].index(tag)
                s_tag = self.tagMap[0][n]
            except KeyError:
                s_tag = []
                print ("Matching NEF tag not found for {}".format(tag))
        return [s_tag]

    def get_nef_loop_tags(self, star_loop_tags):
        aut_tag = []
        nt = []
        for t in star_loop_tags:
            try:
                st = self.get_nef_tag(t)
                if len(st)==1 and st[0] not in nt:
                    nt.append(st[0])
            except ValueError:
                print ("tag not found {}".format(t))
        out_tag = nt
        # if star_loop_tags[0].split(".")[0] == "_nef_chemical_shift":
        #     out_tag.append('_Atom_chem_shift.ID')
        #     out_tag.append('_Atom_chem_shift.Ambiguity_code')
        #     out_tag.append('_Atom_chem_shift.Ambiguity_set_ID')
        #     out_tag.append('_Atom_chem_shift.Assigned_chem_shift_list_ID')
        # if star_loop_tags[0].split(".")[0] == "_nef_distance_restraint":
        #     out_tag.append('_Gen_dist_constraint.Member_logic_code')
        #     out_tag.append('_Gen_dist_constraint.Gen_dist_constraint_list_ID')
        # if star_loop_tags[0].split(".")[0] == "_nef_dihedral_restraint":
        #     out_tag.append('_Torsion_angle_constraint.Torsion_angle_constraint_list_ID')
        # if star_loop_tags[0].split(".")[0] == "_nef_rdc_restraint":
        #     out_tag.append('_RDC_constraint.RDC_constraint_list_ID')
        # if star_loop_tags[0].split(".")[0] == "_nef_peak":
        #     out_tag.append('_Peak_row_format.Spectral_peak_list_ID')
        # if star_loop_tags[0].split(".")[0] == "_nef_sequence":
        #     out_tag.append('_Chem_comp_assembly.Assembly_ID')
        return out_tag

    def get_nmrstar_loop_tags(self, nef_loop_tags):
        aut_tag = []
        nt = []
        for t in nef_loop_tags:
            st = self.get_nmrstar_tag(t)
            if st[0] != st[1]:
                aut_tag.append(st[1])
            nt.append(st[0])
        if len(aut_tag) != 0:
            out_tag = nt + aut_tag
        else:
            out_tag = nt
        if nef_loop_tags[0].split(".")[0] == "_nef_chemical_shift":
            out_tag.append('_Atom_chem_shift.ID')
            out_tag.append('_Atom_chem_shift.Ambiguity_code')
            out_tag.append('_Atom_chem_shift.Ambiguity_set_ID')
            out_tag.append('_Atom_chem_shift.Assigned_chem_shift_list_ID')
        if nef_loop_tags[0].split(".")[0] == "_nef_distance_restraint":
            out_tag.append('_Gen_dist_constraint.Member_logic_code')
            out_tag.append('_Gen_dist_constraint.Gen_dist_constraint_list_ID')
        if nef_loop_tags[0].split(".")[0] == "_nef_dihedral_restraint":
            out_tag.append('_Torsion_angle_constraint.Torsion_angle_constraint_list_ID')
        if nef_loop_tags[0].split(".")[0] == "_nef_rdc_restraint":
            out_tag.append('_RDC_constraint.RDC_constraint_list_ID')
        if nef_loop_tags[0].split(".")[0] == "_nef_peak":
            out_tag.append('_Peak_row_format.Spectral_peak_list_ID')
        if nef_loop_tags[0].split(".")[0] == "_nef_sequence":
            out_tag.append('_Chem_comp_assembly.Assembly_ID')
        return out_tag

    def get_nmrstar_atom(self, res, nef_atom):
        """
        Returns (atom with out wildcard,[IUPAC atom list],ambiguity code)

        """
        ambiguity_code = 1
        atom_type = None
        try:
            atoms = self.atomDict[res]
            atom_list = []
            try:
                ref_atom = re.findall(r'(\S+)([xyXY])([%*])$|(\S+)([%*])$|(\S+)([xyXY]$)', nef_atom)[0]
                atm_set = [ref_atom.index(i) for i in ref_atom if i != ""]
                pattern = None
                if atm_set == [0, 1, 2]:
                    atom_type = ref_atom[0]
                    pattern = re.compile(r'%s\S\d+' % (ref_atom[0]))
                    alist2 = [i for i in atoms if re.search(pattern, i)]
                    xid = sorted(set([int(i[len(ref_atom[0])]) for i in alist2]))
                    if ref_atom[1] == "x" or ref_atom[1] == "X":
                        atom_list = [i for i in alist2 if int(i[len(ref_atom[0])]) == xid[0]]
                    else:
                        atom_list = [i for i in alist2 if int(i[len(ref_atom[0])]) == xid[1]]
                    ambiguity_code = 2
                elif atm_set == [3, 4]:
                    atom_type = ref_atom[3]
                    if ref_atom[4] == "%":
                        pattern = re.compile(r'%s\d+' % (ref_atom[3]))
                    elif ref_atom[4] == "*":
                        pattern = re.compile(r'%s\S+' % (ref_atom[3]))
                    else:
                        logging.critical("Wrong NEF atom {}".format(nef_atom))
                    atom_list = [i for i in atoms if re.search(pattern, i)]
                    ambiguity_code = 1

                elif atm_set == [5, 6]:
                    atom_type = ref_atom[5]
                    pattern = re.compile(r'%s\S+' % (ref_atom[5]))
                    atom_list = [i for i in atoms if re.search(pattern, i)]
                    if len(atom_list) != 2:
                        atom_list = []
                    elif ref_atom[6] == "y" or ref_atom[6] == "Y":
                        # atom_list.reverse()[]
                        atom_list = atom_list[-1:]
                    elif ref_atom[6] == "x" or ref_atom[6] == "X":
                        atom_list = atom_list[:1]
                    else:
                        logging.critical("Wrong NEF atom {}".format(nef_atom))
                    ambiguity_code = 2

                else:
                    logging.critical("Wrong NEF atom {}".format(nef_atom))
            except IndexError:

                # print nef_atom
                pass
                atom_type = nef_atom
            if len(atom_list) == 0:
                if nef_atom in atoms:
                    atom_list.append(nef_atom)
                else:
                    if nef_atom == "H%":  # To handle terminal protons
                        atom_list = ['H1', 'H2', 'H3']
                        atom_type = "H"
        except KeyError:
            # self.logfile.write("%s\tResidue not found,%s,%s\n"%(self.TimeStamp(time.time()),res,nef_atom))
            # print "Residue not found",res,nef_atom
            if res != ".":
                self.logger.critical("Non-standard residue found {}".format(res))
            atom_list = []
            atom_type = nef_atom

            if nef_atom == "H%":
                atom_list = ['H1', 'H2', 'H3']
                atom_type = "H"
        return atom_type, atom_list, ambiguity_code

    def translate_cs_row(self, f_tags, t_tags, row_data, iupac = True):
        """
        Translates row of data in chemical shift loop from NEF into NMR-STAR
        :param f_tags: nef tags
        :type f_tags: list
        :param t_tags: List nmr-star tags
        :type t_tags: list
        :param row_data: List nef data
        :type row_data: list
        :param iupac : convert atom nomenclature into IUPAC
        :type iupac: bool
        :return: List NMR-STAR data
        """
        out_row = []
        new_id = None
        if '_nef_chemical_shift.chain_code' in f_tags and '_nef_chemical_shift.sequence_code' in f_tags:

            cci = f_tags.index('_nef_chemical_shift.chain_code')
            sci = f_tags.index('_nef_chemical_shift.sequence_code')
            try:
                old_id = [i for i in self.seqDict.keys() if i[0] == row_data[cci] and i[1] == row_data[sci]][0]
                new_id = self.seqDict[old_id]
            except AttributeError:
                new_id = (cci, sci)
            except IndexError:
                new_id = (cci, sci)
        if len(f_tags) != len(t_tags):
            atm_index = f_tags.index('_nef_chemical_shift.atom_name')
            res_index = f_tags.index('_nef_chemical_shift.residue_name')
            n_atm = self.get_nmrstar_atom(row_data[res_index], row_data[atm_index])[1]
            ambi = self.get_nmrstar_atom(row_data[res_index], row_data[atm_index])[2]
            if iupac == False:
                n_atm = [row_data[atm_index]]
            for i in n_atm:
                out = [None] * len(t_tags)
                for j in f_tags:
                    stgs = self.get_nmrstar_tag(j)
                    if stgs[0] == stgs[1]:
                        out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                    else:
                        if j == '_nef_chemical_shift.atom_name':
                            out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                            out[t_tags.index(stgs[1])] = i
                        elif j == '_nef_chemical_shift.chain_code':
                            out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                            out[t_tags.index(stgs[1])] = new_id[0]
                        elif j == '_nef_chemical_shift.sequence_code':
                            out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                            out[t_tags.index(stgs[1])] = new_id[1]
                        else:
                            out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                            out[t_tags.index(stgs[1])] = row_data[f_tags.index(j)]
                    out[t_tags.index('_Atom_chem_shift.Ambiguity_code')] = ambi
                    out[t_tags.index('_Atom_chem_shift.Ambiguity_set_ID')] = '.'
                out_row.append(out)
        else:
            out_row.append(row_data)
        return out_row

    @staticmethod
    def get_residue_identifier(tag_list):
        out_list = []
        for j in range(1, 16):
            out = [None] * 2
            chk_string = re.compile(r'\S+.chain_code_{}'.format(j))
            r1 = [chk_string.search(i).group() for i in tag_list if chk_string.search(i)]
            if len(r1) > 0:
                out[0] = r1[0]
            chk_string = re.compile(r'\S+.sequence_code_{}'.format(j))
            r2 = [chk_string.search(i).group() for i in tag_list if chk_string.search(i)]
            if len(r2) > 0:
                out[1] = r2[0]
            if len(r1) > 0 and len(r2) > 0:
                out_list.append(out)
        #             chk_string = re.compile('\S+.residue_name_{}'.format(j))
        #             r=[chk_string.search(i).group() for i in tag_list if chk_string.search(i)]
        #             if len(r)>0: out[2]=r[0]
        #             chk_string = re.compile('\S+.atom_name_{}'.format(j))
        #             r=[chk_string.search(i).group() for i in tag_list if chk_string.search(i)]
        #             if len(r)>0: out[3]=r[0]
        return out_list

    def translate_row(self, f_tags, t_tags, row_data):
        """
        Translates row of data in a loop from NEF into NMR-STAR
        :param f_tags: NEF tags
        :type f_tags: list
        :param t_tags: NMR-STAR tags
        :type t_tags: list
        :param row_data: NEF data
        :type row_data: list
        :return:
        """
        # print (f_tags)
        out_row = []
        res_list = self.get_residue_identifier(f_tags)
        tmp_dict = {}
        for res1 in res_list:
            try:
                tmp_dict[res1[0]] = self.seqDict[(row_data[f_tags.index(res1[0])], row_data[f_tags.index(res1[1])])][0]
            except KeyError:
                tmp_dict[res1[0]] = row_data[f_tags.index(res1[0])]
            try:
                tmp_dict[res1[1]] = self.seqDict[(row_data[f_tags.index(res1[0])], row_data[f_tags.index(res1[1])])][1]
            except KeyError:
                tmp_dict[res1[1]] = row_data[f_tags.index(res1[1])]
        # print (tmp_dict)
        if len(f_tags) != len(t_tags):
            out = [None] * len(t_tags)
            for j in f_tags:
                stgs = self.get_nmrstar_tag(j)
                if stgs[0] == stgs[1]:
                    out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                else:
                    if 'chain_code' in j or 'sequence_code' in j:
                        out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                        out[t_tags.index(stgs[1])] = tmp_dict[j]
                    else:
                        out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                        out[t_tags.index(stgs[1])] = row_data[f_tags.index(j)]

                # else:
                #   print ("ERROR",f_tags)
            out_row.append(out)
        else:
            out_row.append(row_data)
        return out_row

    def translate_seq_row(self, f_tags, t_tags, row_data):
        """
        Translates row of data in sequence  loop from NEF into NMR-STAR
        :param f_tags: NEF tags
        :type f_tags: list
        :param t_tags: NMR-STAR tags
        :type t_tags: list
        :param row_data: NEF data
        :type row_data: list
        :return:
        """

        out_row = []
        if len(f_tags) != len(t_tags):
            out = [None] * len(t_tags)
            for j in f_tags:
                stgs = self.get_nmrstar_tag(j)
                if stgs[0] == stgs[1]:
                    out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                else:
                    if j == '_nef_sequence.chain_code':
                        out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                        out[t_tags.index(stgs[1])] = self.chains.index(row_data[f_tags.index(j)]) + 1
                    elif j == '_nef_sequence.sequence_code':
                        out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                        out[t_tags.index(stgs[1])] = self.cid[
                            self.chains.index(row_data[f_tags.index('_nef_sequence.chain_code')])]
                    else:
                        out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                        out[t_tags.index(stgs[1])] = row_data[f_tags.index(j)]
            out_row.append(out)
        else:
            out_row.append(row_data)
        return out_row

    def translate_seq_row2(self, f_tags, t_tags, row_data):
        """
        Translates row of data in sequence  loop from NEF into NMR-STAR
        :param f_tags: NMR-STAR tags
        :type f_tags: list
        :param t_tags: NEF tags
        :type t_tags: list
        :param row_data: NMR-STAR data
        :type row_data: list
        :return:
        """

        out_row = []
        if len(f_tags) != len(t_tags):
            out = [None] * len(t_tags)
            for j in f_tags:
                ntag = self.get_nef_tag(j)
                stgs = self.get_nmrstar_tag(ntag)
                if stgs[0] == stgs[1]:
                    out[t_tags.index(ntag)] = row_data[f_tags.index(stgs[0])]
                else:
                    out[t_tags.index(ntag)] = row_data[f_tags.index(stgs[1])]
            out_row.append(out)
        else:
            out_row.append(row_data)
        return out_row

    def translate_row2(self, f_tags, t_tags, row_data):
        """
        Translates row of data in sequence  loop from NEF into NMR-STAR
        :param f_tags: NMR-STAR tags
        :type f_tags: list
        :param t_tags: NEF tags
        :type t_tags: list
        :param row_data: NMR-STAR data
        :type row_data: list
        :return:
        """

        out_row = []
        if len(f_tags) != len(t_tags):
            out = [None] * len(t_tags)
            for j in f_tags:
                try:
                    ntag = self.get_nef_tag(j)[0]
                    stgs = self.get_nmrstar_tag(ntag)
                    if stgs[0] == stgs[1]:
                        out[t_tags.index(ntag)] = row_data[f_tags.index(stgs[0])]
                    else:
                        out[t_tags.index(ntag)] = row_data[f_tags.index(stgs[0])]
                except ValueError:
                   # print ("tag not found {}".format(j))
                    pass
            out_row.append(out)
        else:
            out_row.append(row_data)
        return out_row

    def translate_restraint_row(self, f_tags, t_tags, row_data, iupac = True):
        """
        Translates row of data in restraint loop from NEF into NMR-STAR
        :param f_tags: NEF tags
        :type f_tags: list
        :param t_tags: NMR-STAR tags
        :type t_tags: list
        :param row_data: NEF data
        :type row_data: list
        :param iupac: Convert atom nomenclature into IUPAC
        :type iupac: bool
        :return:
        """
        out_row = []
        res_list = self.get_residue_identifier(f_tags)
        # print (res_list,f_tags)
        tmp_dict = {}
        for res1 in res_list:
            try:
                tmp_dict[res1[0]] = self.seqDict[(row_data[f_tags.index(res1[0])], row_data[f_tags.index(res1[1])])][0]
            except KeyError:
                tmp_dict[res1[0]] = row_data[f_tags.index(res1[0])]
            try:
                tmp_dict[res1[1]] = self.seqDict[(row_data[f_tags.index(res1[0])], row_data[f_tags.index(res1[1])])][1]
            except KeyError:
                tmp_dict[res1[1]] = row_data[f_tags.index(res1[1])]
        if len(f_tags) != len(t_tags):
            atm_index1 = f_tags.index('_nef_distance_restraint.atom_name_1')
            res_index1 = f_tags.index('_nef_distance_restraint.residue_name_1')
            atm_index2 = f_tags.index('_nef_distance_restraint.atom_name_2')
            res_index2 = f_tags.index('_nef_distance_restraint.residue_name_2')
            n_atm1 = self.get_nmrstar_atom(row_data[res_index1], row_data[atm_index1])[1]
            n_atm2 = self.get_nmrstar_atom(row_data[res_index2], row_data[atm_index2])[1]

            if iupac == False:
                n_atm1 = [row_data[atm_index1]]
                n_atm2 = [row_data[atm_index2]]

            for i in n_atm1:
                for k in n_atm2:
                    out = [None] * len(t_tags)
                    for j in f_tags:
                        stgs = self.get_nmrstar_tag(j)
                        if stgs[0] == stgs[1]:
                            out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                        else:
                            if j == '_nef_distance_restraint.atom_name_1':
                                out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                                out[t_tags.index(stgs[1])] = i
                            elif 'chain_code_1' in j or 'sequence_code_1' in j:
                                out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                                out[t_tags.index(stgs[1])] = tmp_dict[j]
                            elif j == '_nef_distance_restraint.atom_name_2':
                                out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                                out[t_tags.index(stgs[1])] = k
                            elif 'chain_code_2' in j or 'sequence_code_2' in j:
                                out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                                out[t_tags.index(stgs[1])] = tmp_dict[j]
                            else:
                                out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                                out[t_tags.index(stgs[1])] = row_data[f_tags.index(j)]

                    out_row.append(out)
        else:
            out_row.append(row_data)

        return out_row

    def nef_to_nmrstar(self, nef_file, star_file=None, iupac = True):
        """
        Translates NEF into NMR-STAR
        :param nef_file: input file in NEF format
        :param star_file: output file name
        :param iupac: convert the atom nomenclature to IUPAC True/False default: True
        :return:
        """
        (file_path, file_name) = ntpath.split(os.path.realpath(nef_file))
        is_done = True
        info = []
        warning = []
        error = []
        if star_file is None:
            if iupac== False:
                star_file = file_path + "/" + file_name.split(".")[0] + "noiupac.str"
            else:
                star_file = file_path + "/" + file_name.split(".")[0] + ".str"
        (is_readable, dat_content, nef_data) = self.read_input_file(nef_file)
        try:
            star_data = pynmrstar.Entry.from_scratch(nef_data.entry_id)
        except AttributeError:
            star_data = pynmrstar.Entry.from_scratch(file_name.split(".")[0])
            warning.append('Not a complete Entry')
        if is_readable:
            if dat_content == "Entry":
                self.chains = sorted(list(set(nef_data.get_loops_by_category('nef_sequence')[0].get_tag('chain_code'))))
            elif dat_content == "Saveframe":
                self.chains = sorted(list(set(nef_data[0].get_tag('chain_code'))))
            elif dat_content == "Loop":
                self.chains = sorted(list(set(nef_data.get_tag('chain_code'))))
            else:
                is_done = False
                error.append('File content unknown')
            assembly_id = 0
            cs_list = 0
            rest_list = 0
            ang_list = 0
            rdc_list = 0
            peak_list = 0
            if dat_content == "Entry":
                for saveframe in nef_data:
                    sf = pynmrstar.Saveframe.from_scratch(saveframe.name)

                    for tag in saveframe.tags:
                        if tag[0].lower() == "sf_category":
                            sf.add_tag("Sf_category", self.get_nmrstar_tag(saveframe.category)[0])
                        else:
                            neftag = '{}.{}'.format(saveframe.tag_prefix, tag[0])
                            ntag = self.get_nmrstar_tag(neftag)
                            if len(ntag)>0:
                                sf.add_tag(ntag[0], tag[1])
                            else:
                                sf.add_tag('Details','{}:{}'.format(neftag,tag[1]))
                    if saveframe.category == "nef_nmr_meta_data":
                        sf.add_tag("NMR_STAR_version", "3.2.0.15")
                        # sf.add_tag("Generated_date", self.TimeStamp(time.time()), update=True)
                        try:
                            lp1 = saveframe.get_loop_by_category('_nef_program_script')
                            lp1.add_data(['NEFTranslator', 'NEFtoNMRSTAR', '.'])
                            # print (lp1.tags)
                        except KeyError:
                            pass  # May be better to add audit loop
                    for loop in saveframe:
                        if loop.category == "_nef_sequence":
                            self.cid = []  # Comp_index_ID list
                            for c in self.chains:  # Comp_index_ID initialized with 1
                                self.cid.append(1)
                            self.seqDict = {}

                        if loop.category == '_nef_distance_restraint':
                            r_index_id = 1
                            rest_list+=1
                        if loop.category == '_nef_dihedral_restraint':
                            ang_list+=1
                        if loop.category == '_nef_rdc_restraint':
                            rdc_list+=1
                        if loop.category == "_nef_chemical_shift":
                            cs_list += 1
                        if loop.category == "_nef_peak":
                            peak_list += 1
                        if loop.category == "_nef_sequence":
                            assembly_id += 1
                        lp = pynmrstar.Loop.from_scratch()
                        lp_cols = self.get_nmrstar_loop_tags(loop.get_tag_names())
                        for t in lp_cols:
                            lp.add_tag(t)
                        # print (loop.category,lp.category,lp.get_tag_names(),loop.get_tag_names())
                        cs_id = 1
                        for dat in loop.data:
                            if loop.category == "_nef_sequence":
                                dd = self.translate_seq_row(loop.get_tag_names(), lp.get_tag_names(), dat)
                                self.cid[
                                    self.chains.index(dat[loop.get_tag_names().index('_nef_sequence.chain_code')])] += 1
                                for d in dd:
                                    d[lp.get_tag_names().index('_Chem_comp_assembly.Assembly_ID')] = assembly_id
                                    lp.add_data(d)
                                    self.seqDict[(dat[loop.get_tag_names().index('_nef_sequence.chain_code')],
                                                  dat[loop.get_tag_names().index('_nef_sequence.sequence_code')])] = (
                                        d[lp.get_tag_names().index('_Chem_comp_assembly.Entity_assembly_ID')],
                                        d[lp.get_tag_names().index('_Chem_comp_assembly.Comp_index_ID')])

                            elif loop.category == "_nef_chemical_shift":
                                dd = self.translate_cs_row(loop.get_tag_names(), lp.get_tag_names(), dat, iupac)
                                for d in dd:
                                    d[lp.get_tag_names().index(
                                        '_Atom_chem_shift.Assigned_chem_shift_list_ID')] = cs_list
                                    d[lp.get_tag_names().index(
                                        '_Atom_chem_shift.ID')] = cs_id
                                    lp.add_data(d)
                                    cs_id+=1
                            elif loop.category == '_nef_distance_restraint':
                                dd = self.translate_restraint_row(loop.get_tag_names(), lp.get_tag_names(), dat, iupac)
                                for d in dd:
                                    d[lp.get_tag_names().index('_Gen_dist_constraint.Index_ID')] = r_index_id
                                    d[lp.get_tag_names().index('_Gen_dist_constraint.Gen_dist_constraint_list_ID')] = rest_list
                                    if len(dd) > 1:
                                        d[lp.get_tag_names().index('_Gen_dist_constraint.Member_logic_code')] = "OR"
                                    lp.add_data(d)
                                    r_index_id += 1
                            elif loop.category == '_nef_dihedral_restraint':
                                dd = self.translate_row(loop.get_tag_names(), lp.get_tag_names(), dat)
                                for d in dd:
                                    d[lp.get_tag_names().index(
                                        '_Torsion_angle_constraint.Torsion_angle_constraint_list_ID')] = ang_list
                                    lp.add_data(d)
                            elif loop.category == '_nef_rdc_restraint':
                                dd = self.translate_row(loop.get_tag_names(), lp.get_tag_names(), dat)
                                for d in dd:
                                    d[lp.get_tag_names().index(
                                        '_RDC_constraint.RDC_constraint_list_ID')] = rdc_list
                                    lp.add_data(d)
                            elif loop.category == '_nef_peak':
                                dd = self.translate_row(loop.get_tag_names(), lp.get_tag_names(), dat)
                                for d in dd:
                                    d[lp.get_tag_names().index(
                                        '_Peak_row_format.Spectral_peak_list_ID')] = peak_list
                                    lp.add_data(d)
                            else:
                                dd = self.translate_row(loop.get_tag_names(), lp.get_tag_names(), dat)
                                for d in dd:
                                    lp.add_data(d)

                        # print (loop.data[0])
                        sf.add_loop(lp)
                    if saveframe.category == "nef_chemical_shift_list":
                        sf.add_tag("ID", cs_list)
                    if saveframe.category == "nef_distance_restraint_list":
                        sf.add_tag("ID", rest_list)
                    if saveframe.category == "nef_dihedral_restraint_list":
                        sf.add_tag("ID", ang_list)
                    if saveframe.category == "nef_rdc_restraint_list":
                        sf.add_tag("ID", rest_list)
                    if saveframe.category == "nef_molecular_system":
                        sf.add_tag("ID", assembly_id)
                    star_data.add_saveframe(sf)
                star_data.normalize()
                with open(star_file, 'w') as wstarfile:
                    wstarfile.write(str(star_data))
            elif dat_content == "Saveframe" or dat_content == "Loop":
                if dat_content == "Saveframe":
                    saveframe = nef_data
                    sf = pynmrstar.Saveframe.from_scratch(saveframe.name)
                    for tag in saveframe.tags:
                        if tag[0].lower() == "sf_category":
                            try:

                                sf.add_tag("Sf_category", self.get_nmrstar_tag(saveframe.category)[0])
                            except ValueError:
                                sf.add_tag("Sf_category", self.get_nmrstar_tag(tag[1])[0])
                        else:
                            neftag = '{}.{}'.format(saveframe.tag_prefix, tag[0])
                            sf.add_tag(self.get_nmrstar_tag(neftag)[0], tag[1])
                    if saveframe.category == "nef_nmr_meta_data":
                        sf.add_tag("NMR_STAR_version", "3.2.0.15")

                else:
                    sf = pynmrstar.Saveframe.from_scratch(nef_data.category)
                    if nef_data.category == "_nef_chemical_shift":
                        sf.add_tag("_Assigned_chem_shift_list.Sf_category", 'nef_chemical_shift')
                    saveframe = [nef_data]
                for loop in saveframe:
                    if loop.category == "_nef_sequence":
                        self.cid = []  # Comp_index_ID list
                        for c in self.chains:  # Comp_index_ID initialized with 1
                            self.cid.append(1)
                        self.seqDict = {}

                    if loop.category == '_nef_distance_restraint':
                        r_index_id = 1
                    if loop.category == "_nef_chemical_shift":
                        cs_list += 1
                    lp = pynmrstar.Loop.from_scratch()
                    lp_cols = self.get_nmrstar_loop_tags(loop.get_tag_names())
                    for t in lp_cols:
                        lp.add_tag(t)
                    # print (loop.category,lp.category,lp.get_tag_names(),loop.get_tag_names())
                    for dat in loop.data:
                        if loop.category == "_nef_sequence":
                            dd = self.translate_seq_row(loop.get_tag_names(), lp.get_tag_names(), dat)
                            self.cid[
                                self.chains.index(dat[loop.get_tag_names().index('_nef_sequence.chain_code')])] += 1
                            for d in dd:
                                lp.add_data(d)
                                self.seqDict[(dat[loop.get_tag_names().index('_nef_sequence.chain_code')],
                                              dat[loop.get_tag_names().index('_nef_sequence.sequence_code')])] = (
                                    d[lp.get_tag_names().index('_Chem_comp_assembly.Entity_assembly_ID')],
                                    d[lp.get_tag_names().index('_Chem_comp_assembly.Comp_index_ID')])

                        elif loop.category == "_nef_chemical_shift":
                            dd = self.translate_cs_row(loop.get_tag_names(), lp.get_tag_names(), dat, iupac)
                            for d in dd:
                                d[lp.get_tag_names().index('_Atom_chem_shift.Assigned_chem_shift_list_ID')] = cs_list
                                lp.add_data(d)
                        elif loop.category == '_nef_distance_restraint':
                            dd = self.translate_restraint_row(loop.get_tag_names(), lp.get_tag_names(), dat, iupac)

                            for d in dd:
                                d[lp.get_tag_names().index('_Gen_dist_constraint.Index_ID')] = r_index_id
                                d[lp.get_tag_names().index(
                                    '_Gen_dist_constraint.Gen_dist_constraint_list_ID')] = rest_list
                                if len(dd) > 1:
                                    d[lp.get_tag_names().index('_Gen_dist_constraint.Member_logic_code')] = "OR"
                                lp.add_data(d)
                                r_index_id += 1
                        elif loop.category == '_nef_dihedral_restraint':
                            dd = self.translate_row(loop.get_tag_names(), lp.get_tag_names(), dat)
                            for d in dd:
                                d[lp.get_tag_names().index(
                                    '_Torsion_angle_constraint.Torsion_angle_constraint_list_ID')] = ang_list
                                lp.add_data(d)
                        elif loop.category == '_nef_rdc_restraint':
                            dd = self.translate_row(loop.get_tag_names(), lp.get_tag_names(), dat)
                            for d in dd:
                                d[lp.get_tag_names().index(
                                    '_RDC_constraint.RDC_constraint_list_ID')] = rdc_list
                                lp.add_data(d)
                        elif loop.category == '_nef_peak':
                            dd = self.translate_row(loop.get_tag_names(), lp.get_tag_names(), dat)
                            for d in dd:
                                d[lp.get_tag_names().index(
                                    '_Peak_row_format.Spectral_peak_list_ID')] = peak_list_list
                                lp.add_data(d)
                        else:
                            dd = self.translate_row(loop.get_tag_names(), lp.get_tag_names(), dat)
                            for d in dd:
                                lp.add_data(d)

                    # print (loop.data[0])
                    sf.add_loop(lp)
                star_data.add_saveframe(sf)
            star_data.normalize()
            with open(star_file, 'w') as wstarfile:
                wstarfile.write(str(star_data))
            info.append('File {} successfully written'.format(star_file))

        else:
            is_done = False
            error.append('Input file not readable')
        return is_done, json.dumps({'info': info, 'warning': warning, 'error': error})

    def nef_seq_from_eneity(self,entity):
        nef_seq_sf = pynmrstar.Saveframe.from_scratch(sf_name='nef_molecular_system',tag_prefix='_nef_molecular_system')
        nef_seq_sf.add_tag('sf_category','nef_molecular_system')
        nef_seq_sf.add_tag('sf_sf_framecode','nef_molecular_system')
        seq_loop = pynmrstar.Loop.from_scratch(category='_nef_sequence')
        nef_seq_tags = ['index','chain_code','sequence_code','residue_name','linking','residue_variant']
        for t in nef_seq_tags:
            seq_loop.add_tag(t)
        index = 1
        for e in entity:
            t=e.get_tag_names()
            chain_code = t.index('_Entity_comp_index.Entity_ID')
            seq_code = t.index('_Entity_comp_index.Auth_seq_ID')
            res_name = t.index('_Entity_comp_index.Comp_ID')
            for d in e.data:
                dat = [index,d[chain_code],d[seq_code],d[res_name],'.','.']
                seq_loop.add_data(dat)
        nef_seq_sf.add_loop(seq_loop)
        return nef_seq_sf




    def nmrstar_to_nef(self, star_file, nef_file=None, iupac=True):
        """
        Translates NEF into NMR-STAR
        :param star_file: input file in NMR-STAR format
        :param nef_file: output file name
        :param iupac: convert the atom nomenclature to IUPAC True/False default: True
        :return:
        """
        (file_path, file_name) = ntpath.split(os.path.realpath(star_file))
        is_done = True
        info = []
        warning = []
        error = []
        if nef_file is None:
            if iupac == False:
                nef_file = file_path + "/" + file_name.split(".")[0] + "_2.nef"
            else:
                nef_file = file_path + "/" + file_name.split(".")[0] + "_2.nef"
        (is_readable, dat_content, star_data) = self.read_input_file(star_file)
        try:
            nef_data = pynmrstar.Entry.from_scratch(star_data.entry_id)
        except AttributeError:
            nef_data = pynmrstar.Entry.from_scratch(file_name.split(".")[0])
            warning.append('Not a complete Entry')
        if is_readable:
            if dat_content == "Entry":
                if len(star_data.get_loops_by_category('Chem_comp_assembly'))>0:
                    self.chains = sorted(list(set(star_data.get_loops_by_category('Chem_comp_assembly')[0].
                                                  get_tag('Auth_asym_ID'))))
                elif len(star_data.get_loops_by_category('Entity_comp_index'))>0:
                    nef_seq_sf = self.nef_seq_from_eneity(star_data.get_loops_by_category('Entity_comp_index'))
                else:
                    print ("No Seq information found")
                    exit(1)
            elif dat_content == "Saveframe":
                self.chains = sorted(list(set(star_data[0].get_tag('Auth_asym_ID'))))
            elif dat_content == "Loop":
                self.chains = sorted(list(set(star_data.get_tag('Auth_asym_ID'))))
            else:
                is_done = False
                error.append('File content unknown')
            if dat_content == "Entry":
                for saveframe in star_data:
                    sf = pynmrstar.Saveframe.from_scratch(saveframe.name)

                    for tag in saveframe.tags:
                        if tag[0].lower() == "sf_category":
                            if saveframe.category == 'entity':
                                nef_data.add_saveframe(nef_seq_sf)
                            else:
                                sf.add_tag("sf_category", self.get_nef_tag(saveframe.category)[0])
                        else:
                            neftag = '{}.{}'.format(saveframe.tag_prefix, tag[0])
                            try:
                                sf.add_tag(self.get_nef_tag(neftag)[0], tag[1])
                            except ValueError:
                                print ("Tag not found {}".format(neftag))
                    if saveframe.category == "entry_information":
                        #sf.add_tag("NMR_STAR_version", "3.2.0.15")
                        # sf.add_tag("Generated_date", self.TimeStamp(time.time()), update=True)
                        try:
                            lp1 = saveframe.get_loop_by_category('_Software_applied_methods')
                            lp1.add_data(['NEFTranslator', 'NMRSTARtoNEF', '.'])
                            # print (lp1.tags)
                        except KeyError:
                            pass  # May be better to add audit loop
                    if saveframe.category != 'entity':
                        for loop in saveframe:
                            lp = pynmrstar.Loop.from_scratch()
                            lp_cols = self.get_nef_loop_tags(loop.get_tag_names())
                            for t in lp_cols:
                                lp.add_tag(t)
                            # print (loop.category,lp.category,lp.get_tag_names(),loop.get_tag_names())
                            cs_id = 1

                            for t in lp.get_tag_names():
                                if t.split(".")[1] == "index":
                                    index_id = lp.get_tag_names().index(t)
                            index = 1
                            dd_dummy = []
                            for dat in loop.data:
                                dd = self.translate_row2(loop.get_tag_names(), lp.get_tag_names(), dat)
                                for d in dd:
                                    d[index_id] = 0
                                    if d not in dd_dummy:
                                        dd_dummy.append(d[:])
                                        d[index_id] = index
                                        index+=1
                                        lp.add_data(d)
                            sf.add_loop(lp)
                        nef_data.add_saveframe(sf)
                    with open(nef_file, 'w') as wstarfile:
                        wstarfile.write(str(nef_data))
        else:
            is_done = False
            error.append('Input file not readable')
        return is_done, json.dumps({'info': info, 'warning': warning, 'error': error})


if __name__ == "__main__":
    #fname = sys.argv[1]
    fname = '/Users/kumaran/tmp_val/2jr2/2jr2.nef'
    bt = NEFTranslator()
    print (bt.check_mandatory_tags(fname,file_type='NEF'))
    bt.nef_to_nmrstar(nef_file=fname)
    #bt.nmrstar_to_nef(star_file=fname)
    #print (bt.validate_file('data/2l9r.str','A'))
    #fname = sys.argv[1]
    # f = open('neflist.txt','r').read().split("\n")
    # for fname in f:
    #     print ("Working on {}".format(fname))
    #     bt = NEFTranslator()
    #     bt.nef_to_nmrstar(fname)
