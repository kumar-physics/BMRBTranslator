'''
Created on Jan 18, 2017

@author: root
'''

import ntpath,os,csv,re,time,datetime,string,sys
try:
    import pynmrsta as bmrb
except ImportError as e:
    #print "Using local STAR parser",str(e) 
    (scriptPath,scriptName)=ntpath.split(os.path.realpath(__file__))
    sys.path.append(scriptPath+'/../PyNMRSTAR') #NMR-STAR and NEF-Parser added as a submodule and imported into this project. This is a separate git repository
    try:
        import bmrb
    except ImportError as e:
        print "ERROR: STAR parser from BMRB is not available"
        print str(e)
        exit(1)

class BMRBTranslator(object):
    '''
    Main class for BMRBTranslator
    '''
    __version__=1.0
    (scriptPath,scriptName)=ntpath.split(os.path.realpath(__file__))
    mapFile = scriptPath+'/../lib/NEF_NMRSTAR_equivalence.csv'
    
    
    atomDict = { 'CYS': ['N', 'CA', 'C', 'O', 'CB', 'SG', 'H', 'HA', 'HB2', 'HB3', 'HG'],
                 'ASP': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'OD2', 'H', 'HA', 'HB2', 'HB3', 'HD2'],
                 'SER': ['N', 'CA', 'C', 'O', 'CB', 'OG', 'H', 'HA', 'HB2', 'HB3', 'HG'],
                 'GLN': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'NE2', 'H', 'HA', 'HB2', 'HB3', 'HG2', 'HG3', 'HE21', 'HE22'],
                 'LYS': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'NZ', 'H', 'HA', 'HB2', 'HB3', 'HG2', 'HG3', 'HD2', 'HD3', 'HE2', 'HE3', 'HZ1', 'HZ2', 'HZ3'],
                 'ILE': ['N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'CD1', 'H', 'HA', 'HB', 'HG12', 'HG13', 'HG21', 'HG22', 'HG23', 'HD11', 'HD12', 'HD13'],
                 'PRO': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'HA', 'HB2', 'HB3', 'HG2', 'HG3', 'HD2', 'HD3'],
                 'THR': ['N', 'CA', 'C', 'O', 'CB', 'OG1', 'CG2', 'H', 'HA', 'HB', 'HG1', 'HG21', 'HG22', 'HG23'],
                 'PHE': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'H', 'HA', 'HB2', 'HB3', 'HD1', 'HD2', 'HE1', 'HE2', 'HZ'],
                 'ASN': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'ND2', 'H', 'HA', 'HB2', 'HB3', 'HD21', 'HD22'],
                 'GLY': ['N', 'CA', 'C', 'O', 'H', 'HA2', 'HA3'],
                 'HIS': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2', 'H', 'HA', 'HB2', 'HB3', 'HD1', 'HD2', 'HE1', 'HE2'],
                 'LEU': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'H', 'HA', 'HB2', 'HB3', 'HG', 'HD11', 'HD12', 'HD13', 'HD21', 'HD22', 'HD23'],
                 'ARG': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2', 'H', 'HA', 'HB2', 'HB3', 'HG2', 'HG3', 'HD2', 'HD3', 'HE', 'HH11', 'HH12', 'HH21', 'HH22'],
                 'TRP': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2', 'H', 'HA', 'HB2', 'HB3', 'HD1', 'HE1', 'HE3', 'HZ2', 'HZ3', 'HH2'],
                 'ALA': ['N', 'CA', 'C', 'O', 'CB', 'H', 'HA', 'HB1', 'HB2', 'HB3'],
                 'VAL': ['N', 'CA', 'C', "O'", 'CB', 'CG1', 'CG2', 'H', 'HA', 'HB', 'HG11', 'HG12', 'HG13', 'HG21', 'HG22', 'HG23', "O''"],
                 'GLU': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'OE2', 'H', 'HA', 'HB2', 'HB3', 'HG2', 'HG3', 'HE2'],
                 'TYR': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH', 'H', 'HA', 'HB2', 'HB3', 'HD1', 'HD2', 'HE1', 'HE2', 'HH'],
                 'MET': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'SD', 'CE', 'H', 'HA', 'HB2', 'HB3', 'HG2', 'HG3', 'HE1', 'HE2', 'HE3'],
                 'A': ["C3'", "O3'", "O5'", "C5'", "C4'", "O4'", "C1'", "C2'", 'N9', 'C8', 'N7', 'C5', 'C6', 'N6', 'N1', 'C2', 'N3', 'C4', "O2'", "H5'", "H5''", "H4'", "H3'", "H2'", "HO2'", "H1'", 'H8', 'H61', 'H62', 'H2', "HO5'"],
                 'C': ['P', "C3'", "O3'", 'OP1', 'OP2', "O5'", "C5'", "C4'", "O4'", "C1'", "C2'", 'N1', 'C6', 'C5', 'C4', 'N4', 'N3', 'C2', 'O2', "O2'", "H5'", "H5''", "H4'", "H3'", "H2'", "HO2'", "H1'", 'H41', 'H42', 'H5', 'H6'],
                 'G': ['P', "C3'", "O3'", 'OP1', 'OP2', "O5'", "C5'", "C4'", "O4'", "C1'", "C2'", 'N9', 'C8', 'N7', 'C5', 'C6', 'O6', 'N1', 'C2', 'N2', 'N3', 'C4', "O2'", "H5'", "H5''", "H4'", "H3'", "H2'", "HO2'", "H1'", 'H8', 'H1', 'H21', 'H22'],
                 'U': ['P', "C3'", "O3'", 'OP1', 'OP2', "O5'", "C5'", "C4'", "O4'", "C1'", "C2'", 'N1', 'C6', 'C5', 'C4', 'O4', 'N3', 'C2', 'O2', "O2'", "H5'", "H5''", "H4'", "H3'", "H2'", "HO2'", "H1'", 'H3', 'H5', 'H6', "HO3'"],
                 'DA': ["C3'", "O3'", "O5'", "C5'", "C4'", "O4'", "C1'", "C2'", 'N9', 'C8', 'N7', 'C5', 'C6', 'N6', 'N1', 'C2', 'N3', 'C4', "H5'", "H5''", "H4'", "H3'", "H2'", "H2''", "H1'", 'H8', 'H61', 'H62', 'H2', "HO5'"],
                 'DC': ['P', "C3'", "O3'", 'OP1', 'OP2', "O5'", "C5'", "C4'", "O4'", "C1'", "C2'", 'N1', 'C6', 'C5', 'C4', 'N4', 'N3', 'C2', 'O2', "H5'", "H5''", "H4'", "H3'", "H2'", "H2''", "H1'", 'H41', 'H42', 'H5', 'H6'],
                 'DG': ['P', "C3'", "O3'", 'OP1', 'OP2', "O5'", "C5'", "C4'", "O4'", "C1'", "C2'", 'N9', 'C8', 'N7', 'C5', 'C6', 'O6', 'N1', 'C2', 'N2', 'N3', 'C4', "H5'", "H5''", "H4'", "H3'", "H2'", "H2''", "H1'", 'H8', 'H1', 'H21', 'H22'],
                 'DT': ['P', "C3'", "O3'", 'OP1', 'OP2', "O5'", "C5'", "C4'", "O4'", "C1'", "C2'", 'N1', 'C6', 'C5', 'C4', 'O4', 'N3', 'C2', 'O2', 'C7', "H5'", "H5''", "H4'", "H3'", "H2'", "H2''", "H1'", 'H3', 'H71', 'H72', 'H73', 'H6', "HO3'"]}
    
    
    def TimeStamp(self,ts):
        return datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    
    def __init__(self):
        '''
        Constructor
        '''
        self.ReadMapFile()

        
    def NEFtoStar(self,inFile):
        '''
        Translates NEF file to NMR-STAR
        '''
        self.softwareSpecificSfID=0
        self.nefFile=inFile
        (self.nefFilePath,self.nefFileName) = ntpath.split(inFile)
        self.logFile = inFile.split("."+self.nefFileName.split(".")[-1])[0]+"_.log"
        self.starFile = inFile.split("."+self.nefFileName.split(".")[-1])[0]+"_.str"
        self.log = open(self.logFile,'w')
        self.log.write("\n%s:%s\n"%(string.ljust("Input",25),self.nefFile))
        self.log.write("%s:%s\n"%(string.ljust("Output",25),self.starFile))
        self.log.write("%s:%s\n"%(string.ljust("Log ",25),self.logFile))
        self.log.write("%s:%s\n"%(string.ljust("Translator Version",25),self.__version__))
        self.log.write("%s:%s\n"%(string.ljust("STAR parser Version",25),bmrb._VERSION))
        self.log.write("%s:%s\n\n"%(string.ljust("Date",25),self.TimeStamp(time.time())))
        try:
            self.Log("Reading input file")
            self.nef = bmrb.Entry.from_file(self.nefFile)
        except Exception as e:
            self.Log("Problem with input file",2)
            self.Log(str(e),2)
            print "Error: Can't parse the input NEF file"
            print str(e)
            print "Check the log",self.logFile
            self.log.close()
            exit(1)
        chains=sorted(list(set(self.nef.get_loops_by_category('nef_sequence')[0].get_tag('chain_code'))))
           
        self.star = bmrb.Entry.from_scratch(self.nef.entry_id)
        for saveframe in self.nef:
            
            if saveframe.get_tag("sf_category")[0] in self.tagMap[0]:
                sf = bmrb.Saveframe.from_scratch(saveframe.name)
                self.details = {}
                if saveframe.name == "nef_nmr_meta_data":
                    sf.add_tag("_Entry.NMR_STAR_version","3.2.0.4")
                for tag in saveframe.tags:
                    nef_tag = "%s.%s"%(saveframe.tag_prefix,tag[0])
                    try:
                        star_auth_tag = self.tagMap[1][self.tagMap[0].index(nef_tag)]
                        star_tag = self.tagMap[2][self.tagMap[0].index(nef_tag)]
                    except ValueError:
                        self.Log("Tag %s not found in NEF dictionary"%(nef_tag))
                        self.details[tag[0]] = tag[1]
                        star_auth_tag = ""
                        star_tag = ""
                    if star_auth_tag != "" and star_tag != "":
                        if star_auth_tag == star_tag:
                            if  star_auth_tag.split(".")[-1] == "Sf_category":
                                star_cat = self.tagMap[1][self.tagMap[0].index(tag[1])]
                                sf.add_tag(star_tag,star_cat)
                            else:
                                sf.add_tag(star_tag,tag[1])
                        else:
                            sf.add_tag(star_auth_tag,tag[1])
                            sf.add_tag(star_tag,tag[1])
                    else:
                        self.Log("Equivalent STAR tag for %s is empty"%(nef_tag),2)
                    const_index=1
                for loop in saveframe:
                    lp = bmrb.Loop.from_scratch()
                    missing_col = []
                    auth_col=[]
                    for coln in loop.columns:
                        nef_tag = "%s.%s"%(loop.category,coln)
                        try:
                            star_auth_tag = self.tagMap[1][self.tagMap[0].index(nef_tag)]
                            star_tag = self.tagMap[2][self.tagMap[0].index(nef_tag)]
                        except ValueError:
                            self.Log("Tag %s not found in NEF dictionary"%(nef_tag))
                            self.details[nef_tag] =  str(loop.get_tag(nef_tag))
                            star_auth_tag = ""
                            star_tag = ""
                        #print nef_tag,star_auth_tag,star_tag
                        if star_auth_tag != "" and star_tag != "" :
                            if star_auth_tag != star_tag:
                                
                                auth_col.append(loop.columns.index(coln))
                                lp.add_column(star_auth_tag)
                                lp.add_column(star_tag)
                            else:
                                lp.add_column(star_tag)
                        else:
                            missing_col.append(loop.columns.index(coln))
                    atm_id=[i for i in range(len(lp.columns)) if "Atom_ID" in lp.columns[i]]
                    ent_assembly_id=[i for i in range(len(lp.columns)) if "Entity_assembly_ID" in lp.columns[i]]
                    auth_asym_id=[i for i in range(len(lp.columns)) if "Auth_asym_ID" in lp.columns[i]]
                    if sf.category=="assigned_chemical_shifts":
                        lp.add_column("_Atom_chem_shift.Ambiguity_code")
                    if sf.category=="general_distance_constraints":
                        lp.add_column("_Gen_dist_constraint.Member_logic_code")
                        const_id = 1
                    #if len(auth_col)!=0:
                        #print len(auth_col),len(loop.columns),len(lp.columns),len(loop.data[0]),auth_col
                        #print loop.columns
                        #print lp.columns        
                    for dat in loop.data:
                        if len(auth_col) == 0:
                            if len(missing_col) == 0:
                                lp_data = dat[:]
                            else:
                                tmp_dat = dat[:]
                                for m in sorted(missing_col, reverse = True):
                                    del tmp_dat[m]
                                lp_data = tmp_dat[:]
                                lp.add_data(lp_data[:])
                        else:
                            # If author colum and star col are different we need to map
                            # Mapping is done for atom and chain code 
                            # in addition to that ambiguity code and member logic code is added
                            lp_data_list = [dat[:]]
                            ref_index=0
                            #lp_data_list=[]
                            for k in auth_col:
                                if "chain_code" in loop.columns[k]:
                                    #print ref_index,auth_col
                                    for lp_data in lp_data_list:
                                        try:
                                            lp_data.insert(k+1+ref_index,chains.index(lp_data[k+ref_index])+1)
                                        except ValueError:
                                            lp_data.insert(k+1+ref_index,dat[:][k])
                                    ref_index+=1
                                if "sequence_code" in loop.columns[k]:
                                    for lp_data in lp_data_list:
                                        lp_data.insert(k+1+ref_index,lp_data[k+ref_index])
                                    ref_index+=1
                                if "residue_name" in loop.columns[k]:
                                    for lp_data in lp_data_list:
                                        res_name = dat[:][k]
                                        lp_data.insert(k+1+ref_index,dat[:][k])
                                    ref_index+=1
                                if "atom_name" in loop.columns[k]:
                                    tmp_lp_data=[]
                                    for lp_dat in lp_data_list:
                                        
                                        nef_atm_name = dat[:][k]
                                        #print nef_atm_name, len(lp_data_list),lp_data_list
                                        star_atm_list = self.EquivalentAtom(res_name, nef_atm_name)
                                        if len(star_atm_list) == 0:
                                            
                                            lp_dat.insert(k+1+ref_index,dat[:][k])
                                            #ref_index+=1
                                            self.Log("No matching STAR atom found for %s-%s in loop %s"%(res_name,nef_atm_name,lp.category),1)
                                            if sf.category == "assigned_chemical_shifts":
                                                lp_data.append('1')
                                            if sf.category=="general_distance_constraints" and lp_data[-1]!="OR":
                                                lp_dat.append("OR")
                                        elif len(star_atm_list) == 1:
                                            lp_dat.insert(k+1+ref_index,star_atm_list[0])
                                            #ref_index+=1
                                            if sf.category == "assigned_chemical_shifts":
                                                lp_dat.append('1')
                                            if sf.category=="general_distance_constraints" and lp_data[-1]!="OR":
                                                lp_dat.append("OR")
                                        else:
                                            
                                            for star_atm in star_atm_list:
                                                tmp=lp_dat[:]
                                                tmp.insert(k+1+ref_index,star_atm)
                                                if sf.category == "assigned_chemical_shifts":
                                                    tmp.append('2')
                                                if sf.category=="general_distance_constraints" and tmp[-1]!="OR":
                                                    tmp.append("OR")
                                                
                                                tmp_lp_data.append(tmp)
                                            
                                    if len(tmp_lp_data)>0:
                                         
                                        lp_data_list = [x[:] for x in tmp_lp_data[:]]
                                    ref_index+=1
                                   
                        
                            for lp_data in lp_data_list:
                                if "Index_ID" in lp.columns:
                                    lp_data[lp.columns.index("Index_ID")]=const_index
                                    const_index+=1
                                #print lp.columns
                                #print lp_data
                                if len(missing_col)>0:
                                    for i in missing_col:
                                        del lp_data[-1]
                                #print missing_col
                                #print auth_col
                                lp.add_data(lp_data)
                    if len(lp.columns)>0:
                        sf.add_loop(lp)
                    else:
                        self.details["ss_loop"]=str(loop)
                    #print sf.name,lp.category,lp.columns,missing_col
                if self.details:
                    sf.add_tag("Details","\"%s\""%(str(self.details)))                    
            else:
                self.softwareSpecificSfID+=1
                self.Log("Saveframe category '%s'not found in NEF dictionary"%(saveframe.name))
                if self.softwareSpecificSfID==1:
                    soft_specific_sf = bmrb.Saveframe.from_scratch("software_specific_info")
                    soft_specific_sf.add_tag("_Software_specific_info_list.Sf_category","Software_specific_saveframes")
                    soft_specific_sf.add_tag("_Software_specific_info_list.Sf_framecode","software_specific_info")
                    soft_specific_sf.add_tag("_Software_specific_info_list.ID",1)
                    soft_specific_lp=bmrb.Loop.from_scratch()
                    soft_specific_lp.add_column("_Software_specific_info.Software_saveframe_ID")
                    soft_specific_lp.add_column("_Software_specific_info.Software_saveframe")
                    soft_specific_lp.add_column("_Software_specific_info.Software_specific_info_list_ID")
                    soft_specific_lp_dat=[self.softwareSpecificSfID,"\"%s\""%(str(saveframe)),'1']
                    soft_specific_lp.add_data(soft_specific_lp_dat)
                else:
                    soft_specific_lp_dat=[self.softwareSpecificSfID,"\"%s\""%(str(saveframe)),'1']
                    soft_specific_lp.add_data(soft_specific_lp_dat)
            
            self.star.add_saveframe(sf)
        if self.softwareSpecificSfID:
            soft_specific_sf.add_loop(soft_specific_lp)
            self.star.add_saveframe(soft_specific_sf)
        self.star.normalize()
        #print self.star
        #print soft_specific_sf
        with open(self.starFile,'w') as wstarfile:
            wstarfile.write(str(self.star))
        self.Log("Output file written")
        self.log.close()
    
    
    def Log(self,txt,flag=0):
        if flag==0:
            self.log.write("%s\t%s\n"%(self.TimeStamp(time.time()),txt))
        elif flag==1:
            self.log.write("%s\tWARNING:%s\n"%(self.TimeStamp(time.time()),txt))
        elif flag==2:
            self.log.write("%s\tERROR:%s\n"%(self.TimeStamp(time.time()),txt))
        else:
            print "Logging function takes only three numbers as first argument(0,1,2)"
            exit(1)
            
            
            
            
        
        
        
        
        
    def STARtoNEF(self,inFile):
        '''
        Translates NMR-STAR file NEF
        '''
        self.starFile=inFile
        (self.starFilePath,self.starFileName) = ntpath.split(inFile)
        self.logFile = inFile.split("."+self.starFileName.split(".")[-1])[0]+"_.log"
        self.nefFile = inFile.split("."+self.starFileName.split(".")[-1])[0]+"_.nef"
        self.log = open(self.logFile,'w')
        self.log.write("\n%s:%s\n"%(string.ljust("Input",25),self.starFile))
        self.log.write("%s:%s\n"%(string.ljust("Output",25),self.nefFile))
        self.log.write("%s:%s\n"%(string.ljust("Log ",25),self.logFile))
        self.log.write("%s:%s\n"%(string.ljust("Translator Version",25),self.__version__))
        self.log.write("%s:%s\n"%(string.ljust("STAR parser Version",25),bmrb._VERSION))
        self.log.write("%s:%s\n\n"%(string.ljust("Date",25),self.TimeStamp(time.time())))
        try:
            self.star = bmrb.Entry.from_file(self.starFile)
        except Exception as e:
            self.log.write("%s\t%s"%(self.TimeStamp(time.time()),str(e)))
            print "Error: Can't parse the input NEF file"
            print str(e)
            print "Check the log",self.logFile
            self.log.close()
            exit(1)
            
        
    
        
    def ReadMapFile(self):
        '''Reads the NEF_NMRSTAR_equivalence.csv file and create a mapping as a list'''
        with open(self.mapFile,'rb') as csvfile:
            spamreader = csv.reader(csvfile,delimiter=',')
            map_dat=[]
            for r in spamreader:
                #print r
                if r[0][0]!='#':
                    map_dat.append(r)
        self.tagMap=map(list,zip(*map_dat))
        
    def EquivalentAtom(self,res,nefAtom):
        '''
        Expands regular expression in atom names
        '''
        try:
            atms=self.atomDict[res]
            alist=[]
            try:
                refatm=re.findall(r'(\S+)([xyXY])([%*])$|(\S+)([%*])$|(\S+)([xyXY]$)',nefAtom)[0]  
                set=[refatm.index(i) for i in refatm if i!=""]
                if set==[0,1,2]:
                    pattern=re.compile(r'%s\S\d+'%(refatm[0]))
                    alist=[i for i in atms if re.search(pattern, i)]
                    if refatm[1]=="y" or refatm[1]=="Y":
                        alist.reverse()
                elif set==[3,4]:
                    if refatm[4]=="%":
                        pattern=re.compile(r'%s\d+'%(refatm[3]))
                    elif refatm[4]=="*":
                        pattern=re.compile(r'%s\S+'%(refatm[3]))
                    else:
                        print "something wrong"
                    alist=[i for i in atms if re.search(pattern, i)]
                    
                elif set==[5,6]:
                    pattern=re.compile(r'%s\S+'%(refatm[5]))
                    alist=[i for i in atms if re.search(pattern, i)]
                    if len(alist)!=2:
                        alist=[]
                    elif refatm[6]=="y" or refatm[6]=="Y":
                        #alist.reverse()[]
                        alist=alist[-1:]
                    elif refatm[6]=="x" or refatm[6]=="X":
                        alist=alist[:1]
                    else:
                        print "Something wrong"
                        
                else:
                    print "Wrong regular expression"
            except IndexError:
                #print nefAtom
                pass
            if len(alist)==0:
                if nefAtom in atms:
                    alist.append(nefAtom)
        except KeyError:
            #self.logfile.write("%s\tResidue not found,%s,%s\n"%(self.TimeStamp(time.time()),res,nefAtom))
            #print "Residue not found",res,nefAtom
            self.Log("Non-standard residue found %s"%(res),1)
            alist=[]
        return alist                    
     
        
if __name__=="__main__":
    fname=sys.argv[1]
    p=BMRBTranslator()
    p.NEFtoStar(fname)
    