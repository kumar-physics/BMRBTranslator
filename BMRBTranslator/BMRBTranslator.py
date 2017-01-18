'''
Created on Jan 18, 2017

@author: root
'''

import ntpath,os,csv,re,time,datetime

class BMRBTranslator(object):
    '''
    Main class for BMRBTranslator
    '''
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
        (self.nefFilePath,self.nefFileName) = ntpath.split(inFile)
        self.logFile = inFile.split("."+self.nefFileName.split(".")[-1])[0]+"_.log"
        self.starFile = inFile.split("."+self.nefFileName.split(".")[-1])[0]+"_.str"
        
        
        
        
    def STARtoNEF(self,inFile):
        '''
        Translates NMR-STAR file NEF
        '''
        (self.starFilePath,self.starFileName) = ntpath.split(inFile)
        self.logFile = inFile.split("."+self.starFileName.split(".")[-1])[0]+"_.log"
        self.nefFile = inFile.split("."+self.starFileName.split(".")[-1])[0]+"_.str"
        #print self.logFile
        
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
                refatm=re.findall(r'(\S+)([xy])([%*])$|(\S+)([%*])$|(\S+)([xy]$)',nefAtom)[0]  
                set=[refatm.index(i) for i in refatm if i!=""]
                if set==[0,1,2]:
                    pattern=re.compile(r'%s\S\d+'%(refatm[0]))
                    alist=[i for i in atms if re.search(pattern, i)]
                    if refatm[1]=="y":
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
                    elif refatm[6]=="y":
                        #alist.reverse()[]
                        alist=alist[-1:]
                    elif refatm[6]=="x":
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
            print "Residue not found",res,nefAtom
            alist=[]
        return alist                    
     
        
if __name__=="__main__":
    p=BMRBTranslator()
    print p.EquivalentAtom('ASN', 'H*')