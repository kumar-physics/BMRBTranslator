# NEF TRANSLATOR

This script translates NMR Exchange Format(NEF) file into NMR-STAR file. Both data formats are basically STAR format and
'pynmrstar' package from BMRB is used to read and write both file types. 

###Dependency 
NEFTranslator require 'pynmrstar v 2.6' package which is available via pip installation

### Usage

#### 1. Validate input file

Input file can be either in NEF or NMR-STAR format. It can be a complete file with many saveframes and loops or just a
single saveframe with chemical shift data or just a loop that contain some data. Usually depositor will provide the content
of the file and according to that the file content  parameter is set

- 'A' for all in one file (both chemical shift and restraints)
- 'S' for chemical shift file
- 'R' for restraint file


#### 2. To validate all in one file

```python
In [1]: import NEFTranslator

In [2]: import json

In [3]: nt = NEFTranslator.NEFTranslator()

In [4]: (status_flag,json_data) = nt.validate_file('data/2mqq.nef','A')

In [5]: dat = json.loads(json_data)

In [6]: status_flag
Out[6]: True

In [7]: dat
Out[7]: 
{'info': ['8 saveframes and 8 loops found',
  '8 saveframes and 8 loops found with nef prefix',
  'data/2mqq.nef is a NEF file'],
 'warning': [],
 'error': [],
 'FILE': 'NEF'}

```

status_flag is a boolean which gives the information whether the script is successful or not. If the script fails then
the error tag will have some information. Following example a chemical shift will validated as all in one file. 

```python
In [8]: (status_flag,json_data) = nt.validate_file('data/norest.nef','A')

In [9]: dat = json.loads(json_data)

In [10]: status_flag
Out[10]: False

In [11]: dat
Out[11]: 
{'info': ['11 saveframes and 15 loops found',
  '11 saveframes and 15 loops found with nef prefix',
  'data/norest.nef is a NEF file'],
 'warning': [],
 'error': ['_nef_distance_restraint loop not found'],
 'FILE': 'NEF'}

In [12]: dat['error']
Out[12]: ['_nef_distance_restraint loop not found']
```
If you set file content parameter as 'S', then it will succeed.

#### 3. Extract sequence from chemical shift loop

```python

In [13]: (status_flag,json_data) = nt.get_seq_from_cs_loop('data/norest.nef')

In [14]: status_flag
Out[14]: True

In [15]: dat = json.loads(json_data)

Out[17]: 
[{'A': ['HIS',
   'MET',
   'SER',
   'HIS',
   'THR',
   'GLN',
   'VAL',
   'ILE',
   'GLU',
   'LEU',
   'GLU',
   'ARG',
   'LYS',
   'PHE',
   'SER',
   'HIS',
   'GLN',
   'LYS',
   'TYR',
   'LEU',
   'SER',
   'ALA',
   'PRO',
   'GLU',
   'ARG',
   'ALA',
   'HIS',
   'LEU',
   'ALA',
   'LYS',
   'ASN',
   'LEU',
   'LYS',
   'LEU',
   'THR',
   'GLU',
   'THR',
   'GLN',
   'VAL',
   'LYS',
   'ILE',
   'TRP',
   'PHE',
   'GLN',
   'ASN',
   'ARG',
   'ARG',
   'TYR',
   'LYS',
   'THR',
   'LYS',
   'ARG',
   'LYS',
   'GLN',
   'LEU',
   'SER',
   'SER',
   'GLU',
   'LEU',
   'GLY']}]
```

sequence can be found in 'DATA' tag in json data. It is represented as dictionary with chain id as key and sequence as value.
If the file contains multiple chemical shift loop then  you will get list of dictionaries.

##### NEF to NMR-STAR translation

```python
In [18]: (status_flag,json_data) = nt.nef_to_nmrstar('data/norest.nef')

In [19]: status_flag
Out[19]: True

In [20]: dat = json.loads(json_data)

In [21]: dat
Out[21]: 
{'info': ['File /home/kumaran/git/py-wwpdb_utils_nmr/wwpdb/utils/nmr/NEFTranslator/data/norest.str successfully written'],
 'warning': [],
 'error': []}

```
This will translate input NEF file into NMR-STAR file. If you output file name is not specified, then it will replace .nef with .str.
Output file name can be specified as parameter. 
```python
(status_flag,json_data) = nt.nef_to_nmrstar('data/norest.nef',star_file='data/outfile.str')
```

As mentioned above, this script can be used to for three main purposes. 

- validated input file format; works for both NEF and NMR-STAR `validate_file`
- extract sequence from chemical shift loop 'get_seq_from_cs_loop'
- translate NEF data into NMR-STAR 'nef_to_nmrstar'
