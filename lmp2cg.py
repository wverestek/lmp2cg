#!/usr/bin/env python3  
"""
unfrotunately at the moment unwrapped coordinates in both the data as well as 
the dump file are necessary.
"""

import sys as sys
import os as os
import datetime as datetime

import numpy as np
import networkx as nx
import networkx.algorithms.isomorphism as iso
#import moltemplate as mt
#import matplotlib.pyplot as plt
import code
from inspect import currentframe, getframeinfo


    
# call: interactive(currentframe(),locals())
def interactive(cf,l):
    frameinfo = getframeinfo(cf)
    strvar = '#####DEBUG#####MODE#####DEBUG#####MODE#####\n'\
    'INTERACTIVE DEBUG MODE: \npress CTRL+D to continue script execution\n'\
    +'   FILE: '+frameinfo.filename+'\n'\
    +'   LINE: '+str(frameinfo.lineno)+'\n'\
    +'   FUNC: '+frameinfo.function+'\n'
    code.interact(local=dict(globals(), **l),banner=strvar)


hkeywords = ["atoms", 
             "ellipsoids",
             "lines",
             "triangles",
             "bodies",
             "bonds",
             "angles",
             "dihedrals",
             "impropers",
             "atom types",
             "bond types",
             "angle types",
             "dihedral types",
             "improper types",
             "xlo xhi",
             "ylo yhi",
             "zlo zhi",
             "xy xz yz"]
             
skeywords = [["Masses","atom types"],
             ["Atoms","atoms"],
             ["Ellipsoids","ellipsoids"],
             ["Lines","lines"],
             ["Triangles","triangles"],
             ["Bodies","bodies"],
             ["Bonds","bonds"],
             ["Angles","angles"],
             ["Dihedrals","dihedrals"],
             ["Impropers","impropers"],
             ["Velocities","atoms"],
             ["Pair Coeffs","atom types"],
             ["Bond Coeffs","bond types"],
             ["Angle Coeffs","angle types"],
             ["Dihedral Coeffs","dihedral types"],
             ["Improper Coeffs","improper types"],
             ["BondBond Coeffs","angle types"],
             ["BondAngle Coeffs","angle types"],
             ["MiddleBondTorsion Coeffs","dihedral types"],
             ["EndBondTorsion Coeffs","dihedral types"],
             ["AngleTorsion Coeffs","dihedral types"],
             ["AngleAngleTorsion Coeffs","dihedral types"],
             ["BondBond13 Coeffs","dihedral types"],
             ["AngleAngle Coeffs","improper types"],
             ["Molecules","atoms"]]



class lammps_data:

    # --------------------------------------------------------------------
    def __init__(self,*list):
        self.nselect = 1
    
        if len(list) == 0:
            self.title = "LAMMPS data file"
            self.names = {}
            self.headers = {}
            self.sections = {}
            return

        file = list[0]
        f = open(file)

        self.title = f.readline()
        self.names = {}
        
        headers = {}
        while 1:
            line = f.readline()
            # mod WV 2020-02-18
            #line = line.strip()
            line = line.split('#')[0].strip()
            if len(line) == 0:
                continue
            found = 0
            for keyword in hkeywords:
                if line.find(keyword) >= 0:
                    found = 1
                    words = line.split()
                    if keyword == "xlo xhi" or keyword == "ylo yhi" or \
                        keyword == "zlo zhi":
                        headers[keyword] = (float(words[0]),float(words[1]))
                    elif keyword == "xy xz yz":
                        headers[keyword] = (float(words[0]),float(words[1]),float(words[2]))
                    else:
                        headers[keyword] = int(words[0])
            if not found:
                break

        sections = {}
        while 1:
            found = 0
            for pair in skeywords:
                keyword,length = pair[0],pair[1]
                if keyword == line:
                    found = 1
                    #if not headers.has_key(length):
                    if length not in headers:
                        raise StandardError("data section \'%s\' has no matching header value" % line)
                        # sys.exit("data section %s has no matching header value" % line)
                    f.readline()
                    list = []
                    for i in range(headers[length]): list.append(f.readline().strip())
                    sections[keyword] = list
            if not found:
                #raise StandardError,"invalid section %s in data file" % line
                sys.exit("data section \'%s\' has no matching header value" % line)
            f.readline()
            line = f.readline()
            if not line:
                break
            # mod WV 2020-02-18
            #line = line.strip()
            line = line.split('#')[0].strip()
      
        f.close()
        self.headers = headers
        self.sections = sections
        self.sectionslist = [i for i in self.sections]
    # --------------------------------------------------------------------
    # extract info from data file fields
    def get(self,*list):
        if len(list) == 1:
            field = list[0]
            array = []
            lines = self.sections[field]
            for line in lines:
                words = line.split()
                values = map(float,words)
                array.append(values)
            return array
        elif len(list) == 2:
            field = list[0]
            n = list[1] - 1
            vec = []
            lines = self.sections[field]
            for line in lines:
                words = line.split()
                vec.append(float(words[n]))
            return vec
        else:
            raise StandardError("invalid arguments for data.get()")
    # --------------------------------------------------------------------



class bead_template:

    def __init__(self,*list):
        self.name = []
        self.type = []
        self.style = []
        self.data = []
        
    def readfile(self, file):
        with open(file, 'r') as f:
            # line = 'True'
            # while line:
            for line in f:
                line = line.strip().split()
                # check for empty lines and comments
                if not line:
                    # line = f.readline().strip().split()
                    line = 'True'
                    continue
                if line[0].startswith('#'):
                    # print('start with #:',line)
                    # line = f.readline().strip().split()
                    line = 'True'
                    continue
                
                # check for style
                if line[2].casefold() == 'id':
                    # self.style.append('id')
                    style = 'id'
                elif line[2].casefold() == 'mol':
                    style = 'mol'    
                elif line[2].casefold() == 'del_template':
                    style = 'del_template'    
                elif line[2].casefold() == 'del_id':
                    style = 'del_id'
                elif line[2].casefold() == 'del_type':
                    print('deleting atom type is not supported yet!')
                    style = 'del_type'    
                elif line[2].casefold() == 'del_molid':
                    print('deleting molecule id is not supported yet!')
                    style = 'del_molid'
                else:
                    print('beat template style \''+line[2]+'\' is unknown')
                    sys.exit(0)


                # initially all fields
                idx = len(line)
                # find # for comments
                for index,string in enumerate(line):
                    if '#' in string:
                        idx = index
                        break
                #print(line)
                #print(line[0:idx])

                # build data list and check for ranges, e.g. 1:3 4 5 -> [1, 2, 3, 4, 5]
                data_list = []
                for var in line[3:idx]:
                    range_idx = list(map(int,var.split(':')))
                    if len(range_idx)==1:
                        data_list.append(int(var))
                    elif len(range_idx)==2:
                        for i in list(range(range_idx[0],range_idx[1]+1)):
                            data_list.append(i)
                    else:
                        print('weird stuff happening.... ')
                        sys.exit(0)
                            
                # append to data
                self.name.append(line[0])
                self.type.append(int(line[1]))
                self.style.append(style)
                self.data.append(data_list)
        return



def initialize(arguments):
    # print()
    print('##### Initializing #####')
    # default values
    options = {'indata'  :'lmp.data',\
               'style'   :'full',\
               'template':'cg.template',\
               'outdata' :'cg.data',\
               'pos'     :'com', \
               'remap'   :'False'}

    #          'VOTCA'   :'mapping.xml',\
    #          'MSCG'    :'top.in',\
    
    if len(arguments) == 0:
        print('using default files')
    elif len(arguments)%2 == 0:
        for key,word in zip(arguments[0::2],arguments[1::2]):
            # print('argument list:',key,word)
            if (key in ['--data','-d']):
                options.update({'indata':word})
            elif (key == '--style' or key == '-s'):
                options.update({'style':word})
            elif (key == '--template' or key == '-t'):
                options.update({'template':word})
            elif (key == '--trajectory' or key == '-traj'):
                options.update({'trajectory':word})
            elif (key == '--outdata' or key == '-o'):
                options.update({'outtemplate':word})
            elif (key == '--VOTCA' or key == '-V'):
                options.update({'VOTCA':word})
            elif (key == '--mscg' or key == '--MSCG' or key == '-m'):
                options.update({'MSCG':word})
            elif (key == '--pos' or key == '-pos' or key == '-p'):
                if (word == 'geo' or word == 'GEO'): options.update({'pos':'mean'})
                if (word == 'com' or word == 'COM'): options.update({'pos':'com'})
                if (word == 'charge' or word == 'CHARGE' or word == 'q'): options.update({'pos':'charge'})
            else:
                raise ValueError("unknown keywords:",key,word)
    else:
        raise ValueError('wrong number of arguments!')
               
    for i,j in options.items():
        print('{:10s} {:s}'.format(i,j))

    return(options)

def lammps2networkx(options,d,BT):
    # Masses ----------------------------------------------------------------
    if 'Masses' in d.sectionslist:
        type2mass = {int(i):float(j) for i,j in zip(d.get('Masses',1),d.get('Masses',2))}
    else:
        print('no masses defined in data file, using 1 as default')
        type2mass = {i:1 for i in range(1,d.headers['atom types']+1)}
        
    print('initializing graph... ', end=" ")
    AA = nx.Graph()
    print('done!')
    
    print('AA nodes: extracting data... ', end=" ")
    # Box dimension
    AA.box = list([[d.headers['xlo xhi'][0],d.headers['xlo xhi'][1]],\
                   [d.headers['ylo yhi'][0],d.headers['ylo yhi'][1]],\
                   [d.headers['zlo zhi'][0],d.headers['zlo zhi'][1]]])

    lx = AA.box[0][1] - AA.box[0][0] #d.headers['xlo xhi'][1] - d.headers['xlo xhi'][0]
    ly = AA.box[1][1] - AA.box[1][0] #d.headers['ylo yhi'][1] - d.headers['ylo yhi'][0]
    lz = AA.box[2][1] - AA.box[2][0] #d.headers['zlo zhi'][1] - d.headers['zlo zhi'][0]
    # IDs are always in the first coloumn
    ilist = list(map(int,d.get('Atoms',1)))
    natoms = len(ilist)
        
    if options['style']=='full':
        #interactive(currentframe(),locals())
        # id 
        # ilist = list(map(int,d.get('Atoms',1)))
        # mol
        mlist = list(map(int,d.get('Atoms',2)))
        # type
        tlist = list(map(int,d.get('Atoms',3)))
        # charge
        qlist = list(map(float,d.get('Atoms',4)))
        # position
        xlist = list(map(float,d.get('Atoms',5)))
        ylist = list(map(float,d.get('Atoms',6)))
        zlist = list(map(float,d.get('Atoms',7)))
        # images
        try:
            ix = list(map(int,d.get('Atoms',8)))
            iy = list(map(int,d.get('Atoms',9)))
            iz = list(map(int,d.get('Atoms',10)))
        except:
            ix = [0] * natoms
            iy = [0] * natoms
            iz = [0] * natoms
        
        xlist = [x+i*lx for x,i in zip(xlist,ix)]
        ylist = [y+i*ly for y,i in zip(ylist,iy)]
        zlist = [z+i*lz for z,i in zip(zlist,iz)]
        
        AAlist = [[il,{'mol':ml,\
                      'type':tl,\
                      'q':ql,\
                      'x':xl,\
                      'y':yl,\
                      'z':zl,\
                      'fx':0.0,\
                      'fy':0.0,\
                      'fz':0.0,\
                      'occupied':0,\
                      'bead_type':0,\
                      'bead_name':'',\
                      'bead_id':0}] \
                 for il,ml,tl,ql,xl,yl,zl \
                 in zip(ilist,mlist,tlist,qlist,xlist,ylist,zlist)]
        print('adding to graph... ', end=" ")
        AA.add_nodes_from(AAlist)
        print('done!')
        
        # Edges --------------------------------------------------------------
        print('AA bonds: extracting data... ', end=" ")
        ilist = list(map(int,d.get('Bonds',1)))
        tlist = list(map(int,d.get('Bonds',2)))
        a1list = list(map(int,d.get('Bonds',3)))
        a2list = list(map(int,d.get('Bonds',4)))
        AAlist = [[a1,a2,{'id':il,\
                      'type':tl}] \
                 for a1,a2,il,tl
                 in zip(a1list,a2list,ilist,tlist)]
        print('adding to graph... ', end=" ")
        AA.add_edges_from(AAlist)
        print('done!')

    return AA,type2mass
    
def find_and_mark_subgraph_by_ID(AA,BeadName,BeadType,BeadData):
    # Subgraph that is searched for
    S = AA.subgraph(BeadData)
    print('checking if atoms of bead template \''+str(BeadName)+'\' are connected... ', end=" ", flush=True)
    if nx.algorithms.components.is_connected(S):
        print('done!', flush=True)
    else:
        print("Subgraph is not connected... skipping!", flush=True)
        for i in S.nodes(data=True): print(i)
        print("Subgraph is not connected... skipping!")
        return
    
    print('filtering available atoms by type and occupancy... ', end=" ", flush=True)
    atomtypes = set(nx.get_node_attributes(S,'type').values())

    # searching on a not fully occupied subgraph coudl speed up the search for large systems???
    not_occupied_nodes = [atomid for (atomid,atomtype),occupancy in \
                          zip(nx.get_node_attributes(AA,'type').items(),nx.get_node_attributes(AA,'occupied').values()) \
                          if (occupancy < 1 and atomtype in atomtypes)]
    #not_occupied_nodes = [atomid for (atomid,occupancy) in nx.get_node_attributes(AA,'occupied').items() if (occupancy < 1)]
    print('found '+str(len(not_occupied_nodes))+'... ', end=" ", flush=True)
    H = AA.subgraph(not_occupied_nodes)
    print('done!', flush=True)

    print('setting up the search for isomorphisms... ', end=" ", flush=True)
    # This should contain occupied == 0!!!
    # GM = iso.GraphMatcher(AA,S,node_match=lambda a,b: a['type'] == b['type'])
    # nm = nx.algorithms.isomorphism.numerical_node_match(['type','q','occupied'],['type','q',0],atol=1e-4)
    nm = nx.algorithms.isomorphism.numerical_node_match(['type','q'],['type','q'],atol=1e-4)
    #nm = nx.algorithms.isomorphism.numerical_node_match(['type','occupied'],['type','occupied'],atol=1e-4)
    em = nx.algorithms.isomorphism.numerical_edge_match(['type'],['type'])
    GM = iso.GraphMatcher(H,S,nm,em)
    # ISMAGS seems to be slower...
    # GM = nx.algorithms.isomorphism.ISMAGS(AA,S,nm,em) 
    
    print('searching... ', end=" ", flush=True)
    # get Dict of IDs thta match subgraph
    BeadDict = list(GM.subgraph_isomorphisms_iter())
    # BeadDict = list(GM.match())
    # initilize List and convert Dict to List
    print('done!', flush=True)

    # 'atom_IDs' SHOULD NOT BE SORTED due to the mass thing in the mapping of VOTCA CSG xml file
    # or maybe if VOTCA_OUT==True do the for loops otherwise do the tuple set thingy
    print('sorting and removing doubles... ', end=" ", flush=True)

    # use set() and tuple to remove duplicate entries
    # non hashable -> convert to tuples and use set to remove duplicate netries
    # This does not work if VOTCA mapping file is needed.
    # Due to the sorting the weights are not in sequence anymore
    # BeadList = sorted(set(sorted(map(tuple,[sorted(list(i.keys())) for i in BeadDict]))))

    #interactive(currentframe(),locals())
    # this is not very nice and also slow
    # is there another solution?
    #BeadDict[idx]for idx in idxs
    #n = len(BeadDict)
    #checklist = [1] * n
    #for i in range(n):
    #    if checklist[i]:
    #        for j in range(i+1,n):
    #            if checklist[j]:
    #                if sorted(BeadDict[i].keys()) == sorted(BeadDict[j].keys()):
    #                    checklist[j] = 0
    #BeadList = [BeadDict[i] for i,j in enumerate(checklist) if j]

    # using numpy
    [mins,idxs] = np.unique([min(i.keys()) for i in BeadDict], return_index=True)
    BeadList = [BeadDict[idx] for idx in idxs]
    print('done!', flush=True)
                        
    
    count = 0
    #for i,AtomIDs in enumerate(BeadList[0:10]): print(i,j)
    #for i in range(len(BeadList)):
    for i,AtomIDs in enumerate(BeadList):
        #AtomIDs = BeadList[i]
        S = AA.subgraph(AtomIDs)
        if not (sum(dict(nx.get_node_attributes(S,'occupied')).values())):
            # if 1:
            #print('found:',AtomIDs)
            nx.set_node_attributes(S,1,'occupied')
            nx.set_node_attributes(S,BeadType,'bead_type')
            #nx.set_node_attributes(S,mean(AtomIDs),'bead_id')
            nx.set_node_attributes(S,min(AtomIDs),'bead_id')
            #nx.set_node_attributes(S,i+1,'bead_id')
            nx.set_node_attributes(S,BeadName,'bead_name')
            count+=1
        else:
            print('found, but not set:',AtomIDs)
            print('occupied',list(dict(nx.get_node_attributes(S,'occupied')).values()))
    print('found',len(BeadList),'Clusters and',count,'have been set', flush=True)
    #interactive(currentframe(),locals())
    
def search_bead_templates(AA,BT):
    # for each template
    for i in range(len(BT.data)):
        # search isomorphisms by atom ID ------------------------------------------
        if (BT.style[i] == 'id'):
            print()
            print('Searching: \''+BT.name[i]+'\' with type',BT.type[i],'by',BT.style[i],BT.data[i],'... ')
            find_and_mark_subgraph_by_ID(AA,BT.name[i],BT.type[i],BT.data[i])
            #print('FINDBEADTYPE'+str(BT.type[i])+': There are',len([[atomID,beadType] for atomID,beadType in AA.nodes('bead_type') if beadType==BT.type[i]]),'atoms of beadtype '+str(BT.type[i]))
            #with  open('AA_beadID'+str(i),'w') as f1, open('AA_beadType'+str(i),'w') as f2:
            #    for str1,str2 in AA.nodes('bead_id'): f1.write(str(str1)+'\t'+str(str2)+'\n')
            #    for str1,str2 in AA.nodes('bead_type'): f2.write(str(str1)+'\t'+str(str2)+'\n')
        # delete some nodes if requested ------------------------------------------
        elif (BT.style[i] == 'mol'):
            atomIDs = [atomID for atomID,molID in dict(nx.get_node_attributes(AA,'mol')).items() if molID in BT.data[i]]
            find_and_mark_subgraph_by_ID(AA,BT.name[i],BT.type[i],atomIDs)
        elif (BT.style[i] == 'del_template'):
            print()
            print('Deleting: \''+BT.name[i]+'\' with type',BT.type[i],'by',BT.style[i],len(BT.data[i]),'beads ... ')
            print('old count:',AA.number_of_nodes())
            find_and_mark_subgraph_by_ID(AA,BT.name[i],BT.type[i],BT.data[i])
            AA.remove_nodes_from([atomID for atomID,BeadType in dict(nx.get_node_attributes(AA,'bead_type')).items() if BeadType == BT.type[i]])
            print('new count:',AA.number_of_nodes())
        elif (BT.style[i] == 'del_id'):
            print()
            print('Deleting: \''+BT.name[i]+'\' with type',BT.type[i],'by',BT.style[i],len(BT.data[i]),'beads ... ')
            print('old count:',AA.number_of_nodes(),end=' ')
            AA.remove_nodes_from(BT.data[i])
            print('new count:',AA.number_of_nodes())
        elif (BT.style[i] == 'del_type'):
            print()
            print('Deleting: \''+BT.name[i]+'\' with type',BT.type[i],'by',BT.style[i],len(BT.data[i]),'beads ... ')
            AA.remove_nodes_from([atomID for atomID,atomType in dict(nx.get_node_attributes(AA,'type')).items() if atomType in BT.data[i]])
            print('new count:',AA.number_of_nodes())
        elif (BT.style[i] == 'del_molid'):
            print()
            print('Deleting: \''+BT.name[i]+'\' with type',BT.type[i],'by',BT.style[i],len(BT.data[i]),'beads ... ')
            AA.remove_nodes_from([atomID for atomID,atomType in dict(nx.get_node_attributes(AA,'mol')).items() if atomType in BT.data[i]])
            print('new count:',AA.number_of_nodes())
        #interactive(currentframe(),locals())

    print()
    for i in range(len(BT.data)):
        natoms = len([[atomID,beadType] \
                      for atomID,beadType in AA.nodes('bead_type') \
                      if beadType==BT.type[i]])
        print('SUMMARY'+str(BT.type[i])+': ',str(BT.name[i]),' There are',natoms,'atoms =',natoms/len(BT.data[i]),'beads of beadtype '+str(BT.type[i]))

    #interactive(currentframe(),locals())
    # check for unused atoms
    unoccupied_atoms = [[id,data] for id,data in AA.nodes(data=True) if data['occupied']<1]
    if  unoccupied_atoms:
        print('not all atoms got assigned!')
        print('Please check your bead template definitions')
        print('{:_^10}  {:_^10}  {:_^10}  {:_^10}'.format('ID','mol','type','charge'))
        for i in unoccupied_atoms:
            print('{:10}  {:10}  {:10}  {:10}'.format(i[0],i[1]['mol'],i[1]['type'],i[1]['q']))
    
    # double entries should be already removed in search_subgraphs
    # -> find_and_mark_subgraph_by_ID(AA,BT.name[i],BT.type[i],BT.data[i])
    # but there are probably gaps in the numbering as numbering is based on min(atomID) of each bead
    # therefor we want to renumber the bead IDs 
    atom2bead = dict(nx.get_node_attributes(AA,'bead_id'))
    [atom2bead[i] for i in range(1,50)]
    oldBeadIDlist = sorted(set(atom2bead.values()))
    if 0 in oldBeadIDlist:
        newBeadIDlist = [i for i,j in list(enumerate(oldBeadIDlist,0))]
    else:
        newBeadIDlist = [i for i,j in list(enumerate(oldBeadIDlist,1))]
    oldBead2newBead = dict(zip(oldBeadIDlist,newBeadIDlist))
    nx.set_node_attributes(AA,{i:oldBead2newBead[atom2bead[i]] for i in atom2bead},'bead_id')
    atom2bead = dict(nx.get_node_attributes(AA,'bead_id'))
    return
        
def mean(numbers, weights=None):
    if not weights: weights=[1.0] * len(numbers)
    weights =  list(map(abs,weights))
    return float(sum([i*j for i,j in zip(numbers,weights)])) / sum(weights)

def multiply(a,b):
    return [i*j for i,j in zip(a,b)]

def CGnodes(AA,BT):
    #interactive(currentframe(),locals())

    # dict bead2atom[4] returns a list of node IDs of AA that contribute to bead with ID 4
    atom2bead = dict(nx.get_node_attributes(AA,'bead_id'))

    # remove beads with beadID 0
    # and initilize arrays
    ilist = list(set(atom2bead.values()))
    n =  len(ilist)
    mlist = [0] * n
    tlist = [0] * n
    qlist = [0.0] * n
    nlist = [''] * n
    alist = [''] * n

    # bead to atom dict. first create with bead IDs then fill it from atom2bead
    bead2atom = {beadID:[] for beadID in ilist}
    for atomID,beadID in atom2bead.items():
        #if beadID != 0:
        bead2atom[beadID].append(atomID)

    #interactive(currentframe(),locals())
    type2name = dict({i:j for i,j in zip(BT.type,BT.name)})
    for i,beadID in enumerate(ilist):
        #print('running for',i,'with ID:',beadID)
        mlist[i] = AA.nodes[bead2atom[beadID][0]]['mol']
        tlist[i] = AA.nodes[bead2atom[beadID][0]]['bead_type']
        qlist[i] = 0
        for j in bead2atom[beadID]:
            qlist[i] += AA.nodes[j]['q']
        # tlist[i] = 0 are the atoms that are not assigned to any bead
        #if tlist[i]:
        # this should not be by type
        #nlist[i] = type2name[tlist[i]]
        nlist[i] = AA.nodes[bead2atom[beadID][0]]['bead_name']
        #else: nlist[i] = 0
        alist[i] = bead2atom[beadID]

    # create array for nodes
    CGlist = [[il,{'mol':ml,\
                  'type':tl,\
                  'q':ql,\
                  'x':0.0,\
                  'y':0.0,\
                  'z':0.0,\
                  'fx':0.0,\
                  'fy':0.0,\
                  'fz':0.0,\
                  'name':nl,\
                  'atom_IDs':al}] \
             for il,ml,tl,ql,nl,al in zip(ilist, mlist, tlist, qlist, nlist, alist)]

    # create and fill the graph
    CG = nx.Graph()
    CG.add_nodes_from(CGlist)
        
    # B networkx Graph for CG
    # atom2bead, bead2atom = Dictionaries for mapping
    return CG,atom2bead,bead2atom

def computeCGpositions(AA,atom2bead,bead2atom,CG,weighted_by='geo'):
    #interactive(currentframe(),locals())
    #bead=1
    #atoms = CG.nodes('atom_IDs')[bead]
    beadlist = [int(0)] * CG.number_of_nodes()
    x = [0.0] * CG.number_of_nodes()
    y = [0.0] * CG.number_of_nodes()
    z = [0.0] * CG.number_of_nodes()
    for i,(beadid,atoms) in enumerate(CG.nodes('atom_IDs')):
        # just necessary if unassigned beads
        if beadid > 0:
            beadlist[i] = beadid
            xdata = [AA.nodes('x')[atom] for atom in atoms]
            ydata = [AA.nodes('y')[atom] for atom in atoms]
            zdata = [AA.nodes('z')[atom] for atom in atoms]
            if weighted_by=='geo':
                x[i] = mean(xdata)
                y[i] = mean(ydata)
                z[i] = mean(zdata)

    nx.set_node_attributes(CG,{beadid:{'x':xi,'y':yi,'z':zi} \
                              for beadid,xi,yi,zi \
                              in zip(beadlist,x,y,z) })
    return

def CGbonds(CG,bead2atom,AA):
    bonds = list()
    # get CG bond topology from the FG network
    for beadID,atomIDs in bead2atom.items():
        bead1 = beadID
        boundary_atoms = list(nx.node_boundary(AA,atomIDs))
        #print("bead #",beadID," with atom IDs ", atomIDs," has boundary atoms :",boundary_atoms)
        # undirected Graph each bond show up 2 times
        for atom in boundary_atoms:
            bead2 = AA.nodes[atom]['bead_id']
            if bead1 < bead2:
                if CG.nodes[bead1]['type'] <= CG.nodes[bead2]['type']:
                    bonds.append([bead1,bead2,[str(CG.nodes[bead1]['type'])+'-'+str(CG.nodes[bead2]['type'])]])
                else:
                    bonds.append([bead1,bead2,[str(CG.nodes[bead2]['type'])+'-'+str(CG.nodes[bead1]['type'])]])
    type2comment = dict(enumerate(sorted(set([k[0] for i,j,k in bonds])),1))
    comment2type = {j: i for i,j in type2comment.items()}
    #interactive(currentframe(),locals())
    CGlist = [[b[0],b[1],{'id':n,'type':comment2type[b[2][0]],'comment':b[2][0]}] for n,b in enumerate(bonds,1)]
    CG.add_edges_from(CGlist)
    return

def find_angles(CG):
    angles = list()
    angle_center_atoms = [i for i,j in CG.degree if j >= 2]
    for a_center in angle_center_atoms:
        neighs = list(CG.neighbors(a_center))
        nneighs = len(neighs)
        for i in list(range(nneighs)):
            for j in range(i+1,nneighs):
                if CG.nodes[neighs[i]]['type'] > CG.nodes[neighs[j]]['type']:
                    j,i = i,j
                angles.append([neighs[i],a_center,neighs[j],{'id':0,'type':0,\
                    'comment':str(CG.nodes[neighs[i]]['type'])+'-'+\
                              str(CG.nodes[a_center]['type'])+'-'+\
                              str(CG.nodes[neighs[j]]['type'])}])

    type2comment = dict(enumerate(sorted(set([i[3]['comment'] for i in angles])),1))
    comment2type = {j: i for i,j in type2comment.items()}
    angles = [[i,j,k,{'id':n,'type':comment2type[d['comment']],'comment':d['comment']}] for n,[i,j,k,d] in enumerate(angles,1)]
    return angles

def find_dihedrals(CG):
    dihedrals = list()
    degree = CG.degree
    for i,j in CG.edges:
        if degree[i]>1 and degree[j]>1:
            neighs1 = list(CG.neighbors(i))
            neighs1.remove(j)
            neighs2 = list(CG.neighbors(j))
            neighs2.remove(i)
            #print("Dihedral from",neighs1,"via",i,j,"to",neighs2)
            for h in list(range(len(neighs1))):
                for k in list(range(len(neighs2))):
                    atom1 = neighs1[h]
                    atom2 = i
                    atom3 = j
                    atom4 = neighs2[k]
                    if CG.nodes[atom1]['type'] > CG.nodes[atom4]['type']:
                        atom1,atom2,atom3,atom4 = atom4,atom3,atom2,atom1
                    elif CG.nodes[atom1]['type'] == CG.nodes[atom4]['type'] and CG.nodes[atom2]['type'] > CG.nodes[atom3]['type']:
                        atom1,atom2,atom3,atom4 = atom4,atom3,atom2,atom1
                    dihedrals.append([atom1,atom2,atom3,atom4,\
                                      {'id':0,'type':0,'comment':str(CG.nodes[atom1]['type'])+'-'+\
                                                                 str(CG.nodes[atom2]['type'])+'-'+\
                                                                 str(CG.nodes[atom3]['type'])+'-'+\
                                                                 str(CG.nodes[atom4]['type'])}])

    type2comment = dict(enumerate(sorted(set([i[4]['comment'] for i in dihedrals])),1))
    comment2type = {j: i for i,j in type2comment.items()}

    dihedrals = [[i,j,k,l,{'id':n,'type':comment2type[d['comment']],'comment':d['comment']}] for n,[i,j,k,l,d] in enumerate(dihedrals,1)]
    return dihedrals

def find_impropers(CG):
    # tbd
    return

def graph2lammps(CG, box, outfile):
    #interactive(currentframe(),locals())
    print('generating Atoms section')
    Atoms =  [[node[0],\
               node[1]['mol'],\
               node[1]['type'],\
               node[1]['q'],\
               node[1]['x'],\
               node[1]['y'],\
               node[1]['z'],\
               '# '+str(node[1]['name'])] for node in CG.nodes(data=True)]
    # bond list
    print('generating Bonds section')
    Bonds = [[d['id'],d['type'],i,j,'# '+d['comment']]for i,j,d in CG.edges(data=True)]

    # angle list
    print('generating Angles section')
    Angles = [[d['id'],d['type'],i,j,k,'# '+d['comment']]for i,j,k,d in CG.angles]

    # dihedral list
    print('generating Dihedrals section')
    Dihedrals = [[d['id'],d['type'],i,j,k,l,'# '+d['comment']]for i,j,k,l,d in CG.dihedrals]

    # improper list 
    print('generating Impropers section: tbd')
    # tbd
    
    print('writing LAMMPS data file to',outfile)
    with open(outfile,'w') as f:
        f.write('# This file was created '+datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+\
            ' with MYPROG by '+os.environ['USER']) #+' on '+os.environ['NAME']+'\n')
        f.write('\n')
        f.write(str(len(Atoms))+' atoms\n')
        f.write(str(len(Bonds))+' bonds\n')
        f.write(str(len(Angles))+' angles\n')
        f.write(str(len(Dihedrals))+' dihedrals\n')
        f.write('\n')
        f.write(str(len(set([i[2] for i in Atoms])))+' atom types\n')
        f.write(str(len(set([i[1] for i in Bonds])))+' bond types\n')
        f.write(str(len(set([i[1] for i in Angles])))+' angle types\n')
        f.write(str(len(set([i[1] for i in Dihedrals])))+' dihedral types\n')
        f.write('\n')
        f.write(str(box['xlo xhi'][0])+' '+str(box['xlo xhi'][1])+' xlo xhi\n')
        f.write(str(box['ylo yhi'][0])+' '+str(box['ylo yhi'][1])+' ylo yhi\n')
        f.write(str(box['zlo zhi'][0])+' '+str(box['zlo zhi'][1])+' zlo zhi\n')
        f.write('\n')
        f.write('Atoms # full\n')
        f.write('\n')
        for i in Atoms:
            f.write('{:10d} {:10d} {:10d} {:10f} {:10f} {:10f} {:10f} {:s}\n'.format(*i))
        f.write('\n')
        f.write('Bonds\n')
        f.write('\n')
        for i in Bonds:
            f.write('{:10d} {:10d} {:10d} {:10d} {:s}\n'.format(*i))
        f.write('\n')
        f.write('Angles\n')
        f.write('\n')
        for i in Angles:
            f.write('{:10d} {:10d} {:10d} {:10d} {:10d} {:s}\n'.format(*i))
        f.write('\n')
        f.write('Dihedrals\n')
        f.write('\n')
        for i in Dihedrals:
            f.write('{:10d} {:10d} {:10d} {:10d} {:10d} {:10d} {:s}\n'.format(*i))

def graph2mscg(CG,BT,outfile):
    print('outfile is',outfile)
    # Maybe needs some fancy reordering
    # or should this be done after reading AA in general???
    #
    # At the moment it is required, that atoms appear in
    # an ordered way in the data file
    #
    # what if it is just one single large molecule
    # or not style 'full' with molid
    with open(outfile,'w') as f:
        print('do some fancy stuff here')
        f.write('cgsites '+str(nx.number_of_nodes(CG))+'\n')
        f.write('cgtypes '+str(len(set(nx.get_node_attributes(CG,'type').values())))+'\n')
        for i in BT.type:
            f.write(BT.name[i-1]+'\n')
        molecules = list(set(nx.get_node_attributes(CG,'mol').values()))
        n_molecules = len(molecules)
        bead2mol = nx.get_node_attributes(CG,'mol').items()
        f.write('moltypes '+str(n_molecules)+'\n')
        for molid in molecules:
            mol_atoms = [i for i,j in bead2mol if j==molid]
            S = CG.subgraph(mol_atoms)
            f.write('mol '+str(S.number_of_nodes())+' -1'+'\n')
            f.write('sitetypes\n')
            for i in nx.get_node_attributes(S,'type').values():
                f.write(str(i)+'\n')
            f.write('bonds '+str(S.number_of_edges())+'\n')
            for i,j in S.edges():
                f.write(str(i)+' '+str(j)+'\n')

    print('MSCG not properly done yet, TBD!')
            
    return 
    
def graph2votca(AA, BT, CG, type2mass, outfile):
    print('outfile is topology.xml')
    with open('topology.xml','w') as f:
        f.write('<topology>\n')
        mols = list(set(nx.get_node_attributes(AA,'mol').values()))
        f.write('    <molecules>\n')
        f.write('        <box xx=\"{:f}\" yy=\"{:f}\" zz=\"{:f}\"/>\n'.format((AA.box[0][1]-AA.box[0][0])/10,(AA.box[1][1]-AA.box[1][0])/10,(AA.box[2][1]-AA.box[2][0])/10))
        #interactive(currentframe(),locals())
        f.write('        <molecule name=\"m1AA\" nmols=\"1\" nbeads=\"'+str(AA.number_of_nodes())+'\" >\n')
        #interactive(currentframe(),locals())
        for mol in mols:
            AtomIDs = list(i for i,j in dict(nx.get_node_attributes(AA,'mol')).items() if j==mol)
            SG = AA.subgraph(AtomIDs)
            for atomid,data in SG.nodes(data=True):
                attr = list(data.keys())
                beadstr = '<bead name=\"'+str(atomid)+'\"'
                if 'type' in attr:
                    beadstr += ' type=\"'+str(data['type'])+'\"'
                    beadstr += ' mass=\"'+str(type2mass[data['type']])+'\"'
                if 'q' in attr:
                    beadstr += ' q=\"'+str(data['q'])+'\"'
                if 'mol' in attr:
                    beadstr += ' resid=\"'+str(data['mol'])+'\"'
                beadstr += ' />\n'
                f.write('           '+beadstr)
        f.write('        </molecule>\n\n')
        f.write('    </molecules>\n\n')
        #interactive(currentframe(),locals())    
        
        if hasattr(AA,'edges'):
            f.write('    <bonded>\n')
            bts = list(set(nx.get_edge_attributes(AA,'type').values()))
            for bt in bts:
                atomIDs = [[id1,id2] for id1,id2,data in AA.edges(data=True) if data['type']==bt]
                mystr = 'mol'
                f.write('        <bond>\n')
                f.write('        <name>bt'+str(bt)+'</name>\n')
                f.write('        <!-- atom types '+str(AA.nodes[atomIDs[0][0]]['type'])+'-'\
                                                  +str(AA.nodes[atomIDs[0][1]]['type'])+' -->\n')
                f.write('            <beads>\n')
                for id1,id2 in atomIDs:
                    #f.write('                   m'+str(AA.nodes[i][mystr])+':'+str(i)\
                    #                         +' m'+str(AA.nodes[j][mystr])+':'+str(j)+'\n')
                    f.write('               m1AA:'+str(id1)+' '+'m1AA:'+str(id2)+'\n')
                f.write('            </beads>\n')
                f.write('        </bond>\n\n')
            f.write('    </bonded>\n\n')
        f.write('</topology>\n')

    print('outfile is mapping.xml')
    with open('mapping.xml','w') as f:
        f.write('<cg_molecule>\n')
        f.write('    <name>m1CG</name>\n') # molecule name in cg representation
        f.write('    <ident>m1AA</ident>\n') # molecule name in atomistic topology
        f.write('\n')
        f.write('    <maps>\n')
        #interactive(currentframe(),locals())
        for bs,bn,bt,bd in zip(BT.style,BT.name,BT.type,BT.data):
            if bs.startswith('del'):
                continue
            f.write('        <map>\n')
            f.write('            <name>'+bn+'</name>\n')
            f.write('            <!-- bead type '+str(bt)+' -->\n')
            f.write('            <weights> ')
            mapstr = ''
            for atomID in bd:
                mapstr += ' ' + str(type2mass[AA.nodes[atomID]['type']])
            mapstr += ' </weights>\n'
            f.write(mapstr)
            f.write('        </map>\n')
        f.write('    </maps>\n')
        f.write('\n')

        f.write('    <topology>\n')
        f.write('        <cg_beads>\n')

        #interactive(currentframe(),locals())

        for beadid,data in CG.nodes(data=True):
            f.write('            <cg_bead>\n')
            f.write('                <name>'+str(beadid)+'</name>\n')
            f.write('                <type>'+str(data['type'])+'</type>\n')
            f.write('                <mapping>'+data['name']+'</mapping>\n')
            beadstr = ''
            for atomid in CG.nodes(data=True)[beadid]['atom_IDs']:
                beadstr += ' 1:m1AA:'+str(atomid)
            f.write('                <beads>'+beadstr+'</beads>\n')
            f.write('            </cg_bead>\n')
        f.write('        </cg_beads>\n\n')

        #interactive(currentframe(),locals())
        
        if hasattr(CG,'edges') or hasattr(CG,'angles') \
           or hasattr(CG,'dihedrals')  or hasattr(CG,'impropers'):
            f.write('\n')
            f.write('        <cg_bonded>\n')
            if hasattr(CG,'edges'):
                bond_types = list(set(nx.get_edge_attributes(CG,'type').values()))
                for bond_type in bond_types:
                    atomIDs = [[id1,id2] for id1,id2,data in CG.edges(data=True) if data['type']==bond_type ]
                    f.write('            <bond>\n')
                    f.write('                <name>BondType'+str(bond_type)+'</name>\n')
                    f.write('                <!-- bead types '+str(CG.nodes[atomIDs[0][0]]['type'])+'-'\
                                                              +str(CG.nodes[atomIDs[0][1]]['type'])+' -->\n')
                    f.write('                <beads>\n')
                    for id1,id2 in atomIDs:
                        f.write('                    '+str(id1)+' '+str(id2)+'\n')
                    f.write('                </beads>\n')    
                    f.write('            </bond>\n\n')
                    
            if hasattr(CG,'angles'):
                angle_types = list(set([data['type'] for id1,id2,id3,data in CG.angles]))
                for angle_type in angle_types:
                    atomIDs = [[id1,id2,id3] for id1,id2,id3,data in CG.angles if data['type'] == angle_type]
                    f.write('            <angle>\n')
                    f.write('                <name>AngleType'+str(angle_type)+'</name>\n')
                    f.write('                <!-- bead types '+str(CG.nodes[atomIDs[0][0]]['type'])+'-'\
                                                              +str(CG.nodes[atomIDs[0][1]]['type'])+'-'\
                                                              +str(CG.nodes[atomIDs[0][2]]['type'])+' -->\n')
                    f.write('                <beads>\n')
                    #for id1,id2,id3 in [[id1,id2,id3] \
                    #                    for id1,id2,id3,data in CG.angles \
                    #                    if data['type'] == angle_type]:
                    for id1,id2,id3 in atomIDs:
                        f.write('                    '+str(id1)+' '+str(id2)+' '+str(id3)+'\n')
                    f.write('                </beads>\n')    
                    f.write('            </angle>\n\n')
                    
            if hasattr(CG,'dihedrals'):
                dts = list(set([data['type'] for id1,id2,id3,id4,data in CG.dihedrals]))
                for dt in dts:
                    atomIDs = [[id1,id2,id3,id4] for id1,id2,id3,id4,data in CG.dihedrals if data['type'] == dt]
                    f.write('            <dihedral>\n')
                    f.write('                <name>DihedralType'+str(dt)+'</name>\n')
                    f.write('                <!-- bead types '+str(CG.nodes[atomIDs[0][0]]['type'])+'-'\
                                                              +str(CG.nodes[atomIDs[0][1]]['type'])+'-'\
                                                              +str(CG.nodes[atomIDs[0][2]]['type'])+'-'\
                                                              +str(CG.nodes[atomIDs[0][3]]['type'])+' -->\n')
                    f.write('                <beads>\n')
                    #for id1,id2,id3,id4 in [[id1,id2,id3,id4] \
                    #                        for id1,id2,id3,id4,data in CG.dihedrals \
                    #                        if data['type'] == dt]:
                    for id1,id2,id3,id4 in atomIDs:
                        f.write('                    '+str(id1)+' '+str(id2)+' '+str(id3)+' '+str(id4)+'\n')
                    f.write('                </beads>\n')    
                    f.write('            </dihedral>\n\n')
            if hasattr(CG,'impropers'):
                # do we need a list for different types?
                print('Sorry, this is not included yet!')
                print('please contact the developers.')
                while True:
                    user_inp = input('Do you want to continue? [y/n]: ')
                    if user_inp == 'y':   break
                    elif user_inp == 'n': sys.exit('exiting')
            f.write('        </cg_bonded>\n')

        f.write('    </topology>\n')
        f.write('\n')
        
        f.write('</cg_molecule>\n')
        print('do some fancy stuff here')
    return



def graph2gro(f,NX,wrapped=0):
    #interactive(currentframe(),locals())
    scale = 10
    
    lx = NX.box[0][1]-NX.box[0][0]
    ly = NX.box[1][1]-NX.box[1][0]
    lz = NX.box[2][1]-NX.box[2][0]
    f.write('{:}, t= {:.1f}\n'.format(os.path.basename(__file__),NX.ts))
    f.write('{:5d}\n'.format(NX.number_of_nodes()))
    for atomID,data in NX.nodes(data=True):
        if 'bead_name' in NX.nodes[atomID]: name = data['bead_name']
        else: name = data['name']
            
        #f.write('{:5d}{:.5s}{:5d}{:5d}{:8.3f}{:8.3f}{:8.3f}{:8.4f}{:8.4f}{:8.4f}\n'.format(\
        #    data['mol'], name, data['type'], atomID, \
        #    data['y']/scale,data['x']/scale, data['z']/scale, 0.0, 0.0, 0.0 ))

        x = data['x']
        y = data['y']
        z = data['z']

        if wrapped:
            while x > lx: x-=lx
            while x < 0:  x+=lx
            while y > ly: y-=ly
            while y < 0:  y+=ly
            while z > ly: z-=lz
            while z < 0:  z+=lz
        f.write('{:5d}{:.5s}{:5d}{:5d}{:8.3f}{:8.3f}{:8.3f}{:8.4f}{:8.4f}{:8.4f}\n'.format(\
            data['mol'], name, data['type'], atomID, \
            x/scale,y/scale, z/scale, 0.0, 0.0, 0.0 ))
        # velocities still to come
    f.write('   {:f} {:f} {:f}\n'.format((NX.box[0][1]-NX.box[0][0])/scale, \
                                         (NX.box[1][1]-NX.box[1][0])/scale, \
                                         (NX.box[2][1]-NX.box[2][0])/scale))
    return



def main():

    # check NetworkX version
    if float(nx.__version__) < 2.3:
        print('This programm was tested with NetworkX v2.3 but you have',nx.__version__,'installed.')
        while True:
            user_inp = input('Do you want to continue? [y/n]: ')
            if user_inp == 'y':   break
            elif user_inp == 'n': sys.exit('exiting')
    
    # initialize and read command line arguments ------------------------------
    options = initialize(sys.argv[1:])
    
    # import LAMMPS data file -------------------------------------------------
    print()
    print("##### reading LAMMPS data file:",options['indata'],'#####')
    d = lammps_data(options['indata'])
    print('Sections in data file:')
    for i in d.sections: 
        print(i)

    # read template file ------------------------------------------------------
    print()
    print('##### reading template file:',options['template'],'#####')
    print('reading template file:',options['template'])
    BT = bead_template()
    BT.readfile(options['template'])
    print('{:16} {:>5} {:>5} {:<10}'.format('NAME','TYPE','STYLE','DATA'))
    for i,j,k,l in zip(BT.name,BT.type,BT.style,BT.data):
        print('{:16} {:>5} {:>5}'.format(i,str(j),k),l)

    # networkx ----------------------------------------------------------------
    print()
    print('##### Initiliazing all atom graph and assigning topology #####')
    AA,type2mass = lammps2networkx(options,d,BT)

    # search subgraphs and do BT.style ----------------------------------------
    print()
    print('##### searching for subgraphs #####')
    search_bead_templates(AA,BT)

    # interactive(currentframe(),locals())
    # need to remove unset nodes
    # should this be here?
    AA.remove_nodes_from([i for i,j in dict(nx.get_node_attributes(AA,'occupied')).items() if j == 0])
    
    print()
    print('##### Initiliazing coarse grained graph and assigning topology #####')
    print('CG nodes... ', end=" ")
    # CG = coarse grained networkx
    # AA = all atom networkx
    # BT = BeadTemplate object
    # atom2bead, bead2atom = Dictionaries for mapping
    CG,atom2bead,bead2atom = CGnodes(AA,BT)
    print('positions... ', end=" ")
    computeCGpositions(AA,atom2bead,bead2atom,CG)
    print('done!')
    
    print('CG bonds... ', end=" ")
    # compute the bonds for each bead CGbonds(CG,bead2atom,AA)
    # CG: CG networkx
    # bead2atom: dict with beadID (CG): list of atomIDs (AA)
    # AA: all atom networkx
    CGbonds(CG,bead2atom,AA)
    print('done!')


    print('CG angles... ', end=" ")
    # compute the angles for each bead: find_angles(CG)
    # CG: CG networkx
    # output is of the format
    #     1  ,  2  ,  3  ,  4
    # [[atom1,atom2,atom3,{'id':1,'type':1,'comment':'1-2-3'}],
    #  [atom2,atom3,atom4,{'id':2,'type':2,'comment':'2-3-2'}],...]
    # and can be accessed  with e.g.
    # CG.angles[1][2] to get the third entry of line 2 (everything zero based)
    # or
    # CG.angles[0][3] to get the dictionary of the first angle
    # CG.angles[0][3]['type'] to get the angle type of the first angle
    CG.angles = find_angles(CG)
    print('done!')

    print('CG dihedrals... ', end=" ")
    # compute the dihedrals for each bead: find_dihedrals(CG)
    # CG: CG networkx
    # output is of the format
    #     1  ,  2  ,  3  ,  4   , 5
    # [[atom1,atom2,atom3,atom4,{'id':1,'type':1,'comment':'1-2-3'}],
    #  [atom2,atom3,atom4,atom5,{'id':2,'type':2,'comment':'2-3-2'}],...]
    # and can be accessed  with e.g.
    # CG.dihedrals[1][2] to get the third entry of line 2 (everything zero based)
    # or
    # CG.dihedrals[0][4] to get the dictionary of the first angle
    # CG.dihedrals[0][4]['type'] to get the angle type of the first angle
    CG.dihedrals = find_dihedrals(CG)
    print('done!')

    print('CG impropers... ', end=" ")
    print('tbd ### IF YOU HAVE A USE CASE LET ME KNOW ###')
    
    # write lammps of coarse grained Graph
    if 'outdata' in options.keys():
        print()
        print('##### writing LAMMPS data file',options['outdata'],'#####')
        graph2lammps(CG,d.headers,options['outdata'])

    #if 'MSCG' in options.keys():
    #    print()
    #    print('##### writing MSCG',options['MSCG'],'#####')
    #    graph2mscg(CG,BT,options['MSCG'])

    if 'VOTCA' in options.keys():
        print()
        print('##### writing VOTCA',options['VOTCA'],'#####')
        graph2votca(AA, BT, CG, type2mass, options['VOTCA'])
        
    if 'trajectory' in options.keys():
        print()
        print('##### Trajectory #####')
        print('Trajecory:',options['trajectory'])

        intraj = options['trajectory']
        outtraj = options['trajectory']+'.cg'
        with open(intraj, 'r') as fin, open(outtraj,'w') as fout:
             #, \
             #open('conf.unwrapped.gro','w') as f1,\
             #open('conf.wrapped.gro','w')   as f2:
            natoms = AA.number_of_nodes()
            line = 1
            while line:
                trajheader = [None] * 9
                for i in range(9): 
                    line = fin.readline().strip().split()
                    trajheader[i] = line
                # could be empty line(s)
                if True in [not i for i in trajheader]:
                    break
                timestep = int(trajheader[1][0])
                natoms = int(trajheader[3][0])
                nitems = len(trajheader[8][2:])
                box = [i[:] for i in [[None] * 2] * 3] 
                box[0] = [float(i) for i in trajheader[5]]
                box[1] = [float(i) for i in trajheader[6]]
                box[2] = [float(i) for i in trajheader[7]]
                CG.box = box
                CG.ts = timestep
                AA.box = box
                AA.ts = timestep
                boxl = [j-i for i,j in box]
                boxlhalf = [0.5*i for i in boxl]
                boxm = [0.5*(i+j) for i,j in box]
                bc = trajheader[4][3:]
                print('Reading timestep:',timestep,',',' '.join(trajheader[8][2:]))
                trajdict = {j:i for i,j in enumerate(trajheader[8][2:])}
                trajdata = [i[:] for i in [[None] * nitems] * natoms] 
                for i in range(natoms):
                    line = fin.readline().strip().split()
                    trajdata[i] = line
                
                if 'x' in trajdict.keys(): xstr = 'x'
                if 'xs' in trajdict.keys(): xstr = 'xs'
                if 'xu' in trajdict.keys(): xstr = 'xu'
                if 'y' in trajdict.keys(): ystr = 'y'
                if 'ys' in trajdict.keys(): ystr = 'ys'
                if 'yu' in trajdict.keys(): ystr = 'yu'
                if 'z' in trajdict.keys(): zstr = 'z'
                if 'zs' in trajdict.keys(): zstr = 'zs'
                if 'zu' in trajdict.keys(): zstr = 'zu'
                
                new_pos =  {int(i[trajdict['id']]):\
                              {'x':float(i[trajdict[xstr]]),\
                               'y':float(i[trajdict[ystr]]),\
                               'z':float(i[trajdict[zstr]]),\
                               'fx':float(i[trajdict['fx']]),\
                               'fy':float(i[trajdict['fy']]),\
                               'fz':float(i[trajdict['fz']])} for i in trajdata}
                
                nx.set_node_attributes(AA,new_pos)
                
                for beadID,beadData in CG.nodes(data=True):
                    S = AA.subgraph(beadData['atom_IDs'])
                    if False: #options['remap']=='True':
                        mapx=mapy=mapz=0
                        for i,j in S.edges():
                            dx = S.nodes[j]['x']-S.nodes[i]['x']
                            dy = S.nodes[j]['y']-S.nodes[i]['y']
                            dz = S.nodes[j]['z']-S.nodes[i]['z']
                            if abs(dx)>0.5*boxl[0]: mapx=1
                            if abs(dy)>0.5*boxl[1]: mapy=1
                            if abs(dz)>0.5*boxl[2]: mapz=1
                        if (mapx and xstr=='x')or (mapy and ystr=='y') or (mapz and zstr=='z'):
                            print(mapx,mapy,mapz)
                            [id,x,y,z]=list(zip(*[[i,j['x'],j['y'],j['z']] for i,j in S.nodes(data=True)]))
                            if (mapx and xstr=='x'):
                                nx.set_node_attributes(S,{i:{'x':j['x']-boxl[0]} for i,j in S.nodes(data=True) if j['x']>boxm[0]})
                            if (mapx and xstr=='y'):
                                nx.set_node_attributes(S,{i:{'y':j['y']-boxl[1]} for i,j in S.nodes(data=True) if j['y']>boxm[1]})
                            if (mapx and xstr=='z'):
                                nx.set_node_attributes(S,{i:{'z':j['z']-boxl[2]} for i,j in S.nodes(data=True) if j['z']>boxm[2]})

                                
                    [x,y,z,m,q,fx,fy,fz] = map(list,zip(*[[y['x'],y['y'],y['z'],type2mass[y['type']],y['q'],\
                        y['fx'],y['fy'],y['fz']] for x,y in S.nodes(data=True)]))
                        
                    if options['pos']=='geo':
                        x=mean(x)
                        y=mean(y)
                        z=mean(z)
                    if options['pos']=='com':
                        x=mean(x,m)
                        y=mean(y,m)
                        z=mean(z,m)
                    if options['pos']=='charge':
                        x=mean(x,q)
                        y=mean(y,q)
                        z=mean(z,q)
                    # image flags
                    ix = 0
                    iy = 0
                    iz = 0
                    #interactive(currentframe(),locals())
                    while x < box[0][0]:
                        x += boxl[0]
                        ix -= 1
                    while y < box[1][0]:
                        y += boxl[1]
                        iy -= 1
                    while z < box[2][0]:
                        z += boxl[2]
                        iz -= 1
                    while x > box[0][1]:
                        x -= boxl[0]
                        ix += 1
                    while y > box[1][1]:
                        y -= boxl[1]
                        iy += 1
                    while z > box[2][1]:
                        z -= boxl[2]
                        iz += 1
                        
                    # forces
                    fx = sum(fx)
                    fy = sum(fy)
                    fz = sum(fz)
                    BS = CG.subgraph(beadID)
                    nx.set_node_attributes(BS,{beadID:{'x':x,'y':y,'z':z,\
                                                       'fx':fx,'fy':fy,'fz':fz,\
                                                       'ix':ix,'iy':iy,'iz':iz}})
                    
                for i in range(len(trajheader[:-1])):
                    if i == 3:
                        fout.write(str(CG.number_of_nodes())+'\n')
                    else:
                        fout.write(' '.join(trajheader[i])+'\n')
                fout.write('ITEM: ATOMS id type x y z fx fy fz ix iy iz\n')
                for i,j in CG.nodes(data=True): 
                    # print(i,j['mol'],j['type'],j['q'],j['x'],j['y'],j['z'])
                    fout.write('{:d} {:d} {:10f} {:10f} {:10f} {:10f} {:10f} {:10f} {:d} {:d} {:d}\n'.format(i,j['type'],j['x'],j['y'],j['z'],j['fx'],j['fy'],j['fz'],j['ix'],j['iy'],j['iz']))

                #interactive(currentframe(),locals())
                #graph2gro(f1,AA,0)
                #graph2gro(f2,AA,1)
                    
    return

if __name__ == "__main__":
    main()
