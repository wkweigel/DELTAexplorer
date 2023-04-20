
#The contents of this 


#####################################
###### D E P E N D A N C I E S ######
#####################################

import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import networkx as nx
import matplotlib.pylab as plt
import numpy as np
from pyvis.network import Network
from PIL import Image
from collections import defaultdict
from queue import Queue
import string

ELEMENTS = 'ABCDEFGHIJKLMNOP'
LINKERS = "abcdefghijklmnop"
ELEMENTS_LIST= list(ELEMENTS)
LINKERS_LIST= list(LINKERS)


###########################
###### C L A S S E S ######
###########################

#Class for interpreting and working with DELTA topologies
class SimpleTopology:
    '''
    A class for processing and working with DELTA strings. 
    '''
    def __init__(self, string):
        self.DELTAstring=string
        self.Nodes=self.Get_Nodes(self.DELTAstring, add_dna=True)
        self.Nodes_noDNA=self.Nodes[:-1]
        self.Edges=[]
        self.Get_Edges(self.Nodes)

    def Get_Nodes(self, topology_string, add_dna=True):  
        '''
        Parses a DELTA string into and ordered list of individual topology items (nodes).
        '''
        node_list=[]
        cycle_start_search=True
        cycle_end_search=False
        char_behind=None
        char_ahead=None
        two_char_behind=None
        two_char_ahead=None
        topology_string=list(topology_string)

        for idx, char in enumerate(topology_string):

            #defines peek variables (characters before and after current idx value)
            if idx in range(1, len(topology_string)):
                char_behind=topology_string[idx-1]

            if idx in range(2, len(topology_string)):
                two_char_behind=(topology_string[idx-2])+(char_behind)

            if idx in range(0, len(topology_string)-1):
                char_ahead=topology_string[idx+1]

            if idx in range(0, len(topology_string)-2):
                two_char_ahead=(char_ahead)+(topology_string[idx+2])
            
            # This For loop skips over syntactical topology characters ( "!", "(", and ")" ) and considers them prospectively or retrospectively 
            if char =='!':
                continue
            if char == '(' or char == ')':
                continue
            if char =='-':
                node_list.append('-DNA')
                break
            
            
            # The rest of the For loop identies elements based on leading and trailing syntactical characters
            #Identifies the start of a branched-cyclic element
            if cycle_start_search is True and two_char_behind =='(!':
                cycle_start_search=False
                cycle_end_search=True
                node_list.append('(!' + char + ')')
                continue
            #identifies the start of a cyclic element 
            elif cycle_start_search is True and char_behind == '!':
                cycle_start_search=False
                cycle_end_search=True
                node_list.append('!' + char)
                continue
            #identifies the end of a branched-cyclic element
            elif cycle_end_search is True and two_char_ahead =='!)':
                cycle_end_search=False
                node_list.append('(' + char + '!)')
                continue
            #identifies the end of a cyclic element 
            elif cycle_end_search is True and char_ahead =='!':
                cycle_end_search=False
                node_list.append(char + '!')
                continue
            #identifies a branching element 
            elif char_behind =='(' and char_ahead == ')':
                node_list.append('(' + char + ')')
                continue

            elif char in LINKERS or char in ELEMENTS:
                node_list.append(char)
        if add_dna is True:
            node_list.append("DNA")
        return(node_list)

    def Add_Edge(self, E1, E2): #Adds edge for the two input elements (ie nodes)
        self.Edges.append([E1,E2])

    def Get_Edges(self, node_list): #Finds the edges for a list of ordered topology elements
        """
        Finds the edges for a list of ordered topology elements.
        """
        #This function operates by using a set of nested For loops. 
        #The first loop (item1) iterates through the list of topology items and determines thier type (regular, branched, cyclic).
        #Next, the second For loop (item2) then iterates through the list again and adds edges between item1 and item2 based on the specfic DELTA rules
        edges=[]
        for idxItem1, item1 in enumerate(node_list): #Classifies item1 using an outer For loop. 
            regularItem1=False
            cy_start_item1=False
            cy_end_item1=False
            br_item1=False
            br_cy_start_item1=False
            br_cy_end_item1=False

            if len(item1)==1: #Flags Regular Elements
                regularItem1=True
            if len(item1)==2: #Flags starting Cyclic elements (ie !A)
                if item1[0]=='!':
                    cyc_start=item1
                    cy_start_item1=True
                if item1[1]=="!":
                    cy_end_item1=True
            if len(item1)==3: #Flags Branched elements
                if item1!='DNA':
                    br_item1=True
                    br_point=node_list[idxItem1-1] # Defines the branch point as the element before the Branch element
            if len(item1)==4: #Flags Branched-Cyclic elements.
                if item1[1]=='!': #If Branched-Cyclic element is the start of a cycle
                    br_cyc_start=item1
                    br_cy_start_item1=True
                    br_cyc_point=node_list[idxItem1-1] # Defines the branch point as the element before the Branched-Cyclic element
                if item1[2]== '!': #If Branched-Cyclic element is the end of cycle
                    br_cy_end_item1=True
                    br_cyc_end=item1
                    br_cyc_point=node_list[idxItem1-1] # Defines the branch point as the element before the Branched-Cyclic element
            
            #Classify B
            for idxB, item2 in enumerate(node_list):
                cy_item2=False
                br_cy_item2=False

                if item1=='DNA' and item2=='DNA':
                    if '(' in node_list[(idxB-1)]:
                        self.Add_Edge(node_list[idxB-2],item2)
                        break
                    else:
                        self.Add_Edge(node_list[idxB-1],item2)
                        break

                if item1==item2:
                    continue

                if idxB == idxItem1+1: #If the B is the element following A
                    if regularItem1 is True: #And if A is a regular element (not branched)
                        if node_list[idxItem1+1]=='DNA':#Continue on if the next A element is DNA (since DNA is handled with its own seperate rule above)
                            continue
                        else: #Add edge between consecutive elements A and B
                            self.Add_Edge(item1,item2)
                    if cy_start_item1 or cy_end_item1 is True:
                        if node_list[idxItem1+1]=='DNA':#Continue on if the next A element is DNA (since DNA is handled with its own seperate rule above)
                            continue
                        else: #Add edge between consecutive elements !A and B
                            self.Add_Edge(item1,item2)

                if len(item2)==2: #If B is a cyclic element
                    if item2[1]=='!': #If B is the end point of a cycle
                        cyc_end=item2 #Identify the the end point
                        cy_item2=True #flag the end point as found
                    
                if len(item2)==4: #If B is a branched cyclic element
                    if item2[2] =='!': #If B is the end point of a branched cycle
                        br_cyc_end=item2 #Identify the the end point
                        br_cy_item2=True #flag the end point as found

                if br_item1 is True: #Resolves the branching edge for a regular branching element
                    if idxB == idxItem1+1:
                        self.Add_Edge(br_point,item2)

                if br_cy_start_item1 is True: #Resolves the branching edge for a starting branched cyclic element (the cyclic edge is resolved later down)
                    if idxB == idxItem1+1:
                        self.Add_Edge(br_cyc_point,item2)

                if br_cy_end_item1 is True: #Resolves the branching edge for an ending branched cyclic element (the cyclic edge is resolved later down)
                    if idxB == idxItem1+1:
                        if item2!='DNA':
                            self.Add_Edge(br_cyc_point,item2)
                
                if cy_start_item1 is True: #Resolves cases involving regular cyclic starting points
                    if cy_item2 is True:
                        self.Add_Edge(cyc_start,cyc_end)
                    if br_cy_item2 is True:
                        self.Add_Edge(cyc_start,br_cyc_end)
                
                if br_cy_start_item1 is True: #Resolves cases involving branched cyclic starting points
                    if cy_item2 is True:
                        self.Add_Edge(br_cyc_start,cyc_end)
                    if br_cy_item2 is True:
                        self.Add_Edge(br_cyc_start,br_cyc_end)
        return(self.Edges)


class Topology:
    '''
    A class for processing and working with DELTA strings. 
    '''
    def __init__(self, string):
        self.DELTAstring=string #define the DELTA string
        self.Nodes=self.Get_Nodes(self.DELTAstring, add_dna=True) #Parse string into a list of nodes and adds DNA to end of the list if not already present.
        self.Nodes_noDNA=self.Nodes[:-1] #Creates a node list without DNA
        self.Edges=[]
        self.Get_Edges(self.Nodes) #Produce an edge list from the list of nodes
        self.BondDict=self.Get_BondDict(self.Edges) #Create a bond dictionary with bond counts to each node in the node list

      
        '''
        if '!' not in self.DELTAstring:
            ##self.CycleList=self.Get_Cycles(self.Nodes)
            self.CycleList=self.Get_Cycles(self.Nodes)
        else:
            ##acyclic_list=[item.replace('!', '') for item in self.Nodes]
            acyclic_list=[item.replace('!', '') for item in self.Nodes]
            self.CycleList=self.Get_Cycles(acyclic_list)
        '''
        #self.LinkersList=self.Get_Linkers()
    

    def Prune_All_Syntax(self, list):
        '''
        Removes all syntactical characters and replaces any lowercase letters from a list of DELTA nodes.
        '''
        prune_characters = ['!', '(', ')','-']
        pruned_list=[]
        for item in list:
            for char in prune_characters:
                item=item.replace(char, '')
            pruned_list.append(item)
        return(pruned_list.upper)

    
    def Prune_Syntax(self, list):
        '''
        Removes all syntactical characters from a list of DELTA nodes.
        '''
        prune_characters = ['!', '(', ')','-']
        pruned_list=[]
        for item in list:
            for char in prune_characters:
                item=item.replace(char, '')
            pruned_list.append(item)
        return(pruned_list)

    def Get_Nodes(self, topology_string, add_dna=True):  
        '''
        Parses a DELTA string into and ordered list of individual topology items (nodes).
        '''
        node_list=[]
        cycle_start_search=True
        cycle_end_search=False
        char_behind=None
        char_ahead=None
        two_char_behind=None
        two_char_ahead=None
        topology_string=list(topology_string)

        for idx, char in enumerate(topology_string):

            #defines peek variables (characters before and after current idx value)
            if idx in range(1, len(topology_string)):
                char_behind=topology_string[idx-1]

            if idx in range(2, len(topology_string)):
                two_char_behind=(topology_string[idx-2])+(char_behind)

            if idx in range(0, len(topology_string)-1):
                char_ahead=topology_string[idx+1]

            if idx in range(0, len(topology_string)-2):
                two_char_ahead=(char_ahead)+(topology_string[idx+2])
            
            # This For loop skips over syntactical topology characters ( "!", "(", and ")" ) and considers them prospectively or retrospectively 
            if char =='!':
                continue
            if char == '(' or char == ')':
                continue
            if char =='-':
                node_list.append('-DNA')
                break
            
            # The rest of the For loop identies elements based on leading and trailing syntactical characters
            #Identifies the start of a branched-cyclic element
            if cycle_start_search is True and two_char_behind =='(!':
                cycle_start_search=False
                cycle_end_search=True
                node_list.append('(!' + char + ')')
                continue
            #identifies the start of a cyclic element 
            elif cycle_start_search is True and char_behind == '!':
                cycle_start_search=False
                cycle_end_search=True
                node_list.append('!' + char)
                continue
            #identifies the end of a branched-cyclic element
            elif cycle_end_search is True and two_char_ahead =='!)':
                cycle_end_search=False
                node_list.append('(' + char + '!)')
                continue
            #identifies the end of a cyclic element 
            elif cycle_end_search is True and char_ahead =='!':
                cycle_end_search=False
                node_list.append(char + '!')
                continue
            #identifies a branching element 
            elif char_behind =='(' and char_ahead == ')':
                node_list.append('(' + char + ')')
                continue

            elif char in LINKERS or char in ELEMENTS:
                node_list.append(char)

        if add_dna is True:
            if "DNA" not in self.DELTAstring:
                node_list.append("-DNA")
        return(node_list)

    def Add_Edge(self, E1, E2): #Adds edge for the two input elements (ie nodes)
        self.Edges.append([E1,E2])

    def Get_Edges(self, node_list): #Finds the edges for a list of ordered topology elements
        """
        Finds the edges for a list of ordered topology elements.
        """
        #This function operates by using a set of nested For loops. 
        #The first loop (item1) iterates through the list of topology items and determines thier type (regular, branched, cyclic).
        #Next, the second For loop (item2) then iterates through the list again and adds edges between item1 and item2 based on the specfic DELTA rules
        edges=[]
        for idxItem1, item1 in enumerate(node_list): #Classifies item1 using an outer For loop. 
            regularItem1=False
            cy_start_item1=False
            cy_end_item1=False
            br_item1=False
            br_cy_start_item1=False
            br_cy_end_item1=False

            if len(item1)==1: #Flags Regular Elements
                regularItem1=True
            if len(item1)==2: #Flags starting Cyclic elements (ie !A)
                if item1[0]=='!':
                    cyc_start=item1
                    cy_start_item1=True
                if item1[1]=="!":
                    cy_end_item1=True
            if len(item1)==3: #Flags Branched elements
                if item1!='DNA':
                    br_item1=True
                    br_point=node_list[idxItem1-1] # Defines the branch point as the element before the Branch element
            if len(item1)==4: #Flags Branched-Cyclic elements.
                if item1[1]=='!': #If Branched-Cyclic element is the start of a cycle
                    br_cyc_start=item1
                    br_cy_start_item1=True
                    br_cyc_point=node_list[idxItem1-1] # Defines the branch point as the element before the Branched-Cyclic element
                if item1[2]== '!': #If Branched-Cyclic element is the end of cycle
                    br_cy_end_item1=True
                    br_cyc_end=item1
                    br_cyc_point=node_list[idxItem1-1] # Defines the branch point as the element before the Branched-Cyclic element
            
            #Classify B
            for idxB, item2 in enumerate(node_list):
                cy_item2=False
                br_cy_item2=False

                if item1=='DNA' and item2=='DNA':
                    if '(' in node_list[(idxB-1)]:
                        self.Add_Edge(node_list[idxB-2],item2)
                        break
                    else:
                        self.Add_Edge(node_list[idxB-1],item2)
                        break

                if item1==item2:
                    continue

                if idxB == idxItem1+1: #If the B is the element following A
                    if regularItem1 is True: #And if A is a regular element (not branched)
                        if node_list[idxItem1+1]=='DNA':#Continue on if the next A element is DNA (since DNA is handled with its own seperate rule above)
                            continue
                        else: #Add edge between consecutive elements A and B
                            self.Add_Edge(item1,item2)
                    if cy_start_item1 or cy_end_item1 is True:
                        if node_list[idxItem1+1]=='DNA':#Continue on if the next A element is DNA (since DNA is handled with its own seperate rule above)
                            continue
                        else: #Add edge between consecutive elements !A and B
                            self.Add_Edge(item1,item2)

                if len(item2)==2: #If B is a cyclic element
                    if item2[1]=='!': #If B is the end point of a cycle
                        cyc_end=item2 #Identify the the end point
                        cy_item2=True #flag the end point as found
                    
                if len(item2)==4: #If B is a branched cyclic element
                    if item2[2] =='!': #If B is the end point of a branched cycle
                        br_cyc_end=item2 #Identify the the end point
                        br_cy_item2=True #flag the end point as found

                if br_item1 is True: #Resolves the branching edge for a regular branching element
                    if idxB == idxItem1+1:
                        self.Add_Edge(br_point,item2)

                if br_cy_start_item1 is True: #Resolves the branching edge for a starting branched cyclic element (the cyclic edge is resolved later down)
                    if idxB == idxItem1+1:
                        self.Add_Edge(br_cyc_point,item2)

                if br_cy_end_item1 is True: #Resolves the branching edge for an ending branched cyclic element (the cyclic edge is resolved later down)
                    if idxB == idxItem1+1:
                        if item2!='DNA':
                            self.Add_Edge(br_cyc_point,item2)
                
                if cy_start_item1 is True: #Resolves cases involving regular cyclic starting points
                    if cy_item2 is True:
                        self.Add_Edge(cyc_start,cyc_end)
                    if br_cy_item2 is True:
                        self.Add_Edge(cyc_start,br_cyc_end)
                
                if br_cy_start_item1 is True: #Resolves cases involving branched cyclic starting points
                    if cy_item2 is True:
                        self.Add_Edge(br_cyc_start,cyc_end)
                    if br_cy_item2 is True:
                        self.Add_Edge(br_cyc_start,br_cyc_end)
        return(self.Edges)

    def Construct_Cyclic_Variants(self):
        self.CycleList=self.Get_Cycles(self.Nodes)
        return(self.CycleList)
    def Construct_Cyclic_Variants_with_DNA(self):
        self.CycleList=self.Get_Cycles(self.Nodes)
        return(self.CycleList)

    def Get_Cycles(self, node_list):
        '''
        Finds all valid cyclic permutations of an ACYCLIC topology list or string.
        '''
        cyclic_list=[]
        branch_span=False 
        #Outer for-loop (IdxA) controls all cyclic starting points
        for IdxA, item in enumerate(node_list): 
            temp_listA=node_list.copy() #Copies the toplist as a temp list
            inital_end_point=IdxA+2 #The minimum number of elements needed to complete a cycle is three, therefore the inital endpoint is started at IdxA+2
            if inital_end_point==len(node_list): #Checks to see if the initial endpoint had reached the end of the elements list
                break
            if len(node_list)<3:#Checks if the elements list contains a minimum of 3 terms. By definition a cycle cannot exist for less than 3 terms.
                break
            if '(' in item:#Checks if element is branched and if so, adds a cyclic-branched element in its place
                temp_listA[IdxA]="(!" + item[1] + ")"
            else:
                temp_listA[IdxA]="!" + item #Adds a cyclic element in place of current element
                start_point=IdxA
                for e in range(0, ((len(node_list))-1)):
                    if '(' in node_list[start_point+1]:
                        branch_span=True
            #Inner for-loop (IdxB) controls all cyclic ending points
            for IdxB in range(inital_end_point, len(node_list)): #Considers the range between the inital cycle point and the end of elements list (all cycle endpoints necessarily occur in this range)
                temp_listB=temp_listA.copy()
                temp_cycle=''
                if branch_span is True:
                    branch_span=False
                    continue
                if '(' in node_list[IdxB]: # For handling branch points, ie (X)
                    if branch_span is True:
                        continue
                    for c in node_list[IdxB]:
                        if c == '(' or c == ')':
                            continue
                        else:
                            temp_listB[IdxB]="(" + c + "!)"
                else: # For handling regular elements
                    temp_listB[IdxB]=node_list[IdxB] + "!" 
                for e in temp_listB:
                    temp_cycle=temp_cycle + e
                cyclic_list.append(temp_cycle)
        return(cyclic_list)


    def Get_BondDict(self, edges):
        '''
        Creates a dictionary of nodes and bonds in the form {node:#bonds,...}
        '''
        bond_dict={}

        for A, B in edges:
            if A in bond_dict:
                bond_dict[A]+=1
            if A not in bond_dict:
                bond_dict[A]=1
            if B in bond_dict:
                bond_dict[B]+=1
            if B not in bond_dict:
                bond_dict[B]=1
        return(bond_dict)

    def Get_Linkers(self):
        global bonds
        global linker_list
        
        def Return_String(list):
            string=''
            for element in list:
                if element=='DNA':
                    break
                else:
                    string = string + element
            return(string)

        #returns a variant linker string from a nodelist using the specified index value  
        def Make_Linker_Variant(source_list, replacement_index):
            temp_list=source_list.copy()
            if temp_list[replacement_index] != '-DNA':
                temp_list[replacement_index]=source_list[replacement_index].lower()
                temp_string=Return_String(temp_list)
            return(temp_string, temp_list)

        #Checks to see if the element preceeding the DNA element is branched (looks for a '(' in the last element)
        def Terminal_Variant(elements_list):
            if '(' not in elements_list[-1]: 
                if 'DNA' not in elements_list[-1]:
                            temp_list_a=elements_list.copy() #Create temp copy of the elements list
                            temp_list_a[-1]=elements_list[-1].lower() #Create a linker variant using the final element
                            temp_string_a=Return_String(temp_list_a) #Convert the temp element list to a concatonated string
                            linker_list.append(temp_string_a) #Add the string to the linker list

        # Checks if a given element (and corresonding index) in a given list is a linker
        def Linker_Search(current_list, search_element, search_idx):
            double_skip=False
            single_skip=False
            end_search=False

            # Checks if search element has more than one edge 
            # If it does, replace the element with a lowercase variant
            # Creates string from list and adds string to linkers list
            if bonds.get(search_element)>1: 
                var_str, var_list=Make_Linker_Variant(current_list, search_idx) #creates a variant with the current element replaced with a lowercase char
                linker_list.append(var_str)# Adds variant to list of linkers

                
                #Valency differentiation now performed in a seperate dedicaded function
                """
                if bonds.get(search_element)==2:
                    if bonds.get(search_element)!=3:  
                        divalent_linker_list.append(var_str)
                if bonds.get(search_element)==3:
                    trivalent_linker_list.append(var_str)
                """

            # Determines where the next potential linker postion could be (how many elements should be skipped)
            # if the next element is branched, perform a double skip. If not, perform a single skip.
            # Breaks loop if there is an IndexError (ie the loop is on the last element of the list)  
                try: 
                    if '(' in current_list[search_idx+1]:
                        double_skip=True
                    else:
                        single_skip=True
                except IndexError:
                    end_search=True
                    
            # Sets the approriate subsequent start index
            # If index error is thrown, flag end_search as true
                if double_skip is True: #Defines new index start position or flags the end_search if index does not exist
                    try:
                        start_idx=search_idx+3 #incremented by 3 because linking elements cannot be branched or adjacent to eachother.
                        skip_num=3
                    except IndexError:
                        end_search=True
                        
                if single_skip is True: #Defines new index start position or flags the end_search if index does not exist
                    try:
                        start_idx=search_idx+2 #incremented by 2 because two linking elements cannot be adjacent.
                        skip_num=2
                    except IndexError:
                        end_search=True
                
                if end_search is True: #re-initialize the starting index and the skip number 
                    start_idx=0
                    skip_num=0
                return(var_list, start_idx, end_search, skip_num)

        #List containing all the topologies with linker elements
        linker_list=[]
        #List containing only topologies with divalent linker elements
        
        #Counts bonds (ie connections) where DNA is included as an element
        edges=self.Nodes
        bonds=self.BondDict

        # Main (ie outer) for loop
        # The secondary (ie inner) for loop considers all accetable linker permutations that come after the current iteration in the main loop
        for main_idx, main_element in enumerate(self.Nodes):

            if '(' in main_element: #Skips branch elements. These by definition cannot be linkers.
                continue

            if bonds.get(main_element)==1: #Skips cases where the "A" element has only one connection (ie "A" is not part of a cyclic connection)
                continue
                
            if main_element=='A':
                continue

            # Creates parent_list (containing only a single linker element) that is used in the inner for loop
            #The next secondary for loop identifies additional child linker permutations 
            parent_list, initial_idx, end_search, skip_num = Linker_Search(self.Nodes, main_element, main_idx)

            if end_search is True: # Exits loop if end_search is returned as true
                break

        # Secondary for loop does not have a full iterative structure
        # Iteration of this loop is controlled by initial index manipulation
        # Considers the elements in the parent list starting from the first acceptable index postion after the parent linker's postion
            for parent_idx, parent_element in enumerate(parent_list[initial_idx:], start=initial_idx): 
                try:
                    child_list, child_idx, end_search, skip_num = Linker_Search(parent_list, parent_element, parent_idx)
                except TypeError:
                    pass

                if end_search is True:
                    break
                if initial_idx != len(self.Nodes)-1: #If the initial index (the starting idx for secondary loop) is not the second to last element
                    if parent_idx != len(self.Nodes)-2:
                        try:
                            Terminal_Variant(child_list)
                        except UnboundLocalError:
                            pass

                initial_idx=initial_idx+skip_num
        return(linker_list)
    
    
#Class for performing graphical node, edge, and path operations
class Graph:
    '''
    A class for processing and working with network graphs. 
    '''
    def __init__(self, node_list, edge_list):
        self.edge_list=edge_list
        self.node_list=node_list
        self.adjacency_dict=self.construct_adjacency_dict(self.edge_list)
        self.G=nx.from_edgelist(self.edge_list)

    def construct_adjacency_dict(self, edge_list):
        '''
        Create an adjacency dictionary from an edge list.
        '''
        graph = defaultdict(list)
        for edge in edge_list:
            a, b = edge[0], edge[1]
            graph[a].append(b)
            graph[b].append(a)
        return(graph)
    
    def calculate_distances(self, source):
        '''
        Calculate distances from the specified node to all other nodes in the graph.
        '''
        Q = Queue()
        distance_dict = {k: 999999999 for k in self.adjacency_dict.keys()}
        visited_vertices = list()
        Q.put(source)
        visited_vertices.append(source)
        while not Q.empty():
            vertex = Q.get()
            if vertex == source:
                distance_dict[vertex] = 0
            for u in self.adjacency_dict[vertex]:
                if u not in visited_vertices:
                    # update the distance
                    if distance_dict[u] > distance_dict[vertex] + 1:
                        distance_dict[u] = distance_dict[vertex] + 1
                    Q.put(u)
                    visited_vertices.append(u)
        return distance_dict
        
    def calculate_specific_path(self,a,b):
        paths=[]
        temp_paths=nx.all_simple_paths(self.G, source=a, target=b, cutoff=len(self.node_list))
        paths.extend(temp_paths)
        self.specific_paths=paths
        return(paths)

    def get_longest_specific_path(self,a,b):
        self.calculate_specific_path(a,b)
        self.max_paths=[]
        list_len = [len(i) for i in self.specific_paths]
        for path in self.specific_paths:
            if len(path)==max(list_len):
                self.max_paths.append(path)
                self.longest_path=self.max_paths[0] #only take the first path (if multiple equal-length paths exist)
        return(self.longest_path)
    
    def calculate_all_paths(self):
        paths=[]
        for node in self.node_list:
            temp_paths=nx.all_simple_paths(self.G, source=node, target="DNA", cutoff=len(self.node_list))
            paths.extend(temp_paths)
        self.all_paths=paths
        
    def calculate_all_edge_paths(self):
        edge_paths=[]
        for node in self.node_list:
            temp_edge_paths=nx.all_simple_edge_paths(self.G, source=node, target="DNA", cutoff=len(self.node_list))
            edge_paths.extend(temp_edge_paths)
            self.all_edge_paths=edge_paths

    def reverse_tuple_list(self, tuple_list):
        self.reversed_tuples=[]
        for edge in tuple_list:
            a,b=edge[0],edge[1]
            self.reversed_tuples.append((b,a))
        return(self.reversed_tuples)

    def get_longest_path(self):
        self.calculate_all_paths()
        self.max_paths=[]
        list_len = [len(i) for i in self.all_paths]
        for path in self.all_paths:
            if len(path)==max(list_len):
                self.max_paths.append(path)
                self.longest_path=self.max_paths[0] #only take the first path (if multiple equal-length paths exist)
        return(self.longest_path)
    
    def get_all_longest_paths(self):
        self.calculate_all_paths()
        self.max_paths=[]
        list_len = [len(i) for i in self.all_paths]
        for path in self.all_paths:
            if len(path)==max(list_len):
                self.max_paths.append(path)
        return(self.max_paths)

    def get_all_longest_edge_paths(self):
        self.calculate_all_edge_paths()
        self.max_edge_paths=[]
        list_len = [len(i) for i in self.all_edge_paths]
        for path in self.all_edge_paths:
            if len(path)==max(list_len):
                self.max_edge_paths.append(path)
        return(self.max_edge_paths)
    
    def identify_true_branch_nodes_and_edges(self):
        self.get_longest_path()
        self.true_branch_nodes=[]
        self.true_branch_edges=[]
        for node in self.node_list: #gets a list of the nodes that do not appear in the main path
            if node not in self.longest_path:
                self.true_branch_nodes.append(node)
        if self.true_branch_nodes==None:
            self.true_branch_nodes=None
            self.true_branch_edges=None
        for node in self.true_branch_nodes:  #gets a list of the edges that do not appear in the main path
            for edge in self.edge_list:
                a, b = edge[0], edge[1]
                if node in a or node in b:
                    self.true_branch_edges.append((a,b))
            return(self.true_branch_nodes,self.true_branch_edges)
    
    def perform_longest_path_corrections(self):
        for main_path, main_edge_path in zip(self.max_paths,self.max_edge_paths):
            perform_branch_correction=True

            #identify true branch nodes (nodes that are not part of the longest path)
            #if no branches are found, skip furthur branch correction steps
            try:
                true_branch_nodes,true_branch_edges=self.identify_true_branch_nodes_and_edges()
            except(TypeError):
                perform_branch_correction=False
            
            #Apply
            #this section is skipped if the main path no longer uses any branch points
            if perform_branch_correction is True:
                for branch_node, branch_edge in zip(true_branch_nodes,true_branch_edges):
                    a,b=branch_edge[0],branch_edge[1]
                    for node, edge in zip(main_path, main_edge_path):
                        if a==node:
                            attachment_index=main_path.index(node)
                            main_path.insert(attachment_index+1,'('+ branch_node +')')
                            main_edge_path.insert(attachment_index+1, branch_edge)

            #for use in a complete edge search, calculate the reverse edges for the main path
            reverse_main_edge_path=self.reverse_tuple_list(main_edge_path)

            #convert the edge_list (list of lists --> list of tuples)
            pruned_tuples=[tuple(l) for l in self.edge_list]

            #Create list of tuple edge-pairs (regular edge, reverse edge )
            tuple_pairs=[]
            for regular_tuple, reverse_tuple in zip(main_edge_path, reverse_main_edge_path):
                tuple_pairs.append([regular_tuple,reverse_tuple])
            
            #identify the edge used for cycle closure
            for edge1 in pruned_tuples:
                search_counter=0
                for pair in tuple_pairs:
                    if edge1 not in pair:
                        search_counter+=1
                if search_counter==len(tuple_pairs):
                    cyclic_edge=edge1

            #find the start and end indicies for the ring closure within the main path
            start_point_search=True
            for node in main_path:
                if node in cyclic_edge:
                    if start_point_search is True:
                        start_index=main_path.index(node)
                        main_edge_path.insert(start_index+1,cyclic_edge)
                        start_point_search=False
                        continue
                    if start_point_search is False:
                        end_index=main_path.index(node)
                        start_point_search=True
            

            #insert the cyclic edge into the main path
            try:
                main_edge_path.insert(start_index+1,cyclic_edge)
            except (UnboundLocalError):
                continue
            #reverse_main_edge_path.insert(start_index+1,cyclic_edge)

            #apply ! to the start and end cycle nodes
            try:
                main_path[start_index]='!'+main_path[start_index] 
            except (UnboundLocalError):
                continue
            try:
                main_path[end_index]=main_path[end_index]+'!'
            except (UnboundLocalError):
                continue

            alphabet_map=dict(zip(string.ascii_uppercase, range(0,int(len(main_path)-1))))
            alphabet_map['DNA']='DNA'

            self.corrected_list=[]
            corrected_node=''
            for node_idx,node in enumerate(main_path):
                corrected_node=''
                if node=="DNA":
                    self.corrected_list.append('-DNA')
                    break
                for char_idx, char in enumerate(node):
                    if char.isalpha():
                        corrected_char= list(filter(lambda x: alphabet_map[x] == node_idx, alphabet_map))[0]
                        if char.islower():
                            corrected_char=corrected_char.lower()
                        corrected_node=corrected_node+corrected_char
                    else:
                        corrected_node=corrected_node+char
                self.corrected_list.append(corrected_node)
 

            self.corrected_string=''
            for node in self.corrected_list:
                self.corrected_string=self.corrected_string+node
        try:
            return(self.corrected_list, self.corrected_string)
        except(AttributeError):
            return(['A','A','A','A','A','A','A'], 'NoChange')
        
        
    def perform_longest_path_check(self):
        for main_path, main_edge_path in zip(self.max_paths,self.max_edge_paths):
            perform_branch_correction=True

            #identify true branch nodes (nodes that are not part of the longest path)
            #if no branches are found, skip furthur branch correction steps
            try:
                true_branch_nodes,true_branch_edges=self.identify_true_branch_nodes_and_edges()
            except(TypeError):
                perform_branch_correction=False
            
            #Apply
            #this section is skipped if the main path no longer uses any branch points
            if perform_branch_correction is True:
                for branch_node, branch_edge in zip(true_branch_nodes,true_branch_edges):
                    a,b=branch_edge[0],branch_edge[1]
                    for node, edge in zip(main_path, main_edge_path):
                        if a==node:
                            attachment_index=main_path.index(node)
                            main_path.insert(attachment_index+1,'('+ branch_node +')')
                            main_edge_path.insert(attachment_index+1, branch_edge)

            #for use in a complete edge search, calculate the reverse edges for the main path
            reverse_main_edge_path=self.reverse_tuple_list(main_edge_path)

            #convert the edge_list (list of lists --> list of tuples)
            pruned_tuples=[tuple(l) for l in self.edge_list]

            #Create list of tuple edge-pairs (regular edge, reverse edge )
            tuple_pairs=[]
            for regular_tuple, reverse_tuple in zip(main_edge_path, reverse_main_edge_path):
                tuple_pairs.append([regular_tuple,reverse_tuple])
            
            #identify the edge used for cycle closure
            for edge1 in pruned_tuples:
                search_counter=0
                for pair in tuple_pairs:
                    if edge1 not in pair:
                        search_counter+=1
                if search_counter==len(tuple_pairs):
                    cyclic_edge=edge1

            #find the start and end indicies for the ring closure within the main path
            start_point_search=True
            for node in main_path:
                if node in cyclic_edge:
                    if start_point_search is True:
                        start_index=main_path.index(node)
                        main_edge_path.insert(start_index+1,cyclic_edge)
                        start_point_search=False
                        continue
                    if start_point_search is False:
                        end_index=main_path.index(node)
                        start_point_search=True
            

            #insert the cyclic edge into the main path
            try:
                main_edge_path.insert(start_index+1,cyclic_edge)
            except (UnboundLocalError):
                continue
            #reverse_main_edge_path.insert(start_index+1,cyclic_edge)

            #apply ! to the start and end cycle nodes
            try:
                main_path[start_index]='!'+main_path[start_index] 
            except (UnboundLocalError):
                continue
            try:
                main_path[end_index]=main_path[end_index]+'!'
            except (UnboundLocalError):
                continue

            alphabet_map=dict(zip(string.ascii_uppercase, range(0,int(len(main_path)-1))))
            alphabet_map['DNA']='DNA'


            self.corrected_list=[]
            corrected_node=''
            for node_idx,node in enumerate(main_path):
                corrected_node=''
                if node=="DNA":
                    self.corrected_list.append('-DNA')
                    break
                for char_idx, char in enumerate(node):
                    if char.isalpha():
                        corrected_char= list(filter(lambda x: alphabet_map[x] == node_idx, alphabet_map))[0]
                        if char.islower():
                            corrected_char=corrected_char.lower()
                        corrected_node=corrected_node+corrected_char
                    else:
                        corrected_node=corrected_node+char

                self.corrected_list.append(corrected_node)


            topologycheck=Topology(self.corrected_string)
            bondcheck=True
            for node, bonds in topologycheck.BondDict.items():
                if bonds>3:
                    bondcheck=False
                
            if bondcheck is False:
                return(['A','A','A','A','A','A','A'], 'NoChange')
         
            self.corrected_string=''
            for node in self.corrected_list:
                self.corrected_string=self.corrected_string+node
        try:
            return(self.corrected_list, self.corrected_string)
        except(AttributeError):
            return(['A','A','A','A','A','A','A'], 'NoChange')
            

#Class for separating and holding linkers in diffrent valency groups
class Linker:
    def __init__(self):
        self.linker_counter=0
        self.divalent_counter=0
        self.trivalent_counter=0
    
    divalent_list=[]
    trivalent_list=[]





###############################
###### F U N C T I O N S ######
###############################


def Parse_Topology_String(top_string, add_dna=True):  
    '''
    Parse an input DELTA string into a list of individual elements with retained syntactical characters
    '''
    global elements_list
    global linkers_list
    global top_list

    elements_list= list(ELEMENTS)
    linkers_list= list(LINKERS)
    top_list=[]
    C1_search=True
    C2_search=False
    char_behind=None
    char_ahead=None
    two_char_behind=None
    two_char_ahead=None
    top_string=list(top_string)
    
    for idx, char in enumerate(top_string):

        #defines peek variables (characters before and after current idx value)
        if idx in range(1, len(top_string)):
            char_behind=top_string[idx-1]

        if idx in range(2, len(top_string)):
            two_char_behind=(top_string[idx-2])+(char_behind)

        if idx in range(0, len(top_string)-1):
            char_ahead=top_string[idx+1]

        if idx in range(0, len(top_string)-2):
            two_char_ahead=(char_ahead)+(top_string[idx+2])
        
        # For loop skips over syntactical topology characters ("!", "(", and ")")
        if char =='!':
            continue
        if char == '(' or char == ')':
            continue
        
        # The rest of the For loop identies elements based on leading and trailing syntactical characters
        #Identifies the start of a branched-cyclic element
        if C1_search is True and two_char_behind =='(!':
            C1_search=False
            C2_search=True
            top_list.append('(!' + char + ')')
            continue
        #identifies the start of a cyclic element 
        elif C1_search is True and char_behind == '!':
            C1_search=False
            C2_search=True
            top_list.append('!' + char)
            continue
        #identifies the end of a branched-cyclic element
        elif C2_search is True and two_char_ahead =='!)':
            C2_search=False
            top_list.append('(' + char + '!)')
            continue
        #identifies the end of a cyclic element 
        elif C2_search is True and char_ahead =='!':
            C2_search=False
            top_list.append(char + '!')
            continue
        #identifies a branching element 
        elif char_behind =='(' and char_ahead == ')':
            top_list.append('(' + char + ')')
            continue

        elif char in LINKERS or char in ELEMENTS:
            top_list.append(char)
    if add_dna is True:
        top_list.append("-DNA")
    return(top_list)

def Construct_Main_Tree(base_elements):
    '''
    The external main tree growth function. Generates the main tree nodes that act as parents for the child cyclic and linker permutations.
    '''
    main_tree_node_list=['A']
    growth_control={'A':'inactive', 'A':'active'}
    main_tree_branches_list=[['A', 'AB']]
    ###temp_cyclic_string=''
    ###variant_nodes=[]
    ###variant_branches=[]

    ###selected_nodes=[]
    ###selected_branches=[]

    global cyclic_list
    cyclic_list=[]
    #For each letter in the sequence list, consider all possible branched and cyclic permutations
    #Add the permutations to a tree of nodes according to complexity (number of DEs in the node)
    for letter in base_elements: 
        for current_node,value in list(growth_control.items()): #Iteratively goes through the growth control dictionary
            
            if value == 'active': #Looks for any active keys (ie nodes) in the growth control dict. 
                new_node=str(current_node)+letter #Create a new node by appending the existing active node with an unbranched version of the current letter from the list
                main_tree_node_list.append(new_node)#Add the new node to the node list
                growth_control[new_node]='active'#Add the new node to the growth control dict as active
                growth_control[current_node]='inactive'#Sets the current growth control value to inactive
                if current_node=='A':
                    continue #Skips the branch addition step for the 'A' node since the (A, AB) branch is initialized in the branch list
                
                else:
                    main_tree_branches_list.append([current_node,new_node]) #Add a branch (ie edge) between the current key (ie node) and the new_node

                if '(' not in current_node: #Adds current letter to current node as branch element if the current node has no banching elements
                    new_br_node=str(current_node) + '(' + letter + ')'
                    main_tree_node_list.append(new_br_node)
                    growth_control[new_br_node]='active'
                    main_tree_branches_list.append([current_node,new_br_node])

                if '(' in current_node:#For nodes with existing branches, Adds current letter to current node as branch element if the current node has an unbranched final element
                    templist=Parse_Topology_String(current_node,add_dna=False)
                    if '(' not in templist[-1]:
                        new_br_node=str(current_node) + '(' + letter + ')'
                        main_tree_node_list.append(new_br_node)
                        growth_control[new_br_node]='active'
                        main_tree_branches_list.append([current_node,new_br_node])

    return(main_tree_node_list, main_tree_branches_list)

def Construct_Linker_Groups(MainTreeNodes):
    '''
    Uses the bond order dict to assign linker groups to linker permutations based on the valency of the linker
    '''
    FullLinkerList=[]
    for Node in MainTreeNodes:
        temp_Topology=Topology(Node)
        all_linkers = temp_Topology.Get_Linkers()
        FullLinkerList=FullLinkerList + all_linkers

    for linker in FullLinkerList:
        linker_topology=Topology(linker)
        temp_linker=Linker()

        full_topology=linker_topology.Nodes
        pruned_topology=linker_topology.Prune_Syntax(full_topology)
        for idx, element in enumerate(pruned_topology):
            if element in LINKERS: # apply valencey flags for any linkers in the pruned node list
                if linker_topology.BondDict.get(full_topology[idx])==2:
                    temp_linker.divalent_counter+=1
                if linker_topology.BondDict.get(full_topology[idx])==3:
                    temp_linker.trivalent_counter+=1
        
        #use valencey flags to assign node to appropriate linker group
        if temp_linker.divalent_counter>0 and temp_linker.trivalent_counter==0:
            temp_linker.divalent_list.append(linker)
        if temp_linker.divalent_counter==0 and temp_linker.trivalent_counter>0:
            temp_linker.trivalent_list.append(linker)

    return(FullLinkerList, temp_linker.divalent_list, temp_linker.trivalent_list)

def Map_Child_Main_Nodes(child_list, parent_list):
    '''
    Converts a node to its main node variant and returns the edge between the them.
    '''
    mapped_edges=[]
    for child_node in child_list:
        pre_map_node=child_node.upper()
        if '!' in pre_map_node:
                map_node=pre_map_node.replace('!', '')
        else:
            map_node=pre_map_node
        if map_node in parent_list:
            mapped_edges.append([map_node,child_node])
    return(mapped_edges)

def remove_duplicates(node_list):
    set_nodes = set(node_list) 
    output_list = (list(set_nodes))
    return(output_list)

def Construct_Graph(nodes, edges):
    '''
    Creates lists for the colors and shaped for each element node in the tree
    '''
    node_colors=[]
    for element in nodes:
        if element == 'DNA':
            node_colors.append('#012A4A')
            continue
        for char in element:
            if char not in ELEMENTS: 
                if char not in LINKERS and element != 'DNA':
                    continue
            else: 
                if 'A' in element:
                    node_colors.append('#A9D6E5')
                    break
                if 'B' in element:
                    node_colors.append('#61A5C2')
                    break
                if 'C' in element:
                    node_colors.append('#468FAF')
                    break
                if 'D' in element:
                    node_colors.append('#2A6F97')
                    break
                if 'E' in element:
                    node_colors.append('#014F86')
                    break
                if 'F' in element:
                    node_colors.append('#01497C')
                    break
                if 'G' in element:
                    node_colors.append('#013A63')
                    break
                if 'H' in element:
                    node_colors.append('#013A63')
                    break
            if char in LINKERS and element != 'DNA': 
                node_colors.append('#F88379')

    node_shapes=[]
    for element in nodes:
        if element == 'DNA':
            node_shapes.append('dot')
        for char in element:
            if element=='DNA':
                continue
            if char in LINKERS:
                node_shapes.append('hexagon')
            if char in ELEMENTS:
                node_shapes.append('dot')
            else:
                continue 
    return(node_colors, node_shapes)

def find(string, char):
    for i, lett in enumerate(string):
        if lett == char:
            yield i

def get_matrix(string1):
    topology1=Topology(string1)
    G1=nx.from_edgelist(topology1.Edges)
    array1=nx.to_numpy_array(G1, nodelist=topology1.Nodes)
    return(array1)

def check_isomorphism(string):
    input_array=get_matrix(string)
    corrected_string=correct_cyclic_DELTA_topology(string)
    corrected_array=get_matrix(corrected_string)
    if np.array_equal(input_array,corrected_array) is False:
        return(False,corrected_string)
    else:
        return(True,string)


def check_string(string):
    error_flag=False
    messages=[]
    #Cycle Checks
    if '!' in string:
        cycle_positions=list(find(string,'!'))
        #Multiple cycle check
        if string.count('!')>2:
            messages.append("Cycle Count Error: Delta currently only interprets strings with a single ring closure.")
            error_flag=True
        if len(cycle_positions)==1 :
            messages.append("Cycle Point Error: Both a start point and end point for a ring closure must be specified.")
            error_flag=True
        
        if len(cycle_positions)<1:
            if cycle_positions[1]-cycle_positions[0]<4:
                messages.append("Cycle Size Error: A ring must contain a miniumum of three diversity elements.")
                error_flag=True

        iso_check,corrected_string=check_isomorphism(string)
        if iso_check is False:
            messages.append('A longer main chain was detected.')
            messages.append('Try using: ' + str(corrected_string[:-4]))
            error_flag=True
            
    if '(' in string:
        br_start_positions=list(find(string,'('))
        br_end_positions=list(find(string,')'))
        br_positions=list(zip(br_start_positions, br_end_positions))
        #Adjacent branch element check
        if string.count(')(')>0:
            messages.append("Branching Element Error: At least one non-branching element must seperate branching elements.")
            error_flag=True
        for start, end in br_positions:
            if end-start>2:
                if string[(end-1)] != '!':
                    messages.append("Branching Element Error: Delta currently only interprets branches containing a single diversity element.")
                    error_flag=True
        for start in br_start_positions:
            if string[start+1].islower():
                messages.append("Branching Element Warning: A linker should not occupy a branching element position.")
                error_flag=True
            
    return(error_flag,messages)

def reverse_edge_list(edge_list):
    reversed_edges=[]
    for edge in edge_list:
        a,b=edge[0],edge[1]
        reversed_edges.append([b,a])
    return(reversed_edges)

def reverse_tuple_list(tuple_list):
    reversed_edges=[]
    for edge in tuple_list:
        a,b=edge[0],edge[1]
        reversed_edges.append((b,a))
    return(reversed_edges)

def Prune_Syntax_From_Edges(list):
    '''
    Removes all syntactical characters and replaces any lowercase letters from an DETLA edge list.
    '''
    prune_characters = ['!', '(', ')','-']
    pruned_list=[]
    for edge in list:
        a, b = edge[0], edge[1]
        for char in a:
            if char in prune_characters:
                a=a.replace(char, '')
        for char in b:
            if char in prune_characters:
                b=b.replace(char, '')

        pruned_list.append([a,b])
    return(pruned_list)

def construct_adjacency_dict(edge_list):
    '''
    Create an adjacency dictionary from an edge list.
    '''
    graph = defaultdict(list)
    for edge in edge_list:
        a, b = edge[0], edge[1]
        graph[a].append(b)
        graph[b].append(a)
    return(graph)

def calculate_distances(input_graph, source):
    '''
    Calculate distances from the specified node to all other nodes in the specifed graph.
    '''
    Q = Queue()
    distance_dict = {k: 999999999 for k in input_graph.keys()}
    visited_vertices = list()
    Q.put(source)
    visited_vertices.append(source)
    while not Q.empty():
        vertex = Q.get()
        if vertex == source:
            distance_dict[vertex] = 0
        for u in input_graph[vertex]:
            if u not in visited_vertices:
                # update the distance
                if distance_dict[u] > distance_dict[vertex] + 1:
                    distance_dict[u] = distance_dict[vertex] + 1
                Q.put(u)
                visited_vertices.append(u)
    return distance_dict

def correct_cyclic_DELTA_topology(node):
    cycle_topology=Topology(node)
    pruned_edges=Prune_Syntax_From_Edges(cycle_topology.Edges)
    pruned_nodes=cycle_topology.Prune_Syntax(cycle_topology.Nodes)
    cycle_graph=Graph(pruned_nodes,pruned_edges)
    reverse_edge_list(pruned_edges)
    cycle_graph.get_all_longest_paths()
    cycle_graph.get_all_longest_edge_paths()
    correct_list,correct_string=cycle_graph.perform_longest_path_corrections()
    if correct_string=='NoChange':
        correct_string=node
    return(correct_string)

def Add_DNA(node_list):
    DNA_node_list=[]
    for node in node_list:
        DNA_node_list.append(node+'-DNA')
    return(DNA_node_list)

def Remove_DNA(node_list):
    noDNA_node_list=[]
    for node in node_list:
        if '-DNA' in node:
            noDNA_node_list.append(node.replace('-DNA',''))
        else:
            noDNA_node_list.append(node)
    return(noDNA_node_list)

def Get_Main_Cycles(main_node_list):
    cyclic_variations=[]
    corrected_cycles=[]
    #Get all the cyclic permutations of the Main DNA Nodes
    for main_node in main_node_list:
        topology=Topology(main_node)
        pruned_edges=Prune_Syntax_From_Edges(topology.Edges)
        pruned_nodes=topology.Prune_Syntax(topology.Nodes)
        graph=Graph(pruned_nodes,pruned_edges)
        if (len(topology.Nodes_noDNA))>=3: #remove topologies that are too smal to support ring closures
            possible_cycle_points=[]

            #Get possible cycle points based on bond count
            for node, bonds in topology.BondDict.items():
                if bonds <=2:
                    possible_cycle_points.append(node)
            possible_cycle_points.pop(-1) #DNA is never a cycle point, so remove it.

            #Create syntatically pruned combinatorial pairs of all possible closure points 
            cycle_pairs = Prune_Syntax_From_Edges([(a, b) for idx, a in enumerate(possible_cycle_points) for b in possible_cycle_points[idx + 1:]])
            
            # Create a list only the possible cycle pairs based on ring size .
            # The minimum longest path must be 3 or greater. 
            possible_cycle_pairs=[]
            for pair in cycle_pairs:
                a,b = pair[0], pair[1]
                max_path=graph.get_longest_specific_path(a,b)
                if len(max_path)>=3:
                    possible_cycle_pairs.append((a,b))
            
            # For each pair of cycle points, add the appropriate "!" syntax to the acyclic nodes.
            # Concatonate the nodes into the corresponding cyclic variant string.
            for cycle_pair in possible_cycle_pairs:
                cyclic_variant=''
                start_node, end_node= cycle_pair[0], cycle_pair[1]
                for DE_node in topology.Nodes[:-1]:
                    if start_node in DE_node:
                        cycle_start=DE_node.replace(start_node,'!'+start_node)
                        cyclic_variant=cyclic_variant+cycle_start
                        continue
                    if end_node in DE_node:
                        cycle_end=DE_node.replace(end_node, end_node+'!')
                        cyclic_variant=cyclic_variant+cycle_end
                        continue
                    else:
                        cyclic_variant=cyclic_variant+DE_node
                cyclic_variant=cyclic_variant+'-DNA'      
                cyclic_variations.append(cyclic_variant)

    for cycle in cyclic_variations:
        corrected_cycle_string=correct_cyclic_DELTA_topology(cycle)
        corrected_cycles.append(corrected_cycle_string)

    return(corrected_cycles)

#Count the upper case letters (Diversity Elements) in a topology string
def count_uppercase_letters(string):
    uppercase_counter=0
    for char in string:
        if char.isupper():
            uppercase_counter+=1
    return(uppercase_counter)

def filter_for_n_DEs(DE_num, edge_list):
    PassingNodes=[]
    PassingEdges=[]
    for edge in edge_list:
        node1DEs=count_uppercase_letters(edge[0])
        node2DEs=count_uppercase_letters(edge[1])
        if node1DEs<=DE_num and node2DEs <=DE_num:
            PassingEdges.append(edge)
        if edge[0].isupper() and node2DEs <=DE_num:
            PassingEdges.append(edge)

        if edge[1].isupper() and node1DEs <=DE_num:
            PassingEdges.append(edge)

        if edge[0].isupper() and edge[1].isupper():
            PassingEdges.append(edge)

    for edge in PassingEdges:
        if edge[0] not in PassingNodes:
            PassingNodes.append(edge[0])
        if edge[1] not in PassingNodes:
            PassingNodes.append(edge[1])
    return(PassingNodes, PassingEdges)

def filter_nodes_for_n_DEs(DE_num, node_list):
    PassingNodes=[]
    for node in node_list:
        DE_count=count_uppercase_letters(node[:-4])
        if DE_count<=DE_num:
            PassingNodes.append(node)
    return(PassingNodes)

###############################
###### S T R E A M L I T ######
###############################

st. set_page_config(layout="wide")

# Set header title
st.title('DEL Topology Visualization')

#Define the topologies that populate the node selection dropdown
DefaultNodes=['ABC', 'AbC', '!ABC!', 'Ab(C)', 'A!bCD!', 'A!b(C)DE!', 'Specify Other' ]

#Create the node selection dropdown
node_selection= st.sidebar.selectbox('Choose a node to inspect:', DefaultNodes)

#Conditional text entry field if "Specify Other" is selected
if node_selection=='Specify Other':
    SpecifiedNode=st.sidebar.text_input('Specify Topology (up to 8 elements):', value='ABC')

#Create topology class object for the specified topology
if node_selection=='Specify Other':
    custom_user_node=SimpleTopology(SpecifiedNode)
    Nodes=custom_user_node.Nodes
    Edges=custom_user_node.Edges
else:
    user_node=SimpleTopology(node_selection)
    Nodes=user_node.Nodes
    Edges=user_node.Edges

if node_selection=='Specify Other':
    error=False
    error, error_messages=check_string(SpecifiedNode)
    if error is True:
        for message in error_messages:
            st.sidebar.write(message)  

#create the lists for the colors and shapes of the nodes and edges
colors, shapes= Construct_Graph(Nodes,Edges)

#create the pyvis networkx object to represent the user DEL input
DEL= Network()

#load the nodes, shapes, and colors for the DEL into the neworkx object
DEL.add_nodes(Nodes, color=colors, shape=shapes)

#load the edges for the DEL into the neworkx object
DEL.add_edges(Edges)

# Toggle the permutation viewer on or off (DEFULT index=1)
view_permutation_options = st.sidebar.radio('Show Permutation Explorer Options',('Yes', 'No'), index=1)

#Conditional section only used if the permutation view options are toggled on by the user
if view_permutation_options == 'Yes':
    #Define number of DEL elements from streamlit user input (DEFULT index=3)
    DEL_size= int(st.sidebar.selectbox('Number of structural elements:',('3', '4', '5', '6', '7', '8'), index=3))

    #Define number of diversity elements from streamlit user input (DEFULT index=2)
    Diversity_Elements=int(st.sidebar.selectbox('Number of diversity elements:',('2','3','4','5','6'), index=2))

    # Toggle consideration of cyclic topologies from streamlit user input (DEFULT index=1)
    cycle_check = st.sidebar.radio('Consider cyclic permuations?',('Yes', 'No'), index=0)

    # Toggle display of linker permutations from streamlit user input (DEFULT index=1)
    linker_check = st.sidebar.radio('Show linker permutations?', ('Yes', 'No'), index=0)

    #Control the linker grouping to visualize (DEFULT index=2)
    if linker_check == 'Yes':
        linker_group=st.sidebar.radio('Select linker type to show:', ('All', 'Divalent', 'Trivalent'), index=2)

    # Initialization value for scaling node size by literature precedence
    lit_scale='No'

    # Option for scaling node size by literature precedence appears if cycle check is "Yes"(DEFULT index=1)
    if cycle_check == 'Yes':
        lit_scale= st.sidebar.radio('Scale by literature prevalance?',('Yes', 'No'), index=1)

    # Allows user to switch between Dendridic or Hierarchical tree layout.
    view= st.sidebar.radio ('Tree Layout:', ('Dendridic', 'Hierarchical'))

    #Define base string that gets sliced based on the number of elements the user selects
    DE_string='ABCDEFGHIJK'

    #Slice the string to inculde only the number of DEs from the user input
    DE_selection=DE_string[:DEL_size]
    
    #Slice off the first DE (this is because the first node of the toplogy tree "A", contains no possible edges to branched nodes. "AB" is the only possible extention)
    sequence=DE_selection[1:]

    #convert the sliced sequence to a list of elements
    seq_list=list(sequence)

    #Topology class object creation for the sliced sequence
    #Step is no longer needed?
    #user_sequence=Topology(seq_list)

    #construct the list of nodes and edges for the main topology tree
    MainTreeNodes, MainTreeBranches= Construct_Main_Tree(seq_list)

    #Add the DNA element to all the nodes for determining the linker and cycle permuatations
    DNA_MainTreeNodes=Add_DNA(MainTreeNodes)

    #Get the cyclic permuations for the nodes in the main tree 
    if cycle_check=='Yes':
        CyclicVariants=Get_Main_Cycles(DNA_MainTreeNodes)
        CyclicVariants=remove_duplicates(CyclicVariants)
   
    #Get the linker permuations for the nodes in the main tree
    #Compile the list of child variant nodes (linkers only)
    if linker_check == 'Yes' and cycle_check=='No':
        AllMainLinkers, MainDivalentLinkers, MainTrivalentLinkers=Construct_Linker_Groups(DNA_MainTreeNodes)
        if linker_group=='All':
            linker_selection=AllMainLinkers
        if linker_group=='Divalent':
            linker_selection=MainDivalentLinkers
        if linker_group=='Trivalent':
            linker_selection=MainTrivalentLinkers
        linker_selection=remove_duplicates(linker_selection)
        ChildVariants=linker_selection

    #Get the linker and the cyclic permuations for the nodes in the main tree
    #Compile the list of child variant nodes (cycles and linkers)
    if linker_check == 'Yes' and cycle_check=='Yes':
        AllMainLinkers, MainDivalentLinkers, MainTrivalentLinkers=Construct_Linker_Groups(DNA_MainTreeNodes+CyclicVariants)
        if linker_group=='All':
            linker_selection=AllMainLinkers
        if linker_group=='Divalent':
            linker_selection=MainDivalentLinkers
        if linker_group=='Trivalent':
            linker_selection=MainTrivalentLinkers
        linker_selection=remove_duplicates(linker_selection)
        ChildVariants=CyclicVariants+linker_selection
    
    #Compile the list of child variant nodes (cycles only)
    if linker_check == 'No' and cycle_check=='Yes':
        ChildVariants=CyclicVariants

    #Compile the list of child variant nodes (cycles only)
    if linker_check == 'No' and cycle_check=='No':
        ChildVariants=MainTreeNodes

    #if literature scaling option is selected, load the topology descriptor data compiled from the literature
    if lit_scale == 'Yes':
        df=pd.read_csv('Descriptor_Data.csv')
    #create a dictionary for the literature topologies based on the number of literature occurrances
        temp_dict=df.set_index('Entry').to_dict()
        for k,v in temp_dict.items():
            occurance_dict=v

    #Filter the child variant nodes according to the diversity element setting
    DNA_ChildTreeNodes=filter_nodes_for_n_DEs(Diversity_Elements, ChildVariants)
    
    #Remove the "DNA" from the child node strings
    ChildTreeNodes=Remove_DNA(DNA_ChildTreeNodes)
    
    #Map the child nodes to the main tree (via edgelist)
    ChildEdges=Map_Child_Main_Nodes(ChildTreeNodes,MainTreeNodes)

    #Assemble final nodes and edges for plotting
    FinalTreeNodes=MainTreeNodes+ChildTreeNodes
    FinalTreeEdges=MainTreeBranches+ChildEdges



    TREE= Network(height='800px',width='1000px')

    #set the color and shape for the the tree nodes
    if lit_scale == 'No':
        for node in FinalTreeNodes:
            if node.isupper() is not True:
                if '!' in node:
                    TREE.add_node(node, shape="square", color="#F88379")
                else:
                    TREE.add_node(node, shape="hexagon", color="#F88379")
            if '!' in node:
                TREE.add_node(node, shape="square", color="#fff9ae")
            else:
                TREE.add_node(node, shape="dot")
            
    if lit_scale == 'Yes':
        for node in FinalTreeNodes:
            if node in occurance_dict.keys():
                if node.isupper() is not True:
                    if '!' in node:
                        TREE.add_node(node, shape="square", color="#c6a8e0", size=occurance_dict.get(node)*10)
                    else:
                        TREE.add_node(node, shape="hexagon", color="#c6a8e0", size=occurance_dict.get(node)*10)
                if '!' in node:
                    TREE.add_node(node, shape="square", color="#c6a8e0", size=occurance_dict.get(node)*10)
                else:
                    TREE.add_node(node, shape="dot", color="#c6a8e0", size=occurance_dict.get(node)*10)
            elif node.isupper() is not True:
                if '!' in node:
                    TREE.add_node(node, shape="square", color="#F88379", size=10)
                else:
                    TREE.add_node(node, shape="hexagon", color="#F88379", size=10)
            elif '!' in node:
                    TREE.add_node(node, shape="square", color="#fff9ae", size=10)
            else:
                TREE.add_node(node, shape="dot",size=10)

    #Add the list of edges to the topology tree
    TREE.add_edges(FinalTreeEdges)

    #Uncomment this to generate the HTML file with buttons to adjust viewing options
    #TREE.show_buttons()


    if view=='Hierarchical':
        TREE.set_options("""
        const options = {
        "layout": {
            "hierarchical": {
            "enabled": true,
            "levelSeparation": 300,
            "nodeSpacing": 85,
            "treeSpacing": 150
            }
        },
        "physics": {
            "hierarchicalRepulsion": {
            "centralGravity": 0,
            "avoidOverlap": null
            },
            "minVelocity": 0.75,
            "solver": "hierarchicalRepulsion"
        }
        }
        """)

    #Control toggle for showing the topology tree outside of a streamlit deployment (debug purposes only)
    #TREE.show( 'Tree' + ".html")


    # Save and read graph as HTML file (on Streamlit Sharing)
    try:
        path = '/tmp'
        TREE.save_graph(f'{path}/Tree.html')
        HtmlFile1 = open(f'{path}/Tree.html','r',encoding='utf-8')
    # Save and read graph as HTML file (locally)
    except:
        path = './'
        TREE.save_graph(f'{path}/Tree.html')
        HtmlFile1 = open(f'{path}/Tree.html','r',encoding='utf-8')

MainTreeBranches=[]

#Control toggle for showing the DEL tree outside of a streamlit deployment (debug purposes only)
#DEL.show("del.html")

# Save and read graph as HTML file (on Streamlit Sharing)
try:
    path = '/tmp'
    DEL.save_graph(f'{path}/del.html')
    HtmlFile2 = open(f'{path}/del.html','r',encoding='utf-8')
  
# Save and read graph as HTML file (locally)
except:
    path = './'
    DEL.save_graph(f'{path}/del.html')
    HtmlFile2 = open(f'{path}/del.html','r',encoding='utf-8')
#with col2:
st.header("Topology Explorer")
with st.container():
    components.html(HtmlFile2.read(), width= 800, height=800)
#with col1:
if view_permutation_options == 'Yes':
    st.header("Topology Tree")
    with st.container():
        components.html(HtmlFile1.read(), width= 800, height=800)

GuidePath=Image.open('DELTA_Guide.png')
st.header('About DELTA')
st.image(GuidePath)
