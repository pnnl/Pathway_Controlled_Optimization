

import numpy as np
import re
import cobra






'''  
###################################

Read a .dat file into a dataframe

####################################
'''

def ReadRxnFile(dat_file, reactions, S_matrix,print_output= False):
    
    
    left ='LEFT'
    right = 'RIGHT'
    left_compartment = 'LEFT_COMPARTMENT'
    right_compartment = 'RIGHT_COMPARTMENT'
    enzyme_level = 'ENZYME_LEVEL'
    deltag0 = 'DGZERO'
    deltag0_sigma = 'DGZERO StdDev'
    same_compartment = 'Same Compartment?'
    full_rxn = 'Full Rxn'
    rxn = 'REACTION' 


    fdat = open(dat_file, 'r')
    for line in fdat:
        if (line.startswith('REACTION')):
            rxn_name = line[8:-1].lstrip()
            S_matrix.loc[rxn_name,enzyme_level] = 1.0
            reactions.loc[rxn_name,enzyme_level] = 1.0
            if print_output == True:print(rxn_name)

        if (re.match("^LEFT\s",line)):
            line = substitute_common_metabolite_names(line)
            line = line.upper()
            left_rxn = line[4:-1].lstrip()
            left_rxn = re.sub(r'\s+$', '', left_rxn) #Remove trailing white space
            reactions.loc[rxn_name,left] = left_rxn

        elif (re.match('^RIGHT\s',line)):
            line = substitute_common_metabolite_names(line)            
            line = line.upper()
            right_rxn = line[5:-1].lstrip()
            right_rxn = re.sub(r'\s+$', '', right_rxn) #Remove trailing white space
            reactions.loc[rxn_name,right] = right_rxn
        
        elif (line.startswith(left_compartment)):
            cpt_name = line[16:-1].lstrip()
            reactants = re.split('\s+\+\s+',left_rxn)
            if print_output == True:print(reactants)
            for idx in reactants:
                values = re.split(' ', idx);
                if len(values) == 2:
                    stoichiometry = -np.float64(values[0]);
                    if not re.search(':',molecule):
                        molecule = molecule + ':' + cpt_name
                else:
                    stoichiometry = np.float64(-1.0);
                    molecule = values[0]; 
                    if not re.search(':',molecule):
                        molecule = molecule + ':' + cpt_name
                S_matrix.loc[rxn_name,molecule] = stoichiometry;


        elif (line.startswith(right_compartment)):
            cpt_name = line[17:-1].lstrip()
            products = re.split('\s+\+\s+',right_rxn)
            if print_output == True:print(products)
            for idx in products:
                values = re.split(' ', idx);
                if len(values) == 2:
                    stoichiometry = np.float64(values[0]);
                    molecule = values[1];
                    if not re.search(':',molecule):
                        molecule = molecule + ':' + cpt_name
                else:
                    stoichiometry = np.float64(1.0);
                    molecule = values[0];
                    if not re.search(':',molecule):
                        molecule = molecule + ':' + cpt_name
                S_matrix.loc[rxn_name,molecule] = stoichiometry;

        elif (re.match("^ENZYME_LEVEL\s", line)):
            level = line[12:-1].lstrip()
            reactions.loc[rxn_name,enzyme_level] = float(level)
            S_matrix.loc[rxn_name,enzyme_level] = float(level)
                
        elif re.match('^COMMENT',line):
            continue
        elif re.match(r'//',line):
            continue
        elif re.match('^#',line):
            continue
            
    
    reactions[full_rxn] = reactions[left]+"="+reactions[right]   
    fdat.close()
    return(reactions,S_matrix)

























'''  
###################################

Define Routine to Read a .txt file from a Pathway Tools smart frame

####################################
'''


import re
def ReadRxnFrameFile(rxn_file, reactions, S_matrix,print_output = False):
    
    left ='LEFT'
    right = 'RIGHT'
    left_compartment = 'LEFT_COMPARTMENT'
    right_compartment = 'RIGHT_COMPARTMENT'
    enzyme_level = 'ENZYME_LEVEL'
    deltag0 = 'DGZERO'
    deltag0_sigma = 'DGZERO StdDev'
    same_compartment = 'Same Compartment?'
    full_rxn = 'Full Rxn'
    rxn = 'REACTION'

    fdat = open(rxn_file, 'r')
    for line in fdat:
        if print_output == True: print(line)
        if (line.startswith('#')):
            continue
        line = substitute_common_metabolite_names(line)
        if print_output == True: print(line)
        rxn_frame, rxn, pthwy = line.split('\t')
        if print_output == True: print(rxn_frame, rxn)
        # Make name substitutions:
        rxn = rxn.rstrip('\n')

        rxn = re.sub('coenzyme A','COA',rxn,flags=re.IGNORECASE)
        rxn = re.sub('D-sedoheptulose 7-phosphate','D-sedoheptulose-7-phosphate',rxn,flags=re.IGNORECASE)
        rxn = re.sub('D-glyceraldehyde 3-phosphate','D-glyceraldehyde-3-phosphate',rxn,flags=re.IGNORECASE)
        rxn = re.sub('D-ribose 5-phosphate','D-ribose-5-phosphate',rxn,flags=re.IGNORECASE)
        rxn = re.sub('D-xylulose 5-phosphate','D-xylulose-5-phosphate',rxn,flags=re.IGNORECASE)
        rxn = re.sub('D-erythrose 4-phosphate','D-erythrose-4-phosphate',rxn,flags=re.IGNORECASE)
        rxn = re.sub('beta-D-fructofuranose 6-phosphate','beta-D-fructose-6-phosphate',rxn,flags=re.IGNORECASE)
        rxn = re.sub('beta-D-fructose 1,6-bisphosphate','beta-D-fructose-1,6-bisphosphate',rxn,flags=re.IGNORECASE)
        rxn = re.sub('glycerone phosphate','glycerone_phosphate',rxn,flags=re.IGNORECASE)
        rxn = re.sub('DHAP','glycerone_phosphate',rxn,flags=re.IGNORECASE)
        rxn = re.sub('D-glucopyranose 6-phosphate','Beta-D-glucose-6-phosphate',rxn,flags=re.IGNORECASE)
        rxn = re.sub('D-ribulose 5-phosphate','D-ribulose-5-phosphate',rxn,flags=re.IGNORECASE)
        rxn = re.sub(re.escape('1-(5-phospho-beta-D-ribosyl)-5-[(5-phosphoribosylamino)methylideneamino]imidazole-4-carboxamide'),'5-(5-Phospho-D-ribosyl-aminoformimino)-1-(5-phosphoribosyl)-imidazole-4-carboxamide',rxn,flags=re.IGNORECASE)
        rxn = re.sub(r'\+ [0-9]+ H\+','',rxn) # was re.sub not replace
        rxn = rxn.replace('+ H+','')
        rxn = rxn.upper()
        
        reactions.loc[rxn_frame, full_rxn] = rxn
        reactions.loc[rxn_frame, enzyme_level] = 1.0
        
        reactions.loc[rxn_frame, same_compartment] = 1.0
        if print_output == True:print(rxn_frame,rxn)
        left_rxn, right_rxn = rxn.split('=')
        reactions.loc[rxn_frame,left] = left_rxn
        reactions.loc[rxn_frame,right] = right_rxn
        S_matrix.loc[rxn_frame, enzyme_level] = np.float64(1.0);

        reactants = left_rxn.split(' + ')
        for x in reactants:
            x = x.strip()
            if(x == 'H+'):
                continue
            values = re.split(' ',x);
            if print_output == True:print(values)
            if len(values) == 2:
                stoichiometry = -np.float64(values[0]);
                molecule = values[1];

            else:
                stoichiometry = np.float64(-1.0);
                molecule = values[0]; 

            molecule = molecule+':CYTOPLASM'
            S_matrix.loc[rxn_frame,molecule] = stoichiometry;
        
            
        products = right_rxn.split(' + ')
        for x in products:
            x = x.strip()
            if print_output == True:print(x)
            if(x == 'H+'):
                continue
            if print_output == True:print(x)
            values = re.split(' ',x);
            if len(values) == 2:
                stoichiometry = np.float64(values[0]);
                molecule = values[1];

            else:
                stoichiometry = np.float64(1.0);
                molecule = values[0];
 
            molecule = molecule+':CYTOPLASM'
            S_matrix.loc[rxn_frame,molecule] = stoichiometry;
        if print_output == True: print("\n Next Reaction:")
            
    fdat.close()
    return(reactions,S_matrix)


























'''  
###################################

Read in SBML Models

####################################
'''



def ReadSBMLFile(fsbml, reactions, S_matrix, print_output = False):
    model = cobra.io.read_sbml_model(fsbml)


    left ='LEFT'
    right = 'RIGHT'
    left_compartment = 'LEFT_COMPARTMENT'
    right_compartment = 'RIGHT_COMPARTMENT'
    enzyme_level = 'ENZYME_LEVEL'
    deltag0 = 'DGZERO'
    deltag0_sigma = 'DGZERO StdDev'
    same_compartment = 'Same Compartment?'
    full_rxn = 'Full Rxn'
    rxn = 'REACTION'

    
    for irxn in model.reactions:
    
        chem_eqn = str(irxn.build_reaction_string(use_metabolite_names=True))
        if print_output == True: print(irxn.name,'\n\t', chem_eqn)
        chem_eqn = re.sub('<==>','=',chem_eqn)
        chem_eqn = re.sub('<=>','=',chem_eqn)
        chem_eqn = re.sub(r'\+ [0-9]+ H\+','',chem_eqn) 
        chem_eqn = re.sub(r'\+ [0-9]\.0+ H\+','',chem_eqn) 
        chem_eqn = re.sub(r'\+ H\+','',chem_eqn)
        
        left_chem_eqn, right_chem_eqn = chem_eqn.split('=')
        left_chem_eqn = substitute_common_metabolite_names(left_chem_eqn)
        right_chem_eqn = substitute_common_metabolite_names(right_chem_eqn)
        reactions.loc[irxn.name, left] = left_chem_eqn
        reactions.loc[irxn.name, right] = right_chem_eqn
        reactions.loc[irxn.name, enzyme_level] = np.float64(1.0)
        S_matrix.loc[irxn.name, enzyme_level] = np.float64(1.0)
        reactants = irxn.reactants
        products = irxn.products
        
        left_chem_eqn = ''
        operator = ''
        for rx in reactants:
            if(rx.name == 'H+'):
                continue
            try:
                compart = rubrum_compartment_hash(rx.compartment)
            except ValueError:
                if print_output == True:print("Oops!  That was no valid compartment.")
                compart = 'Error'
                return
            rx_name = str(rx.name)
            rx_name = substitute_common_metabolite_names(rx_name)
            full_metab = rx_name+':'+compart
            full_metab = full_metab.upper()
            if print_output == True: print('\t', rx.name, full_metab, irxn.get_coefficient(rx.id))
            S_matrix.loc[irxn.name, full_metab] = np.float64(irxn.get_coefficient(rx.id))
            coeff = np.abs(irxn.get_coefficient(rx.id))
            if (coeff != 1.0):
                left_chem_eqn = left_chem_eqn + operator + str(coeff) + ' ' +full_metab
            else:
                left_chem_eqn = left_chem_eqn + operator +full_metab
            operator = r' + '
        
        right_chem_eqn = ''
        operator = ''
        for pdt in products:
            if(pdt.name == 'H+'):
                continue
            try:
                compart = rubrum_compartment_hash(pdt.compartment)
            except ValueError:
                if print_output == True: print("Oops!  That was no valid compartment.")
                compart = 'Error'
                return
            pdt_name = str(pdt.name)
            pdt_name = substitute_common_metabolite_names(pdt_name)
            full_metab = pdt_name+':'+ compart
            full_metab = full_metab.upper()
            if print_output == True: print('full_metab =', full_metab)
            
            if print_output == True: print('\t',pdt.name, full_metab, irxn.get_coefficient(pdt.id))
            S_matrix.loc[irxn.name, full_metab] = np.float64(irxn.get_coefficient(pdt.id))
            coeff = irxn.get_coefficient(pdt.id)
            if (coeff != 1.0):
                right_chem_eqn = right_chem_eqn + operator + str(coeff) + ' ' +full_metab
            else:
                right_chem_eqn = right_chem_eqn + operator +full_metab
            operator = r' + '
        
        chem_eqn = left_chem_eqn + ' = ' + right_chem_eqn
        chem_eqn = chem_eqn.upper()
        if print_output == True: print(irxn.name,'\n\t', chem_eqn)
        reactions.loc[irxn.name, full_rxn] = chem_eqn
        if print_output == True: print(irxn.genes)
        for gn in irxn.genes:
            if print_output == True: print('\t',gn)
        if print_output == True: print('\n')
        S_matrix.fillna(0,inplace=True)
    return(reactions,S_matrix)
 




def rubrum_compartment_hash(comparment):
    compartment_dict = {
        'c' : 'CYTOPLASM',
        'CCO__45__IN': 'CYTOPLASM'} # this model uses a summary reaction for electron transport directly 
                                    # to cytoplasmic NAD+/NADH instead of membrane-associated quinone/quinol
    
    return (compartment_dict[comparment])






















'''  
###################################

Substitute for common metabolite names

####################################
'''


def substitute_common_metabolite_names(line,print_output = False):
    line = re.sub('BETA-D-GLUCOSE','D-GLUCOSE',line,flags=re.IGNORECASE )
    line = re.sub(' glutamate','L-GLUTAMATE',line,flags=re.IGNORECASE )
    if ('L-GLUTAMATE' not in line.upper()):
        line = re.sub('glutamate','L-GLUTAMATE',line,flags=re.IGNORECASE )
    line = re.sub("1-\(O-CARBOXYPHENYLAMINO\)-1'-DEOXYRIBULOSE-5'-PHOSPHATE","1-(2-Carboxyphenylamino)-1-deoxy-D-ribulose 5-phosphate",line,flags=re.IGNORECASE )
    line = re.sub("AMMONIUM","NH3",line, flags=re.IGNORECASE)
    line = re.sub("AMMONIA","NH3",line, flags=re.IGNORECASE)
    line = re.sub("^phosphate",r"ORTHOPHOSPHATE",line, flags=re.IGNORECASE)
    line = re.sub("\+ phosphate","+ ORTHOPHOSPHATE",line, flags=re.IGNORECASE)
    line = re.sub('coenzyme A','COA',line,flags=re.IGNORECASE)
    line = re.sub('D-sedoheptulose 7-phosphate','D-sedoheptulose-7-phosphate',line,flags=re.IGNORECASE)
    line = re.sub('D-glyceraldehyde 3-phosphate','D-glyceraldehyde-3-phosphate',line,flags=re.IGNORECASE)
    line = re.sub('D-ribose 5-phosphate','D-ribose-5-phosphate',line,flags=re.IGNORECASE)
    line = re.sub('D-xylulose 5-phosphate','D-xylulose-5-phosphate',line,flags=re.IGNORECASE)
    line = re.sub('D-erythrose 4-phosphate','D-erythrose-4-phosphate',line,flags=re.IGNORECASE)
    line = re.sub('beta-D-fructofuranose 6-phosphate','beta-D-fructose-6-phosphate',line,flags=re.IGNORECASE)
    line = re.sub('beta-D-fructose 1,6-bisphosphate','beta-D-fructose-1,6-bisphosphate',line,flags=re.IGNORECASE)
    line = re.sub('glycerone phosphate','glycerone_phosphate',line,flags=re.IGNORECASE)
    line = re.sub('DHAP','glycerone_phosphate',line,flags=re.IGNORECASE)
    line = re.sub('D-glucopyranose 6-phosphate','Beta-D-glucose-6-phosphate',line,flags=re.IGNORECASE)
    line = re.sub('D-ribulose 5-phosphate','D-ribulose-5-phosphate',line,flags=re.IGNORECASE)
    line = re.sub("a reduced thioredoxin",r"A_REDUCED_THIOREDOXIN",line, flags=re.IGNORECASE)
    line = re.sub("an oxidized thioredoxin",r"AN_OXIDIZED_THIOREDOXIN",line, flags=re.IGNORECASE)
    line = re.sub("hydrogencarbonate",r"H2CO3",line, flags=re.IGNORECASE)
    line = re.sub("L-arginino-succinate",r"L-argininosuccinate",line, flags=re.IGNORECASE)    
    line = re.sub("an N10-formyltetrahydrofolate",r"N10-formyltetrahydrofolate",line, flags=re.IGNORECASE)    
    line = re.sub("a tetrahydrofolate",r"tetrahydrofolate",line, flags=re.IGNORECASE)
    line = re.sub("a 5,10-METHYLENETETRAHYDROFOLATE",r"5,10-METHYLENETETRAHYDROFOLATE",line, flags=re.IGNORECASE)    
    line = re.sub(re.escape('NAD(P)+'),'NADP+',line,flags=re.IGNORECASE)
    line = re.sub(re.escape('NAD(P)H'),'NADPH',line,flags=re.IGNORECASE)    
    line = re.sub(re.escape('1-(5-phospho-beta-D-ribosyl)-5-[(5-phosphoribosylamino)methylideneamino]imidazole-4-carboxamide'),'5-(5-Phospho-D-ribosyl-aminoformimino)-1-(5-phosphoribosyl)-imidazole-4-carboxamide',line,flags=re.IGNORECASE)
    line = re.sub("an oxidized ferredoxin [iron-sulfur] cluster",r"AN_OXIDIZED_FERREDOXIN",line, flags=re.IGNORECASE)
    line = re.sub("a reduced ferredoxin [iron-sulfur] cluster",r"A_REDUCED_FERREDOXIN",line, flags=re.IGNORECASE)
    line = re.sub(re.escape('an oxidized ferredoxin [iron-sulfur] cluster'),'AN_OXIDIZED_FERREDOXIN',line,flags=re.IGNORECASE)
    line = re.sub(re.escape('a reduced ferredoxin [iron-sulfur] cluster'),'A_REDUCED_FERREDOXIN',line,flags=re.IGNORECASE)

    # Modify electron transfer chain to produce summary reaction
    line = re.sub('an electron-transfer quinone','NAD+',line,flags=re.IGNORECASE)
    line = re.sub('an electron-transfer quinol','NADH',line,flags=re.IGNORECASE)

    #rxn = re.sub(r'[a-z] [0-6]','',rxn)
    line = re.sub(r'\+ [0-9]+ H\+','',line) # was re.sub not replace
    line = line.replace('+ H+','')
    if print_output == True: print('Here',line)
    return(line)