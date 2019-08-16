import os

            # name                    ,     Charge
ligands  = [['1-bromohexano.pdb',           0       ],
            ['1_2-dicloro-propano.pdb',     0       ],
            ['1_3-di-bromo-propano.pdb',    0       ],
            ['1,2-diBromo-etano.mol2',      0       ],
            ['1-clorohexano.pdb',           0       ],
            ['1_5-diIodo-pentano.pdb',      0       ],
            ['1_5-diBromo-pentano.pdb',     0       ],
            ['1_3-di-Iodo-propano.pdb',     0       ],
            ['1-Iodo-butano.pdb',           0       ],
            ['1_2-diBromo-propano.pdb',     0       ],
            ['1_5-diCloro-pentano.pdb',     0       ],
            ['1-Iodohexano.pdb',            0       ],
            ['1_2-diIodo-etano.pdb',        0       ],
            ['1_2-diCloro-etano.pdb',       0       ],
            ['1-Iodo-3-cloro-propano.pdb',  0       ],
            ['1_2-diBromo-etano.pdb',       0       ],
            ['1-Bromo-butano.pdb',          0       ],
            ['1_3-di-cloro-propano.pdb',    0       ],
            ['1-Bromo-3-cloro-propano.pdb', 0       ],
            ['1-Cloro-butano.pdb',          0       ]]
            
#['1_2-diBromo-propano.pdb', '1_3-di-cloro-propano.pdb', '1-Bromo-butano.pdb', '1_2-diBromo-etano.mpdb', '1_3-di-Iodo-propano.pdb', '1-Iodohexano.pdb', '1_2-diBromo-etano.pdb', '1_2-diIodo-etano.pdb', '1-Cloro-butano.pdb', '1-Bromo-3-cloro-propano.pdb', '1-Iodo-butano.pdb', '1_3-di-bromo-propano.pdb', '1_5-diBromo-pentano.pdb', '1_5-diCloro-pentano.pdb', '1-bromohexano.pdb', '1-clorohexano.pdb', '1-Iodo-3-cloro-propano.pdb', '1_5-diIodo-pentano.pdb', '1_2-dicloro-propano.pdb', '1_2-diCloro-etano.pdb']


ligands1  = os.listdir('pdb')
ligands   =[]
for ligand in ligands1:
    ligand2 = ligand.split('.')
    if ligand2[-1] != 'pdb':
        pass
    else:
        ligands.append(ligand)


print ligands
print len(ligands)


def run_antechamber (ligand = ligand, charge =  0, antechamber = True, prmchk = True):
    """ Function doc """
    if antechamber:
        string  = ('antechamber -i %s -fi %s -o %s -fo mol2 -c bcc -nc %i' %('pdb/'+ligand, ligand[-3:], 'mol2/'+ligand[:-3]+'mol2', charge))
        print string
        os.system(string)
    
    if prmchk:
        string  = ('parmchk -i %s -f mol2 -o %s' %('mol2/'+ligand[:-3]+'mol2', 'mol2/'+ligand[:-3]+'frcmod'))
        print string
        os.system(string)
    

for ligand in ligands:
    run_antechamber( ligand = ligand, 
                     charge =  0, 
                antechamber = True, 
                     prmchk = True)
