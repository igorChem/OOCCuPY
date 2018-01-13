#/usr/bin/env python
# -*- coding: utf-8 -*-
# fmoinput.py

#=======================================================================

#load modules

import os 
from pdb_class import*

#=======================================================================

'''
project planning 

1. Fragment class to hold de Fragment input information 
2. Input function builder 
	a) main control ans SCF groups 
	b) cube options 
	c) print option 
	d) bases functions text (get from fortran program) 
'''
#=======================================================================
#basis functions


#-----------------------------------------------------------------------
# 3_21G
HO_3_21G =  "  3-21G       9   5 \n    "
HO_3_21G += "1 0 -0.070460   0.164842   0.000000   0.000000   0.502117 \n    "
HO_3_21G += "     0.466939   0.000000   0.000000   0.418037\n    "
HO_3_21G += "0 1 -0.070466   0.164858   0.473392   0.000000  -0.167360\n    "
HO_3_21G += "     0.466968   0.394111   0.000000  -0.139327\n    "
HO_3_21G += "0 1 -0.070466   0.164858  -0.236695  -0.409969  -0.167360\n    "
HO_3_21G += "     0.466968  -0.197056  -0.341311  -0.139327\n    "
HO_3_21G += "0 1 -0.070466   0.164858  -0.236695   0.409969  -0.167360\n    "
HO_3_21G += "     0.466968  -0.197056   0.341311  -0.139327\n    "
HO_3_21G += "0 1  1.006716   0.074826   0.000000   0.000000   0.000000\n    "
HO_3_21G += "    -0.143831   0.000000   0.000000  -0.000004\n"

#-----------------------------------------------------------------------

#sto-g 
sto_g = "  STO-3G       5   5 \n    "
sto_g+="1 0  -0.117784   0.542251   0.000000   0.000000   0.850774 \n    "
sto_g+="0 1  -0.117787   0.542269   0.802107   0.000000  -0.283586 \n    "
sto_g+="0 1  -0.117787   0.542269  -0.401054  -0.694646  -0.283586 \n    "
sto_g+="0 1  -0.117787   0.542269  -0.401054   0.694646  -0.283586 \n    "
sto_g+="0 1   1.003621  -0.015003   0.000000   0.000000   0.000000 \n"
     
#-----------------------------------------------------------------------
     
#MINI
    
mini ="  MINI       5   5 \n    "     
mini+="1 0  -0.104883  0.308874   0.000000   0.000000   0.521806 \n    "
mini+="0 1  -0.104884  0.308875   0.491962   0.000000  -0.173935 \n    "
mini+="0 1  -0.104884  0.308876  -0.245981  -0.426051  -0.173934 \n    "
mini+="0 1  -0.104884  0.308876  -0.245981   0.426051  -0.173934 \n    "
mini+="0 1   0.988209  0.063992   0.000000   0.000000   0.000000 \n"

#-----------------------------------------------------------------------
     
#6-311G*
     
HO_6_311Gpol =" 6-311G*     19   5 \n    "
HO_6_311Gpol+="1 0  -0.026794  -0.086612   0.000000   0.000000   0.252561 \n    "
HO_6_311Gpol+="      0.344233   0.000000   0.000000   0.444258   0.240985 \n    "
HO_6_311Gpol+="      0.000000   0.000000   0.292724  -0.016321  -0.016321 \n    "
HO_6_311Gpol+="      0.071305   0.000000   0.000000   0.000000   \n    "
HO_6_311Gpol+="0 1  -0.026796  -0.086614   0.238115   0.000000  -0.084186 \n    "
HO_6_311Gpol+="      0.344242   0.418849   0.000000  -0.148082   0.240988 \n    "
HO_6_311Gpol+="      0.275979   0.000000  -0.097574   0.061569  -0.016321 \n    "
HO_6_311Gpol+="     -0.006585   0.000000  -0.031795   0.000000   \n    "
HO_6_311Gpol+="0 1  -0.026796  -0.086614  -0.119057  -0.206214  -0.084186 \n    "
HO_6_311Gpol+="      0.344242  -0.209425  -0.362733  -0.148082   0.240989 \n    "
HO_6_311Gpol+="     -0.137989  -0.239005  -0.097574   0.003151   0.042097 \n    "
HO_6_311Gpol+="     -0.006585   0.038944   0.015898   0.027536   \n    "
HO_6_311Gpol+="0 1  -0.026796  -0.086614  -0.119057   0.206214  -0.084186 \n    "
HO_6_311Gpol+="      0.344242  -0.209425   0.362733  -0.148082   0.240989 \n    "
HO_6_311Gpol+="     -0.137989   0.239005  -0.097574   0.003151   0.042097 \n    "
HO_6_311Gpol+="     -0.006585  -0.038944   0.015898  -0.027536   \n    "
HO_6_311Gpol+="0 1   0.571513   0.482733   0.000000   0.000000   0.000000 \n    "
HO_6_311Gpol+="     -0.048217   0.000000   0.000000   0.000000  -0.036258 \n    "
HO_6_311Gpol+="      0.000000   0.000000   0.000000  -0.002854  -0.002854 \n    "
HO_6_311Gpol+="     -0.002854   0.000000   0.000000   0.000000 \n"

#-----------------------------------------------------------------------

#3-21+G
HO_3_21mG = "3-21+G       13   5 \n    "
HO_3_21mG +="1 0   -0.068862   0.158922   0.000000   0.000000   0.488690 \n     "
HO_3_21mG +="       0.470307   0.000000   0.000000   0.405055   0.024784 \n     "
HO_3_21mG +="       0.000000   0.000000   0.034050\n     "
HO_3_21mG +="0 1   -0.068865   0.158929   0.460748   0.000000  -0.162900\n     "
HO_3_21mG +="       0.470302   0.381889   0.000000  -0.135011   0.024780\n     "
HO_3_21mG +="       0.032098   0.000000  -0.011349\n     "
HO_3_21mG +="0 1   -0.068865   0.158929  -0.230374  -0.399020  -0.162900\n     "
HO_3_21mG +="       0.470302  -0.190944  -0.330725  -0.135011   0.024780\n     "
HO_3_21mG +="      -0.016048  -0.027797  -0.011349\n     "
HO_3_21mG +="0 1   -0.068865   0.158929  -0.230374   0.399020  -0.162900\n     "
HO_3_21mG +="       0.470302  -0.190944   0.330725  -0.135011   0.024780\n     "
HO_3_21mG +="      -0.016048   0.027797  -0.011349\n     "
HO_3_21mG +="0 1    1.007087   0.075762   0.000000   0.000000   0.000001\n     "
HO_3_21mG +="      -0.151743   0.000000   0.000000  -0.000005  -0.009474\n     "
HO_3_21mG +="       0.000000   0.000000  -0.000001\n"
     
#-----------------------------------------------------------------------

#6-31G 
HO_6_31G ="  6-31G       9   5 \n    "
HO_6_31G +="1 0   -0.067724  0.300281  0.000000  0.000000  0.606750 \n    "
HO_6_31G +="       0.306535  0.000000  0.000000  0.309793 \n    "
HO_6_31G +="0 1   -0.067730  0.300310  0.572037  0.000000 -0.202234 \n    "
HO_6_31G +="       0.306552  0.292061  0.000000 -0.103255 \n    "
HO_6_31G +="0 1   -0.067730  0.300310 -0.286019 -0.495398 -0.202234 \n    "
HO_6_31G +="       0.306552 -0.146031 -0.252933 -0.103255 \n    "
HO_6_31G +="0 1   -0.067730  0.300310 -0.286019  0.495398 -0.202234 \n    "
HO_6_31G +="       0.306552 -0.146031  0.252933 -0.103255 \n    "
HO_6_31G +="0 1    1.011954 -0.016447  0.000000  0.000000  0.000000 \n    "
HO_6_31G +="      -0.059374  0.000000  0.000000 -0.000001 \n"

#-----------------------------------------------------------------------

#6-31G*
HO_6_31Gpol ="  6-31G@ 15   5 \n    "
HO_6_31Gpol +="1 0    -0.065034   0.288264   0.000000   0.000000   0.604413\n    "
HO_6_31Gpol +="        0.290129   0.000000   0.000000   0.319045  -0.017106\n    "
HO_6_31Gpol +="       -0.017106   0.057935   0.000000   0.000000   0.000000\n    "
HO_6_31Gpol +="0 1    -0.065041   0.288294   0.569833   0.000000  -0.201457\n    "
HO_6_31Gpol +="        0.290147   0.300784   0.000000  -0.106342   0.049599\n    "
HO_6_31Gpol +="       -0.017106  -0.008771   0.000000  -0.027223   0.000000\n    "
HO_6_31Gpol +="0 1    -0.065040   0.288293  -0.284917  -0.493490  -0.201456\n    "
HO_6_31Gpol +="        0.290146  -0.150393  -0.260487  -0.106341  -0.000428\n    "
HO_6_31Gpol +="        0.032923  -0.008771   0.033353   0.013612   0.023577\n    "
HO_6_31Gpol +="0 1    -0.065040   0.288293  -0.284917   0.493490  -0.201456\n    "
HO_6_31Gpol +="        0.290146  -0.150393   0.260487  -0.106341  -0.000428\n    "
HO_6_31Gpol +="        0.032923  -0.008771  -0.033353   0.013612  -0.023577\n    "
HO_6_31Gpol +="0 1     1.010938  -0.011976   0.000000   0.000000   0.000000\n    "
HO_6_31Gpol +="       -0.054085   0.000000   0.000000  -0.000001  -0.003175\n    "
HO_6_31Gpol +="       -0.003175  -0.003175   0.000000   0.000000   0.000000\n"

#-----------------------------------------------------------------------
#6-31G**
HO_6_31Gpol2 ="  6-31G** 15   5 \n    "
HO_6_31Gpol2 += "1 0    -0.068254   0.305270   0.000003   0.000000   0.619132\n    "
HO_6_31Gpol2 += "        0.287030   0.000002   0.000000   0.307201  -0.022701\n    "
HO_6_31Gpol2 += "       -0.022701   0.042170   0.000000   0.000000   0.000000\n    "    
HO_6_31Gpol2 += "0 1    -0.068257   0.305303   0.583705   0.000000  -0.206360\n    "
HO_6_31Gpol2 += "        0.287057   0.289613   0.000000  -0.102393   0.034962\n    "
HO_6_31Gpol2 += "       -0.022700  -0.015495   0.000000  -0.023534   0.000000\n    "
HO_6_31Gpol2 += "0 1    -0.068257   0.305306  -0.291851  -0.505502  -0.206358\n    "
HO_6_31Gpol2 += "        0.287061  -0.144805  -0.250811  -0.102392  -0.008284\n    "
HO_6_31Gpol2 += "        0.020546  -0.015495   0.028830   0.011767   0.020381\n    "
HO_6_31Gpol2 += "0 1    -0.068257   0.305306  -0.291851   0.505502  -0.206358\n    "
HO_6_31Gpol2 += "        0.287061  -0.144805   0.250811  -0.102392  -0.008284\n    "
HO_6_31Gpol2 += "        0.020546  -0.015495  -0.028830   0.011767  -0.020381\n    "
HO_6_31Gpol2 += "0 1     1.010732  -0.013164   0.000000   0.000000   0.000001\n    "
HO_6_31Gpol2 += "       -0.052063   0.000000   0.000000   0.000000  -0.001621\n    "
HO_6_31Gpol2 += "       -0.001621  -0.001620   0.000000   0.000000   0.000000\n"

#-----------------------------------------------------------------------

#6-31++G**
HO_6_31mmGpol2 ="  6-31++G**  19   5 \n    "
HO_6_31mmGpol2 +="1 0     -0.064922  0.305919  0.000000  0.000000  0.622220  \n    "   
HO_6_31mmGpol2 +="         0.288564  0.000000  0.000000  0.309676 -0.024300 \n    "
HO_6_31mmGpol2 +="         0.000000  0.000000  0.008101 -0.022838 -0.022838 \n    "
HO_6_31mmGpol2 +="         0.042186  0.000000  0.000000  0.000000 \n    "
HO_6_31mmGpol2 +="0 1     -0.064927  0.305945  0.586627  0.000000 -0.207401 \n    "
HO_6_31mmGpol2 +="         0.288580  0.291956  0.000000 -0.103223 -0.024312 \n    "
HO_6_31mmGpol2 +="         0.007627  0.000000 -0.002702  0.034961 -0.022836 \n    "
HO_6_31mmGpol2 +="        -0.015616  0.000000 -0.023589  0.000000 \n"
HO_6_31mmGpol2 +="0 1     -0.064927  0.305944 -0.293315 -0.508035 -0.207399 \n    "
HO_6_31mmGpol2 +="         0.288578 -0.145978 -0.252842 -0.103223 -0.024312 \n    "
HO_6_31mmGpol2 +="        -0.003814 -0.006606 -0.002702 -0.008387  0.020512 \n    "
HO_6_31mmGpol2 +="        -0.015616  0.028898  0.011795  0.020428 \n    "
HO_6_31mmGpol2 +="0 1     -0.064927  0.305944 -0.293315  0.508035 -0.207399 \n    "
HO_6_31mmGpol2 +="         0.288578 -0.145978  0.252842 -0.103223 -0.024312 \n    "
HO_6_31mmGpol2 +="        -0.003814  0.006606 -0.002702 -0.008387  0.020512 \n    "
HO_6_31mmGpol2 +="        -0.015616 -0.028898  0.011795 -0.020428 \n    "
HO_6_31mmGpol2 +="0 1      1.011030 -0.014586  0.000000  0.000000  0.000000 \n    "
HO_6_31mmGpol2 +="        -0.054140  0.000000  0.000000  0.000000  0.006441 \n    "
HO_6_31mmGpol2 +="         0.000000  0.000000  0.000000 -0.001572 -0.001572 \n    "
HO_6_31mmGpol2 +="        -0.001572  0.000000  0.000000  0.000000 \n"

DFTB = "HOP_C 4 4 \n"
DFTB+= "1 0 0.562060  0.000000  0.000000  0.827096 \n"
DFTB+= "0 1 0.562060  0.779794  0.000000 -0.275699 \n"
DFTB+= "0 1 0.562060 -0.389897  0.675322 -0.275698 \n"
DFTB+= "0 1 0.562060 -0.389897 -0.675322 -0.275698 \n"

#-----------------------------------------------------------------------
basis_functions = [mini,sto_g,HO_3_21G,HO_6_31Gpol,HO_6_311Gpol,DFTB]
#-----------------------------------------------------------------------
## Data basis text

text_basis = {"MINI":0,"STO-3G":1,"3-21G":2,"6-31G*":3,"6-311G*":4,"dftb":5}

path_to_dftb = "/home/barden/Dropbox/Mestrado_Igor/Revisão/dftb/mio-1-1/"

#=======================================================================

class fragment: 
	
	def __init__(self):
		self.name      = '' 
		self.charge    = 0
		self.scf_typ   = "rhf"
		self.mult      = 1 
		self.atomIndex = []
		self.indat     = []
		self.Cnum      = ''
		self.Oxynum    = ''
		self.Oxynum1p  = ''
		self.Oxynum1l  = ''
		
	def add_res(
			self,
			res
			):
		
		self.name += res.name 
		if len(self.name) >=9:
			self.name = self.name[:8] 
		self.name
		self.charge += res.charge
		self.atomIndex += res.atomsNum	
		self.Cnum = res.alfaC.num 
		self.Oxynum = res.oxygen.num
		self.Oxynum1p = res.carbB.num
		self.Oxynum1l = res.carb.num

class FMO_input:
	
	def __init__(
			self,
			pdbfile,
			nres=1,
			DFT=False,
			MP2=False,
			basis="STO-3G",
			Conv=1,
			elePot=False,
			pcm=True,
			nbody=2,
			sulfur=False
			):

		self.name         = pdbfile[:-4] + '.inp'
		self.pdb_name     = pdbfile
		self.input_text   = '' 
		self.nresFrag     = nres
		self.print_level  = "low"
		self.grid_size    = 0.5
		self.res          = []
		self.basis        = basis
		self.nbody        = nbody
		self.fragm        = []
		self.base         = basis_functions[text_basis[basis]]
		self.scftyp       = "rhf"
		self.dft          = DFT
		self.mp2          = MP2
		self.atoms        = []
		self.conv         = Conv
		self.modprp       = 21
		self.pcm          = pcm 
		self.dfttyp       = "b3lyp"
		self.gbasis       = 'N31'
		self.polar        = '' 
		self.ngauss       = 6
		self.ndfunc       = 0
		self.npfunc       = 0		
		self.diffuseP     = False
		self.diffuseS     = False
		self.sulfur       = sulfur
		self.disp         = "UFF"


		if not self.basis == '':		
			if self.basis == '3-21G':			
				self.gbasis = 'n21'			
				self.ngauss = 3
			elif self.basis == '6-31G':
				self.gbsis = 'n31'
			elif self.basis == '6-31G*':
				self.gbasis = 'n31'
				self.ndfunc = 1
			elif self.basis == '6-311G':
				self.gbasis = 'n311'
			elif self.basis == '6-311G(p,d)':
				self.gbasis = 'n311'
				self.ndfunc = 1
			elif self.basis == '6-311G(2p,2d)':
				self.gbasis = 'n311'
				self.npfunc = 2
				self.ndfunc = 2
			elif self.basis == '6-311+G(2p,2d)':
				self.gbasis = 'n311'
				self.npfunc = 2
				self.ndfunc = 2
				self.diffuseS = True
			elif self.basis == '6-311++G(2p,2d)':	
				self.gbasis = 'n311'
				self.npfunc = 2
				self.ndfunc = 2
				self.diffuseP = True
				self.diffuseS = True
			elif self.basis == 'dftb':
				self.gbasis = 'dftb'
		
		if elePot == True: 
			sefl.modprp = 53		
			
		PL= ''		
		if self.print_level == "low":
			PL="-5"
		
		self.mw = 200
					
		self.input_text += " $contrl runtyp=energy nprint={0} maxit=200 $end \n".format(PL)
		self.input_text += " $system mwords={0} $end \n".format(self.mw)
		self.input_text += " $scf npreo(1)=0,-1,1,9999 $end \n"

		if self.conv == 1:
			self.input_text += " $scf  dirscf=.t. soscf=.t. npunch=0 $end \n"
		elif self.conv == 2:
			self.input_text += " $scf dirscf=.t.  diis=.t. swdiis=1e-3 npunch=0 $end \n"

		self.input_text += " $grid size={0} $end \n".format(self.grid_size)
		self.input_text += " $guess guess=huckel $end\n"

		if not self.basis == "dftb":
			self.input_text += " $basis gbasis = {0} ngauss = {1} $end \n".format(self.gbasis,self.ngauss)
			self.input_text += " $basis ndfunc = {0} npfunc = {1} $end \n".format(self.ndfunc,self.npfunc)
			if not self.polar == '': 
				self.input_text +=  " $basis polar = {0} $end \n".format(self.polar)
		
			if self.diffuseP == True:
				self.input_text +=  " $basis diffp =.true.  $end \n"
				self.input_text +=  " $basis diffs =.true. $end \n"
			elif self.diffuseP == False and self.diffuseS == True:
				self.input_text +=  " $basis diffs =.true. $end \n"
		elif self.basis== "dftb":
				self.input_text += " $basis gbasis = {0}  $end \n".format(self.gbasis)
				self.input_text += " $dftb scc=.true. ndftb=3 modgam=1 modesd=2 disp={0} $end \n".format(self.disp)	
				self.input_text += " $dftbsk \n"
				self.input_text += "   C C "+ path_to_dftb +"/C-C.skf \n"
				self.input_text += "   C H "+ path_to_dftb +"/C-H.skf \n"
				self.input_text += "   C O "+ path_to_dftb +"/C-O.skf \n"
				self.input_text += "   C N "+ path_to_dftb +"/C-N.skf \n"
				self.input_text += "   H C "+ path_to_dftb +"/H-C.skf \n"
				self.input_text += "   H H "+ path_to_dftb +"/H-H.skf \n"
				self.input_text += "   H O "+ path_to_dftb +"/H-O.skf \n"
				self.input_text += "   H N "+ path_to_dftb +"/H-N.skf \n"
				self.input_text += "   O C "+ path_to_dftb +"/O-C.skf \n"
				self.input_text += "   O H "+ path_to_dftb +"/O-H.skf \n"
				self.input_text += "   O N "+ path_to_dftb +"/O-N.skf \n"
				self.input_text += "   O O "+ path_to_dftb +"/O-O.skf \n"
				self.input_text += "   N C "+ path_to_dftb +"/N-C.skf \n"
				self.input_text += "   N H "+ path_to_dftb +"/N-H.skf \n"
				self.input_text += "   N O "+ path_to_dftb +"/N-O.skf \n"
				self.input_text += "   N N "+ path_to_dftb +"/N-N.skf \n"
				self.input_text += "   S S "+ path_to_dftb +"/S-S.skf \n"
				self.input_text += "   S H "+ path_to_dftb +"/S-H.skf \n"
				self.input_text += "   S C "+ path_to_dftb +"/S-C.skf \n"
				self.input_text += "   S O "+ path_to_dftb +"/S-O.skf \n"
				self.input_text += "   S N "+ path_to_dftb +"/S-N.skf \n"
				self.input_text += "   N S "+ path_to_dftb +"/N-S.skf \n"
				self.input_text += "   O S "+ path_to_dftb +"/O-S.skf \n"				
				self.input_text += "   H S "+ path_to_dftb +"/H-S.skf \n" 
				self.input_text += "   C S "+ path_to_dftb +"/C-S.skf \n"
				self.input_text += " $end \n"
				
		if self.dft == True:
			self.input_text += " $dft nrad0=96 nthe0=12 nphi0=24 $end \n"
			self.input_text += " $contrl dfttyp={0} $end \n".format(self.dfttyp)

		if self.pcm == True:
			self.input_text += " $pcm solvnt=water ief=-10 icomp=2 icav=1 idisp=1 ifmo=2 $end \n"
			self.input_text += " $pcmcav radii=suahf $end \n"
			self.input_text += " $tescav ntsall=240 $end  \n"


		self.input_text += " $fmoprp \n      "
		self.input_text += "modorb=3 \n      "
		self.input_text += "maxit=55 \n      "
		self.input_text += "modprp={0} \n      ".format(self.modprp)
		self.input_text += "npcmit=2 \n      "		
		self.input_text += "irest=0 \n"
		self.input_text += "nprint=9 \n"

		if self.conv == 2:
			self.input_text += "      ncvscf=2 mconv(2)=785 \n      "

		self.input_text += " $end \n"
		
	
	def init_frag(self):
		
		obj = protein(self.name)
		obj.pdb_parse(self.pdb_name)
		obj.residue_def()
		obj.charge_res()		
		self.res = obj.chain
		self.atoms = obj.atoms		
			
		if self.nresFrag == 1:
			for i in range(len(self.res)):
				frg = fragment()
				frg.add_res(self.res[i])
				self.fragm.append(frg)			
		elif self.nresFrag == 2: 
			k=0; l=0		
			for i in range(len(self.res)):
				if i % 2 == 0 or i == len(self.res):
					frg = fragment()					
					self.fragm.append(frg)
			for i in range(len(self.fragm)):
				self.fragm[i].add_res(self.res[k])
				if not i == len(self.fragm):
					k+=2
			if not len(self.res) % 2 == 0:
				for i in range(len(self.fragm)-1):
					l+=1
					self.fragm[i].add_res(self.res[l])
					if not i == len(self.fragm):
						l+=1
			else:
				for i in range(len(self.fragm)):
					l+=1
					self.fragm[i].add_res(self.res[l])
					if not i == len(self.fragm):
						l+=1				
	
	def mod_charge(
			self,
			n,
			chg
			):		
		
		init_chg = self.fragm[n].charge
		self.fragm[n].charge = init_chg + chg		

		if not abs(chg) % 2 == 0:			
			self.scftyp = 'rohf'
			self.fragm[n].mult=2

		if self.base == "dftb" and self.scftyp=="rohf":
			self.scftyp="uhf"		
					
	def build_input(
			self,
			n=0
			):				
		
		self.input_text += " $fmo \n      scftyp(1)={0} \n      ".format(self.scftyp)
		self.input_text += "nfrag={0}\n      ".format(len(self.fragm))
		self.input_text += "nbody={0}\n      ".format(self.nbody)
		self.input_text += "nlayer=1\n      "
		
		if self.dft == True:
			self.input_text +="dfttyp(1)={0}\n      ".format(self.dfttyp)
		
		if self.mp2 == True:
			self.input_text +="mplevl(1)=2\n      "			
		
		if self.scftyp == 'rohf':			
			self.input_text += "scffrg({0})=rohf\n      ".format(n)
			self.input_text += "mult({0})=2\n      ".format(n)
		elif self.scftyp =="uhf":
			self.input_text += "scffrg({0})=uhf\n      ".format(n)
			self.input_text += "mult({0})=2\n      ".format(n)
				
		self.input_text += "icharg(1) ="
				
		for i in range(len(self.fragm)):
			if i == len(self.fragm)-1:
				self.input_text += "{0:>2} ".format(self.fragm[i].charge)							
			else:
				self.input_text += "{0:>2}, ".format(self.fragm[i].charge)	
				if not i == 0 and i % 10==0:
					self.input_text +=" \n                "	

		self.input_text +="\n\n     frgnam(1)="
		for i in range(len(self.fragm)):				
			if not i == len(self.fragm)-1:
				self.input_text +="{0:>10}, ".format(self.fragm[i].name)
				if not i == 0 and i % 5==0 or i==4:
					self.input_text +=" \n               "	
			else:
				self.input_text +="{0:>10} ".format(self.fragm[-1].name)

					
		self.input_text +="\n\n    indat(1)=0 \n               "
		
		if self.nresFrag == 1:
			for i in range(len(self.fragm)):
				if i == 0:					
					self.input_text +="1 2 5 -{0:^3} 0 \n              ".format(self.fragm[i].atomIndex[-1])				
				elif i == len(self.fragm)-1:
					self.input_text +="{0:^3} {1:^3} {2:^3} -{3:^3} 0 \n                ".format(self.fragm[i-1].atomIndex[2],self.fragm[i-1].atomIndex[3],self.fragm[i].atomIndex[0],self.fragm[i].atomIndex[-1])					
				else:
					self.input_text +="{0} {1:^3} {2:^3} {3:^3} {4:^3} -{5:^3} 0 \n             ".format(self.fragm[i-1].atomIndex[2],self.fragm[i-1].atomIndex[3],self.fragm[i].atomIndex[0],self.fragm[i].atomIndex[1],self.fragm[i].atomIndex[4],self.fragm[i].atomIndex[-1])
		elif self.nresFrag == 2:
			for i in range(len(self.fragm)):
				if i == 0:
					self.input_text +="1    -{0:^3}      {1:^3}   -{2:^3}    0 \n              ".format(self.fragm[i].Cnum,self.fragm[i].Oxynum1p,self.fragm[i].atomIndex[-1])
				elif i == len(self.fragm)-1:
					self.input_text +="{0:^3}   {1:^3}   {2:^3}   -{3:^3}    0 \n             ".format(self.fragm[i-1].Oxynum1l,self.fragm[i-1].Oxynum,self.fragm[i].atomIndex[0],self.fragm[i].atomIndex[-1])			
				else:
					self.input_text +="{0:^3}   {1:^3}   {2:^3}   -{3:^3}   {4:^3}   -{5:^3}     0 \n              ".format(self.fragm[i-1].Oxynum1l,self.fragm[i-1].Oxynum,self.fragm[i].atomIndex[0],self.fragm[i].Cnum,self.fragm[i].Oxynum1p,self.fragm[i].atomIndex[-1])			
		
		self.input_text +=" $end \n"
		self.input_text +=" $fmohyb \n"
		if self.basis == "dftb":
			self.input_text += self.base
		elif not self.basis == "dftb":
			self.input_text += self.base
			self.input_text += mini

		self.input_text +=" $end \n"
		self.input_text +=" $fmobnd \n"
		
		base_name = self.base.split()[0]		
		
		for i in range(len(self.fragm)-1):
			if i == 0:
				self.input_text +="    -{0} {1} {2}      \n".format(self.fragm[i].Cnum,self.fragm[i].Oxynum1l,base_name)
			else:
				self.input_text +="    -{0} {1} {2}      \n".format(self.fragm[i].Cnum,self.fragm[i].Oxynum1l,base_name)
			
		self.input_text +=" $end \n"	

		self.input_text +=" $data \n"
		self.input_text +=" comentary\n"
		self.input_text +=" c1\n"		
		self.input_text +=" H    1 \n"
		self.input_text +=" C    6 \n"
		self.input_text +=" O    8 \n"
		self.input_text +=" N    7 \n"
		if self.sulfur == True:
			self.input_text +=" S    16 \n"
		self.input_text +=" $end \n"		
		self.input_text +=" $fmoxyz \n"
		
		g = 0 
		for i in range(len(self.atoms)):
			g = i + 1
			self.input_text +="{0:>6} {1:^12} {2:^12} {3:^12} {4:^12} \n".format(g,self.atoms[i].element,self.atoms[i].xcoord,self.atoms[i].ycoord,self.atoms[i].zcoord)     
		
		self.input_text +=" $end"									
			
	def write_input(
			self,
			out_name
			):
			
		inp_file = open(out_name,'w')
		inp_file.write(self.input_text)
		inp_file.close()
	
	def make_input(self):
		
		self.init_frag()
		self.build_input()
		self.write_input(out_name=self.name)
			
#==============================================================================
# End of program #
#==============================================================================
		
		
		
			
		
		
		
		
		
		
		
