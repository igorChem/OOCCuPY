#!/usr/bin/env python
# -*- coding: utf-8 -*-
# gms_in.py

#=======================================================================
# This script was produced by Igor Barden Grillo 
# Computational chemistry Reactions Descriptor project 
# Start date: 17.09.2017
#=======================================================================

# import modules 
from xyz_class import *
from scriptgenerator import *

#=======================================================================
#class and functions to create input fo gamess program 
#=======================================================================


path_to_dftb = "/home/barden/Dropbox/mestrado/Revisão/dftb/mio-1-1"

class gms_inp:
	
#-----------------------------------------------------------------------
	# constructor method
#-----------------------------------------------------------------------
	def __init__(self,xyzname):
		self.name         = xyzname
		self.inp          = None
		self.charge       = 0 
		self.xyzA         = None
		self.xyzB         = None
		self.multiplicity = 1 
		self.ab_initio    = 'HF'
		self.semiEmpi     = 'AM1'
		self.runtyp       = 'energy'
		self.scf_typ      = 'rhf' 
		self.gbasis       = 'N31'
		self.polar        = '' 
		self.ngauss       = 6
		self.ndfunc       = 0
		self.npfunc       = 0
		self.guess        = 'huckel' 
		self.pcm          = False
		self.diffuseP     = False
		self.diffuseS     = False
		self.elsden       = False
		self.elspot       = False
		self.sys_group    = ''
		self.contrl_group = ''
		self.scf_group    = '' 
		self.opt_group    = ''
		self.basis_group  = ''
		self.pcm_group    = ''
		self.elsden_group = ''
		self.grid_cube    = ''
		self.data_group   = ''
		self.elspot_group = ''
		self.guess_group  = '' 
		self.other_groups = ''
		self.dftb_group   = ''
		self.input_text   = '' 
		self.diis         = False
		self.soscf        = True
		self.damp         = False 
		self.shift        = False 
		self.rstrct       = False
		self.opptol       = 0.0001
		self.alg_opt      = 'QA'
		self.dfttyp       = 'b3lyp'
		self.punch		  = 0
		self.basis 		  = ''
		self.dftb         = 0
		self.disp         = "UFF"
		self.morb         = 0

		if not self.basis == '':		
			if self.basis == '3-21G':			
				self.gbasis = 'n21'			
				self.ngauss = 3
			elif self.basis == '6-31G':
				self.gbsis = 'n31'
			elif self.basis == '6-31G*':
				self.gbasis = 'n31'
				self.npfunc = 1
			elif self.basis == '6-311G':
				self.gbasis = 'n311'
			elif self.basis == '6-311G(p,d)':
				self.gbasis = 'n311'
				self.npfunc = 1
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
				
			
#------------------------------------------------------------------------------------------
		# init_groups method
#------------------------------------------------------------------------------------------
	
	def init_groups(self):
		
		self.xyzA = xyz_parser(self.name)		
		self.xyzB = xyz_parser(self.name)

		if not self.multiplicity == 1:
			self.scf_typ = 'uhf'			
		
		#------------------------------------------------------------#
			#contrl group
		#------------------------------------------------------------#			
		self.contrl_group =  ' $contrl runtyp= {0} maxit=150 $end \n'.format(self.runtyp)
		self.contrl_group += ' $contrl  icharg = {0} mult = {1}  $end \n'.format(self.charge,self.multiplicity)
		self.contrl_group += ' $contrl scftyp={0} $end \n'.format(self.scf_typ)
		
		if self.ab_initio == 'DFT':
			self.contrl_group +=' $contrl  dfttyp = b3lyp $end \n'
		elif self.ab_initio == 'MP2':
			self.contrl_group +=' $contrl  mplevl=2 $end \n'
		#------------------------------------------------------------#	
			# sytem grouo
		#------------------------------------------------------------#	
		self.sys_group +=  ' $system mwords = 200 $end \n'
		#------------------------------------------------------------#
			#scf group 
		#------------------------------------------------------------#
		self.scf_group = ' $scf dirscf =.true. npunch={0} $end \n'.format(self.punch)
		self.scf_group += ' $scf npreo(1)=0,-1,1,9999 $end \n'
		
		if self.diis   == False:
			self.scf_group += ' $scf soscf=.true. $end \n'
		elif self.diis == True:
			self.scf_group += ' $scf diis=.true. ethrsh=2.0 swdiis=0.005 $end \n'
			
		if self.shift == True:
			self.scf_group += ' $scf shift=.true. $end \n'
			
		if self.rstrct == True:
			self.scf_group += ' $scf rstrct=.true. $end \n'
		
		if self.damp == True:
			self.scf_group += ' $scf damp=.true. $end \n'
		#------------------------------------------------------------#
			#basis group
		#------------------------------------------------------------#	
		if not self.ab_initio == 'SemiEmpi':
			self.basis_group +=' $basis gbasis = {0} ngauss = {1} $end \n'.format(self.gbasis,self.ngauss)
			self.basis_group +=  ' $basis ndfunc = {0} npfunc = {1} $end \n'.format(self.ndfunc,self.npfunc)
			if not self.polar == '': 
				self.basis_group +=  ' $basis polar = {0} $end \n'.format(self.polar)
	
			if self.diffuseP == True:
				self.basis_group +=  ' $basis diffp =.true.  $end \n'
				self.basis_group +=  ' $basis diffs =.true. $end \n'
			elif self.diffuseP == False and self.diffuseP == True:
				self.basis_group +=  ' $basis diffs =.true. $end \n'
		
		elif self.dftb > 0:
			self.basis_group += ' $basis gbasis =dftb $end \n'
			self.dftb_group +=  ' $dftb scc=.true. ndftb={0} modgam=1 modesd=2 disp={1} itypmx =-1 $end \n'.format(self.dftb,self.disp)
			self.dftb_group +=  ' $dftbsk \n'
			self.dftb_group += "   C C "+ path_to_dftb +"/C-C.skf \n"
			self.dftb_group += "   C H "+ path_to_dftb +"/C-H.skf \n"
			self.dftb_group += "   C O "+ path_to_dftb +"/C-O.skf \n"
			self.dftb_group += "   C N "+ path_to_dftb +"/C-N.skf \n"
			self.dftb_group += "   H C "+ path_to_dftb +"/H-C.skf \n"
			self.dftb_group += "   H H "+ path_to_dftb +"/H-H.skf \n"
			self.dftb_group += "   H O "+ path_to_dftb +"/H-O.skf \n"
			self.dftb_group += "   H N "+ path_to_dftb +"/H-N.skf \n"
			self.dftb_group += "   O C "+ path_to_dftb +"/O-C.skf \n"
			self.dftb_group += "   O H "+ path_to_dftb +"/O-H.skf \n"
			self.dftb_group += "   O N "+ path_to_dftb +"/O-N.skf \n"
			self.dftb_group += "   O O "+ path_to_dftb +"/O-O.skf \n"
			self.dftb_group += "   N C "+ path_to_dftb +"/N-C.skf \n"
			self.dftb_group += "   N H "+ path_to_dftb +"/N-H.skf \n"
			self.dftb_group += "   N O "+ path_to_dftb +"/N-O.skf \n"
			self.dftb_group += "   N N "+ path_to_dftb +"/N-N.skf \n"
			self.dftb_group += "   S S "+ path_to_dftb +"/S-S.skf \n"
			self.dftb_group += "   S H "+ path_to_dftb +"/S-H.skf \n"
			self.dftb_group += "   S C "+ path_to_dftb +"/S-C.skf \n"
			self.dftb_group += "   S O "+ path_to_dftb +"/S-O.skf \n"
			self.dftb_group += "   S N "+ path_to_dftb +"/S-N.skf \n"
			self.dftb_group += "   N S "+ path_to_dftb +"/N-S.skf \n"
			self.dftb_group += "   O S "+ path_to_dftb +"/O-S.skf \n"				
			self.dftb_group += "   H S "+ path_to_dftb +"/H-S.skf \n" 
			self.dftb_group += "   C S "+ path_to_dftb +"/C-S.skf \n"
		else:
			self.basis_group += '$basis gbasis = {0} $end \n'.format(self.semiEmpi)
		
		#-----------------------------------------------------------#
			#opt group
		#-----------------------------------------------------------#	
		
		if self.runtyp == 'optimize':
			self.opt_group += ' $statpt opttol={0} nstep=150 $end \n'.format(self.opttol)
			self.opt_group += ' $statpt method={0} $end \n'.format(self.alg_opt)
		
		#-----------------------------------------------------------#
			#electronic density group
		#-----------------------------------------------------------#
		
		if self.elsden == True:	
			self.elsden_group += ' $eldens  ieden=1 morb={0} where=grid output=punch $end \n'.format(self.morb)
			self.grid_cube += ' $grid modgrd=1 $end \n'
			self.grid_cube += ' $grid size = 0.4 $end \n'
			
			self.xyzA.parse_xyz()
			self.xyzA.get_origin()
			self.grid_cube += self.xyzA.elsden_text
			
		#-----------------------------------------------------------#
			#electrostatic potential group
		#-----------------------------------------------------------#
		
		if self.elspot == True:	
			self.elspot_group += ' $elpot ieden=1 morb=0 where=grid output=punch $end \n'
			self.grid_cube += ' $grid modgrd=1 $end \n'
			self.grid_cube += ' $grid size = 0.25 $end \n'
			
			self.xyzA.parse_xyz()
			self.xyzA.get_origin()
			self.grid_cube += self.xyzA.elsden_text
		
		#-----------------------------------------------------------#
			#implict solvent group (pcm)
		#-----------------------------------------------------------#
		
		if self.pcm == True:
			self.pcm_group += ' $pcm solvnt=h2o $end \n'
			
		#-----------------------------------------------------------#
			#guess group
		#-----------------------------------------------------------#
		
		self.guess += ' $guess guess={0} $end \n'.format(self.guess)
		
			
		#-----------------------------------------------------------#
			#data group
		#-----------------------------------------------------------#
		
		self.data_group =  ' $data \n' 	
		self.data_group += 'Molecule specification \n'
		self.data_group += 'c1 \n'
	
		
		self.xyzB.parse_xyz()
		self.xyzB.get_atomnumber()
		self.data_group += self.xyzB.write_text()		
		self.data_group += ' $end'
			
			
#---------------------------------------------------------------------------------------------
	#join text method
#---------------------------------------------------------------------------------------------
				
	def join_text(self,inp_name='',Script=True):
		
		self.input_text =  self.contrl_group
		self.input_text += self.sys_group
		self.input_text += self.basis_group 
		self.input_text += self.scf_group 
		self.input_text += self.pcm_group 
		self.input_text += self.opt_group
		self.input_text += self.elsden_group
		self.input_text += self.grid_cube
		self.input_text += self.elspot_group
		self.input_text += self.guess_group
		self.input_text += self.dftb_group
		self.input_text += self.other_groups
		self.input_text += self.data_group
		
		
		if  inp_name == '':
			self.inp = self.name
			self.inp = self.inp[:-4]+'.inp'
		else:
			self.inp = inp_name
		
		file_inp = open(self.inp,'w')
		file_inp.write(self.input_text)
		file_inp.close()
		
		if Script == True:
			obj2 = script(name=self.inp[:-4])
			obj2.write()
			
	#----------------------------------------------------------#
		#optimize geometry protocol
	#----------------------------------------------------------#
			
	def optimize_geom(self			   ,
					  chg=0			   ,
					  multi=1		   ,
					  alg='QA'		   ,
					  QMmet= 'HF'      ,
					  basis = '6-311G*',
					  inpnam =''       ,
					  scf = 'rhf'      , 
					  opttol=0.0001   ):
		
		self.basis = basis
		self.scf_typ = scf			  
		self.runtyp = 'optimize'
		self.ab_initio = QMmet
		self.charge = chg
		self.multiplicity = multi
		self.opttol = opttol
		self.alg_opt = alg
		
		
		self.init_groups()
		self.join_text(inp_name=inpnam)
					  
	def SP( self         ,
			chg =0       ,
			QMmet= 'HF'  ,
			basis='3-21G', 
			inpnam= ''   ,
			script = True,
			multi=1     ):

		self.basis = basis
		self.ab_initio = QMmet
		self.charge = chg
		self.multiplicity = multi
			
		self.init_groups()
		self.join_text(inp_name=inpnam,Script=script)	
		
	
	def protein_inp(self          ,
					chg=0         ,
					inpnam= ''    ,
					QMmet = 'HF'  ,
					conv = 0      ,
					diis = False  ,
					elspot = False,
					multi=1	      ):
						
		self.ab_initio = QMmet
		self.charge = chg
		self.multiplicity = multi

		if QMmet == "DFTB2":
			self.dftb = 2
			self.ab_initio = "SemiEmpi"
		elif QMmet == "DFTB3":
			self.dftb = 3
			self.ab_initio = "SemiEmpi"

		
		self.ngauss = 3
		self.gbasis = 'N21'
				
		self.pcm = True
		self.elsden = True
		
		if diis == True:
			self.diis = True
			
		if elspot == True:
			self.elspot= True
			self.elsden = False

		if not conv == 0:
			if conv == 1: 
				self.damp = True
			elif conv == 2: 
				self.shift = True
			elif conv == 3:
				self.rstrct = True 
			elif conv == 4: 
				self.damp = True
				self.shift = True	
			elif conv == 5: 
				self.damp = True
				self.rstrct = True
			elif conv == 6: 
				self.shift = True
				self.rstrct = True
			elif conv == 7:
				self.shift = True
				self.rstrct = True
				self.damp = True
		
		self.init_groups()
		self.join_text(inp_name=inpnam)				
		
	
#========================================================================#

