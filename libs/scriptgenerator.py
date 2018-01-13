#!/usr/bin/env python
# -*- coding: utf-8 -*-
# scriptgenerator.py

# python script for sh script generation for cenapad gamess jobs 

import os

class script:
	
	'''
	'''
	
	def __init__(self,name,nprocs=16,queue="paralela"):
		
		'''
		'''
		self.jobname = name
		self.nprocs = nprocs 
		self.fila = queue
		self.log = name + '.log'
		self.out = name + '.out'
		self.err = name + '.err'
		
		if queue == "paralela":
			pass
		
		
	def write(self):
		
		'''
		'''
		input_text = ''
		input_text += '#!/bin/sh \n'
		input_text += '#PBS -q {0} \n'.format(self.fila)
		input_text += '#PBS -l ncpus={0} \n'.format(self.nprocs)
		input_text += '#PBS -N {0} \n'.format(self.jobname)
		input_text += '#PBS -o {0} \n'.format(self.out)
		input_text += '#PBS -e {0} \n'.format(self.err)
		input_text += 'echo "----------------------------------------" \n'
		input_text += 'echo "Inicio do job:" `date` \n'
		input_text += 'echo "Nodes utilizados para execucao:" \n'
		input_text += 'cat $PBS_NODEFILE \n'
		input_text += 'echo "----------------------------------------" \n'
		input_text += 'cd $PBS_O_WORKDIR \n'
		input_text += 'dirgms=/usr/local/gamess/18_AUG_2016-R1 \n'
		input_text += 'NCPUS={0} \n'.format(self.nprocs)
		input_text += '$dirgms/rungms {0} 00 $NCPUS >& {1} \n'.format(self.jobname,self.log)
		input_text += 'echo "----------------------------------------" \n'
		input_text += 'echo "Final do job:" `date` \n'
		input_text += 'echo "----------------------------------------" \n'
		
		
		scriptfile = open(self.jobname + '.sh','w')
		scriptfile.write(input_text)
		scriptfile.close()
		
		

#os.system('sudo cp /home/igor/Dropbox/Mestrado_Igor/cenapadjobs/scriptgenerator.py /usr/lib/python2.7')

		
		
