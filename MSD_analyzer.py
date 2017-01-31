#!/usr/bin/env python2.7

import os
import sys
import argparse
import numpy as np
import csv
import matplotlib

matplotlib.use('Agg')
from prettytable import PrettyTable

from Pymatgen_Related.Diffusion_error_analyzer import DiffusionErrorAnalyzer
from pymatgen.core import Structure, get_el_sp
from pymatgen.io.vasp.outputs import Vasprun


def diffusion_analyze(args):
	if os.path.exists(os.path.join(args.directory,args.folder_feature+str(args.runs_start),"vasprun.xml")):
		vasprun_name="vasprun.xml"
	elif os.path.exists(os.path.join(args.directory,args.folder_feature+str(args.runs_start),"vasprun.xml.gz")):
		vasprun_name="vasprun.xml.gz"
	else:
		print "no vasprun.xml or vasprun.xml.gz in the folder, please check"
		sys.exit()
	paths=[]
	for i in range(args.runs_start,args.runs_end+1):
		path=os.path.join(args.directory,args.folder_feature+str(i),vasprun_name)
        	paths.append(path)
	start_vasprun=Vasprun(paths[0])
	origin_struct=start_vasprun.structures[0]
	DA=DiffusionErrorAnalyzer.from_files(paths,args.element,mint=args.min_t,maxt=args.max_t)
	plt=DA.get_plot(mode="default")
	element_el_sp=get_el_sp(args.element)
	sum_result_with_msd = DA.get_results(include_msd_t=True)
	if args.msd_data:
		msd_result = {}
		msd_result['msd'] = sum_result_with_msd['msd']
		msd_result['msd_components'] = np.array(sum_result_with_msd['msd_components'])
		msd_result['dt'] = sum_result_with_msd['dt']
	sum_result = sum_result_with_msd.copy()
	del(sum_result['msd'])
	del(sum_result['msd_components'])
	del(sum_result['dt'])
	del(sum_result['msd_true'])
	

	if args.save_msd:
		plt.savefig("msd_pymatgen_unweighted_"+str(args.runs_start)+\
		'_'+str(args.runs_end)+'_'+str(args.element)+".png")
	if args.plot_msd:
		plt.show()
	header_args = ("args","value")
	args_table = PrettyTable(header_args)
	args_dict=vars(args)
	for i in args_dict:
		args_table.add_row([i,args_dict[i]])
	print args_table.get_string(sortby="args")
	header_a = ("Directory", "Formula")
	a = PrettyTable(header_a)
	a.add_row([args.directory,origin_struct.composition.formula])
	print a
	header_b = ("Parameter","Value")
	b = PrettyTable(header_b)
	b.align["Parameter"] = "l"
	b.add_row(['diffusion atom:',args.element])
	b.add_row(['charge of diffusion atom in pymatgen:',element_el_sp.full_electronic_structure[-1][2]])
	b.add_row(['unit of diffusivity','cm^2/s'])
	b.add_row(['unit of conductivity','ms/cm'])
	print b
	header_c = ("Parameter","Value")
	c = PrettyTable(header_c)
	c.align["Parameter"] = "l"
	for i in sum_result:
		c.add_row([i,sum_result[i]])
	c.add_row(['Final sampled time /fs',sum_result_with_msd['dt'][-1]])
	c.add_row(['Total msd /A^2',sum_result_with_msd['msd_true']])
	print c.get_string(sortby="Parameter")
	
	if args.msd_data:
		with open("msd-t_pymatgen_"+str(args.runs_start)+\
			'_'+str(args.runs_end)+'_'+str(args.element)+".csv","w") as f:
				w_csv = csv.writer( f, delimiter = ',')
				data = [msd_result['dt'], msd_result['msd'], msd_result['msd_components'][:,0],\
				msd_result['msd_components'][:,1],msd_result['msd_components'][:,2]]
				w_csv.writerow(['dt','msd','msd_component_0','msd_component_1','msd_component_2'])
				w_csv.writerows(zip(*data))

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="""Analyze the diffusion""")
	parser.add_argument("element", type=str, help="The interested diffusion element"
			"Notice that the poscar has no charge information, so no need to provide the charge information"
			"like O2-, just provide O")
	parser.add_argument("folder_feature", type=str, help="the folders of all MD result start with same letters. "
			" The folders are RUN_0, RUN_1, RUN_2 .... You should provide RUN_. "
			" If folders are run_000, run_001, run_002.. Please change them to run_0,run_1,...")
	parser.add_argument("runs_start", type=int, help="the start number of the folders you are interested in. ")
	parser.add_argument("runs_end", type=int, help="the end number of the folders you are interested in. "
			"The end number is included in the calculation."
			"If folders are run_0, run_1, run_2. The end number is 2")
	parser.add_argument("min_t", type=float, help="min percent of time region ")
	parser.add_argument("max_t", type=float, help="max percent of time region ")
	parser.add_argument("-dir", "--directory", dest="directory", type=str, default='.',
			help="The path of the MD result. "
			"Default is current directory")
	parser.add_argument("--output_msd_data", dest="msd_data", action='store_true', default=False,
			help="print the msd-t. "
			"The msd-t information will be writen into a file msd-t-pymatgen-element-start-num-end-num.csv"
			)
	parser.add_argument("--plot_msd", dest="plot_msd", action='store_true', default=False,
			help="whether plot the msd figure,  Default is not plot"
			)
	parser.add_argument("--not_save_msd", dest="save_msd", action='store_false', default=True, 
			help="save the msa figure as:"
			"msd_pymatgen_weighted-or-not_start-num_end-num_diffusion-element.png."
			"Deafult is save it"
			)
	parser.set_defaults(func=diffusion_analyze)
	args = parser.parse_args()
	args.func(args)
