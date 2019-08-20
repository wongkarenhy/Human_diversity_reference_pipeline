#!/usr/bin/python

import time
import re
import os
import pandas as pd
import collections as cl
import argparse


def find_ovarlap(ctable,btable):


	cchr,cstarts,cends,csize,cname=map(list,zip(*ctable))

	if len(btable)==0:

		return [[x,[]] for x in cname]
	
	bchr,bstarts,bends,bsize,bname=map(list,zip(*btable))
	
	bstarts=[int(x)-5000 for x in bstarts]
	bends=[int(x)+5000 for x in bends]
	
	c_coordinates=cstarts+cends
	b_coordinates=bstarts+bends

	
	lc=len(cstarts)
	lb=len(bstarts)
	l2=2*lc
	
	all_coordinates=c_coordinates+b_coordinates
	
	all_coordinates_sort_index=sorted(range(len(all_coordinates)), key=lambda x:all_coordinates[x])
	
	cspan=set()
	bspan=set()
	overlap={}
	for x in all_coordinates_sort_index:
		
		if x<l2:
			real_index=x%(lc)
			
			if x<lc:
				
				cspan.add(real_index)
				overlap[real_index]=set([x for x in bspan])

			else:
				cspan.remove(real_index)
			
		else:
			
			real_index=(x-l2)%(lb)
			
			if x<lb+l2:
				
				bspan.add(real_index)
				for c in cspan:
					overlap[c].add(real_index)
			
			else:
				bspan.remove(real_index)


	
	data_organize=[[cname[cindex], [bname[bindex] for bindex in bindexes if abs(float(bsize[bindex])-float(csize[cindex]))<min(700,0.2*float(csize[cindex]))]] for cindex ,bindexes in overlap.items()]


	#print [x for x in data_organize if len(x[1])>0]
	
	return data_organize
	
				
	
			
				
			
			
		
	


def run(args):
	
	
	inputfile=args.inputcsv
	bionano=args.bio
	outputfile=args.out



	print "loading sv file\n"	
	csv_table=pd.read_csv(inputfile, sep='\t',keep_default_na=False)
	
	csvdata=csv_table[['ref_chr','ref_start','ref_end','insert_size','INS_id']].values.tolist()
	

	
	print "loading bionano file\n"

	bheader=['SmapEntryID', 'QryContigID', 'RefcontigID1', 'RefcontigID2', 'QryStartPos', 'QryEndPos', 'RefStartPos', 'RefEndPos', 'Confidence', 'Type', 'XmapID1', 'XmapID2', 'LinkID', 'QryStartIdx', 'QryEndIdx', 'RefStartIdx', 'RefEndIdx', 'Zygosity', 'Genotype', 'GenotypeGroup', 'RawConfidence', 'RawConfidenceLeft', 'RawConfidenceRight', 'RawConfidenceCenter', 'SVsize','Enzyme','Sortname']


	bheader=['RefcontigID1','RefStartPos','RefEndPos','Stype','SVsize','Enzyme','Sortname']
	
	bio_table=pd.read_csv(bionano, sep='\t',header=None, names=bheader, comment='#',index_col=False)


	biodata=bio_table[['RefcontigID1','RefStartPos','RefEndPos','SVsize','Sortname']].values.tolist()
	
	allchrs=list(set(csv_table['ref_chr']))
	
	csv_eachchr={x:[] for x in allchrs}
	bio_eachchr={str(x):[] for x in xrange(1,26)}


	print "sorting sv data\n"	
	for data0 in csvdata:
		
		csv_eachchr[data0[0]].append(data0)

	print "sorting bionano data\n"
	for data0 in biodata:
		
		bio_eachchr[str(data0[0])].append(data0)
		

	print "running validation\n"
	alloverlaps=[]
	for chr0,csvdata0 in csv_eachchr.items():
		

		print "running %s\n"%chr0

		if chr0[3:].isdigit():chr0_rename=chr0[3:] 
		elif chr0[3:]=='23': chr0_rename="X" 
		elif chr0[3:]=='24': chr0_rename="Y"
		
		biodata0=bio_eachchr[chr0_rename]
		
		alloverlaps.extend(find_ovarlap(csvdata0,biodata0))

	print "organizing data\n"	
	allnames={x:i for i,x in enumerate(list(csv_table['INS_id']))}
	
	alloverlaps=[[x[0],';'.join(x[1])+';',len(set([a.split('_')[0] for a in x[1]]))] for x in sorted(alloverlaps, key=lambda x:allnames[x[0]])]
	

	alloverlaps=pd.DataFrame.from_records(alloverlaps, columns=['INS_id',"bionano_overlap","concordant_sample_num"])


	output=pd.merge(csv_table, alloverlaps, on='INS_id')

	print "outputing\n"



	output.loc[output["concordant_sample_num"]>0].to_csv(outputfile+'_validated.txt', mode='w',sep='\t', header=None, index=False)


	output.loc[output["concordant_sample_num"]==0].to_csv(outputfile+'_notvalidated.txt', mode='w',sep='\t', header=None, index=False)



def main():

	parser=argparse.ArgumentParser(description="Find the most representitive one")
	#parser.add_argument("-d","--delta",help="delta file" ,dest="delta", type=str, required=True)
	parser.add_argument("-i","--input",help="input csv file" ,dest="inputcsv",type=str,required=True)
	parser.add_argument("-b","--bionano",help="bionano reference file" ,dest="bio",type=str,required=True)
	parser.add_argument("-o","--out",help="output file" ,dest="out",type=str,required=True)
	parser.set_defaults(func=run)
	args=parser.parse_args()
	args.func(args)

if __name__=="__main__":
	main()

