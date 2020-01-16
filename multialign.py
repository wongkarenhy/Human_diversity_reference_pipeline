#!/usr/bin/python
import time
import re
import os
import pandas as pd
import collections as cl
import argparse
import multiprocessing as mul
from wm_tools import metapath

manager = mul.Manager()

lock1 = mul.Lock()



def find_breaks(anchor_sizes,seq):


	coordi_count=0
	for i0,base in enumerate(seq):
		if base!='-':
			if coordi_count==anchor_sizes[0]:
				break
			coordi_count+=1

	front_break=i0
	coordi_count=0
	for i1,base in enumerate(seq[i0:]):
		if base!='-':
			if coordi_count==anchor_sizes[1]:
				break
			coordi_count+=1
	insert_break=i1+front_break
	
	return (front_break,insert_break)


def read_mask(query,ref):
	
	
	findgaps=re.finditer(r'-+', query)
	
	for gap in findgaps:
		gap0=gap.span()
		ref=ref[:gap0[0]]+gap.group()+ref[gap0[0]:]
		

	return ref
	#return ref.replace('t','N').replace('g','N').replace('c','N').replace('a','N')


def gap_score(count):


	last=''
	score=0
	for x in count:
		if x==0:
			y=x
		
		elif x in ['-','n','N']:

			if last in ['-','n','N',0]:
				y=0
			elif last>=1:
				y=-2
			else:
				y=-0.2

		elif last in ['-','n','N']:
			if x>=1:
				y=x-2
			else:
				y=x-0.2

		else:
			y=x


		score+=y
		last=x
	return score

class global_count:
	
	def __init__(self, l):
		
		total=len(l)
		
		self.count=cl.defaultdict(int)
	

		find_masked=0

	
		for x in l:
			self.count[x.upper()]+=1

	
		for x in l:
			#self.count[x]=self.count[x.upper()]*0.1

			if x.islower():

				find_masked=1

				break
		

		Ncount=2*self.count['N']
		gapcount=self.count['-']	


		self.count['N']=-total*0.5
	
		if not find_masked:

			self.count={k:3*v+Ncount+gapcount-total for k,v in self.count.items()}
		

		else:
			
		
			self.count={k:1.0*(3*v+Ncount+gapcount-total)/10 for k,v in self.count.items()}


		self.count['-']=0

		for x in ['a','t','c','g','n']:

			upper_value=self.count.get(x.upper()) 

			if upper_value !=None:
				self.count[x.lower()]=upper_value





class multi_align:
	
	def __init__(self, titles,seqs):
		self.titles=titles
		self.alignments=map(list,seqs)
		self.grouped={}
	
	def findscore(self):
		
		
		
		if len(self.titles)>1:
			
			
			count=[global_count(x) for x in zip(*self.alignments)]
				
			self.score=[sum([gcount.count[x] for x,gcount in zip(alignment,count)]) for alignment in self.alignments]

			
		elif len(self.titles)==1:
			self.score=100			
		
	def pickthebest(self):
		
		if len(self.titles)>1:
		
			score_index=sorted(range(len(self.score)), key=lambda x:self.score[x],reverse=True)
			
			self.highest=self.titles[score_index[0]]
			self.highest_seq=self.alignments[score_index[0]]
			
			self.grouped[len(self.grouped)]=[self.highest]

			self.titles,self.alignments=map(list,zip(*[zip(self.titles,self.alignments)[i] for i in score_index[1:]]))
		
		elif len(self.titles)==1:
			
			self.grouped[len(self.grouped)]=[self.titles[0]]
			self.titles,self.alignments=[],[]
			
			
		
		
	def realign(self, cutoff=0):
		
		if len(self.alignments)>0:
		
			anchor_sizes=map(find_breaks, map(lambda x:map(int,x.split('_')[-2:]), [self.highest]+self.titles),[self.highest_seq]+self.alignments)
			anchor_sizes[1:]=[(x[0]-anchor_sizes[0][0],x[1]-anchor_sizes[0][0]) for x in anchor_sizes[1:]]

			#matches=map(gap_score,zip(*[[0 if (coordi<anchor_sizes[i+1][0] or coordi>=anchor_sizes[i+1][1]) or (bases[0]=='N' or x=='N') else '-' if (bases[0]=='-' and x!='-') or (x=='-' and bases[0]!='-')  else 1 if x==bases[0] else -4 for i,x in enumerate(bases[1:])] for coordi,bases in enumerate(zip(*[self.highest_seq]+self.alignments)[anchor_sizes[0][0]:anchor_sizes[0][1]])]))  
			matches=map(gap_score,zip(*[[0 if (coordi<anchor_sizes[i+1][0] or coordi>=anchor_sizes[i+1][1]) or (bases[0]=='-' and x=='-') else 'N' if (bases[0].upper()=='N' or x.upper()=='N') else '-' if (bases[0]=='-')  else 1 if x==bases[0].upper() else 0.1 if x==bases[0].lower() else -4 if (x.isupper() or (x=='-' and bases[0].isupper()) ) else -0.4 if (x!='-' and bases[0]!='-') else '-'  for i,x in enumerate(bases[1:])] for coordi,bases in enumerate(zip(*[self.highest_seq]+self.alignments)[anchor_sizes[0][0]:anchor_sizes[0][1]])]))  

	
			#matches=[score-sum(4 for _ in re.finditer(r'\-+',''.join(seq))) for score, seq in zip(matches, self.alignments)]
						
			new_title,new_alignments=[],[]

			if len(matches)!=len(self.alignments):
				matches=[0]*len(self.alignments)

			
			for match,title,alignment in zip(matches,self.titles,self.alignments):
				if match<=0:
					new_title.append(title)
					new_alignments.append(alignment)
				else:
					
					self.grouped[len(self.grouped)-1].append(title)
			
			self.titles,self.alignments=new_title,new_alignments
	
	def group(self):
		
		while len(self.titles) >0:
			
			self.findscore()
			self.pickthebest()
			self.realign()
		
		
		return self.grouped
	
def generate_files(insertions,tempfolder,refpath):


	global anchor_size

	insert_name, components, sample,haplo,assem_coordi,starts,ends,ref_chrs,ref_starts,ref_ends,strands,pathes=zip(*insertions)

	ref_starts=map(int,ref_starts)
	ref_ends=map(int,ref_ends)
	starts = map(int, starts)
	ends = map(int, ends)
	
	cutstart=max(0, min(ref_starts)-anchor_size)
	
	cutend=max(ref_ends)+anchor_size
	
	componentfile=tempfolder+str(components[0])+'.fa'

	#print pathes,sample,haplo
	assem_pathes=['%s/%s_pseudohap%s.fasta'%(path,path.split('/')[-1].split('_')[0],haplo0) if '.' in haplo0 else path for path,sample0,haplo0 in zip(pathes,sample,haplo)]
	
	scaffolds=[x.split(':')[0] for x in assem_coordi]
	
	if os.path.isfile(componentfile):
		try:
			os.system('rm %s'%componentfile)
		except:
			pass

	if os.path.isdir(tempfolder+str(components[0]))==False:
		try:
			os.system('mkdir %s'%(tempfolder+str(components[0])))
		except:
			pass
		
	cmds=[]
	skipfiles=[]
	for name, assem_path,scaffold,start,end,ref_chr,ref_start,ref_end,std in zip(insert_name,assem_pathes,scaffolds,starts,ends,ref_chrs,ref_starts,ref_ends,strands):

		if anchor_size>0:
			nextline=''
		else:
			nextline='\n'

		if start==end:

			skipfiles.append('>%s_%d_%d'%(name,abs(ref_start-cutstart),abs(end-start)))
			continue
		

		cmd='printf \">%s_%d_%d\n\" >> %s'%(name,abs(ref_start-cutstart),abs(end-start),componentfile)
	
		cmd0="seq=$(samtools faidx %s %s:%d-%d | sed 1d | tr [a-z] [A-Z]); printf  \"$seq\"  >>%s"%(refpath, ref_chr,cutstart,ref_start-1, componentfile)

		if std=='+':
					
			cmd1='seq=$(samtools faidx %s %s:%d-%d | sed 1d | tr [a-z] [A-Z]); printf  \"$seq%s\" >> %s'%(assem_path,scaffold,min(start,end), max(start,end)-1,nextline,componentfile)
		
		else:
			cmd1='seq=$(samtools faidx -i %s %s:%d-%d | sed 1d | tr [a-z] [A-Z]); printf  \"$seq%s\" >> %s'%(assem_path,scaffold,min(start,end), max(start,end)-1,nextline,componentfile)

		cmd2="seq=$(samtools faidx %s %s:%d-%d | sed 1d | tr [a-z] [A-Z]); printf  \"$seq\n\"  >> %s"%(refpath, ref_chr,ref_end,cutend-1, componentfile)

		if anchor_size>0:
	
			cmds.extend([cmd,cmd0,cmd1,cmd2])
	
		else:
			#print cmd1+'\n'
			cmds.extend([cmd,cmd1])



	#cmds=['samtools faidx %s %s:%d-%d >> %s'%(assem_path,scaffold,min(start,end), max(start,end),componentfile) for assem_path,scaffold,start,end in zip(assem_pathes,scaffolds,starts,ends)]

	
	print 'creating temp sequences file for %s'%componentfile


	if cmds==[]:

		return 'pass',skipfiles

	return componentfile,skipfiles
	map(os.system,cmds)


	print 'masking repetitive region %s'%componentfile

	os.system('RepeatMasker -qq -noint  -pa 1 -species human -xsmall -dir %s %s >&-'%(tempfolder+str(components[0]) ,componentfile))

	maskedfile=[componentfile]+[tempfolder+str(components[0])+'/'+x for x in os.listdir(tempfolder+str(components[0])) if '.masked' in x]
	

	maskedfile=maskedfile[-1]

	print 'running multialignment %s'%maskedfile

	#muscle -in %s -out %s -quiet
	os.system('kalign -i %s  -o %s -f fasta -quiet'%(maskedfile, componentfile+'_out'))

	iter0=0
	while os.path.isfile(componentfile+'_out')==False and iter0<3:

		os.system('kalign -i %s  -o %s -f fasta -quiet'%(maskedfile, componentfile+'_out'))
		time.sleep(100)
		iter0+=1	

	if os.path.isfile(componentfile+'_out')==False:
		return "fail",[]

	print 'integrating both files %s'%(componentfile+'_out')

	if maskedfile!=componentfile or 1==1:

		with open(componentfile+'_out',mode='r') as f:
			reads=f.read().split('>')[1:]
		f.close()


		with open(maskedfile,mode='r') as f:
			reads_ori=f.read().split('>')[1:]
		f.close()

		title_seq=[[x.splitlines()[0],''.join(x.splitlines()[1:])] for x in reads]

		titles=list(zip(*title_seq)[0])

		title_seq_ori=sorted([[x.splitlines()[0],''.join(x.splitlines()[1:])] for x in reads_ori], key=lambda x: titles.index(x[0]))

		title_seq=[[x[0],read_mask(x[1],r[1])] for x,r in zip(title_seq,title_seq_ori)]


		out='\n'.join(['>'+'\n'.join(x) for x in title_seq])

		with open(componentfile+'_out',mode='w') as f:
			f.write(out)

		f.close()
	

	return componentfile,skipfiles
		
	
	
	
	
def compare(insertions):
	
	global tempfolder,outputfile,refpath
	
	component=insertions[0][1]

	#insert_name, components, sample,haplo,assem_coordi,start,end,ref_chr,ref_start,ref_end,pathes=zip(*insertions)

	
	if len(insertions)<=0:
	
		insert_name, components, sample,haplo,assem_coordi,start,end,ref_chr,ref_start,ref_end,strands,pathes=zip(*insertions)
		start = map(int, start)
		end = map(int, end)
	
		insert_name_index=sorted(range(len(insert_name)), key=lambda x:abs(end[x])-abs(start[x]),reverse=True)
		multi=[[insert_name[i] for i in insert_name_index]]
		title_filtered=[]

		
	else:
	
		componentfile,skipfiles=generate_files(insertions,tempfolder,refpath)



		if componentfile=='fail':

			print 'encounting error in component #' + str(component)

			return component

		elif componentfile=='pass':

			print 'passing component #' + str(component)

			multi=[]

		else:

			print 'grouping %s'%componentfile
			if not os.path.isfile(componentfile+'_out'):
				return

			with open(componentfile+'_out',mode='r') as f:
				reads=f.read().split('>')[1:]
			f.close()

			reads_ori=''		

			title_seq=[[x.splitlines()[0],list(''.join(x.splitlines()[1:]))] for x in reads]


			unique_sizes=[len(seq[1])-seq[1].count('N')-seq[1].count('-') for seq in title_seq]

			titles_seqs_unique_sizes=map(list,zip(*sorted([title_seq0+[unique_size] for title_seq0, unique_size in zip(title_seq, unique_sizes) if unique_size>0],key=lambda x:x[-1], reverse=True)))


			if len(titles_seqs_unique_sizes)>0:
				titles,seqs,unique_sizes=titles_seqs_unique_sizes

			else:
				titles,seqs,unique_sizes=[],[],[]

			title_filtered=[x for x in zip(*title_seq)[0] if x not in  titles]


			print 'filtered sequence: ',title_filtered


			del titles_seqs_unique_sizes, title_seq, reads_ori, reads, unique_sizes


			if len(titles)>1:

				multi=multi_align(titles, seqs).group().values()

			else:
				multi=[titles]

			print 'outputing results %s'%componentfile


#	if title_filtered==[]:
		#title_filtered=['NA']

	report=','.join([';'.join(['_'.join(x.split("_")[:-2]) for x in titles]) for titles in multi])+'\t'+';'.join(['_'.join(x.split("_")[:-2]) for x in title_filtered+skipfiles])

	lock1.acquire()
	with open(outputfile, mode='a') as f:
		f.write(str(component)+'\t'+report+'\n')
	f.close()
	lock1.release()

	return ''
		

def run(args):
	
	inputfile=args.inputcsv
	metafile=args.metafile
	testmode=args.testmode
	
	m0=int(args.threads)
	
	global tempfolder,outputfile,refpath,anchor_size
	
	tempfolder=args.tempfolder
	outputfile=args.outputfile
	refpath=args.refpath
	anchor_size=int(args.anchor)
	
	try:
		os.mkdir(tempfolder)
	except:
		pass	

	allpath=metapath(metafile).pathes

	print 'reading database'
	
	allinsertions=pd.read_csv(inputfile,sep='\t',low_memory=False,dtype=str)[['INS_id','component','sample' ,'haplo' ,'assm_coords','adjusted_assm_start','adjusted_assm_end','ref_chr','ref_start','ref_end','strand']].values.tolist()


	allcomponent=cl.defaultdict(list)

	print 'preparing'
	
	for insertion in allinsertions:
		
		if insertion[3]=='unphased':
			pass
		allcomponent[insertion[1]].append(insertion+[allpath[insertion[2]]])
	
	del allinsertions
	allcomponent=[x for x in allcomponent.values() if len(x)>1]


	print 'total components:', len(allcomponent)
	
	if testmode==1:
		for component in allcomponent:
			compare(component)

	
	print 'running multi-alignments'

	p=mul.Pool(processes=m0)
	
	failed=p.map(compare,allcomponent)

	p.close()
	p.join()

	
	return 0	

	

def main():

	parser=argparse.ArgumentParser(description="Find the most representitive one")
	#parser.add_argument("-d","--delta",help="delta file" ,dest="delta", type=str, required=True)
	parser.add_argument("-i","--input",help="input csv file" ,dest="inputcsv",type=str,required=True)
	parser.add_argument("-r","--ref",help="input reference file" ,dest="refpath",type=str,required=True)
	parser.add_argument("-o","--output",help="output sv file" ,dest="outputfile",type=str,default='./multi_results.csv')
	parser.add_argument("-t","--tempfolder",help="use temp folder" ,dest="tempfolder",type=str,default='./temp/')
	parser.add_argument("-m","--metafile",help="path metafile" ,dest="metafile",type=str,default='./TMP_sample_metadata.txt')
	parser.add_argument("-n","--threads",help="number of threads" ,dest="threads",type=str,default='8')
	parser.add_argument("-a","--anchor",help="size of anchor" ,dest="anchor",type=str,default='50')
	parser.add_argument("-T","--test",help="testmode" ,dest="testmode",type=int,default=0)

	parser.set_defaults(func=run)
	args=parser.parse_args()
	args.func(args)

if __name__=="__main__":
	main()
