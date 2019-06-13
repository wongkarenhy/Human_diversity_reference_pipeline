#!/usr/bin/python

import pandas as pd

class metapath:
	
	def __init__(self,metafile):
		
		meta=pd.read_csv(metafile,sep='\t',header=None,usecols=range(10))
		self.pathes={x:path for x,path in zip(list(meta[0]),list(meta[4]))}
		
	


