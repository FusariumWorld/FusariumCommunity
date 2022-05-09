#!/usr/bin/env python3
#!-*-coding:utf-8-*-


import os
from Bio import SearchIO
from Bio import SeqIO
from multiprocessing import Pool

upseq="upstream featured sequence"#e.g. "TCGTCATCGGCCA"
downseq="downstream featured sequence"#e.g. "GAGCTCGG" 

#adapter
Fa0_29=['AAGC','AAGT','AAGG','AAGA','AAAC','AAAT','AAAG','AAAA','AACC','AACT','AACG','AACA','AATC','AATT','AATG','AATA','AGGC','AGGT','AGGG','AGGA','AGAC','AGAT','AGAG','AGAA','AGCC','AGCT','AGCG','AGCA','AGTC','AGTT']
Fa30_59=['AGTG','AGTA','ATGC','ATGT','ATGG','ATGA','ATAC','ATAT','ATAG','ATAA','ATCC','ATCT','ATCG','ATCA','ATTC','ATTT','ATTG','ATTA','ACGT','ACGG','ACGA','ACGC','ACAC','ACAT','ACAG','ACAA','ACCC','ACCT','ACCG','ACCA']
Fa60_89=['ACTC','ACTT','ACTG','ACTA','CAGC','CAGT','CAGG','CAGA','CAAC','CAAT','CAAG','CAAA','CACC','CACT','CACG','CACA','CATC','CATT','CATG','CATA','CGGT','CGGA','CGAC','CGAT','CGAG','CGAA','CGCT','CGCA','CGTC','CGTT']
Fa90_119=['CGTG','CGTA','CTGC','CTGT','CTGG','CTGA','CTAC','CTAT','CTAG','CTAA','CTCC','CTCT','CTCG','CTCA','CTTC','CTTT','CTTG','CTTA','CCGT','CCGA','CCAC','CCAT','CCAG','CCAA','CCCT','CCCA','CCTC','CCTT','CCTG','CCTA']
Fa120_149=['TAGC','TAGT','TAGG','TAGA','TAAC','TAAT','TAAG','TAAA','TACC','TACT','TACG','TACA','TATC','TATT','TATG','TATA','TGGC','TGGT','TGGG','TGGA','TGAC','TGAT','TGAG','TGAA','TGCC','TGCT','TGCG','TGCA','TGTC','TGTT']
Fa150_179=['TGTG','TGTA','TTGC','TTGT','TTGG','TTGA','TTAC','TTAT','TTAG','TTAA','TTCC','TTCT','TTCG','TTCA','TTTC','TTTT','TTTG','TTTA','TCGC','TCGT','TCGG','TCGA','TCAC','TCAT','TCAG','TCAA','TCCC','TCCT','TCCG','TCCA']
Fa180_209=['TCTC','TCTT','TCTG','TCTA','GAGC','GAGT','GAGG','GAGA','GAAC','GAAT','GAAG','GAAA','GACC','GACT','GACG','GACA','GATC','GATT','GATG','GATA','GGGT','GGGA','GGAC','GGAT','GGAG','GGAA','GGCT','GGCA','GGTC','GGTT']
Fa210_239=['GGTG','GGTA','GTGC','GTGT','GTGG','GTGA','GTAC','GTAT','GTAG','GTAA','GTCC','GTCT','GTCG','GTCA','GTTC','GTTT','GTTG','GTTA','GCGT','GCGA','GCAC','GCAT','GCAG','GCAA','GCCT','GCCA','GCTC','GCTT','GCTG','GCTA']
filetree={
	'Shandong':{
			'Lijin':['1',Fa0_29],
			'Jiyang':['1',Fa30_59],
			'Xiajin':['1',Fa60_89],
			'Zhangqiu':['1',Fa90_119],
			'Changqing':['1',Fa120_149],
			'Taian':['1',Fa150_179],
			'Juye':['1',Fa180_209]
			},
	'Henan':{
			'Pingyu':['1',Fa210_239],
			'Nanyang':['2',Fa210_239]
			},
	'Shaanxi':{
			'Pucheng':['2',Fa0_29],
			'Sanyuan':['2',Fa30_59],
			'Linwei':['2',Fa60_89],
			'Huayin':['2',Fa90_119],
			'Chencang':['2',Fa120_149],
			'Jingyang':['2',Fa150_179],
			'Xingping':['2',Fa180_209]
			},
	'Hubei':{
			'Xiangyang':['3',Fa0_29],
			'Jingzhou':['3',Fa30_59]
			},
	'Gansu':{
			'Huixian':['3',Fa60_89]
			},
	'Anhui':{
			'Fengtai':['3',Fa120_149]
			},
	'Jiangsu':{
			'Wujin':['3',Fa150_179],
			'Taicang':['3',Fa180_209],
			'Jiangyan':['3',Fa210_239],
			'Haian':['4',Fa0_29]
			},
	'Xinjiang':{
			'Yili':['4',Fa30_59],
			'Tacheng':['4',Fa60_89]
			}
	}



def findfastq(fastqfile,writetopathwithfilename,upseq,downseq):
	print('Run child task %s.....' %(os.getpid()))
	F=open(writetopathwithfilename,'w')
	for read in SeqIO.parse(fastqfile,'fastq'):
		seq=str(read.seq)
		#recomseq=str(read.reverse_complement().seq)
		if upseq in seq:
			if downseq in seq.split(upseq)[-1]:
				F.write(">"+read.id+"\n"+seq+"\n")
				continue
		read=read.reverse_complement()
		seq=str(read.seq)
		if upseq in seq:
			if downseq in seq.split(upseq)[-1]:
				F.write(">"+read.description+"\n"+seq+"\n")
				continue
	F.close()


	
def fastqclassification_By_province():
	os.system('mkdir tmp/reads')
	for Province, City in filetree.items():
		os.mkdir("./tmp/reads/"+Province)
		for city,label in City.items():
			filepath='./tmp/reads/'+Province+'/'+city
			os.mkdir(filepath)
			#print(label[1])
			p=Pool(16) #start 16 child process number for parallel caculation
			#print('Parent process %s++++++++++++++++.'%os.getpid())
			for x in label[1]:
				writefile=filepath+'/'+x+'.fasta'
				upfeature=x+upseq
				fastqfile=label[0]+'.fastq'
				#print(x)
				p.apply_async(findfastq,args=(fastqfile,writefile,upfeature,downseq,))
				#findfastq(fastqfile,writefile,upfeature,downseq)
			p.close()
			p.join()