#!/usr/bin/env python3
#!-*-coding:utf-8-*-


from Bio import SeqIO
from Bio.Seq import Seq
from Bio import SearchIO
import pandas as pd
import os,re


upseq="upstream featured sequence"#e.g. "TCGTCATCGGCCA"
downseq="downstream featured sequence"#e.g. "GAGCTCGG" 
FusariumstandardFile="EF.fa"

path_of_wholecountryfastAfile="tmp/WholeCountry.fasta"
path_of_wholecountryfastQfile="tmp/all.fastq"
path_of_wholecountryUniquecountfile="tmp/WholeCountryCollasper.fasta"
FusariumstandardFile_makeblastdb="tmp/EFstandard"
wholecountryUniquecountfile_blastnfile="tmp/result.blastn"
Fusariumuniquereads_annotate_file="tmp/AnnotatedWholeCountry.fasta"
FusariumStandAnotateFile="tmp/StandardAnnotatedWholeCountry.fasta"

provincedic={
	'Shandong':'SD',
	'Henan':'HN',
	'Shaanxi':'SAX',
	'Hubei':'HB',
	'Gansu':'GS',
	'Anhui':'AH',
	'Jiangsu':'JS',
	'Xinjiang':'XJ'
	}
citydic={
			'Lijin':'LJ',
			'Jiyang':'JIY',
			'Xiajin':'XJ',
			'Zhangqiu':'ZQ',
			'Changqing':'CQ',
			'Taian':'TA',
			'Juye':'JUY',
			'Pingyu':'PY',
			'Nanyang':'NY',
			'Pucheng':'PC',
			'Sanyuan':'SY',
			'Linwei':'LW',
			'Huayin':'HY',
			'Chencang':'CC',
			'Jingyang':'JY',
			'Xingping':'XP',
			'Xiangyang':'XY',
			'Jingzhou':'JZ',
			'Huixian':'HX',
			'Fengtai':'FT',
			'Wujin':'WJ',
			'Taicang':'TC',
			'Jiangyan':'JY',
			'Haian':'HA',
			'Yili':'YL',
			'Tacheng':'TC'
	}

#merge all fastq file and write them into a new file 
def merge_fastq_file():
	filelist=os.listdir()
	fq=open(path_of_wholecountryfastQfile,'w')
	for file in filelist:
		if 'fastq'in file:
			for read in SeqIO.parse(file,'fastq'):
				SeqIO.write(read,fq,'fastq')
	fq.close()
	print("Merged fastq file is completed!")



#count all unique fasta reads of whole country 
def whole_country_fasta_statistics():
	FastaFile=open(path_of_wholecountryfastAfile,'w')
	i=0
	for read in SeqIO.parse(path_of_wholecountryfastQfile,"fastq"):
		seq=str(read.seq)
		if upseq in seq :
			if downseq in seq.split(upseq)[-1]:
				i=i+1
				FastaFile.write(">"+str(i)+'_'+read.description+"\n"+seq.split(upseq)[1].split(downseq)[0]+"\n")
				continue
		read=read.reverse_complement()
		seq=str(read.seq)
		if upseq in seq :
			if downseq in seq.split(upseq)[-1]:
				i=i+1
				FastaFile.write(">"+str(i)+'_'+read.description+"\n"+seq.split(upseq)[1].split(downseq)[0]+"\n")	
				continue
	print("all reads count：%d\n"%(i))
	FastaFile.close()
	#all unique reads are counted
	FastaCollasper=open(path_of_wholecountryUniquecountfile,"w")
	Dic={}#{'unique read':count(int)}
	for read in SeqIO.parse(path_of_wholecountryfastAfile,"fasta"):
		seq=str(read.seq)
		if seq in Dic.keys():
			Dic[seq]=Dic[seq]+1
			continue
		Dic[seq]=1
	dic=sorted(Dic.items(),key=lambda item:item[1],reverse=True)
	i=0
	for set in dic:
		i=i+1
		FastaCollasper.write(">"+str(i)+'_'+str(set[1])+'\n')
		FastaCollasper.write(str(set[0])+'\n')
		
	print("The number of unique reads：%d\n"%(i))
	FastaCollasper.close()

#conduct blastn for the collasper file
def uniquereads_blastn():
	os.system("makeblastdb -in "+FusariumstandardFile+" -dbtype nucl -out "+FusariumstandardFile_makeblastdb)
	b=os.system("blastn -db "+FusariumstandardFile_makeblastdb+" -query "+path_of_wholecountryUniquecountfile+" -evalue 1e-50 -out "+wholecountryUniquecountfile_blastnfile+" -outfmt 6 -num_threads 4")

#calculating percentage of each otu
def FusariumStatistics(dic): #dic {'40131_1': 'F.graminearum_MN381097.1-PH1_99.533'}
	FusDic={}
	for key, value in dic.items():
		if 'like' in value: #dic {'40124_1': 'OTU19143~like~F.asiaticum_MT113024.1_99.299',}
			value=value.split('~')[0]
			FusDic[value]=int(key.split('_')[1])
			continue
		Fusariumspp=value.split('_')[0]
		if Fusariumspp in FusDic.keys():
			FusDic[Fusariumspp]=FusDic[Fusariumspp]+int(key.split('_')[1])
			continue
		FusDic[Fusariumspp]=int(key.split('_')[1])
	Fuslist=sorted(FusDic.items(),key=lambda item:item[1],reverse=True)
	return Fuslist

#annotating Fusarium unique reads
def Fusariumuniquereads_annotation():
	uniquereads_blastn()
	dic={}
	i=0
	for qresult in SearchIO.parse(wholecountryUniquecountfile_blastnfile,'blast-tab',comments=False):
		#print('Search %s has %i hits'%(qresult.id,len(qresult)))
		if qresult.hsps[0].ident_pct > 99.50:
			dic[qresult.id]=qresult.hsps[0].hit_id+'_'+str(qresult.hsps[0].ident_pct)
			#print(qresult.id+'_'+qresult.hsps[0].hit_id+'_'+str(qresult.hsps[0].ident_pct))
			continue
		i=i+1
		dic[qresult.id]='OTU'+str(i)+'~like~'+qresult.hsps[0].hit_id+'_'+str(qresult.hsps[0].ident_pct)
		#print(qresult.id+'_OTU'+str(i)+'~like~'+qresult.hsps[0].hit_id+'_'+str(qresult.hsps[0].ident_pct))
	#print(dic)	{'40131_1': 'F.graminearum_MN381097.1-PH1_99.533'}
	Fuslist=FusariumStatistics(dic)
	totalreads=sum([i[1] for i in Fuslist])
	print(totalreads)
	FusSta=open('percentage of each otu.txt','w')
	for x in Fuslist:
		FusSta.write(x[0]+"____"+str(x[1])+"~"+str(round(float(x[1]/totalreads),4))+'\n')
	FusSta.close()
	AnnotatedFile=open(Fusariumuniquereads_annotate_file,'w')
	for record in SeqIO.parse(path_of_wholecountryUniquecountfile,'fasta'):
		if record.id in dic.keys():
			AnnotatedFile.write(">"+record.id+"_"+dic[record.id]+"\n")
			AnnotatedFile.write(str(record.seq)+"\n")
	AnnotatedFile.close()

#list all dirs,files from the path of 'dir'
def filename(dir):
	for root, dirs, files in os.walk(dir):
		return dirs

#caculate how many reads each file contains
def readNum(fastafile):
	i=0
	for read in SeqIO.parse(fastafile,'fasta'):
		i=i+1
	return i

#merge all fasta files
def mergefastafile(singlefile,allfile):
	AllF=open(allfile,'a+')
	for read in SeqIO.parse(singlefile,'fasta'):
		AllF.write(">"+read.id+"\n"+str(read.seq)+"\n")
	AllF.close()
	
def readmerge(InitialFas,CollasperFas):
	#calculate unique reads
	FastaCollasper=open(CollasperFas,"w")
	Dic={}
	for read in SeqIO.parse(InitialFas,"fasta"):
		seq=str(read.seq)
		#remove adapter
		seq=seq.split(upseq)[-1].split(downseq)[0]
		
		if seq in Dic.keys():
			Dic[seq]=Dic[seq]+1
			continue
		Dic[seq]=1
	dic=sorted(Dic.items(),key=lambda item:item[1],reverse=True)
	i=0
	#ignore the reads with the frequency low than 5
	for set in dic:
		if int(set[1]) < 5:
			continue
		if int(set[1]) > 50:
			i=i+1
		FastaCollasper.write(">"+str(i)+'_'+str(set[1])+'\n')
		FastaCollasper.write(str(set[0])+'\n')	
	FastaCollasper.close()
	return i

def ParseFusSDFile(FusariumStandSeqFile):
	dic={}
	for read in SeqIO.parse(FusariumStandSeqFile,'fasta'):
		seq=str(read.seq)
		dic[seq]=read.id.split('_')[2]
	return dic
	
def FusariumAnnotation(FusariumStandSeqDic,CollasperFas,FusariumAnnotationFile):
	F=open(FusariumAnnotationFile,'w')
	for read in SeqIO.parse(CollasperFas,'fasta'):
		seq=str(read.seq)
		if seq in FusariumStandSeqDic.keys():
			F.write(">"+read.id+"_"+FusariumStandSeqDic[seq]+"\n"+seq+"\n")
	F.close()
		
	
	
def FusariumStatistics(FusariumStandSeqDic,FusariumAnnotationFile,FusariumStatisticsFile):
	dic={}
	F=open(FusariumStatisticsFile,'w')
	#add keys into dic using keys of FusariumStandSeqDic
	for key, value in FusariumStandSeqDic.items():
		if 'like' in value:
			value=value.split('~')[0]
			dic[value]=0
			continue
		value=value.split('_')[0]
		dic[value]=0
	for read in SeqIO.parse(FusariumAnnotationFile,'fasta'):
		classification=read.id.split('_')[-1]
		if 'like' in classification:
			classification=classification.split('~')[0]
		if classification in dic.keys():
			dic[classification]=dic[classification]+int(read.id.split('_')[1])
			continue
	for key ,value in dic.items():
		F.write(str(key)+'\t'+str(value)+'\n')
	F.close()
	return dic

def generate_standard_annotate_file(Fusariumuniquereads_annotate_file):
	F=open(FusariumStandAnotateFile,'w')
	for read in SeqIO.parse(Fusariumuniquereads_annotate_file,"fasta"):
		if int(read.id.split('_')[1]) > 19:
			SeqIO.write(read,F,'fasta')
	F.close()


def classification_by_city():
	ProvinceList=filename('./tmp/reads')
	#print(ProvinceList)
	FusariumStandSeqDic={}
	generate_standard_annotate_file(Fusariumuniquereads_annotate_file)
	FusariumStandSeqDic=ParseFusSDFile(FusariumStandAnotateFile)
	Allpd=pd.DataFrame([])
	for province in ProvinceList:
		for root,dirs,files in os.walk('./tmp/reads/'+province):
			citylist=dirs
			break
		for city in citylist:
			shortname=provincedic[province]+citydic[city]
			#print(shortname)
			citystatisdic={}
			NewF='./tmp/reads/'+province+'/'+city+'/'+city+'.fasta'
			CollasperFas='./tmp/reads/'+province+'/'+city+'/'+city+'_collasper.fasta'
			AnnotatedFas='./tmp/reads/'+province+'/'+city+'/'+city+'_Annotation.fasta'
			FusariumStatisticsFile='./tmp/reads/'+province+'/'+city+'/'+city+'_Statistics.txt'
			if os.path.exists(NewF):
				os.remove(NewF)
			if os.path.exists(CollasperFas):
				os.remove(CollasperFas)
			if os.path.exists(AnnotatedFas):
				os.remove(AnnotatedFas)
			if os.path.exists(FusariumStatisticsFile):
				os.remove(FusariumStatisticsFile)
			for root, dirs, files in os.walk('./tmp/reads/'+province+'/'+city):
				print(city)
				filelist=files
				break
			for file in filelist:
				filepath='./tmp/reads/'+province+'/'+city+'/'+file
				if readNum(filepath)<40:
					#print('filterout')
					continue
				mergefastafile(filepath,NewF)
			readmerge(NewF,CollasperFas)
			FusariumAnnotation(FusariumStandSeqDic,CollasperFas,AnnotatedFas)#f(IN,IN,OUT)
			FusariumStatistics(FusariumStandSeqDic,AnnotatedFas,FusariumStatisticsFile)#f(IN,IN,OUT)
			FSF=open(FusariumStatisticsFile,'r')
			citystatisdic={}
			for line in FSF.readlines():
				line = line.strip('\n')
				citystatisdic[line.split('\t')[0]]=line.split('\t')[1]
			FSF.close()

			df=pd.DataFrame([citystatisdic],index=[shortname],columns=citystatisdic.keys())
			Allpd=pd.concat([Allpd,df])
	Allpd=Allpd.T
	Allpd.to_csv("allstatistics_by_city.csv")		

def classification_by_head():
	ProvinceList=filename('./tmp/reads')
	generate_standard_annotate_file(Fusariumuniquereads_annotate_file)
	FusariumStandSeqDic=ParseFusSDFile(FusariumStandAnotateFile)
	Allpd=pd.DataFrame([])
	for province in ProvinceList:
		for root,dirs,files in os.walk('./tmp/reads/'+province):
			citylist=dirs
			#print(citylist)
			break
		if citylist==[]:
			continue
		for city in citylist:
			stn=provincedic[province]+citydic[city]
			filepath='./tmp/reads/'+province+'/'+city
			NewF='./tmp/reads/'+province+'/'+city+'/'+city+'.singlehead.fasta'
			CollasperFas='./tmp/reads/'+province+'/'+city+'/'+city+'singlehead.collasper.fasta'
			AnnotatedFas='./tmp/reads/'+province+'/'+city+'/'+city+'singlehead.Annotation.fasta'
			FusariumStatisticsFile='./tmp/reads/'+province+'/'+city+'/'+city+'singlehead.Statistics.txt'
			if os.path.exists(NewF):
				os.remove(NewF)
			if os.path.exists(CollasperFas):
				os.remove(CollasperFas)
			if os.path.exists(AnnotatedFas):
				os.remove(AnnotatedFas)
			if os.path.exists(FusariumStatisticsFile):
				os.remove(FusariumStatisticsFile)
			for root, dirs, files in os.walk(filepath):
				#print(city)
				filelist=files
				#print(filelist)
				break
			for file in filelist:
				if not re.match(r'[ATCG]{4}.fasta',file):
					continue
				if readNum(filepath+'/'+file) < 40:
					print(file)
					continue
				readmerge(filepath+'/'+file,CollasperFas)
				FusariumAnnotation(FusariumStandSeqDic,CollasperFas,AnnotatedFas)
				spstatisdic=FusariumStatistics(FusariumStandSeqDic,AnnotatedFas,FusariumStatisticsFile)
				#print(spstatisdic)
				shortname=stn+'_'+file.split('.')[0]
				df=pd.DataFrame([spstatisdic],index=[shortname],columns=spstatisdic.keys())
				Allpd=pd.concat([Allpd,df])
				#print(Allpd)
	Allpd=Allpd.T
	Allpd.to_csv("allstatistics_singlehead_by_city.csv")
