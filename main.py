#!/usr/bin/env python3
#!-*-coding:utf-8-*-

__version__ = '1.1.4'
__author__ = 'Qiang Wang (qzw0009@nwafu.edu.cn)'


from Fusariumdiversity import phaseone,phasetwo
import os


if __name__=='__main__':
    if 'tmp' not in os.listdir():
        os.system("mkdir tmp")
    phaseone.merge_fastq_file()
    phaseone.whole_country_fasta_statistics()
    phaseone.Fusariumuniquereads_annotation()
    phasetwo.fastqclassification_By_province()
    phaseone.classification_by_city()
    phaseone.classification_by_head()