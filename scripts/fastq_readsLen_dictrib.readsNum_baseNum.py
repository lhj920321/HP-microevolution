#!/usr/local/bin/python
# -*- coding: utf-8 -*-  

import sys
import linecache
fastqF = str(sys.argv[1])
outP = str(sys.argv[2])



def header_ReadsLen():
	readsLen_header = "sampleID" + "\t" + "LenflagL" + "\t" + "LenflagR"+ "\t" + "readsCount" + "\n"
	readsLenFO=open(fastqF_readsLentongji,'a')
	readsLenFO.write(readsLen_header)
	readsLenFO.close()
def header_ReadsBaseNum():
	ReadsBaseNum_header = "sample" + "\t" + "readsNum" + "\t" + "baseNum" + "\n"
	reads_baseNumFO=open(reads_baseNumF,'a')
	reads_baseNumFO.write(ReadsBaseNum_header)
	reads_baseNumFO.close()


def readsLen_tongji(fastqfile,Bin):
	zong_base_number = 0
	hangshu_count = len(open(fastqfile,'rU').readlines()) 
	sample_reads_count = hangshu_count/4
	readLenDic = {} 
	
	for fanwei_zuo in range(0,500000,Bin):
		readLenDic[fanwei_zuo] = 0

	for reads_count in range(0,hangshu_count/4):
		read = linecache.getline(fastqfile,4*reads_count+2)
		#print read
		read_len = len(read)-1
		zong_base_number = zong_base_number + read_len

		readLenFlag = (read_len/Bin)*Bin
		readLenDic[readLenFlag] += 1
	readLenlst = readLenDic.keys()
	readLenlst.sort()
	#print sample_reads_count
	#print zong_base_number
	#print readLenlst
	ReadsBaseNum_out = fastqfile.split("/")[-1] + "\t" + str(sample_reads_count) + "\t" + str(zong_base_number)
	print ReadsBaseNum_out
	reads_baseNumFO=open(reads_baseNumF,'a')
	reads_baseNumFO.write(ReadsBaseNum_out + "\n")
	reads_baseNumFO.close()




	Len_outL = []
	for readLen in readLenlst:
		ReadsLen_outL = []
		#if readLenDic[readLen] != 0:
		readsLen_header = "sampleID" + "\t" + "LenflagL" + "\t" + "LenflagR"+ "\t" + "readsCount" + "\n"
		ReadsLen_outL.append(fastqfile.split("/")[-1])
		ReadsLen_outL.append(str(readLen))
		ReadsLen_outL.append(str(readLen + Bin))
		ReadsLen_outL.append(str(readLenDic[readLen]))
		ReadsLen_out = "\t".join(ReadsLen_outL)
		Len_outL.append(ReadsLen_out)

	Len_out = "\n".join(Len_outL)
	print Len_out
	readsLenFO=open(fastqF_readsLentongji,'a')
	readsLenFO.write(Len_out + "\n")
	readsLenFO.close()












Bin = 1000
#fastqF="/home/liuhj/liuhj_2_44/Nanopore_HP_2/data/pacbio/7083.pacbio.fastq"

out_stat_path = "/home/liuhj/liuhj_2_44/Nanopore_HP_2/" + outP
fastqF_readsLentongji = out_stat_path + "/" + fastqF.split("/")[-1] + ".reads_len_fenbu.txt"
reads_baseNumF = out_stat_path + "/" + fastqF.split("/")[-1]  + ".readsNum_baseNum.txt"

header_ReadsLen()
header_ReadsBaseNum()
readsLen_tongji(fastqF,Bin)


