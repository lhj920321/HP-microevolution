#!/usr/bin/python
# -*- coding: utf-8 -*-

# Complementing DNA 
import sys
import linecache
from optparse import OptionParser


################################################################################################################################
# Parameters
################################################################################################################################
usage = "usage:python  %prog [option]"
parser = OptionParser(usage=usage)
parser.add_option("-v","--common_snv_F",
                    dest = "common_snv_F",
                    default = "",
                    metavar = "file",
                    help = "vcf file [required]")
parser.add_option("-o","--out_P",
                    dest = "out_P",
                    default = "",
                    metavar = "path",
                    help = "output path  [required]")


(options,args) = parser.parse_args()
common_snv_F   = os.path.abspath(options.common_snv_F)
out_P          = os.path.abspath(options.out_P)


common_snv_F_ID = common_snv_F.split("/")[-1]
out_F = out_P + common_snv_F_ID + ".iSNVFreqDistrib.txt"

if (not os.path.exists(out_Path)):
    os.mkdir(out_Path)

if (os.path.exists(out_F)):
	os.remove(out_F)



average_freq_list=[]
common_snv_F_ls = open(common_snv_F,'rU').readlines()
for common_snv_F_l in common_snv_F_ls:
	average_freq = float(common_snv_F_l.split("\t")[-1])
	average_freq_list.append(average_freq)
	
out_F_O=open(out_F,'a')
for freq_index in range(0,102,2):
	fanwei_qian = freq_index
	fanwei_hou = freq_index + 2
	count = 0
	for freq in average_freq_list:
		if freq >= fanwei_qian  and freq < fanwei_hou :
			count = count + 1
	out_line =  str(fanwei_qian) + "\t" + str(count) + "\n"
	print out_line
	out_F_O.write(out_line)	
out_F_O.close()	

		














