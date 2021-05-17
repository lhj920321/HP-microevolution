#!/usr/bin/python
# -*- coding: utf-8 -*-

# Complementing DNA 
import sys
import linecache

person = str(sys.argv[1])

#person_list="7 11 12    21 201 213 13     124 1346 142 75     340 345 144 457"
#person = '75'   



if person == '7':
	ref_samp = '7085'
	samp_list=[7083,7084,7085,7086]

if person == '11':
	ref_samp = '11385'
	samp_list=[11383,11384,11385,11386]
if person == '12':
	ref_samp = '12949'
	samp_list=[12947,12949,12950]	
if person == '21':
	ref_samp = '2170'
	samp_list=[2168,2169,2170]	
	#samp_list=[ 2168]
if person == '201':
	ref_samp = '10202'
	samp_list=[10201,10202]	

if person == '213':
	ref_samp = '10214'
	samp_list=[10213,10214]	
if person == '13':
	ref_samp = '13478'
	samp_list=[13477,13478,13479]	

###third
if person == '124':
	ref_samp = '12488'
	samp_list=[12486,12487,12488]

if person == '1346':
	ref_samp = '13466'
	samp_list=[13465,13466,13467]

if person == '142':
	ref_samp = '14270'
	samp_list=[14268,14269,14270]

if person == '75':
	ref_samp = '7599'
	samp_list=[7596,7599]

#fourth
if person == '340':
	ref_samp = '343'
	samp_list=[340,341,343]

if person == '345':
	ref_samp = '346'
	samp_list=[345,346,347]

if person == '144':
	ref_samp = '14467'
	samp_list=[14465,14467,14468]

if person == '457':
	ref_samp = '459'
	samp_list=[457,459,460]



##check
if person == '345':
	ref_samp = '346'
	#samp_list=[345,346,343]
	samp_list=[343]
if person == '144':
	ref_samp = '14467'
	#samp_list=[14465,14467,347]
	samp_list=[347]





#sample = str(sys.argv[1])
#common_snv = str(sys.argv[2])
#out_F = str(sys.argv[3])	

for sample in samp_list :
	average_freq_list=[]
	common_snv_P="/home/liuhj/liuhj_2_44/Nanopore_HP_2/process/common_snv_indel/BGI_2_PersonReprest_genome/common_snv/"
	common_snv=common_snv_P + str(sample) +  "_2_" + str(ref_samp) +  "_common.varscan.vcf"

	common_snv_ls = open(common_snv,'rU').readlines()
	for common_snv_l in common_snv_ls:
		average_freq = float(common_snv_l.split("\t")[-1])
		average_freq_list.append(average_freq)
	
	out_P="/home/liuhj/liuhj_2_44/Nanopore_HP_2/process/freq_feng_count/BGI_2_PersonReprest_genome/"
	out_F = out_P + str(sample) +  "_2_" + str(ref_samp) +  "_common.freq_feng_count.txt"
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
	print 
	print 
	print 
		
			














