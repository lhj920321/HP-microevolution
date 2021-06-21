
#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys,os
from optparse import OptionParser
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import numpy as np 

def ReadPosi(PosiF):
	PosiList = []
	for PosiFl in open(PosiF).readlines():
		if "Sample" not in PosiFl:
			Posi = PosiFl.split("\t")[2]
			PosiList.append(Posi)
	return PosiList




def main():
	CCS_Samps = ["P18-C_bR-j_2_P18-C_bR-x","P18-C_bR-x_2_P18-C_bR-x",\
		"P19-C_bR-d_2_P19-C_bR-j","P19-C_bR-j_2_P19-C_bR-j",\
		"P19-C_bR-t_2_P19-C_bR-j","P1-E-d_2_P1-E-j",\
		"P1-E-j_2_P1-E-j","P1-E-t_2_P1-E-j","P1-E-x_2_P1-E-j",\
		"P20-C_bR-d_2_P20-C_bR-t","P20-C_bR-j_2_P20-C_bR-t",\
		"P20-C_bR-t_2_P20-C_bR-t","P20-C_bR-x_2_P20-C_bR-t",\
		"P21-C_sL-d_2_P21-C_sL-j","P21-C_sL-j_2_P21-C_sL-j",\
		"P21-C_sL-t_2_P21-C_sL-j","P21-C_sL-x_2_P21-C_sL-j",\
		"P22-C_sL-d_2_P22-C_sL-d","P22-C_sL-t_2_P22-C_sL-d",\
		"P23-C_sL-j_2_P23-C_sL-j","P23-C_sL-x_2_P23-C_sL-j",\
		"P24-C_sC-j_2_P24-C_sC-t","P24-C_sC-t_2_P24-C_sC-t",\
		"P24-C_sC-x_2_P24-C_sC-t","P25-C_sC-j_2_P25-C_sC-t",\
		"P25-C_sC-t_2_P25-C_sC-t","P2-E-t_2_P2-E-t",\
		"P4-E-d_2_P4-E-j","P4-E-j_2_P4-E-j","P4-E-x_2_P4-E-j",\
	"P5-E-d_2_P5-E-j","P6-E-d_2_P6-E-t","P6-E-t_2_P6-E-t",\
	"P7-E-d_2_P7-E-j","P7-E-j_2_P7-E-j"]

	print(len(CCS_Samps))
	 
	CCS_Samps = ["P7-E-d_2_P7-E-j","P7-E-j_2_P7-E-j"]
	No = 1
	for CCS_Samp in CCS_Samps:
		Data_File = CCS_Posi_P + "/" + CCS_Samp + ".SNV.PosiAndFreq.txt"
		print(Data_File)
		CCS_SNVPosi = ReadPosi(Data_File)
		
		Data_File = NGS_Posi_P + "/" + CCS_Samp + ".SNV.PosiAndFreq.txt"
		print(Data_File)
		NGS_SNVPosi = ReadPosi(Data_File)


		out_plotF = out_P + CCS_Samp + ".test.pdf"
		plt.figure(dpi=300,figsize=(8,3))

		g = venn2(subsets = [set(NGS_SNVPosi),set(CCS_SNVPosi)], #绘图数据集
			set_labels = ('NGS', 'CCS'), #设置组名
			set_colors=("#098154","#c72e29"),#设置圈的颜色，中间颜色不能修改
			alpha=0.4,#透明度
			normalize_to=1.0)
		ax1 = plt.subplot(221)

		No += 1

		#plt.show()
with PdfPages('multipage_pdf.pdf') as pdf:
    plt.figure(figsize=(3, 3))
    plt.plot(range(7), [3, 1, 4, 1, 5, 9, 2], 'r-o')
    plt.title('Page One')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()



		plt.savefig(out_plotF)



'''
plt.show()



#plot
plt.rcParams['font.sans-serif']=['SimHei'] # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus']=False # 用来正常显示负号



plt.figure(figsize=(8,8), dpi=80)
plt.figure(1)

ax1 = plt.subplot(221)
ax1.plot(t,s, color="r",linestyle = "--")
ax2 = plt.subplot(222)
ax2.plot(t,s,color="y",linestyle = "-")
ax3 = plt.subplot(223)
ax3.plot(t,s,color="g",linestyle = "-.")
ax4 = plt.subplot(224)
ax4.plot(t,s,color="b",linestyle = ":")

'''

if __name__ == "__main__":
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)

	parser.add_option("-C","--CCS_Posi_P",
					  dest = "CCS_Posi_P",
					  default = "",
					  metavar = "path",
					  help = "path of CCS .  [required]")

	parser.add_option("-N","--NGS_Posi_P",
					  dest = "NGS_Posi_P",
					  default = "",
					  metavar = "path",
					  help = "path of NGS .  [required]")

	parser.add_option("-o","--out_P",
					  dest = "out_P",
					  default = "",
					  metavar = "path",
					  help = "output stat file Path of snv Freq distribution.  [required]")



	(options,args) = parser.parse_args()
	CCS_Posi_P   = os.path.abspath(options.CCS_Posi_P)
	NGS_Posi_P	 = os.path.abspath(options.NGS_Posi_P)
	out_P		 = os.path.abspath(options.out_P)



	if ( not os.path.exists(out_P)):
		os.mkdir(out_P)


	main()





