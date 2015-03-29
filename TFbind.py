import sys
import os
import numpy as np
import string

def gene_interval(gene,interval):
	""" This function creates a bed file temp_interval.bed
	Which is used to count the number of tags in each interval 
	
	Inputs are 
		gene: The name of the gene
		inteval: the binning interval
	The output is the bed file

	"""
	cmd="grep -w "+gene+" refGene.txt > tf_temp_1.txt" # get the gene coordinates
	os.system(cmd)
	file1=open("tf_temp_1.txt","r")
	line1=file1.readline()
	s1=line1.split()
	file1.close()

	if len(s1)==0:
		print "No gene", gene
		return 0

	gene_start=string.atoi(s1[4])
	gene_stop=string.atoi(s1[5])
	
	
	pro_region=np.zeros(3)#promoter region

	
	if s1[3]=='+':#you need strand information to get transcription start site (TSS)
		pro_region[0]=gene_start-interval
		pro_region[1]=gene_start+interval
		pro_region[2]=1
		
	elif s1[3]=='-':
		pro_region[0]=gene_stop-interval
		pro_region[1]=gene_stop+interval
		pro_region[2]=0
		

	return pro_region,s1[2]#return the region and the chromosome


if len(sys.argv) != 6:  
	sys.exit("Usage: run TFbind.py [gene_file] [interval] [num_clust] [bedfile] [p-value]")

'''
This program determines clusters with significant binding of a particular transcription factor in a region around the transcription start site of genes in the cluster. The inputs to the program are

1) gene_file:  A list of all genes in the genome
2) Interval: The region around the TSS to consider in base pairs
3) num_clust: The number of clusters. There must be a corresponding file in Clusters/cluster_*.txt (where * runs from 1 to num_clust) which lists genes in the cluster. 
4) bedfile: A bedfile which contains the TF binding sites
5) p-value: The level of significance say 0.05

'''


file1=open(sys.argv[1],"r")
num_clust=int(sys.argv[3])
interval=int(sys.argv[2])
sig_pvalue=float(sys.argv[5])
genes=[]
#form the bed file for all genes
print 'Generating BED file for all genes...'
outname="LocalBeds/AllGenes_Interval.bed"
file2=open(outname,"w")
clust_size=np.zeros(num_clust)
i=1
while file1:
	line1=file1.readline()
	s1=line1.split()
	if len(s1)==0:
		break
	res=gene_interval(s1[0],interval)
	p_int=res[0]
	

	file2.write(res[1])
	file2.write('\t')
	file2.write(str(int(max(p_int[0],0))))
	file2.write('\t')
	file2.write(str(int(max(p_int[1],0))))
	file2.write('\t')
	file2.write(s1[0])
	file2.write('\n')
	i+=1
	
	
	genes.append(s1[0])
file2.close()
file1.close()
print 'Done generating BED file.'


#form the bed file for individual clusters
print 'Generating BED file for individual clusters...'
for i in range(0,num_clust):

 	file1=open('Clusters/cluster_'+str(i+1)+'.txt',"r")
	outname="LocalBeds/Cluster_"+str(i+1)+"_Interval.bed"
	file2=open(outname,"w")
	
	while file1:
		line1=file1.readline()
		s1=line1.split()
		if len(s1)==0:
			break
		res=gene_interval(s1[0],interval)	
		if res!=0:# zero returned if gene name is not in refGene.txt (the refSeq annotation)
			p_int=res[0]
			file2.write(res[1])
			file2.write('\t')
			file2.write(str(int(max(p_int[0],0))))
			file2.write('\t')
			file2.write(str(int(max(p_int[1],0))))
			file2.write('\t')
			file2.write(s1[0])
			file2.write('\n')

		
	file2.close()

	file1.close()

print 'Done generating cluster BED files.'
	
num_samp=1000000


prop=np.zeros(num_clust)#binding proportion
pval_prop=np.zeros(num_clust)#associated p value

#get the number of genes with binding sites within the desired region of the TSS

cmd='intersectBed -a LocalBeds/AllGenes_Interval.bed -c -b '+sys.argv[4]+' > LocalBeds/AllGenes_Interval2.bed'
os.system(cmd)

file1=open('LocalBeds/AllGenes_Interval2.bed',"r")
M=0
while file1:
	line1=file1.readline()
	s1=line1.split()
	if len(s1)==0:
		break

	if int(s1[4])>0:
		M+=1



file1.close()


#now get the binding proportion for each cluster and determine statistical significance

print 'Determining cluster binding proportions and significance...'
for i in range(0,num_clust):
	
 	file1=open('Clusters/cluster_'+str(i+1)+'.txt',"r")
	
	genes_clust=[]
	while file1:
		line1=file1.readline()
		s1=line1.split()
		if len(s1)==0:
			break
		genes_clust.append(s1[0])

	cmd="intersectBed -a LocalBeds/Cluster_"+str(i+1)+"_Interval.bed -c -b "+sys.argv[4]+" > LocalBeds/Cluster_"+str(i+1)+"_Interval2.bed"
	os.system(cmd)

	file3=open("LocalBeds/Cluster_"+str(i+1)+"_Interval2.bed","r")
	MM=0
	while file1:
		line3=file3.readline()
		s3=line3.split()
		if len(s3)==0:
			break

		if int(s3[4])>0:
			MM+=1



	file3.close()
	s=np.random.hypergeometric(M, len(genes)-M, len(genes_clust),num_samp)
	pvalue=sum(s>=MM)/float(num_samp)
	pval_prop[i]=pvalue
	prop[i]=MM
	clust_size[i]=len(genes_clust)
	print i+1,MM,pvalue

print 'Done!'




np.savetxt('pval_prop_'+str(interval)+'.txt',pval_prop)
np.savetxt('prop_'+str(interval)+'.txt',prop)


