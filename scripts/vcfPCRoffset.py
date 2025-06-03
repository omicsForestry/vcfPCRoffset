#!/usr/bin/env python3

# version 0.1.0

import pysam
import scipy.stats as stats
import argparse
from subprocess import PIPE, Popen
import sys

parser = argparse.ArgumentParser(description='Go through vcf and bam files together, reporting PCR template offsets.')
parser.add_argument('-i','--input', required = True, help = 'Input vcf file')
parser.add_argument('-o','--output', required = True, help = 'Output vcf file')
parser.add_argument('-b','--bam', required = True, help = 'Bam file to check for PCR site discrepancies')
parser.add_argument('-r','--refTag', default='RD', help = 'Reference depth tag in vcf if not the default RD')
parser.add_argument('-v','--varTag', default='AD', help = 'Variant depth tag in vcf if not the default AD')
parser.add_argument('-x','--refVarTag', help = 'Reference/Variant depth tag in vcf if stored in single REF,VAR format')
parser.add_argument('-y','--varRefTag', help = 'Variant/Reference depth tag in vcf if stored in single VAR,REF format')
parser.add_argument('-f','--filter', type=float, help = 'Threshold to filter variant records')
parser.add_argument('-t','--tabix', action='store_true', help = 'Compress input file and make tabix index')

args = parser.parse_args()


def updateVcf(oldVcfFile, newVcfFile, bamFile, refTag, varTag, refVarTag=None, varRefTag=None, filterThresh=None):
	"""
	Iterates through a vcf, checking each mutation in a bam file for PCR bias and writing an output vcf
	
	Args:
		oldVcfFile (str): path to input vcf file
		newVcfFile (str): path to output vcf file
		bamFile (str): path to bam file
		refTag (str): the vcf variant tag denoting number of reference reads - default RD if not specified
		varTag (str): the vcf variant tag denoting number of variant reads - default AD if not specified
		filterThresh: (float): threshold above which to change the filter tag in the output vcf - default no filtering
	"""
	oldVcf=pysam.VariantFile(oldVcfFile, 'r')
	oldVcf.header.formats.add('ROR', '1', 'Float', 'Proportion of offset reference reads')
	oldVcf.header.formats.add('ROV', '1', 'Float', 'Proportion of offset variant reads')
	if filterThresh:
		oldVcf.header.filters.add('PCR_bias', None, None, 'REF/ALT bias detected in PCR templates')
	newVcf=pysam.VariantFile(newVcfFile, 'w', header=oldVcf.header)
	inbam=pysam.Samfile(bamFile, 'rb')
	for variant in oldVcf:
		for sample in variant.samples:
			ref_offset=999
			alt_offest=999
			if refVarTag:
				expAlt=variant.samples[sample][refVarTag][1]/(variant.samples[sample][refVarTag][0] + variant.samples[sample][refVarTag][1])
			elif varRefTag:
				expAlt=variant.samples[sample][varRefTag][0]/(variant.samples[sample][varRefTag][1] + variant.samples[sample][varRefTag][0])
			else:
				expAlt=variant.samples[sample][varTag]/(variant.samples[sample][refTag] + variant.samples[sample][varTag])
			pcrTest=binomProp(inbam, variant.chrom, variant.pos, variant.ref, variant.alts[0], expAlt, 0.05, benferroni=True)
			if pcrTest[2]>0:
				ref_offset=round(pcrTest[3]/pcrTest[2],4)
			if pcrTest[0]>0:
				alt_offset=round(pcrTest[1]/pcrTest[0],4)
			variant.samples[sample]['ROR']=ref_offset
			variant.samples[sample]['ROV']=alt_offset
			if filterThresh:
				if ref_offset > filterThresh or alt_offset > filterThresh:
					variant.filter.add('PCR_bias')
		newVcf.write(variant)
	oldVcf.close()
	newVcf.close()
	inbam.close()

def getSNVfromRead(read, pos, ref, alt):
	"""
	Reports if a read is the reference or mutant allele for an SNV
	
	Args:
		inbam (obj): an open pysam.Samfile object
		pos (int): the position of the mutation
		ref (str): the reference allele
		alt (str): the mutant allele
		
	Returns:
		list of [blank/mut/ref, read start, read end]
	"""
	res='blank'
	seq=read.query_alignment_sequence
	mutPos=(pos-1)-read.pos
	end=read.pos+len(seq)
	if mutPos < len(seq) and read.pos < (pos-2):
		if seq[mutPos]==alt:
			res='mut'
		elif seq[mutPos]==ref:
			res='ref'
		else:
			res='other'
	return([res, read.pos, end])

def getDelfromRead(read, pos, ref, alt):
	"""
	Reports if a read is the reference or mutant allele for a deletion
	
	Args:
		inbam (obj): an open pysam.Samfile object
		pos (int): the position of the mutation
		ref (str): the reference allele
		alt (str): the mutant allele
		
	Returns:
		list of [blank/mut/ref, read start, read end]
	"""

	res='blank'
	mutLen=len(ref)-len(alt)
	seq=read.query_alignment_sequence
	mutPos=(pos)-read.pos
	end=read.pos+len(seq)
	inSize=0
	if mutPos < (len(seq)-mutLen) and read.pos < pos and read.cigartuples:
		upstream=0
		for thisTuple in read.cigartuples:
			if thisTuple[0]==0:
				upstream=upstream+thisTuple[1]
			if thisTuple[0]==2:
				inSize=thisTuple[1]
				break
		if upstream==mutPos and inSize==mutLen:
			res='mut'
		else:
			res='ref'
	return([res, read.pos, end+inSize])

def getInsfromRead(read, pos, ref, alt):
	"""
	Reports if a read is the reference or mutant allele for an insertion
	
	Args:
		inbam (obj): an open pysam.Samfile object
		pos (int): the position of the mutation
		ref (str): the reference allele
		alt (str): the mutant allele
		
	Returns:
		list of [blank/mut/ref, read start, read end]
	"""
	res='blank'
	mutLen=len(alt)-len(ref)
	seq=read.query_alignment_sequence
	mutPos=(pos)-read.pos
	end=read.pos+len(seq)
	inSize=0
	if mutPos < len(seq) and read.pos < pos and read.cigartuples:
		upstream=0
		for thisTuple in read.cigartuples:
			if thisTuple[0]==0:
				upstream=upstream+thisTuple[1]
			if thisTuple[0]==1:
				inSize=thisTuple[1]
				break
		if upstream==mutPos and inSize==mutLen:
			res='mut'
		else:
			res='ref'
	return([res, read.pos, end-inSize])

def collateMutsByReadPosition(inbam, chr, pos, ref, alt):
	"""
	Collects numbers of mutant and reference reads for a mutation
	
	Args:
		inbam (obj): an open pysam.Samfile object
		chr (str): the chromosome of the mmutation
		pos (int): the position of the mutation
		ref (str): the reference allele
		alt (str): the mutant allele
	Returns:
		dictionary of the format {"1000-1101":[1,2], "1002-1102":[10,12]}
		for numbers of ALT and REF reads for each transcript
	"""
	resD={}
	for read in inbam.fetch(chr, pos-2, pos):
		if len(ref)==1 and len(alt)==1:
			thisRes=getSNVfromRead(read, pos, ref, alt)
		elif len(ref)>len(alt):
			thisRes=getDelfromRead(read, pos, ref, alt)
		elif len(alt)>len(ref):
			thisRes=getInsfromRead(read, pos, ref, alt)
		if thisRes[0]!='blank' and thisRes[0]!='other':
			thisPos=str(thisRes[1])+'-'+str(thisRes[2])
			if thisPos not in resD:
				resD[thisPos]=[0,0]
			if thisRes[0]=='mut':
				resD[thisPos][0]+=1
			else:
				resD[thisPos][1]+=1
	return(resD)

def binomProp(inbam, chr, pos, ref, alt, expAlt, pCutoff, benferroni=True):
	"""
	Performs binomial tests on all PCR templates supporting the REF and ALT genotype for a mutation
	
	Args:
		inbam (obj): an open pysam.Samfile object
		chr (str): the chromosome of the mmutation
		pos (int): the position of the mutation
		ref (str): the reference allele
		alt (str): the mutant allele
		expAlt (float): the expected VAF to test against for all templates
		pCutoff (float): p-value cutoff to report a transcript as unbalanced
		benferroni (store_true): whether to carry out multiple testing on each template
	Returns:
		list of [mutant transcripts, unbalanced mutant transcripts, reference transcripts, unbalance reference transcripts]
	"""
	thisD=collateMutsByReadPosition(inbam, chr, pos, ref, alt)
	if thisD=={}:
		return([1,1,1,1])
	mutPsig=0
	mutTot=0
	refPsig=0
	refTot=0
	if benferroni:
		p2use=pCutoff/len(thisD)
	else:
		p2use=pCutoff
	for i in thisD:
		thisP=float(stats.binomtest(k=thisD[i][0], n=sum(thisD[i]), p=expAlt).pvalue)
		mutTot=mutTot+thisD[i][0]
		refTot=refTot+thisD[i][1]
		if thisP < p2use:
			mutPsig=mutPsig+thisD[i][0]
			refPsig=refPsig+thisD[i][1]
	return([mutTot, mutPsig, refTot, refPsig])

def checkIndex(vcf, bam, makeBg=False):
	"""
	Checks if the input vcf file is BGzipped and tabix indexed.
	Then checks if the bam file is index, indexing if not.
	BGzipping will be done with the -t option, tabix indexing will be done if a BGzipped file is available.
	
	Args:
		vcf (str): a vcf file to test
		bam (str): a bam file to test
		makeBg (store_true): whether to BGzip the vcf if not already done - this will mean the unzipped file will no longer be avaiable
	
	Returns:
		store_true: whether a BGzipped file has been made, which changes the file path for the main script
	"""
	p=Popen('ls ' + vcf + '*', shell=True, stdout=PIPE, stderr=PIPE)
	stdout, stderr = p.communicate()
	files=stdout.decode().split('\n')
	bgzip=False
	tabix=False
	for i in files:
		if 'vcf.gz' in i:
			bgzip=True
		if 'vcf.gz.tbi' in i:
			tabix=True
	if not bgzip and not makeBg:
		print('The input vcf is unzipped and unindexed - use the --tabix option to compress and index it.\n***WARNING*** If you need the unzipped veriosn, make a copy.', file=sys.stderr)
		quit()
	elif not bgzip:
		pysam.tabix_index(vcf, preset='vcf', force=True)
	elif not tabix:
		pysam.tabix_index(vcf, preset='vcf', force=True)
	p=Popen('ls ' + bam.replace('.bam','b*'), shell=True, stdout=PIPE, stderr=PIPE)
	stdout, stderr = p.communicate()
	files=stdout.decode().split('\n')
	if bam+'.bai' not in files and bam.replace('.bam','.bai') not in files:
		pysam.index(bam)
	return(bgzip)

def main(theseArgs):
	foundGz=checkIndex(theseArgs.input, theseArgs.bam, theseArgs.tabix)
	if foundGz:
		updateVcf(theseArgs.input, theseArgs.output, theseArgs.bam, theseArgs.refTag, theseArgs.varTag, refVarTag=theseArgs.refVarTag, varRefTag=theseArgs.varRefTag, filterThresh=theseArgs.filter)
	else:
		updateVcf(theseArgs.input+'.gz', theseArgs.output, theseArgs.bam, theseArgs.refTag, theseArgs.varTag, refVarTag=theseArgs.refVarTag, varRefTag=theseArgs.varRefTag, filterThresh=theseArgs.filter)

if __name__ == '__main__':
	main(args)
