import subprocess

class cancer_vcf():
	def __init__(self):
		# start of  command for samtools
		self.faidx_cmd="/usr/bin/samtools faidx /media/MokaNAS2/projects/161012_genome.fa/genome.fa "
		# cosmic vcf file
		self.input_file="/media/MokaNAS2/projects/171127_cancer_vcf/Accel-Amplicon-56G-ClinVar-dbSNP-COSMIC-annotated-target-list.csv"
		# output 4 column file
		self.output_file="/media/MokaNAS2/projects/171127_cancer_vcf/Accel-Amplicon-56G-table.txt"

	def read_script(self):
		"""
		read the (modified) cosmic variants
		for each line produce 2 lines in output:
			one with
			gene, position (chr:pos), variant name (wt) and then the 10bp before the variant, the wt base and the 10bp after
			and one with
			gene, position (chr:pos), variant name (AA_gene_transcriptid) and then the 10bp before the variant, and 11 bases including the alt bases and any wildtype sequences
		this is written to a output file
		"""
		#open the output file
		with open(self.output_file,'w') as output_file:
			#open input file
			with open(self.input_file,'r') as input_list:
			#get to variants in vcf
				file_list=input_list.readlines()
			
			# some counts to describe the success
			bad_count=0
			good_count=0
			
			# for each variant in the list
			for line in file_list:
				#skip the header
				if line.startswith("#"):
					pass
				else:
					#split the line to capture what we need
					splitline=line.split("\t")
					chr=splitline[0]
					pos=int(splitline[1])
					ref=splitline[3]
					alt=splitline[4]
					AA=splitline[7]
					gene=splitline[8]
					transcriptid=splitline[11]
			
					# build variant name
					variant_name=AA+"_"+gene+"_"+transcriptid

					# get wt sequence
					wt_seq=self.get_sequence(chr,pos-10,pos+10)
					
					# for insertions and subsitutions the len(alt) is greater or equal to ref:
					if len(alt) >= len(ref) :
						#get 10bp before alt, not including alt
						seq1=self.get_sequence(chr,pos-10,pos-1)
						#get the 10bp of sequence after (and not including) alt
						seq2=self.get_sequence(chr,pos+len(alt),pos+10)
						# we want the 10bp before a variant, the variant and the 10bp after.
						# in some cases alt is > 1.
						# in these cases we need less of the trailing sequence
						# a maximum of 10bp after the sequence should be unique (total sequence of 21bp)
						# nb this may result in truncating the alt
						combined_alt_after=alt.rstrip()+seq2
						#write seq before alt and  11 bases of the alt and  seq after alt
						sequence_to_write=seq1.rstrip()+combined_alt_after[0:11]
						
					
					#deletions require a different search term
					elif len(alt) < len(ref):
						#get 10bp before ref, not including ref/alt
						seq1=self.get_sequence(chr,pos-10,pos-1)
						#get the 10bp of sequence after (and not including) alt
						seq2=self.get_sequence(chr,pos+len(ref),pos+len(ref)+10)
						# we want the 10bp before a variant, the variant and the 10bp after.
						# in some cases len(alt) is > 1.
						# in these cases we need less of the trailing sequence
						# a maximum of 10bp after the sequence should be unique (total sequence of 21bp)
						# nb this may result in truncating the alt
						combined_alt_after=alt.rstrip()+seq2
						#write seq before alt, alt and seq after alt
						sequence_to_write=seq1.rstrip()+combined_alt_after[0:11]
					
					# sometimes faidx cannot return the sequence. the returned result is the fasta header, containing a ">"
					#print and count the seq where the sequence can't be found
					if ">" in seq1 and ">" in seq2:
						bad_count+=1
						print "no sequence found for:"+str(chr)+":"+str(pos-10)+"-"+str(pos-1) + " and no sequence found for:"+str(chr)+":"+str(pos+len(alt))+"-"+str(pos+10)
					elif ">" in seq1:
						bad_count+=1
						print "no sequence found for:"+str(chr)+":"+str(pos-10)+"-"+str(pos-1)
					elif ">" in seq2:
						bad_count+=1
						print "no sequence found for:"+str(chr)+":"+str(pos+len(alt))+"-"+str(pos+10)
					else:
						# otherwise write the lines to the output file
						output_file.write(gene+"\t"+chr+":"+str(pos)+"\t"+AA+"_"+gene+"_"+transcriptid+"\t"+sequence_to_write+"\n")
						output_file.write(gene+"\t"+chr+":"+str(pos)+"\twt\t"+wt_seq)
						good_count+=1

			print "errors with " + str(bad_count) + " samples"
			print "good samples = " + str(good_count)

	def get_sequence(self,chr,start,stop):
		""" receive genomic range and return the sequence """
		# run faidx and return only the last line
		cmd=self.faidx_cmd+chr+":"+str(start)+"-"+str(stop)+" | tail -1"
		
		# return the sequence
		return subprocess.check_output(cmd,shell=True)
		
def main():
	c=cancer_vcf()
	c.read_script()

if __name__ =="__main__":
	main()
