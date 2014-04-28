# LIVE EDIT FROM WEB, GREMBI'S ACCOUNT
import os
import sys

class ScriptTracker:
    ''' Koshlan Mayer-Blackwell (c) April 2014
    A major challenge for me in bioinformatics is the creation of files that can be interpretted long after the fact. I introduced this to solve that problem. 
    This Class Is For Creating Organization When Running A Script. It works by creating a time stamped folder where outputs and inputs can be reliablly stored. Functions can be reused. Other Useful Blast Functions Are All Included as built in functions'''
    def __init__(self):
        import sys
        import os
        from datetime import datetime
        unique_timestamp = datetime.now().strftime('%Y-%m-%d_%H_%M_%S')
        self.unique_timestamp = unique_timestamp
        self.fn_log = "%s.log" %(unique_timestamp)
        import os 
        os.system("mkdir %s" %(unique_timestamp))
        os.system("mkdir %s/output" %(unique_timestamp))
        os.system("mkdir %s/process" %(unique_timestamp))
        os.system("mkdir %s/input" %(unique_timestamp))
        os.system("mkdir %s/log" %(unique_timestamp))
        os.system('echo \n >> %s/log/%s' %(self.unique_timestamp,self.fn_log))
        self.append_to_log(" ".join(map(str, sys.argv)))
                
    def append_to_file(self,text_to_write, destination_file):
        '''appends text to the end of a file'''
        oh = open(destination_file, 'a') 
        oh.write("%s\t"%(text_to_write))
        oh.close()
    def append_to_log(self, text):
        self.append_to_file(text, self.unique_timestamp + "/log/" + self.fn_log)
    def report_to_std_error(self,text):
        import sys
        sys.stderr.write(text)
    def report_and_append_log(self, text):
        self.report_to_std_error(text)
        self.append_to_log(text)
    def dump_entire_file_to_log(self, file):
        fh = open(self, file)
        oh = open(self.fn_log, 'a')
        for line in fh:
            oh.write(line)
        fh.close()
        oh.close()   
    def copy_files_to_dir(self,file, directory):
        import os
        if file.find("/") != -1:
            filename = file.split("/")[-1]
        else:
            filename = str(file)
        os.system('cp %s %s/%s' %(file, directory, filename))
    def copy_to_inputs(self, file):
        self.copy_files_to_dir(file, self.unique_timestamp+'/input')
    def copy_to_outputs(self, file):
        self.copy_files_to_dir(file, self.unique_timestamp+'/output')
    def copy_to_process(self, file):
        self.copy_files_to_dir(file, self.unique_timestamp+'/process') 
    def make_fasta_dictionary(self, fn_fasta):
        from Bio import SeqIO
        handle = open(fn_fasta, "rU")
        record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
        handle.close()
        return record_dict
    def make_subfasta(self, list, seq_dictionary, destination_file):
        oh = open(destination_file, 'w')
        for x in list:
            oh.write(">%s\n%s\n"%(seq_dictionary[x].description, seq_dictionary[x].seq))
        oh.close()
    def my_makeblastdb(self, input_fasta, type, output_name):
        '''type can be either 'prot' or 'nucl' '''
        import os
        os.system("makeblastdb -in %s -dbtype %s -out %s" %(input_fasta, type, output_name))
    def my_blastp(self, database, search_fasta, output_result):
        import os
        os.system("blastp -db %s -query %s -out %s -outfmt 6" %(database, search_fasta, output_result))
    def sort_blast(self, fn_blast_output_name, fn_blast_best_hit):
         import sys
         import os
         os.system('sort -k1,1 -k12,12gr -k11,11g  %s | sort -u -k1,1 --merge > %s' %(fn_blast_output_name, fn_blast_best_hit))
    def quick_genome_informed_network(self, all_v_all_tabular_blast, lower_PID_bound):
        D = {}
        Dr = {}
        import sys
        import math
        fh = open(all_v_all_tabular_blast, 'r')
        oh = open(all_v_all_tabular_blast + "_" + str(lower_PID_bound) + "_genome_informed_network" , 'w')
        for line in fh:
            if line.find("Warning") != -1:
                continue
            line = line.strip().split()
            try:
                query = line[0]
                query_genome = query.split("_")[0] 
                hit = line[1]
                hit_genome = hit.split("_")[0] 
                pid = math.floor(float(line[2]))
                allign_length = line[3]
            except IndexError or ValueError:
                continue
            if allign_length < 100:
                continue
            if query not in D.keys():
                D[query] = {}
            if hit not in Dr.keys():
                Dr[hit] = {}
           
            token = False
            if query_genome not in Dr[hit].keys():
                Dr[hit][query_genome] = True
                token = True
                token2 = False
            if hit_genome not in D[query].keys():
                D[query][hit_genome] = True
                token2 = True
            if token == True and token2 == True:
                if pid > lower_PID_bound and hit_genome != query_genome:
                    oh.write("%s\t%s\t%s\n"%(query, hit, pid)) 
        fh.close()
        oh.close()
