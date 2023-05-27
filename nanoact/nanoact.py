# %%
from subprocess import Popen, PIPE, run, check_output
import gzip
import hdbscan
import numpy as np
import edlib
from matplotlib import pyplot as plt
from sklearn.manifold import MDS
import re
import shutil
from requests import get, post
import xmltodict
import pandas as pd
import json
import os
from collections import Counter
import urllib.parse
import time

# %%
class NanoAct():
    def __init__(self, TEMP = './temp/'):
        self.TEMP = TEMP
    def _lib_path(self):
        #get the path of the library
        return os.path.dirname(os.path.realpath(__file__))
    def _exec(self, cmd,suppress_output=True):
        if suppress_output:
            with open(os.devnull, 'w') as DEVNULL:
                out = run(cmd, stdout=DEVNULL, stderr=DEVNULL, shell=True)
            return None
        else:
            out = run(cmd, stdout=PIPE, shell=True)
            print (">>", cmd)
            print("Output:")
            print(out.stdout.decode('utf-8'))
            print("Exception:")
            print(out.stderr)
            return out.stdout.decode('utf-8'), out.stderr
    def _clean_temp(self):
        #Clean the temp folder
        try:
            shutil.rmtree(self.TEMP)
            os.mkdir(self.TEMP)
        except FileNotFoundError:
            os.mkdir(self.TEMP)
    def _exec_rt (self, cmd, prefix=""):
        p = Popen(cmd, stdout=PIPE, shell=True)
        for line in iter(p.stdout.readline, b''):
            print ('{}>>> {}'.format(prefix, line.rstrip()))
    def _IUPACde(self, seq):
        seq = str(seq)
        seq = seq.replace('R','[AG]')
        seq = seq.replace('Y','[CT]')
        seq = seq.replace('S','[GC]')
        seq = seq.replace('W','[AT]')
        seq = seq.replace('K','[GT]')
        seq = seq.replace('M','[AC]')
        seq = seq.replace('B','[CGT]')
        seq = seq.replace('D','[AGT]')
        seq = seq.replace('H','[ACT]')
        seq = seq.replace('V','[ACG]')
        seq = seq.replace('N','[ACGT]')
        return seq
    def _extract_degenerate_seq(self, seq):
        return list(self.sre_yield.AllStrings(self.IUPACde(seq.upper())))  # better
    """
    def _align_two_seq(self, seq1,seq2):
        max_score = -99
        max_alignment = None
        seq1 = self.IUPACde(str(seq1).upper())
        seq2 = self.IUPACde(str(seq2).upper())
        for s1 in self.extract_degenerate_seq(seq1):
            for s2 in self.extract_degenerate_seq(seq2):
                #Align two sequence
                aligner = Align.PairwiseAligner()
                aligner.mode = 'local'
                aligner.match_score = 1
                aligner.mismatch_score = -1
                aligner.open_gap_score = -1
                aligner.extend_gap_score = -1
                alignments = aligner.align(s1, s2)
                #Get the best alignment
                for alignment in alignments:
                    if alignment.score > max_score:
                        max_score = alignment.score
                        max_alignment = alignment
                    break
        return max_alignment,max_score
    """
    def _fastq_reader(self, handle):
        #Custom fastq reader, which can handle the case when the quality score is inconsistent with the sequence length
        while True:           
            line = handle.readline()
            if not line:
                return
            if line[0] != "@":
                raise ValueError("Records in Fastq files should start with '@'")
            title = line[1:].rstrip()
            seq = handle.readline().rstrip()
            handle.readline() # skip the line starting with "+"
            qual = handle.readline().rstrip()
            if len(qual) != len(seq):
                print("Inconsistency found", title)
                #Temporary workaround (return fake qual)
                yield {"title": title, "seq": seq, "qual": "A"*len(seq)}
                pass
            else:
                yield {"title": title, "seq": seq, "qual": qual}
    def _fastq_writer(self,title,seq,qual,handle):
        handle.write("@{}\n{}\n+\n{}\n".format(title,seq,qual))
    def _fasta_reader(self, handle):
        #Custom fasta reader
        line = handle.readline()
        while True:
            if not line:
                return
            if line[0] == ">":
                title = line[1:].rstrip()
                seq = ""
                while True:
                    line = handle.readline()
                    if not line or line[0] == ">":
                        yield {"title": title, "seq": seq}
                        break
                    else:
                        seq += line.rstrip()
    def _count_seq_num (self,fastq):
        with open(fastq,'r') as count_len:
            a = count_len.readlines()    
        return len(a)/4
    def _fastq_rename_title(self, src,des): #modified seq title in fastq to number
        with open(src, "r") as f:
            with open(des, "w") as f2:
                counter = 0
                output = ""
                for line in f:
                    if line.startswith("@"):
                        output += "@seq" + str(counter) + "\n"
                        counter += 1
                    else:
                        output += line
                f2.write(output)
    def _pairwise_distance(self, s1,s2):
        #Calculate the pairwise distance between two sequences
        #This is a wrapper for edlib.align
        return edlib.align(str(s1).upper(), str(s2).upper())['editDistance']/(min(len(s1), len(s2))+0.1)*100
    def _reverse_complement(self, s):
        s = list(s)
        for pos, nuc in enumerate(s):
            if nuc == 'A':
                s[pos] = 'T'
            elif nuc == 'T':
                s[pos] = 'A'
            elif nuc == 'C':
                s[pos] = 'G'
            elif nuc == 'G':
                s[pos] = 'C'
            elif nuc == 'a':
                s[pos] = 't'
            elif nuc == 't':
                s[pos] = 'a'
            elif nuc == 'c':
                s[pos] = 'g'
            elif nuc == 'g':
                s[pos] = 'c'
        return ''.join(s[::-1])
    def NCBIblast(self, seqs = ">a\nTTGTCTCCAAGATTAAGCCATGCATGTCTAAGTATAAGCAATTATACCGCGGGGGCACGAATGGCTCATTATATAAGTTATCGTTTATTTGATAGCACATTACTACATGGATAACTGTGG\n>b\nTAATACATGCTAAAAATCCCGACTTCGGAAGGGATGTATTTATTGGGTCGCTTAACGCCCTTCAGGCTTCCTGGTGATT\n" ):
        program = "blastn&MEGABLAST=on"
        database = "nt"
        encoded_queries = urllib.parse.quote(seqs)
        WORD_SIZE = 32
        EXPECT = 0.001
        # build the request
        args = "CMD=Put&PROGRAM=" + program + "&DATABASE=" + database + "&QUERY=" + encoded_queries + "&WORD_SIZE=" + str(WORD_SIZE) + "&EXPECT=" + str(EXPECT)
        url = 'https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi'
        response = post(url, data=args)
        #print("BLASTING {} sequences".format(len(seqs.split(">"))-1))
        # parse out the request id
        rid = ""
        for line in response.text.split('\n'):
            if line.startswith('    RID = '):
                rid = line.split()[2]
        #Search submitted
        print("Query", rid, "submitted.")
        print("You can check the status at https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=" + rid + "")
        print("And results here: https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&RID=" + rid + "")
        # poll for results
        retry = 30
        while True:
            time.sleep(30)
            retry -= 1
            if retry == 0:
                print("Search", rid, "timed out")
                return None
            url = 'https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=' + rid
            response = get(url)
            if 'Status=WAITING' in response.text:
                print("Searching...")
                continue
            if 'Status=FAILED' in response.text:
                print("Search", rid, "failed; please report to blast-help@ncbi.nlm.nih.gov.")
            if 'Status=UNKNOWN' in response.text:
                print("Search", rid, "expired.")
            if 'Status=READY' in response.text:
                if 'ThereAreHits=yes' in response.text:
                    print("Search complete, retrieving results...")
                    break
                else:
                    print("No hits found.")
        # retrieve and display results
        url = f'https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&RID={rid}&FORMAT_TYPE=XML'
        print("Retrieving results from", url)
        response = get(url)
        #Convert the XML to a dictionary
        blast_dict = xmltodict.parse(response.text)
        #Get the first hit of each query
        pool = {}

        #If there is only one query, the xml format is different
        #blast_dict['BlastOutput']['BlastOutput_iterations']['Iteration'] will be a dict instead of a list
        #So we need to convert it to a list
        if type(blast_dict['BlastOutput']['BlastOutput_iterations']['Iteration']) == dict:
            blast_dict['BlastOutput']['BlastOutput_iterations']['Iteration'] = [blast_dict['BlastOutput']['BlastOutput_iterations']['Iteration']]

        for rec in blast_dict['BlastOutput']['BlastOutput_iterations']['Iteration']:
            try:
                seq_name = rec['Iteration_query-def']
            except Exception as e:
                print(e)
                continue
            try: 
                hit = rec['Iteration_hits']['Hit'][0]
            except:
                hit = None
            if hit:
                acc = hit['Hit_accession']
                if type(hit['Hit_hsps']['Hsp']) == list:
                    hit_hsp = hit['Hit_hsps']['Hsp'][0]
                else:
                    hit_hsp = hit['Hit_hsps']['Hsp']
                hit_seq = hit_hsp['Hsp_hseq'].replace('-', '')
                hit_def = hit['Hit_def']
                similarity = round(int(hit_hsp['Hsp_identity'])/int(hit_hsp['Hsp_align-len']),2)
                
                pool[seq_name] = {'acc': acc, 'hit_seq': hit_seq, 'hit_def': hit_def, 'similarity': similarity, 'org': ""}
                #Get taxon info
                taxon_info_URI = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&retmode=xml&id={}'.format(acc)
                try:
                    r = get(taxon_info_URI)
                    #xml to json
                    r = xmltodict.parse(r.text)
                    org = r['GBSet']['GBSeq']['GBSeq_organism']
                    #taxa = r['GBSet']['GBSeq']['GBSeq_taxonomy']
                    #Get taxid in db_xref
                    #Get db_xref
                    taxid=""
                    for i in r['GBSet']['GBSeq']['GBSeq_feature-table']['GBFeature']:
                        if i['GBFeature_key'] == 'source':
                            for j in i['GBFeature_quals']['GBQualifier']:
                                if j['GBQualifier_name'] == 'db_xref':
                                    taxid = j['GBQualifier_value'].split(':')[-1]
                    pool[seq_name].update({'org': org, 'taxid': taxid})
                except Exception as e:
                    print(e)
                    pass
                #Get all ranks
                ranks = {"kingdom":"incertae sedis", "phylum":"incertae sedis", "class":"incertae sedis", "order":"incertae sedis", "family":"incertae sedis", "genus":"incertae sedis"}
                try:
                    #Get details from taxid
                    if taxid == "":
                        continue
                    taxid_info_URI = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&rettype=xml&id={}'.format(taxid)
                    r = get(taxid_info_URI)
                    r = xmltodict.parse(r.text)
                    for i in r['TaxaSet']['Taxon']['LineageEx']['Taxon']:
                        if i['Rank'] in ranks.keys():
                            ranks[i['Rank']] = i['ScientificName']
                    pool[seq_name].update(ranks)
                except Exception as e:
                    print(e)
                    pass
        return pool
    def blast_2 (self, src, des, name="blast.csv", funguild = True, startswith="con_", max_query_length=500, batch = 5):
        #Collect all sequences
        pool_df = pd.DataFrame() 
        query_seqs=[] 
        for f in os.scandir(src):
            if f.name.startswith(startswith) and f.name.endswith(".fas"):
                with open(f.path, 'r') as handle:
                    seqs = list(self._fasta_reader(handle))
                    for s in seqs:
                        pool_df = pd.concat([pool_df, pd.DataFrame([s])], ignore_index=True)
                        if len(s['seq']) >= max_query_length:
                            #If sequence is too long, preserve only max_query_length in middle
                            diff = len(s['seq']) - max_query_length
                            s['seq'] = s['seq'][int(diff/2):int(diff/2)+max_query_length]
                        query_seqs.append(f">{s['title']}\n{s['seq']}")       
        #set title as index
        pool_df.set_index('title', inplace=True)
        for index, row in pool_df.iterrows():
            pool_df.loc[index, 'length'] = str(len(row['seq']))
            try:
                #2110_cluster_-1_r2154.fas	
                #{sample}_cluster_{cluster_no}_r{reads_count}.fas
                sample, cluster_no, reads_count = re.search("(.*)_cluster_([-0-9]+)_r(\d+).fas", index).groups()
                pool_df.loc[index, 'sample'] = sample
                pool_df.loc[index, 'cluster_no'] = cluster_no
                pool_df.loc[index, 'reads_count'] = reads_count
            except Exception as e:
                print(e)
                pass    
        #Blast all sequences
        i = 0
        blast_result_pool = {}
        while i < len(query_seqs):
            print("Blasting", i, "to", i+batch, "of", len(query_seqs))
            query = "\n".join(query_seqs[i:i+batch])
            retry = 3
            while retry >= 0:
                try:
                    blast_result = self.NCBIblast(query)
                    if blast_result != None:
                        blast_result_pool.update(blast_result)
                        break
                    else:
                        retry -= 1
                except Exception as e:
                    print(e)
                    retry -= 1
            i+=batch
        #print(blast_result_pool)
        for sample in blast_result_pool.keys():
            for key in blast_result_pool[sample].keys():
                pool_df.loc[sample,key] = blast_result_pool[sample][key]

            #Check funguild

            if funguild and blast_result_pool[sample]['org'] != "":
                funguild_des = []
                try:
                    funguild_des = json.loads(get(f"https://www.mycoportal.org/funguild/services/api/db_return.php?qDB=funguild_db&qField=taxon&qText={blast_result_pool[sample]['org']}").text)
                except:
                    pass
                if funguild_des != []:
                    #print(funguild_des)
                    try:
                        pool_df.loc[sample,'funguild'] = funguild_des[0]['guild']
                    except:
                        pass
                    try:
                        pool_df.loc[sample,'funguild_notes'] = funguild_des[0]['notes']
                    except:
                        pass
        pool_df.to_csv(f"{des}/{name}", encoding ='utf-8-sig')
        return f"{des}/{name}"

    def _mafft (self, src, des):
        mafft_bin = self._lib_path() + "/bin/mafft.bat"
        #./mafft.bat --genafpair --maxiterate 1000 2110_cluster_1_r442.fas > output.fas
        cmd = f"{mafft_bin} {src} > {des}"
        self._exec(cmd,suppress_output=True)
    def orientation(self, src, des, tsv,search_region=200):
        #Input: a folder containing all the fas files, fas_file should be named as {sample_id}.fas
        #Input2: a barcode index file, containing following columns: SampleID, FwPrimer, RvPrimer
        #Search_region: number of bases in the beginning of raw read to search for primer
        #Output: a folder containing sequences with the right orientation
        try:
            os.makedirs(des, exist_ok=True)
        except:
            pass
        bar_idx = pd.read_csv(tsv, sep='\t')
        #Read fasta file
        for f in os.scandir(src):
            if f.is_file() and f.name.endswith(".fas"):
                sample_id = f.name[:-4]
                try:
                    F = bar_idx[bar_idx['SampleID'].astype(str)==sample_id]['FwPrimer'].values[0]
                    R = bar_idx[bar_idx['SampleID'].astype(str)==sample_id]['RvPrimer'].values[0]
                except IndexError:
                    print(f"Sample {sample_id} not found in the barcode file, skipping")
                    continue
                with open(f.path) as handle:
                    with open(f"{des}/{f.name}", "w") as output:
                        for record in self._fasta_reader(handle):
                            #Check if the sequence is in the right orientation
                            aln_f = edlib.align(F.upper(), record['seq'].upper()[:search_region], mode="HW", task="locations")
                            aln_r = edlib.align(R.upper(), record['seq'].upper()[:search_region], mode="HW", task="locations")
                            if aln_f['editDistance'] > aln_r['editDistance']:
                                record['seq'] = self._reverse_complement(record['seq'])
                            output.write(f">{record['title']}\n{record['seq']}\n")
        return des
    def mafft_consensus (self, src, des):
        try:
            os.makedirs(des)
        except:
            pass
        abs_des = os.path.abspath(des)
        for f in os.scandir(src):
            if f.is_file() and f.name.endswith(".fas"):
                #Align sequences
                print("Working on", f.name, "...")
                self._mafft(f.path, f"{abs_des}/aln_{f.name}")
                #naive consensus
                with open(f"{abs_des}/aln_{f.name}") as handle:
                    records = list(self._fasta_reader(handle))
                    consensus = ""
                    for i in range(len(records[0]['seq'])):
                        col = [r['seq'][i] for r in records]
                        #get_most_common
                        com = Counter(col).most_common(1)[0][0]
                        if com != "-":
                            consensus += com
                    with open(f"{abs_des}/con_{f.name}", "w") as out:
                        out.write(f">{f.name}\n{consensus}")

        return abs_des
 
    def _get_sample_id_single (self, seq, barcode_hash_table, mismatch_ratio_f = 0.15,mismatch_ratio_r = 0.15):
        # Define a helper function to identify the sample ID of a sequence read based on its barcode
        ids = []
        integrity = []
        seqs = []
        # Convert the input sequence read to uppercase and extract the beginning and end segments
        seq = seq.upper()
        seqF = seq[:150]
        seqR = seq[-150:]
        
        for id in barcode_hash_table:
            # Iterate through the hash table of barcodes and their corresponding index sequences
            FwIndex = barcode_hash_table[id]["FwIndex"].upper()
            RvAnchor = barcode_hash_table[id]["RvAnchor"].upper()
            # Extract the forward and reverse index sequences of the current barcode and convert to uppercase
            FwIndex_check = edlib.align(FwIndex,seqF, mode="HW", k=int(len(FwIndex) * mismatch_ratio_f), task="locations")
            RvAnchor_check = edlib.align(RvAnchor,seqR, mode="HW", k=int(len(RvAnchor) * mismatch_ratio_r), task="locations")
            # Align the forward and reverse index sequences with the beginning and end of the input sequence read, allowing for a certain number of mismatches defined by the mismatch ratio
            if FwIndex_check["editDistance"] != -1:
                #mark founded region to lower case
                seqF = seqF[:FwIndex_check["locations"][0][0]] + seqF[FwIndex_check["locations"][0][0]:FwIndex_check["locations"][0][1]].lower() + seqF[FwIndex_check["locations"][0][1]:]
                if RvAnchor_check["editDistance"] != -1:
                    #mark founded region to lower case
                    seqR = seqR[:RvAnchor_check["locations"][0][0]] + seqR[RvAnchor_check["locations"][0][0]:RvAnchor_check["locations"][0][1]].lower() + seqR[RvAnchor_check["locations"][0][1]:]
                    ids.append(id)
                    integrity.append(True)
                else:
                    ids.append(id)
                    integrity.append(False)
                seqs.append(seqF +seq[150:-150] + seqR)
        # If a barcode is identified, return the corresponding sample ID and a boolean indicating whether the barcode was matched with sufficient integrity
        return ids, integrity, seqs
    def _get_sample_id_dual (self, seq, barcode_hash_table, mismatch_ratio_f = 0.15, mismatch_ratio_r = 0.15):
        ids = []
        seqs = []
        # Convert the input sequence read to uppercase and extract the beginning and end segments
        seq = seq.upper()
        seqF = seq[:150]
        seqR = seq[-150:]
        seq_REV = self._reverse_complement(seq)
        seq_REV_F = seq_REV[:150]
        seq_REV_R = seq_REV[-150:]
        
        for id in barcode_hash_table:
            # Iterate through the hash table of barcodes and their corresponding index sequences
            FwIndex = barcode_hash_table[id]["FwIndex"].upper()
            RvIndex = barcode_hash_table[id]["RvIndex"].upper()
            #Reverse complement the reverse index sequence
            RvIndex = self._reverse_complement(RvIndex)
            # Extract the forward and reverse index sequences of the current barcode and convert to uppercase
            FwIndex_check = edlib.align(FwIndex,seqF, mode="HW", k=int(len(FwIndex) * mismatch_ratio_f), task="locations")
            RvIndex_check = edlib.align(RvIndex,seqR, mode="HW", k=int(len(RvIndex) * mismatch_ratio_r), task="locations")
            #Check read in reverse complement
            Fwindex_check_R = edlib.align(FwIndex,seq_REV_F, mode="HW", k=int(len(FwIndex) * mismatch_ratio_f), task="locations")
            Rvindex_check_R = edlib.align(RvIndex,seq_REV_R, mode="HW", k=int(len(RvIndex) * mismatch_ratio_r), task="locations")
            #print(FwIndex,seqF,RvIndex,seqR)
            #print("FwIndex_check",FwIndex_check, "RvIndex_check",RvIndex_check)
            if FwIndex_check["editDistance"] != -1 and RvIndex_check["editDistance"] != -1:
                #mark founded region to lower case

                seqF = seqF[:FwIndex_check["locations"][0][0]] + seqF[FwIndex_check["locations"][0][0]:FwIndex_check["locations"][0][1]].lower() + seqF[FwIndex_check["locations"][0][1]:]
                seqR = seqR[:RvIndex_check["locations"][0][0]] + seqR[RvIndex_check["locations"][0][0]:RvIndex_check["locations"][0][1]].lower() + seqR[RvIndex_check["locations"][0][1]:]
                
                ids.append(id)
                seqs.append(seqF +seq[150:-150] + seqR)
            elif Fwindex_check_R["editDistance"] != -1 and Rvindex_check_R["editDistance"] != -1:
                #mark founded region to lower case
                seq_REV_F = seq_REV_F[:Fwindex_check_R["locations"][0][0]] + seq_REV_F[Fwindex_check_R["locations"][0][0]:Fwindex_check_R["locations"][0][1]].lower() + seq_REV_F[Fwindex_check_R["locations"][0][1]:]
                seq_REV_R = seq_REV_R[:Rvindex_check_R["locations"][0][0]] + seq_REV_R[Rvindex_check_R["locations"][0][0]:Rvindex_check_R["locations"][0][1]].lower() + seq_REV_R[Rvindex_check_R["locations"][0][1]:]
                ids.append(id)
                seqs.append(seq_REV_F +seq[150:-150] + seq_REV_R)
        return ids, seqs

    def singlebar(self, src, des, BARCODE_INDEX_FILE, mismatch_ratio_f = 0.15,mismatch_ratio_r = 0.15, expected_length_variation = 0.3):
        """
        1.Process the raw sequencing file using a fastq_reader function which reads in four lines at a time representing one sequencing read.
        2.For each read, call the _get_sample_id_single function to identify its corresponding sample ID and barcode integrity.
        3.If a sample ID is identified, append the read to the corresponding sample's output file.
        4.If a barcode is truncated, append it to the "TRUNCATED" dictionary.
        5.If multiple barcodes are identified, append it to the "MULTIPLE" dictionary.
        6.If no barcode is identified, append it to the "UNKNOWN" dictionary.
        7.Finally, output the demultiplexed reads into different output files based on their barcode.
        """
        # Define the main function for demultiplexing
        if BARCODE_INDEX_FILE.endswith("tsv"):
            sep = "\t"
        elif BARCODE_INDEX_FILE.endswith("csv"):
            sep = ","
        else:
            raise ValueError("BARCODE_INDEX_FILE must be a tsv or csv file")
        BARCODE_IDX_DF = pd.read_csv(BARCODE_INDEX_FILE, sep=sep)
        # Read in the barcode index file as a pandas DataFrame
        if not all([x in BARCODE_IDX_DF.columns for x in ["SampleID", "FwIndex", "RvAnchor", "ExpectedLength"]]):
            raise ValueError("BARCODE_INDEX_FILE must have SampleID, FwIndex, RvAnchor columns")
        print ("BARCODE_INDEX_FILE loaded")
        # Check whether the barcode index file has the required columns
        barcode_hash_table = {}
        for index, row in BARCODE_IDX_DF.iterrows():
            barcode_hash_table[row["SampleID"]] = {"FwIndex": str(row["FwIndex"]), "RvAnchor": str(row["RvAnchor"]), "ExpectedLength": row["ExpectedLength"]}
        # Store the barcode index file as a hash table of barcode-sample ID pairs
        pool = {}
        counter = 0
        pool["UNKNOWN"] = []
        pool["MULTIPLE"] = []
        pool["TRUNCATED"] = []
        pool["IncorrectLength"] = []
        for id in barcode_hash_table:
            pool[id] = []
        # Initialize dictionaries for output files for each sample and additional dictionaries for reads with unknown, multiple, or truncated barcodes
        with open(src, "r") as handle:
            for record in self._fastq_reader(handle):
                ids, integrity, seqs= self._get_sample_id_single(record["seq"], barcode_hash_table, mismatch_ratio_f, mismatch_ratio_r)
                if len(ids) == 1:
                    #if only one barcode is identified, append the read to the corresponding sample's output file
                    record["seq"] = seqs[0]
                    if integrity[0] == False:
                        pool["TRUNCATED"].append(record)
                    else:
                        #Check if seq in ExpectedLength
                        if (len(seq) - barcode_hash_table[ids[0]]["ExpectedLength"])**2 < (barcode_hash_table[ids[0]]["ExpectedLength"] * expected_length_variation)**2:
                            pool[ids[0]].append(record)
                        else:
                            pool["IncorrectLength"].append(record)
                elif len(ids) > 1:
                    pool["MULTIPLE"].append(record)
                else:
                    pool["UNKNOWN"].append(record)
                counter += 1
                if counter % 10000 == 0:
                    print(counter)
        #Save to separate fastq file
        try:
            os.makedirs(f"{des}/", exist_ok=True)
            os.makedirs(f"{des}/trash/", exist_ok=True)
        except:
            pass
        stat_df = pd.DataFrame(columns=["SampleID", "Count"])
        for bin in pool:
            stat_df = pd.concat([stat_df, pd.DataFrame({"SampleID": [bin], "Count": [len(pool[bin])]})])
            if bin in ['MULTIPLE', 'UNKNOWN', 'TRUNCATED','IncorrectLength']:
                path = f"{des}/trash/{bin}"
            else:
                path = f"{des}/{bin}"
            #Save to separate fastq file
            with open(path+".fastq", "w") as handle:
                for record in pool[bin]:
                    handle.write(record["title"] + "\n")
                    handle.write(record["seq"] + "\n")
                    handle.write("+" + "\n")
                    handle.write(record["qual"] + "\n")
            #Save to separate fasta file
            with open(path+".fas", "w") as handle:
                for record in pool[bin]:
                    handle.write(">"+ record["title"] + "\n")
                    handle.write(record["seq"] + "\n")
        stat_df.to_csv(f"{des}/2_Singlebar_stat.csv", index=False)
        #Print out the number of reads discarded due to unknown, multiple, or truncated barcodes
        FAILED_NUM = len(pool['UNKNOWN']) + len(pool['MULTIPLE']) + len(pool['TRUNCATED'])
        print (f"{counter-FAILED_NUM}/{counter} ({(counter-FAILED_NUM)/counter*100:.2f}%) reads were demultiplexed successfully")
        return des

    def dualbar(self, src,des, BARCODE_INDEX_FILE, mismatch_ratio_f = 0.15,mismatch_ratio_r = 0.15, expected_length_variation = 0.3):
        """
        1.Process the raw sequencing file using a fastq_reader function which reads in four lines at a time representing one sequencing read.
        2.For each read, call the _get_sample_id_dual function to identify its corresponding sample ID and barcode integrity.
        3.If a sample ID is identified, append the read to the corresponding sample's output file.
        4.If only one barcode is identified, append it to the "SINGLE" dictionary.
        5.If multiple barcodes are identified, append it to the "MULTIPLE" dictionary.
        6.If no barcode is identified, append it to the "UNKNOWN" dictionary.
        7.Finally, output the demultiplexed reads into different output files based on their barcode.
        """
        if BARCODE_INDEX_FILE.endswith(".tsv"):
            sep = "\t"
        elif BARCODE_INDEX_FILE.endswith(".csv"):
            sep = ","
        else:
            raise ValueError("Barcode index file must be either a .tsv or .csv file")
        BARCODE_IDX_DF = pd.read_csv(BARCODE_INDEX_FILE, sep = sep)
        #Check all columns are present
        if not all(x in BARCODE_IDX_DF.columns for x in ["SampleID","FwIndex","RvIndex", "ExpectedLength"]):
            raise ValueError("Barcode index file must contain columns: SampleID, FwIndex, RvIndex, ExpectedLength")
        print("BARCODE_INDEX_FILE loaded")
        # Create a hash table of barcodes and their corresponding index sequences
        barcode_hash_table = {}
        for index, row in BARCODE_IDX_DF.iterrows():
            barcode_hash_table[row["SampleID"]] = {"FwIndex": str(row["FwIndex"]), "RvIndex": str(row["RvIndex"]), "ExpectedLength": row["ExpectedLength"]}
        # Create a dictionary to store the number of reads that are demultiplexed into each sample
        pool = {}
        counter = 0
        # Create a dictionary to store the number of reads that are demultiplexed into each sample
        pool['Unknown'] = []
        pool['Multiple'] = []
        pool['IncorrectLength'] = []
        # Initialize the output files for each sample
        for id in barcode_hash_table.keys():
            pool[id] = []

        with open(src,"r") as handle:
            for record in self._fastq_reader(handle):
                #Get the sequence and quality score
                seq = record["seq"]
                #Get the barcode sequence
                ids, seqs = self._get_sample_id_dual(seq, barcode_hash_table, mismatch_ratio_f, mismatch_ratio_r)
                #Check if the barcode is identified
                if len(ids) == 1:
                    #check if seq in ExpectedLength
                    if (len(seqs[0]) - barcode_hash_table[ids[0]]["ExpectedLength"])**2 < (barcode_hash_table[ids[0]]["ExpectedLength"] * expected_length_variation)**2:
                        pool[ids[0]].append(record)
                    else:
                        pool["IncorrectLength"].append(record)
                elif len(ids) > 1:
                    pool["Multiple"].append(record)
                else:
                    pool["Unknown"].append(record)
                counter += 1
                if counter % 10000 == 0:
                    print(counter)

        #Save to separate fastq file
        try:
            os.makedirs(f"{des}/", exist_ok=True)
            os.makedirs(f"{des}/trash/", exist_ok=True)
        except:
            pass
        stat_df = pd.DataFrame(columns=["SampleID", "Count"])
        for bin in pool:
            stat_df = pd.concat([stat_df, pd.DataFrame({"SampleID": [bin], "Count": [len(pool[bin])]})])
            if bin in ['Multiple', 'Unknown','IncorrectLength']:
                path = f"{des}/trash/{bin}"
            else:
                path = f"{des}/{bin}"
            #Save to separate fastq file
            with open(path+".fastq", "w") as handle:
                for record in pool[bin]:
                    handle.write(record["title"] + "\n")
                    handle.write(record["seq"] + "\n")
                    handle.write("+" + "\n")
                    handle.write(record["qual"] + "\n")
            #Save to separate fasta file
            with open(path+".fas", "w") as handle:
                for record in pool[bin]:
                    handle.write(">"+ record["title"] + "\n")
                    handle.write(record["seq"] + "\n")
        stat_df.to_csv(f"{des}/2_Dualbar_stat.csv", index=False)
        #Print out the number of reads discarded due to unknown, multiple, or truncated barcodes
        FAILED_NUM = len(pool['Unknown']) + len(pool['Multiple']) + len(pool['IncorrectLength'])
        print (f"{counter-FAILED_NUM}/{counter} ({(counter-FAILED_NUM)/counter*100:.2f}%) reads were demultiplexed successfully")
        return des  

    def combine_fastq(self, src, des, name = "all.fastq"):
        try:
            os.makedirs(des, exist_ok=True)
        except Exception as e:
            print(e)
            pass
        with open(f'{des}/{name}', 'w') as outfile:
            for f in os.scandir(src):
                if f.name.endswith(".fastq.gz"):
                    print("Found fastq file: {}".format(f.name))
                    with gzip.open(f.path, 'rt') as infile:
                        for line in infile:
                            outfile.write(line)
        return f'{des}/{name}'.format(self.TEMP)
    def nanoflit(self, src, des, name = "all.fastq", NANOFILT_QSCORE = 8,  NANOFILT_MIN_LEN = 400, NANOFILT_MAX_LEN = 8000):
        try:
            os.makedirs(des, exist_ok=True)
        except Exception as e:
            pass
        print("Start Nanoflit...")
        des += f"/{name}"
        self._exec(f"NanoFilt -q {NANOFILT_QSCORE} --length {NANOFILT_MIN_LEN} --maxlength {NANOFILT_MAX_LEN} {src} > {des}")
        raw_fastq_lines = sum(1 for line in open(src)) /4
        filtered_fastq_line = sum(1 for line in open(des)) /4
        print("Raw reads: {}, Passed: {}({}%)".format(raw_fastq_lines, filtered_fastq_line, int(filtered_fastq_line/raw_fastq_lines*100)))
        return des
    def _average_quality(self, quality_string):
        """
        Calculate the average quality score of a given quality string.

        Args:
            quality_string (str): quality string in FASTQ format.

        Returns:
            average_quality (float): the average quality score.
        """
        quality_scores = [ord(char) - 33 for char in quality_string]
        average_quality = sum(quality_scores) / len(quality_scores)
        return average_quality
    def qualityfilt(self, src, des, name="all.fastq", QSCORE = 8, MIN_LEN = 400, MAX_LEN = 8000):
        try:
            os.makedirs(des, exist_ok=True)
        except Exception as e:
            pass
        print("Start Qualityfilt...")
        des += f"/{name}"
        total = 0
        passed = 0
        with open(src, 'r') as infile:
            with open(des, 'w') as outfile:
                while True:
                    title = infile.readline().strip()
                    if not title:
                        break
                    seq = infile.readline().strip()
                    plus = infile.readline().strip()
                    qual = infile.readline().strip()
                    total += 1
                    if self._average_quality(qual) >= QSCORE and len(seq) >= MIN_LEN and len(seq) <= MAX_LEN:
                        outfile.write(title + "\n")
                        outfile.write(seq + "\n")
                        outfile.write(plus + "\n")
                        outfile.write(qual + "\n")
                        passed += 1
        print(f"{passed}/{total} ({passed/total*100:.2f}%) reads were passed quality filter")
        return des
    def minibar(self, src, des, BARCODE_INDEX_FILE, MINIBAR_INDEX_DIS):
        src = src
        #Check if the barcode index file is valid
        print("Checking barcode index file...")
        _, err = self._exec(f"minibar.py {BARCODE_INDEX_FILE} -info cols")
        if err:
            raise Exception("Invalid barcode index file")
        else:
            print(f"{BARCODE_INDEX_FILE} is valid")
        cwd = os.getcwd()
        os.chdir(des)
        out, err = self._exec(f"minibar.py -F -C -e {MINIBAR_INDEX_DIS} {BARCODE_INDEX_FILE} {src} 2>&1")
        os.chdir(cwd)
        return des
    def batch_to_fasta(self, src, des):
        #Convert all fastq files in a folder to fasta
        print("Start converting fastq to fasta...")
        for f in os.scandir(src):
            if f.name.endswith(".fastq"):
                print("Converting {}".format(f.name))
                with open(f"{src}/{f.name}", 'r') as infile:
                    with open(f"{des}/{f.name[:-2]}",'w') as outfile:
                        for s in self._fastq_reader(infile):
                            outfile.write(">{}\n{}\n".format(s['title'],s['seq']))
        return des
    def distance_matrix(self, handle, TRUNCATE_HEAD_TAIL = True):
        raw = list(self._fasta_reader(handle))
        if TRUNCATE_HEAD_TAIL:
            for s in raw:
                try:
                    s['seq'] = self._truncate_head_and_tail(s['seq'])
                except Exception as e:
                    1
                    #print("Labeled HEAD not found in ", s['title'])
        print("Number of records:", len(raw))
        #Create distance matrix
        #fill with -1
        dm = np.full((len(raw), len(raw)), 0)
        for i in range(0, len(raw)):
            for j in range(i+1, len(raw)):
                d = min(self._pairwise_distance(raw[i]["seq"], raw[j]["seq"]), self._pairwise_distance(raw[i]["seq"], self._reverse_complement(raw[j]["seq"])))
                dm[i][j] = d
                dm[j][i] = dm[i][j]
        return dm
    def hdbscan(self, dm, min_cluster_size = 6, min_samples = 1):
        #HDBSCAN clustering
        clusterer = hdbscan.HDBSCAN(min_cluster_size = min_cluster_size, min_samples = min_samples)
        clusterer.fit(dm)
        return clusterer.labels_
    def clusters_fastas(self, src, des, min_cluster_size = 0.3, mds = True):
        try:
            os.makedirs(des, exist_ok=True)
        except Exception as e:
            print(e)
            pass
        for f in os.scandir(src):
            if f.name.endswith(".fas"):
                clustered_seq = {}
                print("Clustering {}".format(f.name))
                #Read, calculate distance matrix, cluster
                with open(f.path, 'r') as infile:
                    dm = self.distance_matrix(infile)
                    abs_cluster_size = max(2,int(dm.shape[0]*min_cluster_size))
                    print("abs_cluster_size: ", abs_cluster_size)
                    try:
                        labels = self.hdbscan(dm, abs_cluster_size) #Use relative_cluster_size
                        #print(f.name, labels)
                    except Exception as e:
                        print(e)
                        print("Clustering failed")
                        continue
                infile.close()
                #Organize sequences by cluster
                with open(f.path, 'r') as infile:
                    seqs = list(self._fasta_reader(infile))
                for i, l in enumerate(labels):
                    if l not in clustered_seq:
                        clustered_seq[l] = []
                    clustered_seq[l].append(seqs[i])
                infile.close()
                print ("Number of clusters:", len(clustered_seq))
                #Write to file
                for l in clustered_seq:
                    with open(f"{des}/{f.name[:-4]}_cluster_{l}_r{len(clustered_seq[l])}.fas", 'w') as outfile:
                        for s in clustered_seq[l]:
                            outfile.write(">{}\n{}\n".format(s['title'],s['seq']))
                #Visualize cluster result with mds
                if mds:
                    dm_norm = dm / dm.max()
                    mds = MDS(n_components=2,random_state=5566, dissimilarity='precomputed', normalized_stress="auto")
                    mds_results = mds.fit_transform(dm_norm)
                    fig, ax = plt.subplots(figsize=(15,15))
                    #ax.scatter(df['PC1'], df['PC2'], c=cluster_labels, cmap='rainbow', s=18)
                    ax.scatter(mds_results[:,0], mds_results[:,1], c=labels, cmap='rainbow', s=18)  
                    #Save the plot
                    plt.savefig(f"{des}/{f.name[:-4]}_MDS.jpg", dpi=56)
                    #Do not show the plot
                    plt.close()
    def lamassemble (self, src, des, mat = "/content/lamassemble/train/promethion.mat"):
        for f in os.scandir(src):
            if f.name.endswith(".fas"):
                self._exec(f'lamassemble /content/lamassemble/train/promethion.mat -a {f.path} > {des}/aln_{f.name}')
                self._exec(f'lamassemble /content/lamassemble/train/promethion.mat -c -n {f.name[:-4]} {des}/aln_{f.name} > {des}/con_{f.name}')
        return des              
    def deHead (self, src, des, start_offset = 0 , end_offset = 0):
        #Note: this function can only be applied to fasta files generated by singlebar which labels reads as follows:
        #Head regions are labeled by minibar as HEAD(uppercase)+barcode(lowercase)+SEQ WE NEED(uppercase)+barcode(lowercase)+TAIL(uppercase)
        #Start_offset, end_offset:  adjust the position of cut off point
        #                           which is useful when  we also want to remove the primer region
        try:
            os.makedirs(des, exist_ok=True)
        except Exception as e:
            print(e)
        for f in os.scandir(src):
            if f.name.endswith(".fas"):   
                with open (f.path, 'r') as infile:
                    with open (f"{des}/{f.name}", 'w') as outfile:
                        for s in self._fasta_reader(infile):
                            try:
                                #Head regions are labeled by minibar as HEAD(uppercase)+barcode(lowercase)+SEQ WE NEED(uppercase)+barcode(lowercase)+TAIL(uppercase) 
                                r = re.search("([A-Z]+)([a-z]+)([A-Z]+)([a-z]+)([A-Z]+)", s['seq'])
                                s['seq'] = s['seq'][r.start(3)+start_offset:r.end(3)-end_offset]
                                outfile.write(">{}\n{}\n".format(s['title'],s['seq']))
                            except Exception as e:
                                pass
                                #print("Labeled HEAD not found in ", s['title'])
        return des

    def trim_reads (self, src, des,  
                    BARCODE_INDEX_FILE,fw_col = "FwPrimer",rv_col = "RvPrimer",
                    fw_offset = 0, rv_offset = 0,
                    mismatch_ratio_f = 0.15,mismatch_ratio_r = 0.15,
                    discard_no_match = True,
                    check_both_directions = True,
                    reverse_complement_rv = True,
                    ):
        # This function removes primers from fasta files based on the provided barcode index file.
        # src: folder containing fasta files
        # des: folder to save the output fasta files
        # BARCODE_INDEX_FILE: tsv or csv file containing SampleID, FwPrimer, and RvPrimer columns
        # fw_col: column name for the forward primer sequence in the barcode index file
        # rv_col: column name for the reverse primer sequence in the barcode index file
        # mismatch_ratio_f: maximum allowed mismatch ratio for the forward primer sequence
        # mismatch_ratio_r: maximum allowed mismatch ratio for the reverse primer sequence
        # fw_offset: offset for the forward primer sequence
        # rv_offset: offset for the reverse primer sequence
        # discard_no_match: whether to discard sequences that do not match the primers
        # check_both_directions: whether to check both forward and reverse complement sequences, if primers are found in the reverse complement sequence, the output sequence will be reverse complemented
        # reverse_complement_rv: whether to reverse complement the reverse primer sequence
        try:
            os.makedirs(des, exist_ok=True)
        except Exception as e:
            print(e)
        try:
            if BARCODE_INDEX_FILE.endswith(".tsv"):
                df = pd.read_csv(BARCODE_INDEX_FILE, sep="\t")
            elif BARCODE_INDEX_FILE.endswith(".csv"):
                df = pd.read_csv(BARCODE_INDEX_FILE)
            else:
                raise Exception("BARCODE_INDEX_FILE should be a tsv or csv file")
            #Check if SampleID, FwPrimer, RvPrimer columns exist
            if not all (x in df.columns for x in ["SampleID", fw_col, rv_col]):
                raise Exception(f"BARCODE_INDEX_FILE should contain SampleID, {fw_col}, {rv_col} columns")
        except Exception as e:
            print(e)
        for f in os.scandir(src):
            if f.name.endswith(".fas"):   
                print("Processing", f.name)
                SampleID = f.name.replace(".fas","")
                try:
                    fw_trim = df.loc[df['SampleID'] == SampleID][fw_col].values[0]
                    rv_trim = df.loc[df['SampleID'] == SampleID][rv_col].values[0]
                    #Upper case
                    fw_trim = fw_trim.upper()
                    rv_trim = rv_trim.upper()
                    if reverse_complement_rv:
                        rv_trim = self._reverse_complement(rv_trim)
                except Exception as e:
                    print(e)
                    print("SampleID not found in BARCODE_INDEX_FILE")
                    continue
                with open (f.path, 'r') as infile:
                    with open (f"{des}/{f.name}", 'w') as outfile:
                        #Record trimmed, no-match, and total reads
                        trimmed_F = 0
                        trimmed_R = 0
                        total = 0
                        for s in self._fasta_reader(infile):
                            total += 1
                            seq_upper = s['seq'].upper()
                            try:
                                #use edlib to find fw_trim and rv_trim in the sequence
                                fw_loc = edlib.align(fw_trim, seq_upper, mode="HW", task="locations", k=int(len(fw_trim)*mismatch_ratio_f))
                                rv_loc = edlib.align(rv_trim, seq_upper, mode="HW", task="locations", k=int(len(rv_trim)*mismatch_ratio_r))
                                #Check if both fw_loc and rv_loc are found
                                if fw_loc['locations'] and rv_loc['locations'] and (rv_loc['locations'][0][0]-fw_loc['locations'][0][1]) > 0:
                                    #write trimmed sequence
                                    #print("F", fw_loc, rv_loc)
                                    outfile.write(">{}\n{}\n".format(s['title'],s['seq'][fw_loc['locations'][0][1]+fw_offset:rv_loc['locations'][0][0]-rv_offset]))
                                    trimmed_F += 1
                                elif check_both_directions:
                                    #Reverse complement the sequence and try again
                                    seq_upper = self._reverse_complement(seq_upper)
                                    fw_loc = edlib.align(fw_trim, seq_upper, mode="HW", task="locations", k=int(len(fw_trim)*mismatch_ratio_f))
                                    rv_loc = edlib.align(rv_trim, seq_upper, mode="HW", task="locations", k=int(len(rv_trim)*mismatch_ratio_r))
                                    if fw_loc['locations'] and rv_loc['locations'] and (rv_loc['locations'][0][0]-fw_loc['locations'][0][1]) > 0:
                                        #reverse complement the sequence
                                        reversed_seq = self._reverse_complement(s['seq'])
                                        #write trimmed sequence
                                        #print("R", fw_loc, rv_loc)
                                        outfile.write(">{}\n{}\n".format(s['title'],reversed_seq[fw_loc['locations'][0][1]+fw_offset:rv_loc['locations'][0][0]-rv_offset]))
                                        trimmed_R += 1
                                    elif discard_no_match:
                                        pass
                                    else:
                                        #write original sequence
                                        outfile.write(">Not_found_{}\n{}\n".format(s['title'],s['seq']))
                                else:
                                    if discard_no_match != True:
                                        #write original sequence
                                        outfile.write(">Not_found_{}\n{}\n".format(s['title'],s['seq']))
                            except Exception as e:
                                print(e)
                                continue
                print("Total reads:", total, "Trimmed(foward):", trimmed_F, "Trimmed(reverse):", trimmed_R, "No match:", total-trimmed_F-trimmed_R)








    """
    def medaka (self, src, des, startswith="con_"):
        for f in os.scandir(src):
            if f.name.startswith(startswith) and f.name.endswith(".fas"):
                aln_file = f"{assemble}/aln_{f.name[4:]}"
                shutil.rmtree("./medaka_seq/", ignore_errors=True)
                shutil.rmtree("./medaka/", ignore_errors=True)
                os.makedirs("./medaka_seq/", exist_ok=True)
                shutil.copyfile(f.path, './medaka_seq/con_temp.fas')
                shutil.copyfile(aln_file, './medaka_seq/aln_temp.fas')
                self._exec('medaka_consensus -d "./medaka_seq/con_temp.fas" -i "./medaka_seq/aln_temp.fas" -o "./medaka" 2>&1')
                try:
                    shutil.copyfile('./medaka/consensus.fasta', f'{des}/{f.name}')
                except:
                    pass
    """
    def blast(self, src, des, name="blast.csv", funguild = True, startswith="con_"):
        pool = []
        for f in os.scandir(src):
            row = {}
            if f.name.startswith(startswith) and f.name.endswith(".fas"):
                print("Blasting", f.name)
                with open(f.path, 'r') as handle:
                    raw = list(self._fasta_reader(handle))
                    for s in raw:
                        row['name'] = s['title']
                        row['seq'] = s['seq']
                        row['length'] = len(s['seq'])
                        try:
                            #con_sample_B11TUSc50_rDNA_cluster_1_r6
                            info = s['title'].split('cluster')[1]
                            row['cluster'] = info.split("_")[1]
                            row['reads'] = info.split("_")[2].replace("r","")
                        except:
                            pass
                try:
                    blast = self._blast(row['seq'])
                    row['organism'] = blast['org']
                    row['taxa'] = blast['taxa']
                    row['BLAST_simil'] = blast['sim']
                    row['BLAST_acc'] = blast['acc']
                    row['BLAST_seq'] = blast['seq']
                except:
                    pass

                #return {'acc':acc, 'org':org, 'taxa':taxo, 'sim':sim, 'seq':seq}
                #Check funguild
                if funguild:
                    funguild_des = []
                    try:
                        funguild_des = json.loads(get(f"https://www.mycoportal.org/funguild/services/api/db_return.php?qDB=funguild_db&qField=taxon&qText={row['organism']}").text)
                    except:
                        pass
                    if funguild_des != []:
                        #print(funguild_des)
                        try:
                            row['funguild'] = funguild_des[0]['guild']
                        except:
                            pass
                        try:
                            row['funguild_notes'] = funguild_des[0]['notes']
                        except:
                            pass
                print(row)
                pool.append(row)
        
        #pd.DataFrame(pool)[['name','cluster','reads', 'organism','taxa','seq', 'BLAST_simil','BLAST_acc','BLAST_seq', 'funguild', 'funguild_notes']].to_csv(f"{des}/blast.csv", index=False)
        pd.DataFrame(pool).to_csv(f"{des}/{name}", index=False)
        return f"{des}/{name}"
    def mmseqs_cluster(self, src, des, mmseqs="/nanoact/bin/mmseqs", min_seq_id=0.5, cov_mode=0, k=14, threads=8, s=7.5, cluster_mode=0, min_read_num = 0):
        #Get current library file  path
        lib = os.path.dirname(os.path.realpath(__file__))
        mmseqs = f"{lib}/bin/mmseqs"
        try:
            os.makedirs(des)
        except:
            pass
        abs_des = os.path.abspath(des)
        for f in os.scandir(src):
            if f.is_file() and f.name.endswith(".fas"):
                print("Clustering", f.name)
                sample = f.name.split(".fas")[0]
                #clean up temp folder
                self._clean_temp()
                #build db
                #print("Creating db")
                self._exec(f"{mmseqs} createdb {f.path} {self.TEMP}/db")
                #cluster
                #print("Clustering")
                self._exec(f"{mmseqs} cluster {self.TEMP}/db {self.TEMP}/cluster {self.TEMP}/tmp --min-seq-id {min_seq_id} --cov-mode {cov_mode} -k {k} --threads {threads} -s {s} --cluster-mode {cluster_mode}")
                #export tsv
                #self._exec(f"{mmseqs} createtsv {self.TEMP}/db {self.TEMP}/db {self.TEMP}/cluster {self.TEMP}/cluster.tsv")
                #export fasta
                #print("Parsing result")
                self._exec(f"{mmseqs} createseqfiledb {self.TEMP}/db {self.TEMP}/cluster {self.TEMP}/cluster.seq")
                self._exec(f"{mmseqs} result2flat {self.TEMP}/db {self.TEMP}/db {self.TEMP}/cluster.seq {self.TEMP}/cluster.fas")
                try:
                    with open(f"{self.TEMP}/cluster.fas", 'r') as handle:
                        bin = {}
                        cluster_no = -1
                        for rec in self._fasta_reader(handle):
                            if rec['seq'] == "":
                                cluster_no +=1
                                bin[cluster_no] = []
                                continue
                            else:
                                bin[cluster_no].append(rec)
                except Exception as e:
                    print("Error reading output file", e)
                    continue
                #save each cluster to file
                print(f"Number of clusters", len(bin))
                for cluster_no in bin:
                    if len(bin[cluster_no]) < min_read_num:
                        continue
                    with open(f"{abs_des}/{sample}_cluster_{cluster_no}_r{len(bin[cluster_no])}.fas", 'w') as handle:
                        for rec in bin[cluster_no]:
                            handle.write(f">{rec['title']}\n{rec['seq']}\n")
        return des

    def vsearch_OTUs(self, src, des, vsearch="/nanoact/bin/vsearch", id=0.9):
        #Get current library file  path
        lib = os.path.dirname(os.path.realpath(__file__))
        vsearch = f"{lib}/bin/vsearch"
        try:
            os.makedirs(des)
        except:
            pass
        abs_des = os.path.abspath(des)
        for f in os.scandir(src):
            if f.is_file() and f.name.endswith(".fas"):
                print("Clustering", f.name)   
                sample = f.name.split(".fas")[0]
                #clean up temp folder
                self._clean_temp()
                self._exec(f"{vsearch} --cluster_size {f.path} --id {id} --strand plus --sizein --sizeout --fasta_width 0 --uc {self.TEMP}/all.clustered.uc --relabel OTU_ --centroids {self.TEMP}/all.otus.fasta --otutabout {self.TEMP}/all.otutab.txt --clusters {self.TEMP}/cluster ",
                           suppress_output=True
                           )
                #read cluster uc file
                try:
                    uc = pd.read_csv(f"{self.TEMP}/all.clustered.uc", sep="\t", header=None)
                except:
                    print("Error reading output file", sample)
                    continue
                #Writing each cluster to file
                uc = uc[uc[0].isin(["S","H"])]
                uc.sort_values(by=8,ascending=False,inplace=True)
                seqs = list(self._fasta_reader(open(f.path,"r")))
                seqs = sorted(seqs,key=lambda d: d['title'])
                #export row 8 and row 1 as a list with {key:8 and value:1}
                seq_name_clust = uc[[8,1]].to_dict(orient="records")
                #separate each cluster to bin
                bin = {}
                for name_clust in seq_name_clust:
                    for seq in seqs:
                        if name_clust[8] in seq['title']:
                            seq_fas = f">{seq['title']}\n{seq['seq']}\n"
                            try:
                                bin[name_clust[1]].append(seq_fas)
                            except KeyError:
                                bin[name_clust[1]] = [seq_fas]
                            #remove seq from seqs to speed up next search
                            seqs.remove(seq)
                            break
                for cluster_no in bin:
                    with open(f"{abs_des}/{sample}_cluster_{cluster_no}_r{len(bin[cluster_no])}.fas", 'w') as handle:
                        for seq in bin[cluster_no]:
                            handle.write(seq)
                #Copy otu table to destination
                shutil.copy(f"{self.TEMP}/all.otutab.txt", f"{abs_des}/{sample}_otu_table.txt")

    def cd_hit_est(self, src, des, cd_hit_est="./nanoact/bin/cd-hit-est", id=0.9, n=10):
        #Get current library file  path
        lib = os.path.dirname(os.path.realpath(__file__))
        cd_hit_est = f"{lib}/bin/cd-hit-est"
        try:
            os.makedirs(des)
        except:
            pass
        abs_des = os.path.abspath(des)
        for f in os.scandir(src):
            if f.is_file() and f.name.endswith(".fas"):
                print("Clustering", f.name)   
                sample = f.name.split(".fas")[0]
                #clean up temp folder
                self._clean_temp()
                #-i input file
                #-o output file
                #-c sequence identity threshold
                #-n Suggested word size: 8,9,10 for thresholds 0.90 ~ 1.0 7 for thresholds 0.88 ~ 0.9 6 for thresholds 0.85 ~ 0.88 5 for thresholds 0.80 ~ 0.85 4 for thresholds 0.75 ~ 0.8
                #-d length of description in .clstr file, default 20
                self._exec(f"{cd_hit_est} -i {f.path} -o {self.TEMP}/cdhit.fas -c {id} -n {n} -d 0")
                #Try read clstr file
                try:
                    with open(f"{self.TEMP}/cdhit.fas.clstr", 'r') as handle:
                        bin = {}
                        for line in handle:
                            if line.startswith(">Cluster"):
                                cluster_no = line.split()[1]
                                bin[cluster_no] = []
                            else:
                                title = line.split(">")[1].split("...")[0]
                                bin[cluster_no].append(title)
                except:
                    print("Error reading output file", sample)
                    continue

                #Reading original sequence file
                with open(f.path, 'r') as handle:
                    seqs = list(self._fasta_reader(handle))
                #Write each cluster to a file
                for cluster in bin:
                    with open(f"{abs_des}/{sample}_cluster_{cluster}_r{len(bin[cluster])}.fas", 'w') as handle:
                        for seq_title in bin[cluster]:
                            for read in seqs:
                                if seq_title in read['title']:
                                    handle.write(f">{read['title']}\n{read['seq']}\n")
                                    #remove used read from the list, so that it won't be used again
                                    seqs.remove(read)
                                    break





    def _calculate_5mer_frequency(self, sequence):
        frequency = {}
        for i in range(len(sequence) - 4):
            kmer = sequence[i:i+5]
            frequency[kmer] = frequency.get(kmer, 0) + 1
        return frequency
    
    def fas_to_5mer(self, fas_path):
        frequency_vectors = []
        sequences = []
        for rec in self._fasta_reader(open(fas_path,"r")):
            frequency_vectors.append(self._calculate_5mer_frequency(rec['seq']))
            sequences.append(rec['title'])
        # Create a DataFrame to store the frequency vectors
        df = pd.DataFrame(frequency_vectors)
        # Add a column for the sequence identifiers
        df['Sequence'] = ['Sequence {}'.format(i+1) for i in range(len(sequences))]
        # Reorder the columns to have 'Sequence' as the first column
        df = df[['Sequence'] + list(df.columns[:-1])]
        #fill nan with 0
        df.fillna(0, inplace=True)
        return df

