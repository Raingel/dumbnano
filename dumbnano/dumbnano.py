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
class NanoAmpliParser():
    def __init__(self, TEMP = './temp/'):
        self.TEMP = TEMP
    def _lib_path(self):
        #get the path of the library
        return os.path.dirname(os.path.realpath(__file__))
    def _exec(self, cmd,suppress_output=False):
        if suppress_output:
            with open(os.devnull, 'w') as DEVNULL:
                out = run(cmd, stdout=DEVNULL, stderr=DEVNULL, shell=True)
            return None
        else:
            out = run(cmd, stdout=PIPE, shell=True)
            #print (">>", cmd)
            #print("Output:")
            #print(out.stdout.decode('utf-8'))
            #print("Exception:")
            #print(out.stderr)
            return out.stdout.decode('utf-8'), out.stderr
    
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
        #Reverse complement a sequence
        s = s.replace("A","T")
        s = s.replace("T","A")
        s = s.replace("C","G")
        s = s.replace("G","C")
        #for lower case
        s = s.replace("a","t")
        s = s.replace("t","a")
        s = s.replace("c","g")
        s = s.replace("g","c")
        return s[::-1]
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
        response = get(url)
        #Convert the XML to a dictionary
        blast_dict = xmltodict.parse(response.text)
        #Get the first hit of each query
        pool = {}
        for rec in blast_dict['BlastOutput']['BlastOutput_iterations']['Iteration']:
            seq_name = rec['Iteration_query-def']
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
                similarity = hit_hsp['Hsp_identity']
                #Get taxon info
                r = get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&retmode=xml&id={}'.format(acc))
                #xml to json
                r = xmltodict.parse(r.text)
                org = r['GBSet']['GBSeq']['GBSeq_organism']
                taxa = r['GBSet']['GBSeq']['GBSeq_taxonomy']
                pool[seq_name] = {'acc': acc, 'hit_seq': hit_seq, 'hit_def': hit_def, 'org': org, 'taxa': taxa}
        return pool
    def _get_sample_id (self, seq, barcode_hash_table, mismatch_ratio_f = 0.15,mismatch_ratio_r = 0.15):
        # Define a helper function to identify the sample ID of a sequence read based on its barcode
        ids = []
        integrity = []
        seq = seq.upper()
        seqF = seq[:100]
        seqR = seq[-100:]
        # Convert the input sequence read to uppercase and extract the beginning and end segments
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
        #Paste back modified seq
        seq = seqF +seq[100:-100] + seqR
        # If a barcode is identified, return the corresponding sample ID and a boolean indicating whether the barcode was matched with sufficient integrity
        return ids, integrity, seq
    def _mafft (self, src, des):
        mafft_bin = self._lib_path() + "/bin/mafft.bat"
        #./mafft.bat --genafpair --maxiterate 1000 2110_cluster_1_r442.fas > output.fas
        cmd = f"{mafft_bin} {src} > {des}"
        self._exec(cmd,suppress_output=True)
    def orientation(self, src, des, tsv):
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
                            aln_f = edlib.align(F.upper(), record['seq'].upper()[:100], mode="HW", task="locations")
                            aln_r = edlib.align(R.upper(), record['seq'].upper()[:100], mode="HW", task="locations")
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
    def singlebar(self, src, des, BARCODE_INDEX_FILE, mismatch_ratio_f = 0.15,mismatch_ratio_r = 0.15):
        """
        1.Process the raw sequencing file using a fastq_reader function which reads in four lines at a time representing one sequencing read.
        2.For each read, call the _get_sample_id function to identify its corresponding sample ID and barcode integrity.
        3.If a sample ID is identified, append the read to the corresponding sample's output file.
        4.If a barcode is truncated, append it to the "TRUNCATED" dictionary.
        5.If multiple barcodes are identified, append it to the "MULTIPLE" dictionary.
        6.If no barcode is identified, append it to the "UNKNOWN" dictionary.
        7.Finally, output the demultiplexed reads into different output files based on their barcode.
        """
        # Define the main function for demultiplexing
        BARCODE_IDX_DF = pd.read_csv(BARCODE_INDEX_FILE, sep="\\t")
        # Read in the barcode index file as a pandas DataFrame
        if not all([x in BARCODE_IDX_DF.columns for x in ["SampleID", "FwIndex", "RvAnchor"]]):
            raise ValueError("BARCODE_INDEX_FILE must have SampleID, FwIndex, RvAnchor columns")
        # Check whether the barcode index file has the required columns
        barcode_hash_table = {}
        for index, row in BARCODE_IDX_DF.iterrows():
            barcode_hash_table[row["SampleID"]] = {"FwIndex": row["FwIndex"], "RvAnchor": row["RvAnchor"]}
        # Store the barcode index file as a hash table of barcode-sample ID pairs
        pool = {}
        counter = 1
        pool["UNKNOWN"] = []
        pool["MULTIPLE"] = []
        pool["TRUNCATED"] = []
        for id in barcode_hash_table:
            pool[id] = []
        # Initialize dictionaries for output files for each sample and additional dictionaries for reads with unknown, multiple, or truncated barcodes
        with open(src, "r") as handle:
            for record in self._fastq_reader(handle):
                ids, integrity, seq= self._get_sample_id(record["seq"], barcode_hash_table, mismatch_ratio_f, mismatch_ratio_r)
                record["seq"] = seq
                if len(ids) == 1:
                    if integrity[0] == False:
                        pool["TRUNCATED"].append(record)
                    else:
                        pool[ids[0]].append(record)
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
            if bin in ['MULTIPLE', 'UNKNOWN', 'TRUNCATED']:
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
        for f in os.scandir(src):
            if f.name.endswith(".fas"):   
                try:
                    os.makedirs(des, exist_ok=True)
                except Exception as e:
                    print(e)
                    continue
                with open (f.path, 'r') as infile:
                    with open (f"{des}/{f.name}", 'w') as outfile:
                        for s in self._fasta_reader(infile):
                            try:
                                #Head regions are labeled by minibar as HEAD(uppercase)+barcode(lowercase)+SEQ WE NEED(uppercase)+barcode(lowercase)+TAIL(uppercase) 
                                r = re.search("([A-Z]+)([a-z]+)([A-Z]+)([a-z]+)([A-Z]+)", s['seq'])
                                s['seq'] = s['seq'][r.start(3)+start_offset:r.end(3)-end_offset]
                                outfile.write(">{}\n{}\n".format(s['title'],s['seq']))
                            except Exception as e:
                                1
                                #print("Labeled HEAD not found in ", s['title'])
        return des
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
                        if len(s['seq']) <= max_query_length:
                            #If sequence is too long, preserve only max_query_length in middle
                            s['seq'] = s['seq'][:int(max_query_length//2)] + s['seq'][-int(max_query_length//2):]
                        query_seqs.append(f">{s['title']}\n{s['seq']}")       
        #set title as index
        pool_df.set_index('title', inplace=True)
        for index, row in pool_df.iterrows():
            pool_df.loc[index, 'length'] = str(len(row['seq']))
            try:
                #2110_cluster_-1_r2154.fas	
                #{sample}_cluster_{cluster_no}_r{reads_count}.fas
                sample, cluster_no, reads_count = re.search("(.*)_cluster_([-0-9]+)_r(\d+).fas", row['title']).groups()
                pool_df.loc[index, 'sample'] = sample
                pool_df.loc[index, 'cluster_no'] = cluster_no
                pool_df.loc[index, 'reads_count'] = reads_count
            except:
                pass    
        #Blast all sequences
        i = 0
        blast_result_pool = {}
        while i < len(query_seqs):
            print("Blasting", i, "to", i+batch, "of", len(query_seqs))
            query = "\n".join(query_seqs[i:i+batch])
            retry = 3
            while True:
                try:
                    blast_result = self.NCBIblast(query)
                    if blast_result != None:
                        blast_result_pool.update(blast_result)
                        break
                    else:
                        retry -= 1
                except:
                    retry -= 1
            i+=batch
        for sample in blast_result_pool.keys():
            pool_df.loc[sample, 'acc'] = blast_result_pool[sample]['acc']
            pool_df.loc[sample, 'hit_seq'] = blast_result_pool[sample]['hit_seq']
            pool_df.loc[sample, 'hit_def'] = blast_result_pool[sample]['hit_def']
            pool_df.loc[sample, 'org'] = blast_result_pool[sample]['org']
            pool_df.loc[sample, 'taxa'] = blast_result_pool[sample]['taxa']
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
        return pool_df
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




