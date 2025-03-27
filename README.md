<h2 id="overview">Overview</h2>
This document describes the analysis pipeline for building and utilizing the **InsectTFDB** database, designed to catalog transcription factors (TFs) across a wide range of insect species. The pipeline consists of data collection, transcription factor identification, functional annotation, clustering analysis, and database construction.



<h2 id="requirements">Requirements</h2>
<h3 id="software">Software</h3>
1. gffread (v0.12.7)
2. HMMER (v3.4)
3. DIAMOND (v2.1.10)
4. eggNOG-mapper (v2.1.9)
5. MMseqs2 (version 13.45111)

<h3 id="data-sources">Data Sources</h3>
InsectBase 2.0 and NCBI for genome assemblies and corresponding official gene sets.
AnimalTFDB 4.0 and InterPro for Hidden Markov Model (HMM) profiles of transcription factor families.

<h2 id="analysis-pipeline">Analysis Pipeline</h2>
<h3 id="64dfde5a">Step 1: Data Collection</h3>
<h4 id="wZImj">1.1 Genome Assembly Retrieval</h4>
Genomes and official gene sets are retrieved from **InsectBase 2.0** and **NCBI**, filtered by BUSCO completeness scores greater than 80%. 

<h4 id="qWJma">1.2 GFF collection</h4>
 For genomes not derived from InsectBase 2.0, annotations were uniformly completed following the pipeline at [https://github.com/meiyang12/Genome-annotation-pipeline](https://github.com/meiyang12/Genome-annotation-pipeline).

<h3 id="ofkVn">Step 2: Transcription Factor Identification</h3>
<h4 id="jZZcW">2.1  Extract Protein Sequences  </h4>
Use `gffread` to extract CDS from the GFF file and translate them into protein sequences.

```bash
gffread genome.gff -g genome.fa -y proteins.fa
```

Where:

+ `genome.gff` is the annotation file;
+ `genome.fa` is the corresponding genome sequence;
+ `proteins.fa` is the output protein FASTA file.

---

<h4 id="jQ52m">2.2 Use HMMER to Identify Transcription Factors</h4>
Download the TF HMM models (from AnimalTFDB or InterPro), and use `hmmscan` for scanning.

```bash
hmmscan --domtblout tf.domtblout --cpu 40 TF.hmm proteins.fa > tf.hmmscan.txt
```

Use the following script to filter matches with E-value ≤ 1e-4:

```python
import pandas as pd

with open("tf.domtblout") as f:
    lines = [l for l in f if not l.startswith("#")]

rows = []
for line in lines:
    parts = line.strip().split()
    evalue = float(parts[6])
    if evalue <= 1e-4:
        rows.append((parts[0], parts[3], evalue))

df = pd.DataFrame(rows, columns=["Query", "Domain", "E-value"])
df.to_csv("filtered_tf_hits.csv", index=False)
```

---

<h3 id="wvEyM">Step 3: Functional Annotation</h3>
<h4 id="hZnGc">3.1 Sequence Alignment with UniProt (DIAMOND)</h4>
```bash
diamond blastp -d uniprot.dmnd -q proteins.fa \
  -o uniprot_matches.tsv \
  --sensitive -k 1 --query-cover 50 -e 1e-5 -p 20
```

+ `uniprot.dmnd` is the pre-built DIAMOND database for UniProt.

---

<h4 id="SJkYi">3.2 Functional Annotation using eggNOG-mapper (GO / KEGG) (<font style="color:rgb(123, 126, 138);">eggNOG 6 database)</font></h4>
```bash
emapper.py --cpu 20 --mp_start_method forkserver \
--data_dir /dev/shm/ -o $query_name \
--output_dir /emapper_out \
--temp_dir /emapper_tmp \
--override -m diamond \
--dmnd_ignore_warnings \
-i /queries.fasta \
--evalue 0.001 --score 60 --pident 40 --query_cover 20 --subject_cover 20 \
--itype proteins --tax_scope auto --target_orthologs all \
--go_evidence non-electronic --pfam_realign none --report_orthologs \
--decorate_gff yes \
--excel > /emapper.out 2> /emapper.err
```

The output includes `.annotations` files, containing functional annotations, GO terms, and KEGG pathway information. <font style="color:rgb(123, 126, 138);"></font>

---

<h3 id="LUs5g">Step 4: Clustering Analysis</h3>
Use `MMseqs2` for clustering the identified TF proteins to analyze evolutionary relationships.

```bash
mmseqs createdb proteins.fa tf_db

mmseqs cluster tf_db tf_cluster tmp_dir --min-seq-id 0.8 -c 0.8 --cov-mode 1 --threads 8

mmseqs createtsv tf_db tf_db tf_cluster tf_cluster.tsv
```





<h4 id="ee4ih">Table S2. DNA binding domains of transcription factors </h4>
| <font style="color:black;">Family</font> | <font style="color:black;">Domain</font> | <font style="color:black;">PFAM</font> | <font style="color:black;">Group</font> | <font style="color:black;">Cutoff (</font>_<font style="color:black;">E</font>_<font style="color:black;">-value)</font> |
| --- | --- | --- | --- | --- |
| <font style="color:black;">AP-2</font> | <font style="color:black;">TF_AP-2 domain</font> | [<u><font style="color:#0563c1;">PF03299</font></u>](http://pfam.xfam.org/family/PF03299) | <font style="color:black;">Basic Domians group</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">bHLH</font> | <font style="color:black;">HLH domain</font> | [<u><font style="color:#0563c1;">PF00010</font></u>](http://pfam.xfam.org/family/PF00010) | <font style="color:black;">Basic Domians group</font> | <font style="color:black;">1.0E-02</font> |
| <font style="color:black;">Nrf1</font> | <font style="color:black;">Nrf1_DNA-bind domain</font> | [<u><font style="color:#0563c1;">PF10491</font></u>](http://pfam.xfam.org/family/PF10491) | <font style="color:black;">Basic Domians group</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">RFX</font> | <font style="color:black;">RFX domain</font> | [<u><font style="color:#0563c1;">PF02257</font></u>](http://pfam.xfam.org/family/PF02257) | <font style="color:black;">Basic Domians group</font> | <font style="color:black;">1.0E-10</font> |
| <font style="color:black;">TF_bZIP</font> | <font style="color:black;">bZIP domain</font> | <font style="color:black;">AnimalTFDB</font> | <font style="color:black;">Basic Domians group</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">TSC22</font> | <font style="color:black;">TSC22 domain</font> | [<u><font style="color:#0563c1;">PF01166</font></u>](http://pfam.xfam.org/family/PF01166) | <font style="color:black;">Basic Domians group</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">CBF</font> | <font style="color:black;">CBF_beta domain</font> | [<u><font style="color:#0563c1;">PF02312</font></u>](http://pfam.xfam.org/family/PF02312) | <font style="color:black;">Beta-Scaffold Factors</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">CP2</font> | <font style="color:black;">CP2 domain</font> | [<u><font style="color:#0563c1;">PF04516</font></u>](http://pfam.xfam.org/family/PF04516) | <font style="color:black;">Beta-Scaffold Factors</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">CSD</font> | <font style="color:black;">CSD domain</font> | [<u><font style="color:#0563c1;">PF00313</font></u>](http://pfam.xfam.org/family/PF00313) | <font style="color:black;">Beta-Scaffold Factors</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">CSL</font> | <font style="color:black;">BTD domain</font> | [<u><font style="color:#0563c1;">PF09270</font></u>](http://pfam.xfam.org/family/PF09270) | <font style="color:black;">Beta-Scaffold Factors</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">GCM</font> | <font style="color:black;">GCM domain</font> | [<u><font style="color:#0563c1;">PF03615</font></u>](http://pfam.xfam.org/family/PF03615) | <font style="color:black;">Beta-Scaffold Factors</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">P53</font> | <font style="color:black;">P53 domain</font> | [<u><font style="color:#0563c1;">PF00870</font></u>](http://pfam.xfam.org/family/PF00870) | <font style="color:black;">Beta-Scaffold Factors</font> | <font style="color:black;">1.0E-03</font> |
| <font style="color:black;">RHD</font> | <font style="color:black;">RHD domain</font> | [<u><font style="color:#0563c1;">PF00554</font></u>](http://pfam.xfam.org/family/PF00554) | <font style="color:black;">Beta-Scaffold Factors</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">Runt</font> | <font style="color:black;">Runt domain</font> | [<u><font style="color:#0563c1;">PF00853</font></u>](http://pfam.xfam.org/family/PF00853) | <font style="color:black;">Beta-Scaffold Factors</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">STAT</font> | <font style="color:black;">STAT_bind domain</font> | [<u><font style="color:#0563c1;">PF02864</font></u>](http://pfam.xfam.org/family/PF02864) | <font style="color:black;">Beta-Scaffold Factors</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">ARID</font> | <font style="color:black;">ARID domain</font> | [<u><font style="color:#0563c1;">PF01388</font></u>](http://pfam.xfam.org/family/PF01388) | <font style="color:black;">Helix-turn-helix</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">COE</font> | <font style="color:black;">COE domain</font> | <font style="color:black;">AnimalTFDB</font> | <font style="color:black;">Helix-turn-helix</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">CUT</font> | <font style="color:black;">Homeobox|CUT</font> | [<u><font style="color:#0563c1;">PF02376</font></u>](http://pfam.xfam.org/family/PF02376) | <font style="color:black;">Helix-turn-helix</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">E2F</font> | <font style="color:black;">E2F_TDP domain</font> | [<u><font style="color:#0563c1;">PF02319</font></u>](http://pfam.xfam.org/family/PF02319) | <font style="color:black;">Helix-turn-helix</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">ETS</font> | <font style="color:black;">Ets domain</font> | [<u><font style="color:#0563c1;">PF00178</font></u>](http://pfam.xfam.org/family/PF00178) | <font style="color:black;">Helix-turn-helix</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">Fork_head</font> | <font style="color:black;">Fork_head domain</font> | [<u><font style="color:#0563c1;">PF00250</font></u>](http://pfam.xfam.org/family/PF00250) | <font style="color:black;">Helix-turn-helix</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">Homeobox</font> | <font style="color:black;">Homeobox</font> | [<u><font style="color:#0563c1;">PF00046</font></u>](http://pfam.xfam.org/family/PF00046) | <font style="color:black;">Helix-turn-helix</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">HPD</font> | <font style="color:black;">HPD domain</font> | [<u><font style="color:#0563c1;">PF05044</font></u>](http://pfam.xfam.org/family/PF05044) | <font style="color:black;">Helix-turn-helix</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">HSF</font> | <font style="color:black;">HSF_DNA-bind domain</font> | [<u><font style="color:#0563c1;">PF00447</font></u>](http://pfam.xfam.org/family/PF00447) | <font style="color:black;">Helix-turn-helix</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">HTH</font> | <font style="color:black;">HTH_psq domain</font> | [<u><font style="color:#0563c1;">PF05225</font></u>](http://pfam.xfam.org/family/PF05225) | <font style="color:black;">Helix-turn-helix</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">IRF</font> | <font style="color:black;">IRF domain</font> | [<u><font style="color:#0563c1;">PF00605</font></u>](http://pfam.xfam.org/family/PF00605) | <font style="color:black;">Helix-turn-helix</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">MYB</font> | <font style="color:black;">Myb_DNA-binding domain</font> | [<u><font style="color:#0563c1;">PF00249</font></u>](http://pfam.xfam.org/family/PF00249) | <font style="color:black;">Helix-turn-helix</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">PAX</font> | <font style="color:black;">PAX domain</font> | [<u><font style="color:#0563c1;">PF00292</font></u>](http://pfam.xfam.org/family/PF00292) | <font style="color:black;">Helix-turn-helix</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">Pou</font> | <font style="color:black;">Homeobox|Pou</font> | [<u><font style="color:#0563c1;">PF00157</font></u>](http://pfam.xfam.org/family/PF00157) | <font style="color:black;">Helix-turn-helix</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">SRF</font> | <font style="color:black;">SRF domain</font> | [<u><font style="color:#0563c1;">PF00319</font></u>](http://pfam.xfam.org/family/PF00319) | <font style="color:black;">Helix-turn-helix</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">TEA</font> | <font style="color:black;">TEA domain</font> | [<u><font style="color:#0563c1;">PF01285</font></u>](http://pfam.xfam.org/family/PF01285) | <font style="color:black;">Helix-turn-helix</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">GTF2I</font> | <font style="color:black;">GTF2I domain</font> | [<u><font style="color:#0563c1;">PF02946</font></u>](http://pfam.xfam.org/family/PF02946) | <font style="color:black;">Other Alpha-Helix Group</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">HMG</font> | <font style="color:black;">HMG_box domain</font> | [<u><font style="color:#0563c1;">PF00505</font></u>](http://pfam.xfam.org/family/PF00505) | <font style="color:black;">Other Alpha-Helix Group</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">NF-YA</font> | <font style="color:black;">CBFB_NFYA domain</font> | [<u><font style="color:#0563c1;">PF02045</font></u>](http://pfam.xfam.org/family/PF02045) | <font style="color:black;">Other Alpha-Helix Group</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">NF-YB</font> | <font style="color:black;">NF-YB domain</font> | <font style="color:black;">AnimalTFDB</font> | <font style="color:black;">Other Alpha-Helix Group</font> | <font style="color:black;">1.0E-18</font> |
| <font style="color:black;">NF-YC</font> | <font style="color:black;">NF-YC domain</font> | <font style="color:black;">AnimalTFDB</font> | <font style="color:black;">Other Alpha-Helix Group</font> | <font style="color:black;">1.0E-18</font> |
| <font style="color:black;">SAND</font> | <font style="color:black;">SAND domain</font> | [<u><font style="color:#0563c1;">PF01342</font></u>](http://pfam.xfam.org/family/PF01342) | <font style="color:black;">Other Alpha-Helix Group</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">AF-4</font> | <font style="color:black;">AF-4 domain</font> | [<u><font style="color:#0563c1;">PF05110</font></u>](http://pfam.xfam.org/family/PF05110) | <font style="color:black;">Unclassified Structure</font> | <font style="color:black;">1.0E-03</font> |
| <font style="color:black;">CG-1</font> | <font style="color:black;">CG-1 domain</font> | [<u><font style="color:#0563c1;">PF03859</font></u>](http://pfam.xfam.org/family/PF03859) | <font style="color:black;">Unclassified Structure</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">CSRNP_N</font> | <font style="color:black;">CSRNP_N domain</font> | [<u><font style="color:#0563c1;">PF16019</font></u>](http://pfam.xfam.org/family/PF16019) | <font style="color:black;">Unclassified Structure</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">CTF_NFI</font> | <font style="color:black;">CTF/NFI and MH1 domain</font> | [<u><font style="color:#0563c1;">PF00859</font></u>](http://pfam.xfam.org/family/PF00859) | <font style="color:black;">Unclassified Structure</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">DACH</font> | <font style="color:black;">DACH domain</font> | <font style="color:black;">AnimalTFDB</font> | <font style="color:black;">Unclassified Structure</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">GCFC</font> | <font style="color:black;">GCFC domain</font> | [<u><font style="color:#0563c1;">PF07842</font></u>](http://pfam.xfam.org/family/PF07842) | <font style="color:black;">Unclassified Structure</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">HMGA</font> | <font style="color:black;">HMGA domain</font> | <font style="color:black;">AnimalTFDB</font> | <font style="color:black;">Unclassified Structure</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">LRRFIP</font> | <font style="color:black;">LRRFIP domain</font> | [<u><font style="color:#0563c1;">PF09738</font></u>](http://pfam.xfam.org/family/PF09738) | <font style="color:black;">Unclassified Structure</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">MBD</font> | <font style="color:black;">MBD domain</font> | [<u><font style="color:#0563c1;">PF01429</font></u>](http://pfam.xfam.org/family/PF01429) | <font style="color:black;">Unclassified Structure</font> | <font style="color:black;">1.0E-05</font> |
| <font style="color:black;">MH1</font> | <font style="color:black;">MH1 domain</font> | [<u><font style="color:#0563c1;">PF03165</font></u>](http://pfam.xfam.org/family/PF03165) | <font style="color:black;">Unclassified Structure</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">NCU-G1</font> | <font style="color:black;">NCU-G1 domain</font> | [<u><font style="color:#0563c1;">PF15065</font></u>](http://pfam.xfam.org/family/PF15065) | <font style="color:black;">Unclassified Structure</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">NDT80_PhoG</font> | <font style="color:black;">NDT80_PhoG domain</font> | [<u><font style="color:#0563c1;">PF05224</font></u>](http://pfam.xfam.org/family/PF05224) | <font style="color:black;">Unclassified Structure</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">PC4</font> | <font style="color:black;">PC4 domain</font> | [<u><font style="color:#0563c1;">PF02229</font></u>](http://pfam.xfam.org/family/PF02229) | <font style="color:black;">Unclassified Structure</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">T-box</font> | <font style="color:black;">T-box domain</font> | [<u><font style="color:#0563c1;">PF00907</font></u>](http://pfam.xfam.org/family/PF00907) | <font style="color:black;">Unclassified Structure</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">Tub</font> | <font style="color:black;">Tub domain</font> | [<u><font style="color:#0563c1;">PF01167</font></u>](http://pfam.xfam.org/family/PF01167) | <font style="color:black;">Unclassified Structure</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">DM</font> | <font style="color:black;">DM domain</font> | [<u><font style="color:#0563c1;">PF00751</font></u>](http://pfam.xfam.org/family/PF00751) | <font style="color:black;">Zinc-Coordinating Group</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">ESR-like</font> | <font style="color:black;">zf-C4|ESR-like</font> | <font style="color:black;">AnimalTFDB</font> | <font style="color:black;">Zinc-Coordinating Group</font> | <font style="color:black;">1.0E-34</font> |
| <font style="color:black;">NGFIB-like</font> | <font style="color:black;">zf-C4|NGFIB-like</font> | <font style="color:black;">AnimalTFDB</font> | <font style="color:black;">Zinc-Coordinating Group</font> | <font style="color:black;">1.0E-24</font> |
| <font style="color:black;">RXR-like</font> | <font style="color:black;">zf-C4|RXR-like</font> | <font style="color:black;">AnimalTFDB</font> | <font style="color:black;">Zinc-Coordinating Group</font> | <font style="color:black;">1.0E-30</font> |
| <font style="color:black;">SF-like</font> | <font style="color:black;">zf-C4|SF-like</font> | <font style="color:black;">AnimalTFDB</font> | <font style="color:black;">Zinc-Coordinating Group</font> | <font style="color:black;">1.0E-40</font> |
| <font style="color:black;">THAP</font> | <font style="color:black;">THAP domain</font> | [<u><font style="color:#0563c1;">PF05485</font></u>](http://pfam.xfam.org/family/PF05485) | <font style="color:black;">Zinc-Coordinating Group</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">THR-like</font> | <font style="color:black;">zf-C4|THR-like</font> | <font style="color:black;">AnimalTFDB</font> | <font style="color:black;">Zinc-Coordinating Group</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">ZBTB</font> | <font style="color:black;">zf-C2H2|ZBTB</font> | [<u><font style="color:#0563c1;">PF00651</font></u>](http://pfam.xfam.org/family/PF00651) | <font style="color:black;">Zinc-Coordinating Group</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">zf-BED</font> | <font style="color:black;">zf-BED domain</font> | [<u><font style="color:#0563c1;">PF02892</font></u>](http://pfam.xfam.org/family/PF02892) | <font style="color:black;">Zinc-Coordinating Group</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">zf-C2H2</font> | <font style="color:black;">zf-C2H2 domain</font> | [<u><font style="color:#0563c1;">PF00096</font></u>](http://pfam.xfam.org/family/PF00096) | <font style="color:black;">Zinc-Coordinating Group</font> | <font style="color:black;">1.0E-05</font> |
| <font style="color:black;">zf-C2HC</font> | <font style="color:black;">zf-C2HC domain</font> | [<u><font style="color:#0563c1;">PF01530</font></u>](http://pfam.xfam.org/family/PF01530) | <font style="color:black;">Zinc-Coordinating Group</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">zf-CCCH</font> | <font style="color:black;">zf-CCCH domain</font> | [<u><font style="color:#0563c1;">PF00642</font></u>](http://pfam.xfam.org/family/PF00642) | <font style="color:black;">Zinc-Coordinating Group</font> | <font style="color:black;">1.0E-20</font> |
| <font style="color:black;">zf-GAGA</font> | <font style="color:black;">zf-GAGA domain</font> | [<u><font style="color:#0563c1;">PF09237</font></u>](http://pfam.xfam.org/family/PF09237) | <font style="color:black;">Zinc-Coordinating Group</font> | <font style="color:black;">1.0E-05</font> |
| <font style="color:black;">zf-GATA</font> | <font style="color:black;">zf-GATA domain</font> | [<u><font style="color:#0563c1;">PF00320</font></u>](http://pfam.xfam.org/family/PF00320) | <font style="color:black;">Zinc-Coordinating Group</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">zf-LITAF-like</font> | <font style="color:black;">zf-LITAF-like domain</font> | [<u><font style="color:#0563c1;">PF10601</font></u>](http://pfam.xfam.org/family/PF10601) | <font style="color:black;">Zinc-Coordinating Group</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">zf-MIZ</font> | <font style="color:black;">zf-MIZ domain</font> | [<u><font style="color:#0563c1;">PF02891</font></u>](http://pfam.xfam.org/family/PF02891) | <font style="color:black;">Zinc-Coordinating Group</font> | <font style="color:black;">1.0E-04</font> |
| <font style="color:black;">zf-NF-X1</font> | <font style="color:black;">zf-NF-X1 domain</font> | [<u><font style="color:#0563c1;">PF01422</font></u>](http://pfam.xfam.org/family/PF01422) | <font style="color:black;">Zinc-Coordinating Group</font> | <font style="color:black;">1.0E-04</font> |


