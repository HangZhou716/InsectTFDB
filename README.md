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
Genomes and official gene sets are retrieved from "InsectBase 2.0" and "NCBI", filtered by BUSCO completeness scores greater than 80%. 

<h4 id="qWJma">1.2 GFF collection</h4>
 For genomes not derived from InsectBase 2.0, annotations were uniformly completed following the pipeline at https://github.com/meiyang12/Genome-annotation-pipeline.

<h3 id="ofkVn">Step 2: Transcription Factor Identification</h3>
<h4 id="jZZcW">2.1  Extract Protein Sequences  </h4>
Use "gffread" to extract CDS from the GFF file and translate them into protein sequences.

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
  --sensitive -k 1 --id 30 --query-cover 50 -e 1e-5 -p 20
```

"uniprot.dmnd" is the pre-built DIAMOND database for UniProt.

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

The output includes ".annotations" files, containing functional annotations, GO terms, and KEGG pathway information. <font style="color:rgb(123, 126, 138);"></font>

---

<h3 id="LUs5g">Step 4: Clustering Analysis</h3>
Use "MMseqs2" for clustering the identified TF proteins to analyze evolutionary relationships.

```bash
mmseqs createdb proteins.fa tf_db

mmseqs cluster tf_db tf_cluster tmp_dir --min-seq-id 0.8 -c 0.8 --cov-mode 1 --threads 8

mmseqs createtsv tf_db tf_db tf_cluster tf_cluster.tsv
```





<h4 id="ee4ih">Table S2. DNA binding domains of transcription factors </h4>

| Family     | Domain                        | PFAM      | Group                    | Cutoff (E-value) |
|------------|-------------------------------|-----------|--------------------------|------------------|
| AP-2       | TF_AP-2 domain                 | PF03299   | Basic Domians group       | 1.0E-04          |
| bHLH       | HLH domain                     | PF00010   | Basic Domians group       | 1.0E-02          |
| Nrf1       | Nrf1_DNA-bind domain           | PF10491   | Basic Domians group       | 1.0E-04          |
| RFX        | RFX domain                     | PF02257   | Basic Domians group       | 1.0E-10          |
| TF_bZIP    | bZIP domain                    | AnimalTFDB| Basic Domians group       | 1.0E-04          |
| TSC22      | TSC22 domain                   | PF01166   | Basic Domians group       | 1.0E-04          |
| CBF        | CBF_beta domain                | PF02312   | Beta-Scaffold Factors     | 1.0E-04          |
| CP2        | CP2 domain                     | PF04516   | Beta-Scaffold Factors     | 1.0E-04          |
| CSD        | CSD domain                     | PF00313   | Beta-Scaffold Factors     | 1.0E-04          |
| CSL        | BTD domain                     | PF09270   | Beta-Scaffold Factors     | 1.0E-04          |
| GCM        | GCM domain                     | PF03615   | Beta-Scaffold Factors     | 1.0E-04          |
| P53        | P53 domain                     | PF00870   | Beta-Scaffold Factors     | 1.0E-03          |
| RHD        | RHD domain                     | PF00554   | Beta-Scaffold Factors     | 1.0E-04          |
| Runt       | Runt domain                    | PF00853   | Beta-Scaffold Factors     | 1.0E-04          |
| STAT       | STAT_bind domain               | PF02864   | Beta-Scaffold Factors     | 1.0E-04          |
| ARID       | ARID domain                    | PF01388   | Helix-turn-helix          | 1.0E-04          |
| COE        | COE domain                     | AnimalTFDB| Helix-turn-helix          | 1.0E-04          |
| CUT        | Homeobox_CUT                  | PF02376   | Helix-turn-helix          | 1.0E-04          |
| E2F        | E2F_TDP domain                 | PF02319   | Helix-turn-helix          | 1.0E-04          |
| ETS        | Ets domain                     | PF00178   | Helix-turn-helix          | 1.0E-04          |
| Fork_head  | Fork_head domain               | PF00250   | Helix-turn-helix          | 1.0E-04          |
| Homeobox   | Homeobox                       | PF00046   | Helix-turn-helix          | 1.0E-04          |
| HPD        | HPD domain                     | PF05044   | Helix-turn-helix          | 1.0E-04          |
| HSF        | HSF_DNA-bind domain            | PF00447   | Helix-turn-helix          | 1.0E-04          |
| HTH        | HTH_psq domain                 | PF05225   | Helix-turn-helix          | 1.0E-04          |
| IRF        | IRF domain                     | PF00605   | Helix-turn-helix          | 1.0E-04          |
| MYB        | Myb_DNA-binding domain         | PF00249   | Helix-turn-helix          | 1.0E-04          |
| PAX        | PAX domain                     | PF00292   | Helix-turn-helix          | 1.0E-04          |
| Pou        | Homeobox_Pou                  | PF00157   | Helix-turn-helix          | 1.0E-04          |
| SRF        | SRF domain                     | PF00319   | Helix-turn-helix          | 1.0E-04          |
| TEA        | TEA domain                     | PF01285   | Helix-turn-helix          | 1.0E-04          |
| GTF2I      | GTF2I domain                   | PF02946   | Other Alpha-Helix Group   | 1.0E-04          |
| HMG        | HMG_box domain                 | PF00505   | Other Alpha-Helix Group   | 1.0E-04          |
| NF-YA      | CBFB_NFYA domain               | PF02045   | Other Alpha-Helix Group   | 1.0E-04          |
| NF-YB      | NF-YB domain                   | AnimalTFDB| Other Alpha-Helix Group   | 1.0E-18          |
| NF-YC      | NF-YC domain                   | AnimalTFDB| Other Alpha-Helix Group   | 1.0E-18          |
| SAND       | SAND domain                    | PF01342   | Other Alpha-Helix Group   | 1.0E-04          |
| AF-4       | AF-4 domain                    | PF05110   | Unclassified Structure    | 1.0E-03          |
| CG-1       | CG-1 domain                    | PF03859   | Unclassified Structure    | 1.0E-04          |
| CSRNP_N    | CSRNP_N domain                 | PF16019   | Unclassified Structure    | 1.0E-04          |
| CTF_NFI    | CTF/NFI and MH1 domain         | PF00859   | Unclassified Structure    | 1.0E-04          |
| DACH       | DACH domain                    | AnimalTFDB| Unclassified Structure    | 1.0E-04          |
| GCFC       | GCFC domain                    | PF07842   | Unclassified Structure    | 1.0E-04          |
| HMGA       | HMGA domain                    | AnimalTFDB| Unclassified Structure    | 1.0E-04          |
| LRRFIP     | LRRFIP domain                  | PF09738   | Unclassified Structure    | 1.0E-04          |
| MBD        | MBD domain                     | PF01429   | Unclassified Structure    | 1.0E-05          |
| MH1        | MH1 domain                     | PF03165   | Unclassified Structure    | 1.0E-04          |
| NCU-G1     | NCU-G1 domain                  | PF15065   | Unclassified Structure    | 1.0E-04          |
| NDT80_PhoG | NDT80_PhoG domain              | PF05224   | Unclassified Structure    | 1.0E-04          |
| PC4        | PC4 domain                     | PF02229   | Unclassified Structure    | 1.0E-04          |
| T-box      | T-box domain                   | PF00907   | Unclassified Structure    | 1.0E-04          |
| Tub        | Tub domain                     | PF01167   | Unclassified Structure    | 1.0E-04          |
| DM         | DM domain                      | PF00751   | Zinc-Coordinating Group   | 1.0E-04          |
| ESR-like   | zf-C4_ESR-like                 | AnimalTFDB| Zinc-Coordinating Group   | 1.0E-34          |
| NGFIB-like | zf-C4_NGFIB-like               | AnimalTFDB| Zinc-Coordinating Group   | 1.0E-24          |
| RXR-like   | zf-C4_RXR-like                 | AnimalTFDB| Zinc-Coordinating Group   | 1.0E-30          |
| SF-like    | zf-C4_SF-like                  | AnimalTFDB| Zinc-Coordinating Group   | 1.0E-40          |
| THAP       | THAP domain                    | PF05485   | Zinc-Coordinating Group   | 1.0E-04          |
| THR-like   | zf-C4_THR-like                 | AnimalTFDB| Zinc-Coordinating Group   | 1.0E-04          |
| ZBTB       | zf-C2H2_ZBTB                   | PF00651   | Zinc-Coordinating Group   | 1.0E-04          |
| zf-BED     | zf-BED domain                  | PF02892   | Zinc-Coordinating Group   | 1.0E-04          |
| zf-C2H2    | zf-C2H2 domain                 | PF00096   | Zinc-Coordinating Group   | 1.0E-05          |
| zf-C2HC    | zf-C2HC domain                 | PF01530   | Zinc-Coordinating Group   | 1.0E-04          |
| zf-CCCH    | zf-CCCH domain                 | PF00642   | Zinc-Coordinating Group   | 1.0E-20          |
| zf-GAGA    | zf-GAGA domain                 | PF09237   | Zinc-Coordinating Group   | 1.0E-05          |
| zf-GATA    | zf-GATA domain                 | PF00320   | Zinc-Coordinating Group   | 1.0E-04          |
| zf-LITAF-like | zf-LITAF-like domain        | PF10601   | Zinc-Coordinating Group   | 1.0E-04          |
| zf-MIZ     | zf-MIZ domain                  | PF02891   | Zinc-Coordinating Group   | 1.0E-04          |
| zf-NF-X1   | zf-NF-X1 domain                | PF01422   | Zinc-Coordinating Group   | 1.0E-04          |


