#!/usr/local/bioinf/bin/python3
#Modules
import os
import sys
import inspect
from time import sleep
import subprocess
import pandas as pd
import numpy as np
import argparse

print('''
----------------------------------------------------
-----------------------Start------------------------
''')
print('''
-----------------Building Pipeline------------------
''')
#User
user=subprocess.check_output('echo $USER', shell=True)
user=user.decode('utf-8')
user=user.replace("\n","")

#Args
parser = argparse.ArgumentParser()
parser.add_argument('--Input','-I',required=True)
parser.add_argument('--outputDir','-O', required=True, help='Path of the Directory where all output should be written to')
parser.add_argument('--GTFannotation','-GTF', required=True)
parser.add_argument('--GATKReferece','-GatkR', required=True)
parser.add_argument('--KnownIndels','-KI',required=True)
parser.add_argument('--FastaReferece','-Fasta',required=True)
parser.add_argument('--STARindex','-SI',required=True)
parser.add_argument('--LocalTmpDir','-LTD', help='Path to local tmp dir on the compute knode of HPC clusters')
parser.add_argument('--GlobalTmpDir','-GTD', help='Directory in wich temporary result file schould be written to, if not set a temp directory in your home folder will be created')
parser.add_argument('--Threads','-T', help= 'Number between 8 and ∞. It the max number of threads used for computing (it is onely used for Mapping) if not set onley one thread is used per process')
parser.add_argument('--GATK','-GATK', required=True)
parser.add_argument('--picard','-picard',required=True)
parser.add_argument('--VEP','-VEP',required=True)
parser.add_argument('--KallsitoIndex','-KID',required=True)
parser.add_argument('--VEPcacheDir','-VC',required=True)
parser.add_argument('--STARFusion','-SF',required=True)
parser.add_argument('--CTATresourceLib','-CTATlib',required=True)

args=parser.parse_args()
inputDirRNA=args.Input
outputDir=args.outputDir
GTFann=args.GTFannotation
GATKref=args.GATKReferece
KnownIndels=args.KnownIndels
STARindex=args.STARindex
arg_tmp_global=args.GlobalTmpDir
Threads=args.Threads
Threads=float(Threads)
ThreadsDiv1=int(round(Threads-2))
ThreadsDiv2=int(round(Threads/10*3))
ThreadsElse=1 if Threads is None else 4
ThreadsGATK=int(round(ThreadsElse/2))
ThreadsElse=str(ThreadsElse)
ThreadsBWA=1 if Threads is None else ThreadsDiv1
ThreadsSamtools=1 if Threads is None else ThreadsDiv2
Threads=args.Threads
Threads=1 if Threads is None else Threads
Threads=int(Threads)
ThreadsAll=str(Threads)
KallsitoIndex=args.KallsitoIndex
GATK=args.GATK
picard=args.picard
VEP=args.VEP
VEPcacheDir=args.VEPcacheDir
STARFusion=args.STARFusion
CTATlib=args.CTATresourceLib

#IO and directories
Input = [os.path.basename(x) for x in os.listdir(inputDirRNA) if x.endswith('R1.fastq.gz')]
Basename = ' '.join(Input).replace('R1.fastq.gz','').split()
CWD = (sys.path[0])
if arg_tmp_global is None:
    os.system('mkdir -p /home/%s/EXOMEtmp'%user)
    tmp_global=('/home/%s/EXOMEtmp'%user)
else:
    tmp_global=arg_tmp_global
arg_tmp_local=args.LocalTmpDir
if arg_tmp_local is None:
    tmp_local=tmp_global
else:
    tmp_local=arg_tmp_local
tmp_local=('%s/RNApipe'%tmp_local)
NodeList = subprocess.check_output('qconf -sh', shell=True)
NodeList=NodeList.decode('utf-8')
NodeList=NodeList.replace("\n",",")
NodeList=NodeList.split(',')
NodeList.pop()
for item in NodeList:
    os.system('''ssh %s 'mkdir -p %s' '''%(item,tmp_local))
os.system('mkdir -p %s/Jobscipts'%tmp_global)
Jobscripts=('%s/Jobscipts'%tmp_global)

#Jobscript_template SGE
class SGE_Jobscript:
    def __init__(self,Threads,Jobname,Error,Log,Command):
        self.Threads=Threads
        self.Jobname=Jobname
        self.Error=Error
        self.Log=Log
        self.Command=Command

    def SGE_Jobscript_template(self):
        print('''#!/bin/sh
#$ -S /bin/sh
        
#$ -pe smp '''+self.Threads,'''
#$ -V 
        
#### Jobdescription at qstat
#$ -N '''+self.Jobname,'''

#### Error Outputfile
#$ -e '''+self.Error,'''
#$ -o '''+self.Log,'''

#### Resubmit
#$ -r y

'''+self.Command)

#Mapping RNA
os.system('mkdir -p %s/Resuts/Mapping/RNA'%outputDir)
Results_mapping_RNA=('%s/Resuts/Mapping/RNA'%outputDir)
os.system('mkdir -p %s/Log/Mapping/RNA'%outputDir)
LogsSTAR=('%s/Log/Mapping/RNA'%outputDir)
for item in Basename:
    Threads=ThreadsAll
    Name=('Mapping_RNA%s'%item)
    Error=('%s/%s_error.dat'%(LogsSTAR,item))
    Log=('%s/%s_log.dat'%(LogsSTAR,item))
    Command=('STAR --genomeDir %s --readFilesIn %s/%sR1.fastq.gz %s/%sR2.fastq.gz --outFileNamePrefix %s/%s --runThreadN %s \
--outReadsUnmapped None --twopassMode Basic --readFilesCommand "gunzip -c" --outSAMstrandField intronMotif --outSAMunmapped Within \
--outSAMtype BAM SortedByCoordinate --chimSegmentMin 12 --chimJunctionOverhangMin 12 \
--chimOutJunctionFormat 1 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5 \
--chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG -4 --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 --outSAMattrRGline ID:ExomeSeq SM:%s/%s_RNA.dedup.sorted.bam PL:Illumina'
%(STARindex,inputDirRNA,item,inputDirRNA,item,Results_mapping_RNA,item,Threads,Results_mapping_RNA,item,))
    Mapping_RNA = SGE_Jobscript(Threads,Name,Error,Log,Command)
    sys.stdout = open('%s/Mapping_RNA_%s.sh'%(Jobscripts,item),'w')
    Mapping_RNA.SGE_Jobscript_template()
    sys.stdout.close()

#Transcript Qiantification
os.system('mkdir -p %s/Resuts/Kallisto'%outputDir)
Results_Kallisto=('%s/Resuts/Kallisto'%outputDir)
os.system('mkdir -p %s/Log/Kallisto'%outputDir)
LogsKallisto=('%s/Log/Kallisto'%outputDir)
for item in Basename:
    Threads=ThreadsElse
    Name=('Kallisto_%s'%item)
    Error=('%s/%s_error.dat'%(LogsKallisto,item))
    Log=('%s/%s_log.dat'%(LogsKallisto,item))
    Command=('kallisto quant -i %s  -o %s/%s  -b 100 -l 200 -s 30 --threads=%s %s/%sR1.fastq.gz %s/%sR2.fastq.gz'%(KallsitoIndex,Results_Kallisto,item,ThreadsElse,inputDirRNA,item,inputDirRNA,item))
    Kallisto = SGE_Jobscript(Threads,Name,Error,Log,Command)
    sys.stdout = open('%s/Kallisto_%s.sh'%(Jobscripts,item),'w')
    Kallisto.SGE_Jobscript_template()
    sys.stdout.close()

#RNA_Sorting_dublicateRemouval_and_Indexing
os.system('mkdir -p %s/Log/RNA_SortDedupIndex'%outputDir)
LogsIndexing=('%s/Log/RNA_SortDedupIndex'%outputDir)
for item in Basename:
    Threads=str(int(ThreadsSamtools*3))
    Name=('RNA_SortDedupIndex_%s'%item)
    Error=('%s/%s_error.dat'%(LogsIndexing,item))
    Log=('%s/%s_log.dat'%(LogsIndexing,item))
    Command=('mkfifo %s/%s_RNA.NameSorted.fifo &&\n\
mkfifo %s/%s_RNA.Fixmate.fifo &&\n\
mkfifo %s/%s_RNA.Sorted.fifo &&\n\
chmod 777 %s/*.fifo\n\
sleep 5 &&\n\
samtools sort -@ %s -n -o %s/%s_RNA.NameSorted.fifo %s/%sAligned.sortedByCoord.out.bam  &\n\
samtools fixmate -@ %s -m %s/%s_RNA.NameSorted.fifo %s/%s_RNA.Fixmate.fifo &\n\
samtools sort -@ %s -o %s/%s_RNA.Sorted.fifo %s/%s_RNA.Fixmate.fifo  &\n\
samtools markdup -r -@ %s %s/%s_RNA.Sorted.fifo %s/%s_RNA.dedup.sorted.bam &&\n\
samtools index %s/%s_RNA.dedup.sorted.bam &&\n\
rm -f %s/%s_RNA.NameSorted.fifo &&\n\
rm -f %s/%s_RNA.Fixmate.fifo &&\n\
rm -f %s/%s_RNA.Sorted.fifo'%(tmp_local,item,
tmp_local,item,
tmp_local,item,
tmp_local,
ThreadsSamtools,tmp_local,item,Results_mapping_RNA,item,
ThreadsSamtools,tmp_local,item,tmp_local,item,
ThreadsSamtools,tmp_local,item,tmp_local,item,
ThreadsSamtools,tmp_local,item,Results_mapping_RNA,item,
Results_mapping_RNA,item,
tmp_local,item,
tmp_local,item,
tmp_local,item))
    Indexing = SGE_Jobscript(Threads,Name,Error,Log,Command)
    sys.stdout = open('%s/RNA_SortDedupIndex_%s.sh'%(Jobscripts,item),'w')
    Indexing.SGE_Jobscript_template()
    sys.stdout.close()

#SplitNcigarReads and BQSR RNA Data
os.system('mkdir -p %s/Log/RNA_SplitNCigarReadsAndBQSR'%outputDir)
LogsIndexing=('%s/Log/RNA_SplitNCigarReadsAndBQSR'%outputDir)
for item in Basename:
    Threads=ThreadsElse
    Name=('RNA_SplitNCigarReadsAndBQSR_%s'%item)
    Error=('%s/%s_error.dat'%(LogsIndexing,item))
    Log=('%s/%s_log.dat'%(LogsIndexing,item))
    Command=('java -jar -Djava.io.tmpdir=%s/ -XX:ConcGCThreads=%s -XX:ParallelGCThreads=%s %s SplitNCigarReads -R %s -I %s/%s_RNA.dedup.sorted.bam -O %s/%s_RNA.SplitNCigarReads.bam &&\n\
java -jar -Djava.io.tmpdir=%s/ -XX:ConcGCThreads=%s -XX:ParallelGCThreads=%s %s BaseRecalibrator -R %s -I %s/%s_RNA.SplitNCigarReads.bam --known-sites %s -O %s/%s_RNA.recalDataTable.tab &&\n\
java -jar -Djava.io.tmpdir=%s/ -XX:ConcGCThreads=%s -XX:ParallelGCThreads=%s %s ApplyBQSR -R %s -I %s/%s_RNA.SplitNCigarReads.bam --bqsr-recal-file %s/%s_RNA.recalDataTable.tab -O %s/%s_RNA.dedup.sorted.BQSR.bam &&\n\
samtools index %s/%s_RNA.dedup.sorted.BQSR.bam &&\n\
rm -f %s/%s_RNA.recalDataTable.tab &&\n\
rm -f %s/%s_RNA.SplitNCigarReads.bam'%(tmp_local,ThreadsGATK,ThreadsGATK,GATK,GATKref,Results_mapping_RNA,item,tmp_local,item,
tmp_local,ThreadsGATK,ThreadsGATK,GATK,GATKref,tmp_local,item,KnownIndels,tmp_local,item,
tmp_local,ThreadsGATK,ThreadsGATK,GATK,GATKref,tmp_local,item,tmp_local,item,Results_mapping_RNA,item,
Results_mapping_RNA,item,
tmp_local,item,
tmp_local,item))
    RNA_SplitNCigarReadsAndBQSR = SGE_Jobscript(Threads,Name,Error,Log,Command)
    sys.stdout = open('%s/RNA_SplitNCigarReadsAndBQSR_%s.sh'%(Jobscripts,item),'w')
    RNA_SplitNCigarReadsAndBQSR.SGE_Jobscript_template()
    sys.stdout.close()

#Fusion Calling
os.system('mkdir -p %s/Resuts/Fusions'%outputDir)
Fusions=('%s/Resuts/Fusions'%outputDir)
os.system('mkdir -p %s/Log/STARFusion'%outputDir)
LogsSTARFusion=('%s/Log/STARFusion'%outputDir)
for item in Basename:
    Threads=str(1)
    Name=('STARFusion_%s'%item)
    Error=('%s/%s_error.dat'%(LogsSTARFusion,item))
    Log=('%s/%s_log.dat'%(LogsSTARFusion,item))
    Command=('mkdir -p %s/%s_Fusions &&\n\
%s --genome_lib_dir %s -J %s/%sChimeric.out.junction --output_dir %s/%s_Fusions'
%(Fusions,item,
STARFusion,CTATlib,Results_mapping_RNA,item,Fusions,item))
    STARFusionJob = SGE_Jobscript(Threads,Name,Error,Log,Command)
    sys.stdout = open('%s/STARFusion_%s.sh'%(Jobscripts,item),'w')
    STARFusionJob.SGE_Jobscript_template()
    sys.stdout.close()

#Haplotype caller 
os.system('mkdir -p %s/Resuts/Variants/Unfiltered/HaplotypeCaller'%outputDir)
VariantsHaplotypeCaller=('%s/Resuts/Variants/Unfiltered/HaplotypeCaller'%outputDir)
os.system('mkdir -p %s/Log/VariantCalling/HaplotypeCaller'%outputDir)
LogsHaplotypeCaller=('%s/Log/VariantCalling/HaplotypeCaller'%outputDir)
for item in Basename:
    Threads=ThreadsElse
    Name=('VC_HapplotypeCaller_%s'%item)
    Error=('%s/%s_error.dat'%(LogsHaplotypeCaller,item))
    Log=('%s/%s_log.dat'%(LogsHaplotypeCaller,item))
    Command=('java -jar -Djava.io.tmpdir=%s/ -XX:ConcGCThreads=%s -XX:ParallelGCThreads=%s %s HaplotypeCaller -R %s -I %s/%s_RNA.dedup.sorted.BQSR.bam --stand-call-conf 30.0 -O %s/%s_RNA.vcf '
    %(tmp_local,ThreadsElse,ThreadsElse,GATK,GATKref,Results_mapping_RNA,item,VariantsHaplotypeCaller,item))
    VCHapplotypeCaller = SGE_Jobscript(Threads,Name,Error,Log,Command)
    sys.stdout = open('%s/VC_HapplotypeCaller_%s.sh'%(Jobscripts,item),'w')
    VCHapplotypeCaller.SGE_Jobscript_template()
    sys.stdout.close()

#Variant Filtering
os.system('mkdir -p %s/Variants/Filtered'%tmp_global)
FilterdVariants=('%s/Variants/Filtered'%tmp_global)
os.system('mkdir -p %s/Log/VariantFiltering'%outputDir)
LogsVariantFiltering=('%s/Log/VariantFiltering'%outputDir)
for item in Basename:
    Threads=ThreadsElse
    Name=('VF_%s'%item)
    Error=('%s/%s_error.dat'%(LogsVariantFiltering,item))
    Log=('%s/%s_log.dat'%(LogsVariantFiltering,item))
    Command=('java -jar -Djava.io.tmpdir=%s -XX:ConcGCThreads=%s -XX:ParallelGCThreads=%s %s FilterVcf INPUT=%s/%s_diseased.vcf OUTPUT=%s/%s_FilteredVariants.vcf MIN_DP=10'
    %(tmp_local,ThreadsElse,ThreadsElse,picard,VariantsHaplotypeCaller,item,FilterdVariants,item))
    VariantFiltering = SGE_Jobscript(Threads,Name,Error,Log,Command)
    sys.stdout = open('%s/VariantFiltering_%s.sh'%(Jobscripts,item),'w')
    VariantFiltering.SGE_Jobscript_template()
    sys.stdout.close()

#Annotation
os.system('mkdir -p %s/Resuts/Annotation'%outputDir)
AnnotationDir=('%s/Resuts/Annotation'%outputDir)
os.system('mkdir -p %s/Log/Annotation'%outputDir)
LogsAnnotation=('%s/Log/Annotation'%outputDir)
for item in Basename:
    Threads=ThreadsElse
    Name=('Annotation_%s'%item)
    Error=('%s/%s_error.dat'%(LogsAnnotation,item))
    Log=('%s/%s_log.dat'%(LogsAnnotation,item))
    Command=('%s -i %s/%s_FilteredVariants.vcf -o %s/%s_VEPann.txt --cache --dir /data/databases/vep/ --assembly GRCh38 --offline --sift b --variant_class --polyphen b --humdiv --gene_phenotype --regulatory --symbol --force_overwrite'
    %(VEP,FilterdVariants,item,AnnotationDir,item))
    Annotation = SGE_Jobscript(Threads,Name,Error,Log,Command)
    sys.stdout = open('%s/Annotation_%s.sh'%(Jobscripts,item),'w')
    Annotation.SGE_Jobscript_template()
    sys.stdout.close()

sys.stdout = sys.__stdout__

#for item in Basename:
#    os.system('qsub -q long.q %s/Kallisto_%s.sh'%(Jobscripts,item))
#    os.system('qsub -q long.q %s/Mapping_RNA_%s.sh'%(Jobscripts,item))

StateMapping = subprocess.check_output("qstat", shell=True)
while b'Mapping' in StateMapping:
    sleep (1)
    StateMapping = subprocess.check_output("qstat", shell=True)
    if b'Mapping' not in StateMapping:
        continue



#Counting gene Features
Star_output = [os.path.basename(x) for x in os.listdir(Results_mapping_RNA) if x.endswith(".bam")]
Star_output_fullPath = [Results_mapping_RNA +('/')+ item for item in Star_output] 
FeatureCounts_input = ' '.join(Star_output_fullPath)

os.system('mkdir -p %s/Resuts/FeatureCounts'%outputDir)
Results_FeatureCounts=('%s/Resuts/FeatureCounts'%outputDir)
os.system('mkdir -p %s/Log/FeatureCounts'%outputDir)
LogsFeatureCounts=('%s/Log/FeatureCounts'%outputDir)
Threads=ThreadsElse
Name=('FeatureCounts')
Error=('%s/%s_error.dat'%(LogsFeatureCounts,item))
Log=('%s/%s_log.dat'%(LogsFeatureCounts,item))
Command=('featureCounts -t exon -g gene_id --primary -a %s -o %s/Raw.featureCounts.tab %s '%(GTFann, Results_FeatureCounts, FeatureCounts_input))
FeatureCounts = SGE_Jobscript(Threads,Name,Error,Log,Command)
sys.stdout = open('%s/FeatureCounts.sh'%(Jobscripts),'w')
FeatureCounts.SGE_Jobscript_template()
sys.stdout.close()

sys.stdout = sys.__stdout__

for item in Basename:
    os.system('qsub -q long.q %s/STARFusion_%s.sh'%(Jobscripts,item))
    os.system('qsub -q long.q %s/RNA_SortDedupIndex_%s.sh'%(Jobscripts,item))

os.system('qsub -q long.q %s/FeatureCounts.sh'%(Jobscripts))

print('''
-----------------Indexing Bam Files-----------------
''')

StateIndexingRNA = subprocess.check_output("qstat", shell=True)
while b'RNA' in StateIndexingRNA:
    sleep (1)
    StateIndexingRNA = subprocess.check_output("qstat", shell=True)
    if b'RNA' not in StateIndexingRNA:
        continue

for item in Basename:
    os.system('qsub -q long.q %s/RNA_SplitNCigarReadsAndBQSR_%s.sh'%(Jobscripts,item))

StateRNA_Split = subprocess.check_output("qstat", shell=True)
while b'RNA_Split' in StateRNA_Split:
    sleep (1)
    StateRNA_Split = subprocess.check_output("qstat", shell=True)
    if b'RNA_Split' not in StateRNA_Split:
        continue

for item in Basename:
    os.system('qsub -q long.q %s/VC_HapplotypeCaller_%s.sh'%(Jobscripts,item))

State = subprocess.check_output("qstat", shell=True)
while b'FeatureCou' in State:
    sleep (1)
    State = subprocess.check_output("qstat", shell=True)
    if b'FeatureCou' not in State:
        continue

print('''
------------------Calculating TPMs------------------
''')

df_feature_counts_raw=pd.read_table('%s/Raw.featureCounts.tab'%(Results_FeatureCounts), skiprows=1)
df_feature_counts_for_TPM=df_feature_counts_raw.drop(columns=['Chr','Start','End','Strand'])
df_feature_counts_for_DE=df_feature_counts_raw.drop(columns=['Chr','Start','End','Strand','Length'])

for item in Basename:
    df_feature_counts_for_DE.rename(columns={'%s/%sAligned.sortedByCoord.out.bam'%(Results_mapping_RNA,item):'%s'%item}, inplace=True)
    df_feature_counts_for_TPM.rename(columns={'%s/%sAligned.sortedByCoord.out.bam'%(Results_mapping_RNA,item):'%s'%item}, inplace=True)

df_feature_counts_for_DE.to_csv('%s/featureCounts_for_DEseq2.tab'%(Results_FeatureCounts), sep='\t', index=None)

for item in Basename:
    df_feature_counts_for_TPM['%sTPM_tmp'%item]=df_feature_counts_for_TPM.Length/df_feature_counts_for_TPM['%s'%item]
    
df_feature_counts_for_TPM=df_feature_counts_for_TPM.replace(np.inf, 0 ) 
df_TPM=df_feature_counts_for_DE['Geneid']

df_TPM= pd.DataFrame()
df_TPM['Geneid']=df_feature_counts_for_DE['Geneid']
for item in Basename:
    column='%sTPM_tmp'%item
    Sum_for_tpm=df_feature_counts_for_TPM[column].sum()
    Sumunder1=1/Sum_for_tpm
    TPM=(Sumunder1*df_feature_counts_for_TPM[column])*1000000
    df_TPM['%sTPM'%item]=TPM

df_TPM.to_csv('%s/TPMs.tab'%(Results_FeatureCounts), sep='\t', index=None)

StateVC = subprocess.check_output("qstat", shell=True)
while b'VC' in StateVC:
    sleep (1)
    StateVC = subprocess.check_output("qstat", shell=True)
    if b'VC' not in StateVC:
        continue


for item in Basename:
    os.system('qsub -q long.q %s/VariantFiltering_%s.sh'%(Jobscripts,item))

print('''
-----------------Filtering Variants-----------------
''')

StateVF = subprocess.check_output("qstat", shell=True)
while b'VF' in StateVF:
    sleep (1)
    StateVF = subprocess.check_output("qstat", shell=True)
    if b'VF' not in StateVF:
        continue

for item in Basename:
    os.system('qsub -q long.q %s/Annotation_%s.sh'%(Jobscripts,item))

print('''
----------------Annotating Variants-----------------
''')

StateAnnotation = subprocess.check_output("qstat", shell=True)
while b'Annotation' in StateAnnotation:
    sleep (1)
    StateAnnotation = subprocess.check_output("qstat", shell=True)
    if b'Annotation' not in StateAnnotation:
        continue    

print('''
------------------------Done------------------------
----------------------------------------------------
''')
