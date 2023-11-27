#Modules
import os
import sys
import inspect
from time import sleep
import subprocess
import pandas as pd
import numpy as np
import argparse

CWD = (sys.path[0])

print('''
----------------------------------------------------
-----------------------Start------------------------
''')
#Args
parser = argparse.ArgumentParser()
parser.add_argument('--InputDirNormal', '-IN', required=True, help='Path of the Directory with your input fasq files')
parser.add_argument('--InputDirDiseased', '-ID', required=True, help='Path of the Directory with your input fasq files')
parser.add_argument('--outputDir','-O', required=True, help='Path of the Directory where all output should be written to')



args=parser.parse_args()
inputDirNormal=args.InputDirNormal 
inputDirDiseased=args.InputDirDiseased
outputDir=args.outputDir


#Path to resources
Path_to_BWA_reference=('/data/genomes/hg38/index/bwa/gdc/GRCh38.d1.vd1/GRCh38.d1.vd1.fa')
Path_to_hg38refGene_gtf=('/data/genomes/hg38/annotation/ucsc/hg38refGene.gtf')
Path_to_GATK_reference=('/data/genomes/hg38/fasta/gdc/GRCh38.d1.vd1/GRCh38.d1.vd1.fa')
Path_to_known_indel_sites=('/data/databases/GATKresourceBundle/Mills_and_1000G_gold_standard.indels.hg38.vcf')
Path_to_fasta_reference=('/data/genomes/hg38/fasta/gdc/GRCh38.d1.vd1/GRCh38.d1.vd1.fa')
Path_to_input_files_cancer=inputDirDiseased
Path_to_input_files_normal=inputDirNormal

#temp Directorys
local_tmp=('/local/scratch/gronauer')
os.system('mkdir -p %s/tmp'%CWD)
global_tmp=('%s/tmp'%CWD)
os.system('mkdir -p %s/JobScripts'%CWD)
Job_Scripts =('%s/JobScripts'%CWD) 

#Log and error Directorys
os.system('mkdir -p %s/Log'%CWD)
log=('%s/Log'%CWD)
os.system('mkdir -p %s/Error'%CWD)
error=('%s/Error'%CWD)

#Path To output Directorys
Path_to_Fastq_files_cancer=Path_to_input_files_cancer
Path_to_Fastq_files_normal=Path_to_input_files_normal

os.system('mkdir -p /data/projects/2020/OvarianCancerHH/Exome_Analysis/Results/Mapped_normals')
Path_to_Mapped_normals=('/data/projects/2020/OvarianCancerHH/Exome_Analysis/Results/Mapped_normals')
os.system('mkdir -p /data/projects/2020/OvarianCancerHH/Exome_Analysis/Results/Mapped_tumor')
Path_to_Mapped_tumor=('/data/projects/2020/OvarianCancerHH/Exome_Analysis/Results/Mapped_tumor')

os.system('mkdir -p %s/Variants/HaplotypeCaller'%global_tmp)
Path_to_HaplotypeCaller_variant_calls=('%s/Variants/HaplotypeCaller'%global_tmp)
os.system('mkdir -p %s/Variants/Mutect2'%global_tmp)
Path_to_Mutect2_variant_calls=('%s/Variants/Mutect2'%global_tmp)
os.system('mkdir -p %s/Variants/SomaticSniper'%global_tmp)
Path_to_SomaticSniper_variant_calls=('%s/Variants/SomaticSniper'%global_tmp)
os.system('mkdir -p %s/Variants/Varscan2'%global_tmp)
Path_to_Varscan2_variant_calls=('%s/Variants/Varscan2'%global_tmp)
os.system('mkdir -p %s/Variants/Strelka'%global_tmp)
Path_to_Strelka_dir=('%s/Variants/Strelka'%global_tmp)
os.system('mkdir -p /data/projects/2020/OvarianCancerHH/Exome_Analysis/Results/Variants/unfiltered')
Path_to_all_vcfs_unfiltered=('/data/projects/2020/OvarianCancerHH/Exome_Analysis/Results/Variants/unfiltered')
os. system('mkdir -p %s/Variants/filtering_tmp'%global_tmp)
Variant_filtering_tmp = ('%s/Variants/filtering_tmp'%global_tmp)
os.system('mkdir -p /data/projects/2020/OvarianCancerHH/Exome_Analysis/Results/Variants/filtered')
Path_to_all_vcfs_filtered=('/data/projects/2020/OvarianCancerHH/Exome_Analysis/Results/Variants/filtered')

files_path1 = [os.path.basename(x) for x in os.listdir(Path_to_input_files_cancer) if x.endswith('R1.fastq.gz')]
Basename = ' '.join(files_path1).replace('_R1.fastq.gz','').split()


#Mapping 
os.system('mkdir -p %s/Mapping'%log)
os.system('mkdir -p %s/Mapping'%error)
Mapping_log=('%s/Mapping'%log)
Mapping_error=('%s/Mapping'%error)

#Mapping normal 
for item in Basename:
    with open(('%s/Mapping_normal%s.sh'%(Job_Scripts, item)), 'w') as Mapping_normal_Jobscript:
        Mapping_normal_Jobscript.write('''
#!/bin/sh
#$ -S /bin/sh

#$ -pe smp 10
#$ -cwd
#$ -V

#### Jobdescription at qstat
#$ -N Mapping_normal%s 

#### Error Outputfile
#$ -e %s/%s_normal.dat
#$ -o %s/%s_normal.dat

#### Resubmit
#$ -r y

ulimit -n 2000

mkdir -p %s/mapping_normal
bwa mem -t 10 -M %s %s/%s_normal_R1.fastq.gz %s/%s_normal_R2.fastq.gz > %s/mapping_normal/%s_normal.bam
'''%(item, 
Mapping_error,item,
Mapping_log,item,
global_tmp,
Path_to_BWA_reference,Path_to_Fastq_files_normal, item, Path_to_Fastq_files_normal, item, global_tmp, item))
Mapping_normal_Jobscript.close

#Mapping tumor
for item in Basename:
    with open(('%s/Mapping_tumor%s.sh'%(Job_Scripts, item)), 'w') as Mapping_tumor_Jobscript:
        Mapping_tumor_Jobscript.write('''
#!/bin/sh
#$ -S /bin/sh

#$ -pe smp 10
#$ -cwd
#$ -V

#### Jobdescription at qstat
#$ -N Mapping_tumor%s 

#### Error Outputfile
#$ -e %s/%s_tumor.dat
#$ -o %s/%s_tumor.dat

#### Resubmit
#$ -r y

ulimit -n 2000

mkdir -p %s/mapping_tumor
bwa mem -t 10 -M %s %s/%s_tumor_R1.fastq.gz %s/%s_tumor_R2.fastq.gz > %s/mapping_tumor/%s_tumor.bam
'''%(item,
Mapping_error,item,
Mapping_log,item,
global_tmp,
Path_to_BWA_reference, Path_to_Fastq_files_cancer, item, Path_to_Fastq_files_cancer, item, global_tmp, item))
Mapping_tumor_Jobscript.close

for item in Basename:
    os.system('qsub -q long.q  %s/Mapping_tumor%s.sh'%(item,Job_Scripts,item))
    os.system('qsub -q long.q %s/Mapping_normal%s.sh'%(item,Job_Scripts,item))

#Duplicate Removual
os.system('mkdir -p %s/Duplicate_remouval'%log)
os.system('mkdir -p %s/Duplicate_remouval'%error)
Duplicate_remouval_log=('%s/Duplicate_remouval'%log)
Duplicate_remouval_error=('%s/Duplicate_remouval'%error)

for item in Basename:
    with open ('%s/Duplicate_remouval%s_normal.sh'%(Job_Scripts, item), 'w') as Duplicate_remouval_jobscript_normal:
        Duplicate_remouval_jobscript_normal.write('''
#!/bin/sh
#$ -S /bin/sh

#$ -pe smp 2
#$ -cwd
#$ -V

#### Jobdescription at qstat
#$ -N Duplicate_remouval_normal%s 

#### Error Outputfile
#$ -e %s/%s_normal.dat
#$ -o %s/%s_normal.dat

#### Resubmit
#$ -r y

ulimit -n 2000      
samtools sort -@2 -m 20G %s/mapping_normal/%s_normal.bam >%s/mapping_normal/%s_sorted_normal.bam
mkdir -p %s/AddOrReplaceReadGroups_tmp%s;
mkdir -p %s/AddOrReplaceReadGroups_tmp/%s;
java -jar -Djava.io.tmpdir=%s/AddOrReplaceReadGroups_tmp%s/ -XX:ConcGCThreads=2 -XX:ParallelGCThreads=2 /usr/local/bioinf/picard/picard-tools-2.21.4/picard.jar AddOrReplaceReadGroups I=%s/mapping_normal/%s_sorted_normal.bam O=%s/AddOrReplaceReadGroups_tmp/%s/%s_normal_with_readgroup.bam RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=%s/%s_normal_ddup.bam;
mkdir -p %s/MarkDuplicates_tmp%s;
java -jar -Djava.io.tmpdir=%s/MarkDuplicates_tmp%s/ -XX:ConcGCThreads=2 -XX:ParallelGCThreads=2 /usr/local/bioinf/picard/picard-tools-2.21.4/picard.jar MarkDuplicates I=%s/AddOrReplaceReadGroups_tmp/%s/%s_normal_with_readgroup.bam  O=%s/%s_normal_ddup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT QUIET=true M=output.metrics && 
rm -f %s/mapping_tumor/%s_sorted_normal.bam &&
rm -r -f %s/AddOrReplaceReadGroups_tmp%s &&
rm -r -f %s/AddOrReplaceReadGroups_tmp/%s &&
rm -f %s/mapping_normal/%s_normal.bam &&
rm -r -f%s/MarkDuplicates_tmp%s
'''%(item,
Duplicate_remouval_error,item,
Duplicate_remouval_log,item,
global_tmp,item,global_tmp,item,
local_tmp,item,
global_tmp,item,
local_tmp,item,global_tmp,item,global_tmp,item,item,Path_to_Mapped_normals,item,
local_tmp,item,
local_tmp,item,global_tmp,item,item,Path_to_Mapped_normals,item,
global_tmp,item,
local_tmp,item,
global_tmp,item,
global_tmp,item,
local_tmp,item
))
Duplicate_remouval_jobscript_normal.close

for item in Basename:
    with open ('%s/Duplicate_remouval%s_tumor.sh'%(Job_Scripts, item), 'w') as Duplicate_remouval_jobscript_tumor:
        Duplicate_remouval_jobscript_tumor.write('''
#!/bin/sh
#$ -S /bin/sh

#$ -pe smp 2
#$ -cwd
#$ -V

#### Jobdescription at qstat
#$ -N Dublicate_remouval_tumor%s 

#### Error Outputfile
#$ -e %s/%s_tumor.dat
#$ -o %s/%s_tumor.dat

#### Resubmit
#$ -r y

ulimit -n 2000        
samtools sort -@2 -m 20G %s/mapping_tumor/%s_tumor.bam > %s/mapping_tumor/%s_sorted_tumor.bam
mkdir -p %s/AddOrReplaceReadGroups_tmp%s;
mkdir -p %s/AddOrReplaceReadGroups_tmp/%s;
java -jar -Djava.io.tmpdir=%s/AddOrReplaceReadGroups_tmp%s/ -XX:ConcGCThreads=2 -XX:ParallelGCThreads=2 /usr/local/bioinf/picard/picard-tools-2.21.4/picard.jar AddOrReplaceReadGroups I=%s/mapping_tumor/%s_sorted_tumor.bam O=%s/AddOrReplaceReadGroups_tmp/%s/%s_tumor_with_readgroup.bam RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=%s/%s_tumor_ddup.bam;
mkdir -p %s/MarkDuplicates_tmp%s;
java -jar -Djava.io.tmpdir=%s/MarkDuplicates_tmp%s/ -XX:ConcGCThreads=2 -XX:ParallelGCThreads=2 /usr/local/bioinf/picard/picard-tools-2.21.4/picard.jar MarkDuplicates I=%s/AddOrReplaceReadGroups_tmp/%s/%s_tumor_with_readgroup.bam  O=%s/%s_tumor_ddup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT QUIET=true M=output.metrics && 
rm -f %s/mapping_tumor/%s_sorted_tumor.bam &&
rm -r -f %s/AddOrReplaceReadGroups_tmp%s &&
rm -r -f %s/AddOrReplaceReadGroups_tmp/%s &&
rm -f %s/mapping_tumor/%s_tumor.bam &&
rm -r -f%s/MarkDuplicates_tmp%s
'''%(item,
Duplicate_remouval_error,item,
Duplicate_remouval_log,item,
global_tmp,item,global_tmp,item,
local_tmp,item,
global_tmp,item,
local_tmp,item,global_tmp,item,global_tmp,item,item,Path_to_Mapped_tumor,item,
local_tmp,item,
local_tmp,item,global_tmp,item,item,Path_to_Mapped_tumor,item,
global_tmp,item,
local_tmp,item,
global_tmp,item,
global_tmp,item,
local_tmp,item))
Duplicate_remouval_jobscript_tumor.close

for item in Basename:
    os.system('qsub -q long.q -hold_jid Mapping_tumor%s %s/Duplicate_remouval%s_tumor.sh'%(item, Job_Scripts, item))
    os.system('qsub -q long.q -hold_jid Mapping_normal%s %s/Duplicate_remouval%s_normal.sh'%(item, Job_Scripts, item))

#Indexing
os.system('mkdir -p %s/Indexing'%log)
os.system('mkdir -p %s/Indexing'%error)
Indexing_log=('%s/Indexing'%log)
Indexing_error=('%s/Indexing'%error)

for item in Basename:
    with open ('%s/Indexing%s.sh'%(Job_Scripts, item), 'w') as Indexing_jobscript:
        Indexing_jobscript.write('''
#!/bin/sh
#$ -S /bin/sh

#$ -pe smp 1
#$ -cwd
#$ -V

#### Jobdescription at qstat
#$ -N Indexing%s 

#### Error Outputfile
#$ -e %s/%s.dat
#$ -o %s/%s.dat

#### Resubmit
#$ -r y

ulimit -n 2000        
samtools index %s/%s_normal_ddup.bam %s/%s_normal_ddup.bam.bai ;
samtools index %s/%s_tumor_ddup.bam %s/%s_tumor_ddup.bam.bai ;
cp %s/%s_normal_ddup.bam.bai %s/%s_normal_ddup.bai ;
cp %s/%s_tumor_ddup.bam.bai %s/%s_tumor_ddup.bai 
'''%(item,
Indexing_error,item,
Indexing_log,item,
Path_to_Mapped_normals,item,Path_to_Mapped_normals,item,
Path_to_Mapped_tumor,item,Path_to_Mapped_tumor,item,
Path_to_Mapped_normals,item,Path_to_Mapped_normals,item,
Path_to_Mapped_tumor,item,Path_to_Mapped_tumor,item
))
Indexing_jobscript.close

print('''
--------------!!Processing Aligments!!--------------
-----------------------Bam2fq-----------------------
-----------------------Mapping----------------------
------------------Removing Dublicates---------------
''')
StateBam2fq = subprocess.check_output("qstat", shell=True)
while b'Bam2fq' in StateBam2fq:
    sleep (1)
    StateBam2fq = subprocess.check_output("qstat", shell=True)
    if b'Bam2fq' not in StateBam2fq:
        continue

StateMapping = subprocess.check_output("qstat", shell=True)
while b'Mapping' in StateMapping:
    sleep (1)
    StateMapping = subprocess.check_output("qstat", shell=True)
    if b'Mapping' not in StateMapping:
        continue

StateDuplicate_remouval = subprocess.check_output("qstat", shell=True)
while b'Duplicate' in StateDuplicate_remouval:
    sleep (1)
    StateDuplicate_remouval = subprocess.check_output("qstat", shell=True)
    if b'Duplicate' not in StateDuplicate_remouval:
        continue

StateBams_normal = subprocess.check_output("ls /data/projects/2020/OvarianCancerHH/Exome_Analysis/Results/Mapped_normals/", shell=True)
StateBams_tumor = subprocess.check_output("ls /data/projects/2020/OvarianCancerHH/Exome_Analysis/Results/Mapped_tumor/", shell=True)

for item in Basename:
    Bam_file_tumor='%s_tumor_ddup.bam'%item
    bt = bytes(Bam_file_tumor, 'utf-8')
    while bt not in StateBams_tumor:
            sleep(1)
            StateBams_tumor = subprocess.check_output("ls /data/projects/2020/OvarianCancerHH/Exome_Analysis/Results/Mapped_tumor/", shell=True)
            if bt in StateBams_tumor:
                continue

for item in Basename:
    Bam_file_normal='%s_normal_ddup.bam'%item
    bn = bytes(Bam_file_normal, 'utf-8')
    while bn not in StateBams_normal:
            sleep(1)
            StateBams_normal = subprocess.check_output("ls /data/projects/2020/OvarianCancerHH/Exome_Analysis/Results/Mapped_normals/", shell=True)
            if bn in StateBams_normal:
                continue

os.system('rm -r %s/mapping_normal'%global_tmp)
os.system('rm -r %s/mapping_tumor'%global_tmp)
os.system('rm -r %s/AddOrReplaceReadGroups_tmp'%global_tmp)

for item in Basename:
    os.system('qsub -q long.q %s/Indexing%s.sh'%(Job_Scripts, item))

print('''
--------------------!!Indexing!!--------------------
''')

StateIndexing = subprocess.check_output("qstat", shell=True)
while b'Indexing' in StateIndexing:
    sleep (1)
    StateIndexing = subprocess.check_output("qstat", shell=True)
    if b'Indexing' not in StateIndexing:
        continue    

for item in Basename:
    Bam_file_tumor='%s_tumor_ddup.bai'%item
    StateBams_tumor = subprocess.check_output("ls /data/projects/2020/OvarianCancerHH/Exome_Analysis/Results/Mapped_tumor/", shell=True)
    bti = bytes(Bam_file_tumor, 'utf-8')
    while bti not in StateBams_tumor:
            sleep(1)
            StateBams_tumor = subprocess.check_output("ls /data/projects/2020/OvarianCancerHH/Exome_Analysis/Results/Mapped_tumor/", shell=True)
            if bti in StateBams_tumor:
                continue

for item in Basename:
    Bam_file_normal='%s_normal_ddup.bai'%item
    StateBams_normal = subprocess.check_output("ls /data/projects/2020/OvarianCancerHH/Exome_Analysis/Results/Mapped_normals/", shell=True)
    bni = bytes(Bam_file_normal, 'utf-8')
    while bni not in StateBams_normal:
            sleep(1)
            StateBams_normal = subprocess.check_output("ls /data/projects/2020/OvarianCancerHH/Exome_Analysis/Results/Mapped_normals/", shell=True)
            if bni in StateBams_normal:
                continue

os.system('rm -r %s/mapping_tumor'%global_tmp)
os.system('rm -r %s/mapping_normal'%global_tmp)


#Variant calling
os.system('mkdir -p %s/Varinat_calling'%log)
os.system('mkdir -p %s/Varinat_calling'%error)
Varinat_calling_log=('%s/Varinat_calling'%log)
Varinat_calling_error=('%s/Varinat_calling'%error)

    #Germline Mutations with Haplotype caller
for item in Basename:
    with open ('%s/VC_Haplotype_caller%s.sh'%(Job_Scripts, item), 'w') as HaplotypeCaller_jobscript:
        HaplotypeCaller_jobscript.write('''
#!/bin/sh
#$ -S /bin/sh

#$ -pe smp 2
#$ -cwd
#$ -V

#### Jobdescription at qstat
#$ -N VC

#### Error Outputfile
#$ -e %s/%s_HaplotypeCaller.dat
#$ -o %s/%s_HaplotypeCaller.dat

#### Resubmit
#$ -r y

ulimit -n 2000        
mkdir -p %s/HaplotypeCaller_tmp%s;
java -jar -Djava.io.tmpdir=%s/HaplotypeCaller_tmp%s/ -XX:ConcGCThreads=2 -XX:ParallelGCThreads=2 /usr/local/bioinf/gatk/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar HaplotypeCaller -R %s -I %s/%s_tumor_ddup.bam --stand-call-conf 30.0 --base-O %s/%s_tumor_Happlotype_caller.vcf 
java -jar -Djava.io.tmpdir=%s/HaplotypeCaller_tmp%s/ -XX:ConcGCThreads=2 -XX:ParallelGCThreads=2 /usr/local/bioinf/gatk/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar HaplotypeCaller -R %s -I %s/%s_normal_ddup.bam --stand-call-conf 30.0 -O %s/%s_normal_Happlotype_caller.vcf;
cp %s/%s_tumor_Happlotype_caller.vcf %s/%s_tumor_Happlotype_caller.vcf
cp %s/%s_normal_Happlotype_caller.vcf %s/%s_normal_Happlotype_caller.vcf
'''%(Varinat_calling_error,item,
Varinat_calling_log,item,
local_tmp,item,
local_tmp,item,Path_to_GATK_reference,Path_to_Mapped_tumor,item,Path_to_HaplotypeCaller_variant_calls,item,
local_tmp,item,Path_to_GATK_reference,Path_to_Mapped_normals,item,Path_to_HaplotypeCaller_variant_calls,item,
Path_to_HaplotypeCaller_variant_calls,item,Path_to_all_vcfs_unfiltered,item,
Path_to_HaplotypeCaller_variant_calls,item,Path_to_all_vcfs_unfiltered,item))
HaplotypeCaller_jobscript.close

    #Somatic variants with Mutect2 
for item in Basename:
    with open ('%s/VC_Mutect2%s.sh'%(Job_Scripts, item), 'w') as Mutect2_jobscript:
        Mutect2_jobscript.write('''
#!/bin/sh
#$ -S /bin/sh

#$ -pe smp 1
#$ -cwd
#$ -V

#### Jobdescription at qstat
#$ -N VC

#### Error Outputfile
#$ -e %s/%s_Mutect2.dat
#$ -o %s/%s_Mutect2.dat

#### Resubmit
#$ -r y


java -jar /usr/local/bioinf/gatk/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar Mutect2 -R %s -I %s/%s_tumor_ddup.bam -I %s/%s_normal_ddup.bam -normal %s/%s_normal_ddup.bam -O %s/%s_somatic_Mutect2.vcf; 
cp %s/%s_somatic_Mutect2.vcf %s/%s_somatic_Mutect2.vcf
'''%(Varinat_calling_error,item,
Varinat_calling_log,item,
Path_to_GATK_reference,Path_to_Mapped_tumor,item,Path_to_Mapped_normals,item,Path_to_Mapped_normals,item,Path_to_Mutect2_variant_calls,item,
Path_to_Mutect2_variant_calls,item,Path_to_all_vcfs_unfiltered,item))
Mutect2_jobscript.close

    #Somatic variants with SomaticSniper
for item in Basename:
    with open ('%s/VC_SomaticSniper%s.sh'%(Job_Scripts, item), 'w') as SomaticSniper_jobscript:
        SomaticSniper_jobscript.write('''
#!/bin/sh
#$ -S /bin/sh

#$ -pe smp 1
#$ -cwd
#$ -V

#### Jobdescription at qstat
#$ -N VC

#### Error Outputfile
#$ -e %s/%s_SomaticSniper.dat
#$ -o %s/%s_SomaticSniper.dat

#### Resubmit
#$ -r y

ulimit -n 2000
bam-somaticsniper -Q 40 -G -L -F vcf -f %s %s/%s_tumor_ddup.bam %s/%s_normal_ddup.bam %s/%s_somatic_SomaticSniper.vcf;
cp %s/%s_somatic_SomaticSniper.vcf %s/%s_somatic_SomaticSniper.vcf
'''%(Varinat_calling_error,item,
Varinat_calling_log,item,
Path_to_fasta_reference,Path_to_Mapped_tumor,item,Path_to_Mapped_normals,item,Path_to_SomaticSniper_variant_calls,item,
Path_to_SomaticSniper_variant_calls,item,Path_to_all_vcfs_unfiltered,item
))
SomaticSniper_jobscript.close

    #Somatic and Germline Variant calling with Varsacan2 
for item in Basename:
    with open ('%s/VC_Varscan2%s.sh'%(Job_Scripts, item), 'w') as Varscan2_jobscript:
        Varscan2_jobscript.write('''
#!/bin/sh
#$ -S /bin/sh

#$ -pe smp 2
#$ -cwd
#$ -V

#### Jobdescription at qstat
#$ -N VC 

#### Error Outputfile
#$ -e %s/%s_Varsacan2.dat
#$ -o %s/%s_Varsacan2.dat

#### Resubmit
#$ -r y

ulimit -n 2000
mkdir -p %s/mpileups/%s;
samtools mpileup -f %s %s/%s_tumor_ddup.bam -o %s/mpileups/%s/%s_tumor.pileup
samtools mpileup -f %s %s/%s_normal_ddup.bam -o %s/mpileups/%s/%s_normal.pileup;
java -jar /usr/local/bioinf/varscan/varscan-2.4.3/VarScan.v2.4.3.jar somatic %s/mpileups/%s/%s_normal.pileup %s/mpileups/%s/%s_tumor.pileup %s/%s_Varscan2_somatic --output-vcf 1;
cp %s/%s_Varscan2_somatic.snp.vcf %s/%s_Varscan2_somatic.snp.vcf
cp %s/%s_Varscan2_somatic.indel.vcf %s/%s_Varscan2_somatic.indel.vcf
'''%(Varinat_calling_error,item,
Varinat_calling_log,item,
global_tmp,item,
Path_to_fasta_reference,Path_to_Mapped_tumor,item,global_tmp,item,item,
Path_to_fasta_reference,Path_to_Mapped_normals,item,global_tmp,item,item,
global_tmp,item,item, global_tmp,item,item,Path_to_Varscan2_variant_calls,item,
Path_to_Varscan2_variant_calls,item,Path_to_all_vcfs_unfiltered,item,
Path_to_Varscan2_variant_calls,item,Path_to_all_vcfs_unfiltered,item))
Varscan2_jobscript.close

    #Somatic and Germline Variant calling with Strelka
for item in Basename:
    with open ('%s/VC_Strelka2%s.sh'%(Job_Scripts, item), 'w') as Strelka2_jobscript:
        Strelka2_jobscript.write('''
#!/bin/sh
#$ -S /bin/sh

#$ -pe smp 3
#$ -cwd
#$ -V

#### Jobdescription at qstat
#$ -N VC 

#### Error Outputfile
#$ -e %s/%s_Strelka2.dat
#$ -o %s/%s_Strelka2.dat

#### Resubmit
#$ -r y


mkdir -p %s/Strelka2_normal/%s
mkdir -p %s/Strelka2_tumor/%s
mkdir -p %s/Strelka2_somatic/%s;
/usr/local/bioinf/strelka2/strelka-2.8.4/bin/configureStrelkaGermlineWorkflow.py --bam %s/%s_tumor_ddup.bam --referenceFasta %s --runDir %s/Strelka2_tumor/%s 
/usr/local/bioinf/strelka2/strelka-2.8.4/bin/configureStrelkaGermlineWorkflow.py  --bam %s/%s_normal_ddup.bam --referenceFasta %s --runDir %s/Strelka2_normal/%s 
/usr/local/bioinf/strelka2/strelka-2.8.4/bin/configureStrelkaSomaticWorkflow.py --normalBam %s/%s_normal_ddup.bam --tumorBam %s/%s_tumor_ddup.bam --referenceFasta %s --runDir %s/Strelka2_somatic/%s; 
%s/Strelka2_tumor/%s/runWorkflow.py -m local -j 4
%s/Strelka2_normal/%s/runWorkflow.py -m local -j 4
%s/Strelka2_somatic/%s/runWorkflow.py -m local -j 4;
cp %s/Strelka2_somatic/%s/results/variants/somatic.snvs.vcf.gz %s/%s_strelka2_somatic.vcf.gz
cp %s/Strelka2_somatic/%s/results/variants/somatic.indels.vcf.gz %s/%s_strelka2_somatic_indels.vcf.gz
'''%(Varinat_calling_error,item,
Varinat_calling_log,item,
Path_to_Strelka_dir,item,
Path_to_Strelka_dir,item,
Path_to_Strelka_dir,item,
Path_to_Mapped_tumor,item,Path_to_fasta_reference,Path_to_Strelka_dir,item,
Path_to_Mapped_normals,item,Path_to_fasta_reference,Path_to_Strelka_dir,item,
Path_to_Mapped_normals,item,Path_to_Mapped_tumor,item,Path_to_fasta_reference,Path_to_Strelka_dir,item,
Path_to_Strelka_dir,item,
Path_to_Strelka_dir,item,
Path_to_Strelka_dir,item,
Path_to_Strelka_dir,item,Path_to_all_vcfs_unfiltered,item,
Path_to_Strelka_dir,item,Path_to_all_vcfs_unfiltered,item))
Strelka2_jobscript.close

for item in Basename:
    os.system('qsub -q long.q %s/VC_Haplotype_caller%s.sh'%(Job_Scripts, item))
    os.system('qsub -q long.q %s/VC_Mutect2%s.sh'%(Job_Scripts, item))
    os.system('qsub -q long.q %s/VC_SomaticSniper%s.sh'%(Job_Scripts, item))
    os.system('qsub -q long.q %s/VC_Varscan2%s.sh'%(Job_Scripts, item))
    os.system('qsub -q long.q %s/VC_Strelka2%s.sh'%(Job_Scripts, item))

print('''
-----------------!!Variant Calling!!----------------
''')


StateVC = subprocess.check_output("qstat", shell=True)
while b'VC' in StateVC:
    sleep (1)
    StateVC = subprocess.check_output("qstat", shell=True)
    if b'VC' not in StateVC:
        continue

os.system('gunzip -f %s/*.gz'%Path_to_all_vcfs_unfiltered)

#VariantFiltration 
os.system('mkdir -p %s/Varinat_Filtering'%log)
os.system('mkdir -p %s/Varinat_Filtering'%error)
Varinat_Filtering_log=('%s/Varinat_Filtering'%log)
Varinat_Filtering_error=('%s/Varinat_Filtering'%error)

for item in Basename:
    with open('%s/%s_Mutect2_fitration_script.sh'%(Job_Scripts, item),'w') as Mutect2_filtration_Jobscript:
        Mutect2_filtration_Jobscript.write('''
        #!/bin/sh
#$ -S /bin/sh

#$ -pe smp 1
#$ -cwd
#$ -V

#### Jobdescription at qstat
#$ -N FMut2 

#### Error Outputfile
#$ -e %s/%s_FilterMutectCalls.dat
#$ -o %s/%s_FilterMutectCalls.dat

#### Resubmit
#$ -r y

java -jar /usr/local/bioinf/gatk/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar FilterMutectCalls -R %s -V %s/%s_somatic_Mutect2.vcf -O %s/%s_FilterMutectCalls_somatic_Mutect2.vcf;
bcftools norm -m - %s/%s_FilterMutectCalls_somatic_Mutect2.vcf > %s/%s_biallelic_somatic_Mutect2.vcf
'''%(Varinat_Filtering_error,item,
Varinat_Filtering_log,item,
Path_to_GATK_reference,Path_to_Mutect2_variant_calls,item,Variant_filtering_tmp,item,
Variant_filtering_tmp,item,Variant_filtering_tmp,item
))

for item in Basename:
    os.system('qsub %s/%s_Mutect2_fitration_script.sh'%(Job_Scripts, item))

stateFMut2 = subprocess.check_output(' qstat', shell=True)
while b'FMut2' in stateFMut2:
    sleep (1)
    stateFMut2 = subprocess.check_output(' qstat', shell=True)
    if b'FMut2' not in stateFMut2:
        continue

for item in Basename:
    os.system('''grep '#' %s/%s_biallelic_somatic_Mutect2.vcf > %s/%s_Filtered_somatic_Mutect2.vcf'''%(Variant_filtering_tmp,item,Variant_filtering_tmp,item))
    os.system('''grep 'PASS' %s/%s_biallelic_somatic_Mutect2.vcf >> %s/%s_Filtered_somatic_Mutect2.vcf'''%(Variant_filtering_tmp,item,Variant_filtering_tmp,item))

for item in Basename:
    with open('%s/%s_Variant_fitration_script.sh'%(Job_Scripts, item),'w') as Variant_filtration_Jobscript:
        Variant_filtration_Jobscript.write('''
#!/bin/sh
#$ -S /bin/sh

#$ -pe smp 2
#$ -cwd
#$ -V

#### Jobdescription at qstat
#$ -N VF 

#### Error Outputfile
#$ -e %s/%s_FilterVcf.dat
#$ -o %s/%s_FilterVcf.dat

#### Resubmit
#$ -r y

mkdir -p %s/FilterVcf%s;
java -jar -Djava.io.tmpdir=%s/FilterVcf%s/ -XX:ConcGCThreads=2 -XX:ParallelGCThreads=2 /usr/local/bioinf/picard/picard-tools-2.21.4/picard.jar FilterVcf INPUT=%s/%s_tumor_Happlotype_caller.vcf OUTPUT=%s/%s_FilterVcf_tumor_Happlotype_caller.vcf MIN_DP=10
java -jar -Djava.io.tmpdir=%s/FilterVcf%s/ -XX:ConcGCThreads=2 -XX:ParallelGCThreads=2 /usr/local/bioinf/picard/picard-tools-2.21.4/picard.jar FilterVcf INPUT=%s/%s_normal_Happlotype_caller.vcf OUTPUT=%s/%s_FilterVcf_normal_Happlotype_caller.vcf MIN_DP=10
java -jar -Djava.io.tmpdir=%s/FilterVcf%s/ -XX:ConcGCThreads=2 -XX:ParallelGCThreads=2 /usr/local/bioinf/picard/picard-tools-2.21.4/picard.jar FilterVcf INPUT=%s/%s_Filtered_somatic_Mutect2.vcf OUTPUT=%s/%s_FilterVcf_somatic_Mutect2.vcf MIN_DP=10
java -jar -Djava.io.tmpdir=%s/FilterVcf%s/ -XX:ConcGCThreads=2 -XX:ParallelGCThreads=2 /usr/local/bioinf/picard/picard-tools-2.21.4/picard.jar FilterVcf INPUT=%s/%s_somatic_SomaticSniper.vcf OUTPUT=%s/%s_FilterVcf_somatic_SomaticSniper.vcf MIN_DP=10
java -jar -Djava.io.tmpdir=%s/FilterVcf%s/ -XX:ConcGCThreads=2 -XX:ParallelGCThreads=2 /usr/local/bioinf/picard/picard-tools-2.21.4/picard.jar FilterVcf INPUT=%s/%s_strelka2_tumor.vcf OUTPUT=%s/%sFilterVcf_strelka2_tumor.vcf MIN_DP=10
java -jar -Djava.io.tmpdir=%s/FilterVcf%s/ -XX:ConcGCThreads=2 -XX:ParallelGCThreads=2 /usr/local/bioinf/picard/picard-tools-2.21.4/picard.jar FilterVcf INPUT=%s/%s_strelka2_normal.vcf OUTPUT=%s/%s_FilterVcf_strelka2_normal.vcf MIN_DP=10
java -jar -Djava.io.tmpdir=%s/FilterVcf%s/ -XX:ConcGCThreads=2 -XX:ParallelGCThreads=2 /usr/local/bioinf/picard/picard-tools-2.21.4/picard.jar FilterVcf INPUT=%s/%s_strelka2_somatic.vcf OUTPUT=%s/%s_FilterVcf_strelka2_somatic.vcf MIN_DP=10

'''%(Varinat_Filtering_error,item,
Varinat_Filtering_log,item,
local_tmp,item,
local_tmp,item,Path_to_all_vcfs_unfiltered,item,Variant_filtering_tmp,item,
local_tmp,item,Path_to_all_vcfs_unfiltered,item,Variant_filtering_tmp,item,
local_tmp,item,Variant_filtering_tmp,item,Variant_filtering_tmp,item,
local_tmp,item,Path_to_all_vcfs_unfiltered,item,Variant_filtering_tmp,item,
local_tmp,item,Path_to_all_vcfs_unfiltered,item,Variant_filtering_tmp,item,
local_tmp,item,Path_to_all_vcfs_unfiltered,item,Variant_filtering_tmp,item,
local_tmp,item,Path_to_all_vcfs_unfiltered,item,Variant_filtering_tmp,item,
))

for item in Basename:
    os.system('qsub %s/%s_Variant_fitration_script.sh'%(Job_Scripts, item))

print('''
-----------------Fitering Variants------------------
''')

StateVF = subprocess.check_output("qstat", shell=True)
while b'VF' in StateVF:
    sleep (1)
    StateVF = subprocess.check_output("qstat", shell=True)
    if b'VF' not in StateVF:
        continue

#Custom_filter_and_format_Strelka_snvs
for item in Basename:
    exclude = [i for i, line in enumerate(open('%s/%s_FilterVcf_strelka2_somatic.vcf'%(Variant_filtering_tmp,item))) if line.startswith('#')]
    df=pd.read_csv('%s/%s_FilterVcf_strelka2_somatic.vcf'%(Variant_filtering_tmp,item), sep=':|\t', engine='python', skiprows=exclude, dtype=str, names=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMA1','FORMA2','FORMAT3','FORMAT4','FORMAT5','FORMAT6','FORMAT7','FORMAT8','FORMAT9','normalAU','normalCU','normalDP','normalFDP','normalFT','normalGU','normalSDP','normalSUBDP','normalTU','tumorAU','tumorCU','tumorDP','tumorFDP','tumorFT','tumorGU','tumorSDP','tumorSUBDP','tumorTU'])
    df = df.stack().str.replace(',','.').unstack()
    

    df['normalAU'] = pd.to_numeric(df['normalAU'])
    df['normalTU'] = pd.to_numeric(df['normalTU'])
    df['normalCU'] = pd.to_numeric(df['normalCU'])
    df['normalGU'] = pd.to_numeric(df['normalGU'])
    df['tumorAU'] = pd.to_numeric(df['tumorAU'])
    df['tumorTU'] = pd.to_numeric(df['tumorTU'])
    df['tumorCU'] = pd.to_numeric(df['tumorCU'])
    df['tumorGU'] = pd.to_numeric(df['tumorGU'])

    def refCounts_normal(row):
        if row['REF'] == 'A':
            val = row['normalAU']
        elif row['REF'] == 'T':
            val = row['normalTU']
        elif row['REF'] == 'G':
            val = row['normalGU']
        elif row['REF'] == 'C':
            val = row['normalCU']
        else:
            val = 'Na'
        return val

    df['refCounts_normal']=df.apply(refCounts_normal,axis=1)

    def altCounts_normal(row):
        if row['ALT'] == 'A':
            val = row['normalAU']
        elif row['ALT'] == 'T':
            val = row['normalTU']
        elif row['ALT'] == 'G':
            val = row['normalGU']
        elif row['ALT'] == 'C':
            val = row['normalCU']
        else:
            val = 'Na'
        return val

    df['altCounts_normal']=df.apply(altCounts_normal,axis=1)

    df['allelFreq_normal']=df['altCounts_normal']/(df['altCounts_normal']+df['refCounts_normal'])

    def refCounts_tumor(row):
        if row['REF'] == 'A':
            val = row['tumorAU']
        elif row['REF'] == 'T':
            val = row['tumorTU']
        elif row['REF'] == 'G':
            val = row['tumorGU']
        elif row['REF'] == 'C':
            val = row['tumorCU']
        else:
            pass
        return val

    df['refCounts_tumor']=df.apply(refCounts_tumor,axis=1)

    def altCounts_tumor(row):
        if row['ALT'] == 'A':
            val = row['tumorAU']
        elif row['ALT'] == 'T':
            val = row['tumorTU']
        elif row['ALT'] == 'G':
            val = row['tumorGU']
        elif row['ALT'] == 'C':
            val = row['tumorCU']
        else:
            pass
        return val

    df['altCounts_tumor']=df.apply(altCounts_tumor,axis=1)

    df['allelFreq_tumor']=df['altCounts_tumor']/(df['altCounts_tumor']+df['refCounts_tumor'])

    def somatic_strelka(row):
        if row['allelFreq_tumor'] >= 0.05 and row['allelFreq_normal'] < 0.05:
            val = 'yes'
        else:
            val= 'no'
        return val

    df['Somatic']=df.apply(somatic_strelka, axis=1)

    df_out=df.drop(columns=['QUAL','FILTER','INFO','FORMA1','FORMA2','FORMAT3','FORMAT4','FORMAT5','FORMAT6','FORMAT7','FORMAT8','FORMAT9','normalAU','normalCU','normalDP','normalFDP','normalFT','normalGU','normalSDP','normalSUBDP','normalTU','tumorAU','tumorCU','tumorDP','tumorFDP','tumorFT','tumorGU','tumorSDP','tumorSUBDP','tumorTU','refCounts_normal','altCounts_normal','refCounts_tumor','altCounts_tumor'])
    df_out.to_csv('%s/%sStrelka2_somatic_VCF_with_allelFreq.csv'%(Variant_filtering_tmp, item))
    exclude = [i for i, line in enumerate(open('%s/%sStrelka2_somatic_VCF_with_allelFreq.csv'%(Variant_filtering_tmp,item))) if 'no' in line]
    df2=pd.read_csv('%s/%sStrelka2_somatic_VCF_with_allelFreq.csv'%(Variant_filtering_tmp, item), sep=',', engine='python', skiprows=exclude, names=['1','#CHROM','POS','ID','REF','ALT','ALTfreqNormal','ALTfreqTumor','SOMATIC'])
    df2=df2.drop(columns=['1'])
    df2.to_csv('%s/%s_all_somatic_calls_Strelka.csv'%(Variant_filtering_tmp,item))

#Custom_filter_and_format_Strelka_indels

for item in Basename:
    exclude = [i for i, line in enumerate(open('%s/%s_strelka2_somatic_indels.vcf'%(Path_to_all_vcfs_unfiltered,item,))) if line.startswith('#')]
    df=pd.read_csv('%s/%s_strelka2_somatic_indels.vcf'%(Path_to_all_vcfs_unfiltered,item,), sep=':|\t', engine='python', skiprows=exclude, dtype=str, names=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT_DP','FORMAT_DP2','FORMAT_TAR','FORMAT_TIR','FORMAT_TOR','FORMAT_DP50','FORMATFDP_50','FORMAT_SUBDP50','FORMAT_BCN50','normal_DP','normal_DP2','normal_TAR','normal_TIR','normal_TOR','normal_DP50','normal_FDP50','normal_SUBDP50','normal_BCN50','tumor_DP','tumor_DP2','tumor_TAR','tumor_TIR','tumor_TOR','tumor_DP50','tumor_FDP50','tumor_SUBDP50','tumor_BCN50'])
    df = df.stack().str.replace(',','.').unstack()

    df['normal_TIR'] = pd.to_numeric(df['normal_TIR'])
    df['normal_TAR'] = pd.to_numeric(df['normal_TAR'])
    df['tumor_TIR'] = pd.to_numeric(df['tumor_TIR'])
    df['tumor_TAR'] = pd.to_numeric(df['tumor_TAR'])
    df['normal_DP'] = pd.to_numeric(df['normal_DP'])
    df['tumor_DP'] = pd.to_numeric(df['tumor_DP'])

    df['allelFreq_normal']=df['normal_TIR']/(df['normal_TIR']+df['normal_TAR'])
    df['allelFreq_tumor']=df['tumor_TIR']/(df['tumor_TIR']+df['tumor_TAR'])

    def somatic_strelka_indels(row):
        if row['allelFreq_tumor'] >= 0.05 and row['allelFreq_normal'] < 0.05:
            val = 'yes'
        else:
            val= 'no'
        return val

    df['Somatic']=df.apply(somatic_strelka, axis=1)

    def sequencing_Depth_StrelkaIndel(row):
        if row['normal_DP'] >=10 and row['tumor_DP']>=10:
            val = 'Pass'
        else:
            val = 'low'
        return val

    df['Depth']=df.apply(sequencing_Depth_StrelkaIndel, axis=1)

    df_out=df.drop(columns=['QUAL','FILTER','INFO','FORMAT_DP','FORMAT_DP2','FORMAT_TAR','FORMAT_TIR','FORMAT_TOR','FORMAT_DP50','FORMATFDP_50','FORMAT_SUBDP50','FORMAT_BCN50','normal_DP','normal_DP2','normal_TAR','normal_TIR','normal_TOR','normal_DP50','normal_FDP50','normal_SUBDP50','normal_BCN50','tumor_DP','tumor_DP2','tumor_TAR','tumor_TIR','tumor_TOR','tumor_DP50','tumor_FDP50','tumor_SUBDP50','tumor_BCN50'])
    df_out.to_csv('%s/%sStrelka2_somatic_indel_VCF_with_allelFreq.csv'%(Variant_filtering_tmp, item))
    exclude = [i for i, line in enumerate(open('%s/%sStrelka2_somatic_indel_VCF_with_allelFreq.csv'%(Variant_filtering_tmp,item))) if 'no' in line]
    df2=pd.read_csv('%s/%sStrelka2_somatic_indel_VCF_with_allelFreq.csv'%(Variant_filtering_tmp, item), sep=',', engine='python', skiprows=exclude, names=['1','#CHROM','POS','ID','REF','ALT','ALTfreqNormal','ALTfreqTumor','SOMATIC','Depth'])
    df2=df2.drop(columns=['1'])
    df2.to_csv('%s/%s_Strelka_indels_DP_Filter.csv'%(Variant_filtering_tmp,item))
    exclude = [i for i, line in enumerate(open('%s/%s_Strelka_indels_DP_Filter.csv'%(Variant_filtering_tmp,item))) if 'low' in line]
    df3=pd.read_csv('%s/%s_Strelka_indels_DP_Filter.csv'%(Variant_filtering_tmp, item), sep=',', engine='python', skiprows=exclude)
    df3=df3.drop(columns=['Unnamed: 0','Depth'])
    print(df3)
    df3.to_csv('%s/%s_all_somatic_calls_Strelka_indels.csv'%(Variant_filtering_tmp,item))

    os.system('''grep -v '#CHROM' %s/%s_all_somatic_calls_Strelka_indels.csv >>%s/%s_all_somatic_calls_Strelka.csv'''%(Variant_filtering_tmp,item,Variant_filtering_tmp,item))


#Custom_filter_and_format_Mutect2

for item in Basename:
    os.system('''grep 'PGT' %s/%s_FilterVcf_somatic_Mutect2.vcf > %s/%s_filtered_somatic_Mutect2_PGT.vcf'''%(Variant_filtering_tmp,item,Variant_filtering_tmp,item))
    os.system('''grep -v 'PGT' %s/%s_FilterVcf_somatic_Mutect2.vcf > %s/%s_filtered_somatic_Mutect2_nonPGT.vcf'''%(Variant_filtering_tmp,item,Variant_filtering_tmp,item))

    #Non PGT
    exclude = [i for i, line in enumerate(open('%s/%s_filtered_somatic_Mutect2_nonPGT.vcf'%(Variant_filtering_tmp,item))) if line.startswith('#')]
    df=pd.read_csv('%s/%s_filtered_somatic_Mutect2_nonPGT.vcf'%(Variant_filtering_tmp,item), sep=':|\t', engine='python', skiprows=exclude, dtype=str, names=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','GTformat','ADformat','AFformat','DPformat','F1R2format','F2R1format','FTformat','SBformat','GTnormal','ADnormal','AFnormal','DPnormal','F1R2normal','F2R1normal','FTnormal','SBnormal','GTtumor','ADtumor','AFtumor','DPtumor','F1R2tumor','F2R1tumor','FTtumor','SBtumor'])
  
    df['AFtumor'] = pd.to_numeric(df['AFtumor'])
    df['AFnormal'] = pd.to_numeric(df['AFnormal'])

    def somatic_mutect2(row):
        if row['AFtumor'] >= 0.05 and row['AFnormal'] < 0.05:
            val = 'yes'
        else:
            val= 'no'
        return val

    df['Somatic']=df.apply(somatic_mutect2, axis=1)

    df_out=df.drop(columns=['QUAL','FILTER','INFO','GTformat','ADformat','AFformat','DPformat','F1R2format','F2R1format','FTformat','SBformat','GTnormal','ADnormal','DPnormal','F1R2normal','F2R1normal','FTnormal','SBnormal','GTtumor','ADtumor','DPtumor','F1R2tumor','F2R1tumor','FTtumor','SBtumor'])
    df_out.to_csv('%s/%ssomatic_Mutect2_VCF_with_allelFreq.csv'%(Variant_filtering_tmp,item))
    exclude = [i for i, line in enumerate(open('%s/%ssomatic_Mutect2_VCF_with_allelFreq.csv'%(Variant_filtering_tmp, item))) if 'no' in line]
    df2=pd.read_csv('%s/%ssomatic_Mutect2_VCF_with_allelFreq.csv'%(Variant_filtering_tmp, item), sep=',', engine='python', skiprows=exclude, names=['1','#CHROM','POS','ID','REF','ALT','ALTfreqNormal','ALTfreqTumor','SOMATIC'])
    df2=df2.drop(columns=['1'])
    df2.to_csv('%s/%s_all_somatic_calls_Mutect2.csv'%(Variant_filtering_tmp,item))
   
    #PGT
    exclude = [i for i, line in enumerate(open('%s/%s_filtered_somatic_Mutect2_PGT.vcf'%(Variant_filtering_tmp,item))) if line.startswith('#')]
    df=pd.read_csv('%s/%s_filtered_somatic_Mutect2_PGT.vcf'%(Variant_filtering_tmp,item), sep=':|\t', engine='python', skiprows=exclude, dtype=str, names=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','GTformat','ADformat','AFformat','DPformat','F1R2format','F2R1format','FTformat','PGTformat','PIDformat','PSformat','SBformat','GTnormal','ADnormal','AFnormal','DPnormal','F1R2normal','F2R1normal','FTnormal','PGTnormal','PIDnormal','PSnormal','SBnormal','GTtumor','ADtumor','AFtumor','DPtumor','F1R2tumor','F2R1tumor','FTtumor','PGTtumor','PIDtumor','PStumor','SBtumor'])
  
    df['AFtumor'] = pd.to_numeric(df['AFtumor'])
    df['AFnormal'] = pd.to_numeric(df['AFnormal'])

    def somatic_mutect2_pgt(row):
        if row['AFtumor'] >= 0.05 and row['AFnormal'] < 0.05:
            val = 'yes'
        else:
            val= 'no'
        return val

    df['Somatic']=df.apply(somatic_mutect2_pgt, axis=1)

    df_out=df.drop(columns=['QUAL','FILTER','INFO','GTformat','ADformat','AFformat','DPformat','F1R2format','F2R1format','FTformat','PGTformat','PIDformat','PSformat','SBformat','GTnormal','ADnormal','DPnormal','F1R2normal','F2R1normal','FTnormal','PGTnormal','PIDnormal','PSnormal','SBnormal','GTtumor','ADtumor','DPtumor','F1R2tumor','F2R1tumor','FTtumor','PGTtumor','PIDtumor','PStumor','SBtumor'])
    df_out.to_csv('%s/%ssomatic_Mutect2_VCF_with_allelFreq_PGT.csv'%(Variant_filtering_tmp,item))
    exclude = [i for i, line in enumerate(open('%s/%ssomatic_Mutect2_VCF_with_allelFreq_PGT.csv'%(Variant_filtering_tmp, item))) if 'no' in line]
    df2=pd.read_csv('%s/%ssomatic_Mutect2_VCF_with_allelFreq_PGT.csv'%(Variant_filtering_tmp, item), sep=',', engine='python', skiprows=exclude, names=['1','#CHROM','POS','ID','REF','ALT','ALTfreqNormal','ALTfreqTumor','SOMATIC'])
    df2=df2.drop(columns=['1'])
    df2.to_csv('%s/%s_all_somatic_calls_Mutect2_PGT.csv'%(Variant_filtering_tmp,item))

    os.system('''grep 'yes' %s/%s_all_somatic_calls_Mutect2_PGT.csv >> %s/%s_all_somatic_calls_Mutect2.csv'''%(Variant_filtering_tmp,item,Variant_filtering_tmp,item))

    
    
#Varscan2 somatic Filter

for item in Basename:
    os.system("sed 's/%//g' {}/{}_Varscan2_somatic.snp.vcf > {}/{}_Varscan2_somatic_forFiltering.snp.vcf".format(Path_to_Varscan2_variant_calls,item,Variant_filtering_tmp,item))
    os.system("sed 's/%//g' {}/{}_Varscan2_somatic.indel.vcf > {}/{}_Varscan2_somatic_forFiltering.indel.vcf".format(Path_to_Varscan2_variant_calls,item,Variant_filtering_tmp,item))
    exclude = [i for i, line in enumerate(open('%s/%s_Varscan2_somatic_forFiltering.snp.vcf'%(Variant_filtering_tmp,item))) if line.startswith('#')]
    exclude_indel = [i for i, line in enumerate(open('%s/%s_Varscan2_somatic_forFiltering.snp.vcf'%(Variant_filtering_tmp,item))) if line.startswith('#')]
    df=pd.read_csv('%s/%s_Varscan2_somatic_forFiltering.snp.vcf'%(Variant_filtering_tmp,item), sep=':|\t', engine='python', skiprows=exclude, dtype=str, names=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','GTformat','GQformat','DPformat','RDformat','ADformat','FREQformat','DP4format','GTnormal','GQnormal','DPnormal','RDnormal','ADnormal','FREQnormal','DP4normal','GTtumor','GQtumor','DPtumor','RDtumor','ADtumor','FREQtumor','DP4tumor'])
    df_indel=pd.read_csv('%s/%s_Varscan2_somatic_forFiltering.indel.vcf'%(Variant_filtering_tmp,item), sep=':|\t', engine='python', skiprows=exclude_indel, dtype=str, names=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','GTformat','GQformat','DPformat','RDformat','ADformat','FREQformat','DP4format','GTnormal','GQnormal','DPnormal','RDnormal','ADnormal','FREQnormal','DP4normal','GTtumor','GQtumor','DPtumor','RDtumor','ADtumor','FREQtumor','DP4tumor'])

    df['FREQnormal'] = pd.to_numeric(df['FREQnormal'])
    df['FREQtumor'] = pd.to_numeric(df['FREQtumor'])
    df_indel['FREQnormal'] = pd.to_numeric(df_indel['FREQnormal'])
    df_indel['FREQtumor'] = pd.to_numeric(df_indel['FREQtumor'])
    df['DPnormal'] = pd.to_numeric(df['DPnormal'])
    df_indel['DPnormal'] = pd.to_numeric(df_indel['DPnormal'])
    df['DPtumor'] = pd.to_numeric(df['DPtumor'])
    df_indel['DPtumor'] = pd.to_numeric(df_indel['DPtumor'])

    def somatic_varscan(row):
        if row['FREQtumor'] >= 5 and row['FREQnormal'] < 5:
            val = 'yes'
        else:
            val = 'no'
        return val
    
    def sequencing_Depth(row):
        if row['DPnormal'] >=10 and row['DPtumor']>=10:
            val = 'Pass'
        else:
            val = 'low'
        return val

    df['Somatic']=df.apply(somatic_varscan, axis=1)
    df['Depth']=df.apply(sequencing_Depth, axis=1)
    df_indel['Somatic']=df_indel.apply(somatic_varscan, axis=1)
    df_indel['Depth']=df_indel.apply(sequencing_Depth, axis=1)
    df_out=df.drop(columns=['QUAL','FILTER','INFO','GTformat','GQformat','DPformat','RDformat','ADformat','FREQformat','DP4format','GTnormal','GQnormal','DPnormal','RDnormal','ADnormal','DP4normal','GTtumor','GQtumor','DPtumor','RDtumor','ADtumor','DP4tumor']) 
    df_out.to_csv('%s/%ssomatic_Varscan2_VCF_with_allelFreq.csv'%(Variant_filtering_tmp,item))
    exclude = [i for i, line in enumerate(open('%s/%ssomatic_Varscan2_VCF_with_allelFreq.csv'%(Variant_filtering_tmp, item))) if 'no' in line] 
    df2=pd.read_csv('%s/%ssomatic_Varscan2_VCF_with_allelFreq.csv'%(Variant_filtering_tmp, item), sep=',', engine='python', skiprows=exclude, names=['#CHROM','POS','ID','REF','ALT','ALTfreqNormal','ALTfreqTumor','SOMATIC','Depth'])
    df2.to_csv('%s/%s_all_somatic_calls_forDPfilter_Varscan2.csv'%(Variant_filtering_tmp,item))
    exclude = [i for i, line in enumerate(open('%s/%s_all_somatic_calls_forDPfilter_Varscan2.csv'%(Variant_filtering_tmp, item))) if 'low' in line] 
    df3=pd.read_csv('%s/%s_all_somatic_calls_forDPfilter_Varscan2.csv'%(Variant_filtering_tmp, item), sep=',', engine='python', skiprows=exclude)
    df4=df3.drop(columns=['Unnamed: 0','Depth'])

    df_indel_out=df_indel.drop(columns=['QUAL','FILTER','INFO','GTformat','GQformat','DPformat','RDformat','ADformat','FREQformat','DP4format','GTnormal','GQnormal','DPnormal','RDnormal','ADnormal','DP4normal','GTtumor','GQtumor','DPtumor','RDtumor','ADtumor','DP4tumor']) 
    df_indel_out.to_csv('%s/%ssomatic_Varscan2_indel_with_allelFreq.csv'%(Variant_filtering_tmp,item))
    exclude_indel = [i for i, line in enumerate(open('%s/%ssomatic_Varscan2_indel_with_allelFreq.csv'%(Variant_filtering_tmp, item))) if 'no' in line] 
    df_indel2=pd.read_csv('%s/%ssomatic_Varscan2_indel_with_allelFreq.csv'%(Variant_filtering_tmp, item), sep=',', engine='python', skiprows=exclude_indel,names=['#CHROM','POS','ID','REF','ALT','ALTfreqNormal','ALTfreqTumor','SOMATIC','Depth'])
    df_indel2.to_csv('%s/%s_all_somatic_indelcalls_forDPfilter_Varscan2.csv'%(Variant_filtering_tmp,item))
    exclude = [i for i, line in enumerate(open('%s/%s_all_somatic_indelcalls_forDPfilter_Varscan2.csv'%(Variant_filtering_tmp, item))) if 'low' in line] 
    df_indel3=pd.read_csv('%s/%s_all_somatic_indelcalls_forDPfilter_Varscan2.csv'%(Variant_filtering_tmp, item), sep=',', engine='python', skiprows=exclude)
    #print(df_indel3)
    df_indel4=df_indel3.drop(columns=['Unnamed: 0','Depth'])
    #print(df_indel4)

    df_all_varscan_calls=df4.merge(df_indel4, how='outer')
    df_all_varscan_calls.to_csv('%s/%s_all_somatic_calls_Varscan2.csv'%(Variant_filtering_tmp,item))


#SomaticSniper somatic Filter

for item in Basename:
    os.system('bcftools norm -m - %s/%s_somatic_SomaticSniper.vcf > %s/%s_biallelic_somatic_SomaticSniper.vcf'%(Path_to_all_vcfs_unfiltered ,item,Variant_filtering_tmp,item))
    exclude = [i for i, line in enumerate(open('%s/%s_biallelic_somatic_SomaticSniper.vcf'%(Variant_filtering_tmp,item))) if line.startswith('#')]
    df=pd.read_csv('%s/%s_biallelic_somatic_SomaticSniper.vcf'%(Variant_filtering_tmp,item), sep=':|\t', engine='python', skiprows=exclude, dtype=str, names=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','GT','IGT','DP','DP4','BCOUNT','GQ','JGQ','VAQ','BQ','MQ','AMQ','SS','SSC','GTnormal','IGTnormal','DPnormal','DP4normal','BCOUNTnormal','GQnormal','JGQnormal','VAQnormal','BQnormal','MQnormal','AMQnormal','SSnormal','SSCnormal','GTtumor','IGTtumor','DPtumor','DP4tumor','BCOUNTtumor','GQtumor','JGQtumor','VAQtumor','BQtumor','MQtumor','AMQtumor','SStumor','SSCtumor'])
    df.to_csv('%s/%s_somaticSniper_test.csv'%(Variant_filtering_tmp,item))
    df=df.drop(columns=['QUAL','FILTER','INFO','GT','IGT','DP','DP4','BCOUNT','GQ','JGQ','VAQ','BQ','MQ','AMQ','SS','SSC','GTnormal','IGTnormal','DP4normal','GQnormal','JGQnormal','VAQnormal','BQnormal','MQnormal','AMQnormal','SSnormal','SSCnormal','GTtumor','IGTtumor','DP4tumor','GQtumor','JGQtumor','VAQtumor','BQtumor','MQtumor','AMQtumor','SStumor','SSCtumor'])
    df.to_csv('%s/%s_somaticSniper_filter.csv'%(Variant_filtering_tmp,item))
    os.system('''sed -i 's/"//g' {}/{}_somaticSniper_filter.csv'''.format(Variant_filtering_tmp,item))
    df=pd.read_csv('%s/%s_somaticSniper_filter.csv'%(Variant_filtering_tmp,item), sep=',',skiprows=1,  names=['1','CHROM','POS','ID','REF','ALT','DPnormal','normalCountsA','normalCountsC','normalCountsG','normalCountsT','DPtumor','tumorCountsA','tumorCountsC','tumorCountsG','tumorCountsT'])
    df.to_csv('%s/%s_somaticSniper_filter.csv'%(Variant_filtering_tmp,item))
    df=pd.read_csv('%s/%s_somaticSniper_filter.csv'%(Variant_filtering_tmp,item), sep=',')

    def ALT_baseCount_tumor(row):
        if row['ALT'] == 'A':
            val = row['tumorCountsA']
        elif row['ALT'] == 'T':
            val = row['tumorCountsT']
        elif row['ALT'] == 'G':
            val = row['tumorCountsG']
        elif row['ALT'] == 'C':
            val = row['tumorCountsC']
        else:
            val = 'Na'
        return val

    def ALT_baseCount_normal (row):
        if row['ALT'] == 'A':
            val = row['normalCountsA']
        elif row['ALT'] == 'T':
            val = row['normalCountsT']
        elif row['ALT'] == 'G':
            val = row['normalCountsG']
        elif row['ALT'] == 'C':
            val = row['normalCountsC']
        else:
            val = 'Na'
        return val

    df['ALT_baseCount_normal']=df.apply(ALT_baseCount_normal,axis=1)
    df['ALT_baseCount_tumor']=df.apply(ALT_baseCount_tumor,axis=1)

    df['ALT_baseCount_normal'] = pd.to_numeric(df['ALT_baseCount_normal'])
    df['ALT_baseCount_tumor'] = pd.to_numeric(df['ALT_baseCount_tumor'])
    df['DPnormal'] = pd.to_numeric(df['DPnormal'])
    df['DPtumor'] = pd.to_numeric(df['DPtumor'])

    df['Allel_Freq_normal']=(100*df['ALT_baseCount_normal'])/df['DPnormal']
    df['Allel_Freq_tumor']=(100*df['ALT_baseCount_tumor'])/df['DPtumor']
    df=df.drop(columns=['Unnamed: 0','1'])

    def somatic_SomaticSniper(row):
        if row['Allel_Freq_tumor'] >= 5 and row['Allel_Freq_normal'] < 5:
            val = 'yes'
        else:
            val = 'no'
        return val
    
    def sequencing_Depth_SomaticSniper(row):
        if row['DPnormal'] >=10 and row['DPtumor']>=10:
            val = 'Pass'
        else:
            val = 'low'
        return val

    df['Somatic']=df.apply(somatic_SomaticSniper, axis=1)
    df['Depth']=df.apply(sequencing_Depth_SomaticSniper, axis=1)
    
    df_out=df.drop(columns=['DPnormal','normalCountsA','normalCountsC','normalCountsG','normalCountsT','DPtumor','tumorCountsA','tumorCountsC','tumorCountsG','tumorCountsT','ALT_baseCount_normal','ALT_baseCount_tumor']) 
    df_out.to_csv('%s/%ssomatic_somaticSniper_VCF_with_allelFreq.csv'%(Variant_filtering_tmp,item))
    exclude = [i for i, line in enumerate(open('%s/%ssomatic_somaticSniper_VCF_with_allelFreq.csv'%(Variant_filtering_tmp, item))) if 'no' in line] 
    df2=pd.read_csv('%s/%ssomatic_somaticSniper_VCF_with_allelFreq.csv'%(Variant_filtering_tmp, item), sep=',', engine='python', skiprows=exclude, names=['#CHROM','POS','ID','REF','ALT','allel_freqNormal','allel_freqTumor','SOMATIC','Depth'])
    df2.to_csv('%s/%s_all_somatic_calls_forDPfilter_somaticSniper.csv'%(Variant_filtering_tmp,item))
    exclude = [i for i, line in enumerate(open('%s/%s_all_somatic_calls_forDPfilter_somaticSniper.csv'%(Variant_filtering_tmp, item))) if 'low' in line] 
    df3=pd.read_csv('%s/%s_all_somatic_calls_forDPfilter_somaticSniper.csv'%(Variant_filtering_tmp, item), sep=',', engine='python', skiprows=exclude)
    df4=df3.drop(columns=['Depth','Unnamed: 0'])
    df4.to_csv('%s/%s_all_somatic_calls_somaticSniper.csv'%(Variant_filtering_tmp,item))
    

#Merging and Validation of calls
for item in Basename:
    df_mutect2=pd.read_csv('%s/%s_all_somatic_calls_Mutect2.csv'%(Variant_filtering_tmp,item))
    df_mutect2=df_mutect2.drop(columns=['Unnamed: 0','ID','ALTfreqNormal','ALTfreqTumor','SOMATIC'])
    df_strelka=pd.read_csv('%s/%s_all_somatic_calls_Strelka.csv'%(Variant_filtering_tmp,item))
    df_strelka=df_strelka.drop(columns=['Unnamed: 0','ID','ALTfreqNormal','ALTfreqTumor','SOMATIC'])
    df_somaticSniper=pd.read_csv('%s/%s_all_somatic_calls_somaticSniper.csv'%(Variant_filtering_tmp,item))
    df_somaticSniper=df_somaticSniper.drop(columns=['Unnamed: 0','ID','allel_freqNormal','allel_freqTumor','SOMATIC'])
    df_Varscan2=pd.read_csv('%s/%s_all_somatic_calls_Varscan2.csv'%(Variant_filtering_tmp,item))
    df_Varscan2=df_Varscan2.drop(columns=['Unnamed: 0','ID','ALTfreqNormal','ALTfreqTumor','SOMATIC'])


    
    df_Mutect2_and_strelka=df_mutect2.merge(df_strelka,how='inner') #AB
    df_Mutect2_and_SomaticSniper=df_mutect2.merge(df_somaticSniper,how='inner') #AC
    df_Mutect2_and_Varscan2=df_mutect2.merge(df_Varscan2,how='inner') #AD
    df_Strelka_and_SomaticSniper=df_strelka.merge(df_somaticSniper,how='inner') #BC
    df_Strelka_and_Varscan2=df_strelka.merge(df_Varscan2,how='inner') #BD
    df_SomaticSniper_and_Varscan2=df_somaticSniper.merge(df_Varscan2,how='inner') #CD

    df_AB_CD=df_Mutect2_and_strelka.merge(df_SomaticSniper_and_Varscan2, how='outer', indicator='AB_CD') #AB_CD
    df_AB_CD=df_AB_CD.rename(columns={'_merge':'VariantCaller'})
    df_AB_CD=df_AB_CD.replace({'both':'Mutect2;Strelka2;SomaticSniper;Varscan2','left_only':'Mutect2;Strelka2','right_only':'SomaticSniper;Varscan2'})
    df_AD_BC=df_Mutect2_and_Varscan2.merge(df_Strelka_and_SomaticSniper, how='outer', indicator='AD_BC') #AD_BC
    df_AD_BC=df_AD_BC.rename(columns={'_merge':'VariantCaller'})
    df_AD_BC=df_AD_BC.replace({'both':'Mutect2;Strelka2;SomaticSniper;Varscan2','left_only':'Mutect2;Varscan2','right_only':'Strelka2;SomaticSniper'})
    df_BD_AC=df_Strelka_and_Varscan2.merge(df_Mutect2_and_SomaticSniper, how='outer', indicator='BD_AC') #BD_AC
    df_BD_AC=df_BD_AC.rename(columns={'_merge':'VariantCaller'})
    df_BD_AC=df_BD_AC.replace({'both':'Mutect2;Strelka2;SomaticSniper;Varscan2','left_only':'Strelka2;Varscan2','right_only':'Mutect2;SomaticSniper'})

    df_AB_CD_AD_BC=df_AB_CD.merge(df_AD_BC, how='outer', indicator='ABCD_ADBC')
    df_AB_CD_AD_BC.to_csv('%s/ABCD_ADBC.csv'%Variant_filtering_tmp)

    def Caller(row):
        if row['ABCD_ADBC'] == 'left_only':
            val1 = row['AB_CD']
        elif row['ABCD_ADBC'] == 'right_only':
            val1 = row['AD_BC']
        elif row['ABCD_ADBC'] == 'both':
            val1 =  row['AB_CD']+';'+row['AB_CD']
        return val1

    df_AB_CD_AD_BC['Caller']=df_AB_CD_AD_BC.apply(Caller, axis=1)
    df_AB_CD_AD_BC=df_AB_CD_AD_BC.drop(columns=['ABCD_ADBC','AB_CD','AD_BC','AB_CD'])
    df_AB_CD_AD_BC['Caller']=df_AB_CD_AD_BC['Caller'].str.split(';').apply(set).str.join(';')
    df_AB_CD_AD_BC_BD_AC=df_AB_CD_AD_BC.merge(df_BD_AC, how='outer', indicator='ABCDADBC_BDAC')

    def Caller2(row):
        if row['ABCDADBC_BDAC'] == 'left_only':
            val1 = row['Caller']
        elif row['ABCDADBC_BDAC'] == 'right_only':
            val1 = row['BD_AC']
        elif row['ABCDADBC_BDAC'] == 'both':
            val1 =  row['Caller']+';'+row['BD_AC']
        return val1
    
    df_AB_CD_AD_BC_BD_AC['Format']=df_AB_CD_AD_BC_BD_AC.apply(Caller2, axis=1)
    df_AB_CD_AD_BC_BD_AC=df_AB_CD_AD_BC_BD_AC.drop(columns=['Caller','BD_AC','ABCDADBC_BDAC'])
    df_AB_CD_AD_BC_BD_AC['Format']=df_AB_CD_AD_BC_BD_AC['Format'].str.split(';').apply(set).str.join(',')
    df_AB_CD_AD_BC_BD_AC['Format1']='VariantCaller='
    df_AB_CD_AD_BC_BD_AC['Format']=df_AB_CD_AD_BC_BD_AC['Format1']+df_AB_CD_AD_BC_BD_AC['Format']
    df_AB_CD_AD_BC_BD_AC=df_AB_CD_AD_BC_BD_AC.drop(columns=['Format1'])
    
    df_all_Calls_2_of_4=df_AB_CD_AD_BC_BD_AC[['#CHROM']].copy()
    df_all_Calls_2_of_4['POS']=df_AB_CD_AD_BC_BD_AC['POS']
    df_all_Calls_2_of_4['ID']='.'
    df_all_Calls_2_of_4['REF']=df_AB_CD_AD_BC_BD_AC['REF']
    df_all_Calls_2_of_4['ALT']=df_AB_CD_AD_BC_BD_AC['ALT']
    df_all_Calls_2_of_4['QUAL']='.'
    df_all_Calls_2_of_4['FILTER']='.'
    df_all_Calls_2_of_4['INFO']=df_AB_CD_AD_BC_BD_AC['Format']
    df_all_Calls_2_of_4['FORMAT']='.'
    df_all_Calls_2_of_4=df_all_Calls_2_of_4.sort_values(['#CHROM','POS'], ascending=True)
    df_all_Calls_2_of_4.to_csv('%s/%s_Somatic_calls_2of4.vcf'%(Variant_filtering_tmp,item),sep='\t',index=False)

    #Called by 3 of 4 Callers !!!!Work!!!
  
    df_Mutect2_and_strelka_and_SomaticSniper=df_Mutect2_and_strelka.merge(df_somaticSniper) #ABC
    df_Mutect2_and_strelka_and_Varscan2=df_Mutect2_and_strelka.merge(df_Varscan2) #ABD
    df_Mutect2_and_SomaticSniper_and_Varscan=df_Mutect2_and_SomaticSniper.merge(df_Varscan2) #ACD
    df_Strelka_and_SomaticSniper_and_Varscan2=df_Strelka_and_SomaticSniper.merge(df_Varscan2) #BCD

    df_somatic_Variants_ABC_ABD=df_Mutect2_and_strelka_and_SomaticSniper.merge(df_Mutect2_and_strelka_and_Varscan2, how='outer', indicator='ABC_ABD')
    df_somatic_Variants_ABC_ABD=df_somatic_Variants_ABC_ABD.replace({'both':'Mutect2;Strelka2;SomaticSniper;Varscan2','right_only':'Mutect2;Strelka2;Varscan2','left_only':'Mutec2;Strelka2;SomaticSniper'})
    df_somatic_Variants_ACD_BCD=df_Mutect2_and_SomaticSniper_and_Varscan.merge(df_Strelka_and_SomaticSniper_and_Varscan2, how='outer', indicator='ACD_BCD')
    df_somatic_Variants_ACD_BCD=df_somatic_Variants_ACD_BCD.replace({'both':'Mutect2;Strelka2;SomaticSniper;Varscan2','right_only':'Mutect2;SomaticSniper;Varscan2','left_only':'Mutec2;SomaticSniper;Varscan2'})
    df_all_somatic_variants_Validated_3of4=df_somatic_Variants_ABC_ABD.merge(df_somatic_Variants_ACD_BCD, how='outer', indicator=True)
    
    def Caller3(row):
        if row['_merge'] == 'left_only':
            val1 = row['ABC_ABD']
        elif row['_merge'] == 'right_only':
            val1 = row['ACD_BCD']
        elif row['_merge'] == 'both':
            val1 =  row['ABC_ABD']+';'+row['ACD_BCD']
        return val1
   
    df_all_somatic_variants_Validated_3of4['INFO']=df_all_somatic_variants_Validated_3of4.apply(Caller3, axis=1)
    df_all_somatic_variants_Validated_3of4=df_all_somatic_variants_Validated_3of4.drop(columns=['ABC_ABD','ACD_BCD','_merge'])
    df_all_somatic_variants_Validated_3of4_vcf=df_all_somatic_variants_Validated_3of4[['#CHROM']].copy()
    df_all_somatic_variants_Validated_3of4_vcf['POS']=df_all_somatic_variants_Validated_3of4['POS']
    df_all_somatic_variants_Validated_3of4_vcf['ID']='.'
    df_all_somatic_variants_Validated_3of4_vcf['REF']=df_all_somatic_variants_Validated_3of4['REF']
    df_all_somatic_variants_Validated_3of4_vcf['ALT']=df_all_somatic_variants_Validated_3of4['ALT']
    df_all_somatic_variants_Validated_3of4_vcf['QUAL']='.'
    df_all_somatic_variants_Validated_3of4_vcf['FILTER']='.'
    df_all_somatic_variants_Validated_3of4_vcf['INFO']=df_all_somatic_variants_Validated_3of4['INFO']
    df_all_somatic_variants_Validated_3of4_vcf['FORMAT']='.'
    df_all_somatic_variants_Validated_3of4_vcf=df_all_somatic_variants_Validated_3of4_vcf.sort_values(['#CHROM','POS'], ascending=True)
    df_all_somatic_variants_Validated_3of4_vcf.to_csv('%s/%s_Somatic_SNP_and_indels_3of4_no_header.vcf'%(Variant_filtering_tmp,item),sep='\t',index=False)
    

    with open('%s/%s_ALL_somatic_Variants_2of4.vcf'%(Path_to_all_vcfs_filtered,item),'w') as ALL_2of4:
        ALL_2of4.write('''##fileformat=VCFv4.2
##FORMAT=<ID=VariantCaller,Number=1,Type=String,Description="VariantCaller_which_called_Varant">
''')
    ALL_2of4
    os.system('''cat %s/%s_Somatic_calls_2of4.vcf >> %s/%s_ALL_somatic_Variants_2of4.vcf'''%(Variant_filtering_tmp,item,Path_to_all_vcfs_filtered,item))
    

    with open('%s/%s_ALL_somatic_Variants_3of4.vcf'%(Path_to_all_vcfs_filtered,item),'w') as ALL_3of4:
        ALL_3of4.write('''##fileformat=VCFv4.2
##FORMAT=<ID=VariantCaller,Number=1,Type=String,Description="VariantCaller_which_called_Varant">
''')
    ALL_3of4.close    
    os.system('''cat %s/%s_Somatic_SNP_and_indels_3of4_no_header.vcf >> %s/%s_ALL_somatic_Variants_3of4.vcf'''%(Variant_filtering_tmp,item,Path_to_all_vcfs_filtered,item))

os.system('rm -r %s/mpileups'%global_tmp)
os.system('rm -r %s/Variants'%global_tmp)

#Variant Anotation and effectprediction

os.system('mkdir -p %s/Annotation'%log)
os.system('mkdir -p %s/Annotation'%error)
Annotation_log=('%s/Annotation'%log)
Annotation_error=('%s/Annotation'%error)

for item in Basename:
    with open('%s/%s_VEP_annotation.sh'%(Job_Scripts, item),'w') as VEP_annotation_Jobscript:
        VEP_annotation_Jobscript.write('''
#!/bin/sh
#$ -S /bin/sh

#$ -pe smp 1
#$ -cwd
#$ -V

#### Jobdescription at qstat
#$ -N VEP_ann 

#### Error Outputfile
#$ -e %s/%s_VEPann.dat
#$ -o %s/%s_VEPann.dat

#### Resubmit
#$ -r y

/usr/local/bioinf/vep/ensembl-vep/vep -i %s/%s_ALL_somatic_Variants_2of4.vcf -o %s/%s_VEPann.txt --cache --dir /data/databases/vep/ --assembly GRCh38 --offline --sift b --variant_class --polyphen b --humdiv --gene_phenotype --regulatory --symbol --force_overwrite
'''%(Annotation_error,item,
Annotation_log,item,
Path_to_all_vcfs_filtered,item,
Path_to_all_vcfs_filtered,item))
VEP_annotation_Jobscript.close

for item in Basename:
    with open('%s/%s_snpEff_annotation.sh'%(Job_Scripts, item),'w') as snpEff_annotation_Jobscript:
        snpEff_annotation_Jobscript.write('''
#!/bin/sh
#$ -S /bin/sh

#$ -pe smp 1
#$ -cwd
#$ -V

#### Jobdescription at qstat
#$ -N snpEff_ann 

#### Error Outputfile
#$ -e %s/%s_snpEffann.dat
#$ -o %s/%s_snpEffann.dat

#### Resubmit
#$ -r y

java -Xmx8g -jar /usr/local/bioinf/snpEff/snpEff_4.3T/snpEff/snpEff.jar /data/databases/snpEff/GRCh38.86 %s/%s_ALL_somatic_Variants_2of4.vcf > %s/%s_snpEff.txt
'''%(Annotation_error,item,
Annotation_log,item,
Path_to_all_vcfs_filtered,item,
Path_to_all_vcfs_filtered,item))
snpEff_annotation_Jobscript.close

for item in Basename:
    os.system('qsub  %s/%s_snpEff_annotation.sh'%(Job_Scripts, item))
    os.system('qsub -q long.q %s/%s_VEP_annotation.sh'%(Job_Scripts, item))

print('''
----------------------!!Done!!-----------------------
''')
