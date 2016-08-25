import os
import sys
import glob
import logging
import argparse
import subprocess

def main(trim_path, bowtie_path, picard_path, gatk_path,
        input_path, index_path, dbsnp_path, adapter_path,
        ref_path, out_path):

    #Get complete path
    trim_path = os.path.abspath(trim_path)
    bowtie_path = os.path.abspath(bowtie_path)
    picard_path = os.path.abspath(picard_path)
    gatk_path = os.path.abspath(gatk_path)
    input_path = [os.path.abspath(files) for files in input_path]
    index_path = os.path.abspath(index_path)
    dbsnp_path = os.path.abspath(dbsnp_path)
    ref_path = os.path.abspath(ref_path)
    out_path = os.path.abspath(out_path)
    adpater_path = os.path.abspath(adapter_path)
    #Check if paths exist
    if not os.path.exists(trim_path):
        raise FileNotFoundError('Trimmomatic not found at {0}'.format(trim_path))
    
    if not os.path.exists(bowtie_path):
        raise FileNotFoundError('Bowtie not found at {0}'.format(bowtie_path))

    if not os.path.exists(picard_path):
        raise FileNotFoundError('Picard not found at {0}'.format(picard_path))

    if not os.path.exists(gatk_path):
        raise FileNotFoundError('Gatk not found at {0}'.format(gatk_path))

    for files in input_path:
        if not os.path.exists(files):
            raise FileNotFoundError('Fastq files not found at {0}'.format(files))
    
    indicies = glob.glob('{0}.*.bt2'.format(index_path))
    print(indicies)
    if len(indicies) == 0:
        raise FileNotFoundError('Bowtie index not found at {0}'.format(index_path))
    
    if not os.path.exists(ref_path):
        raise FileNotFoundError('Reference fasta file not found at {0}'.format(ref_path))

    if not os.path.exists(dbsnp_path):
        raise FileNotFoundError('dbSNP file not found at {0}'.format(dbsnp_path))
    
    if not os.path.exists(adapter_path):
        raise FileNotFoundError('Adapter file not found at {0}'.format(adapter_path))


    #Creat output directory
    if not os.path.exists(out_path):
        os.mkdir(out_path)


    #Trim fastq files
    read1 = input_path[0]
    read2 = input_path[1]
    tread1 = '{1}_trimmed.fq'.format(out_path, os.path.splitext(read1)[0])
    tread2 = '{1}_trimmed.fq'.format(out_path, os.path.splitext(read2)[0])
    sread1 = '{1}_unused.fq'.format(out_path, os.path.splitext(read1)[0])
    sread2 = '{1}_unused.fq'.format(out_path, os.path.splitext(read2)[0])

    tcmd = ['java', '-jar', trim_path, 'PE', '-phred33', read1, read2, tread1,
            sread1, tread2, sread2, 'ILLUMINACLIP:{0}:2:30:10'.format(adapter_path),
            'LEADING:0', 'TRAILING:0', 'SLIDINGWINDOW:4:15', 'MINLEN:36']
            
    trun = subprocess.Popen(tcmd, shell=False)
    trun.wait() 
    
    if trun.returncode != 0:
        print('Fastq trimming failed; Exiting program')
        sys.exit()
         

    #Align the reads using bowtie
    sam_path = '{1}.sam'.format(out_path, os.path.splitext(tread1)[0])
    bcmd = [ bowtie_path, '-x', index_path, '-S', sam_path, '-p', '1' , '-1',
            tread1, '-2', tread2]
    
    brun = subprocess.Popen(bcmd, shell=False)
    brun.wait()
    
    if brun.returncode != 0:
        print('Bowtie failed; Exiting program')
        sys.exit()

    #Add read group information
    add_path = '{0}/{1}_RG.bam'.format(out_path, os.path.splitext(os.path.basename(sam_path))[0])
    acmd = ['java', '-Xmx1g', '-jar', picard_path, 'AddOrReplaceReadGroups',
        'I='+sam_path , 'O='+add_path, 'SORT_ORDER=coordinate', 'RGID=Test', 
        'RGLB=ExomeSeq', 'RGPL=Illumina', 'RGPU=HiSeq2500', 'RGSM=Test', 
        'RGCN=AtlantaGenomeCenter', 'RGDS=ExomeSeq', 'RGDT=2016-08-24', 'RGPI=null', 
        'RGPG=Test', 'RGPM=Test', 'CREATE_INDEX=true']
    
    arun = subprocess.Popen(acmd, shell=False)
    arun.wait()
    
    if arun.returncode != 0:
        print('Picard add read groups failed; Exiting program')
        sys.exit()

    #Mark PCR duplicates
    dup_path = '{0}/{1}_MD.bam'.format(out_path, os.path.splitext(os.path.basename(sam_path))[0])
    met_path = '{0}/{1}_MD.metrics'.format(out_path, os.path.splitext(os.path.basename(sam_path))[0])
    mdcmd = ['java', '-Xmx1g', '-jar', picard_path, 'MarkDuplicates', 'I='+add_path, 
        'O='+dup_path, 'METRICS_FILE='+met_path, 'REMOVE_DUPLICATES=false', 
        'ASSUME_SORTED=true', 'CREATE_INDEX=true']
    
    mdrun = subprocess.Popen(mdcmd, shell=False)
    mdrun.wait()
    if mdrun.returncode != 0:
        print('Picard mark duplicate failed; Exiting program')
        sys.exit()

    #Fix mate information
    fix_path = '{0}/{1}_FM.bam'.format(out_path, os.path.splitext(os.path.basename(sam_path))[0])
    fcmd = ['java', '-Xmx1g', '-jar', picard_path, 'FixMateInformation',
        'I='+dup_path, 'O='+fix_path, 'ASSUME_SORTED=true', 'ADD_MATE_CIGAR=true',
        'CREATE_INDEX=true']

    frun = subprocess.Popen(fcmd, shell=False)
    frun.wait()
    
    if frun.returncode != 0:
        print('Picard fix mate information failed; Exiting program')
        sys.exit()
   
    #Run realigner target creator
    interval_path = '{0}/{1}.intervals'.format(out_path, os.path.splitext(os.path.basename(sam_path))[0]) 
    
    trcmd = ['java', '-jar', gatk_path, '-T', 'RealignerTargetCreator', '-o',
        interval_path, '-nt', '1', '-I', fix_path, '-R', ref_path, '-known',
        dbsnp_path]
    
    trrun = subprocess.Popen(trcmd, shell=False)
    trrun.wait()
    
    if trrun.returncode != 0:
        print('Realigner Target creator failed; Exiting program')
        sys.exit()
     

    #Run indel realigner
    ral_path = '{0}/{1}_IR.bam'.format(out_path, os.path.splitext(os.path.basename(sam_path))[0])
    
    recmd = ['java', '-jar', gatk_path, '-T', 'IndelRealigner',
        '--targetIntervals', interval_path, '-o', ral_path,
        '-I', fix_path, '-R', ref_path]

    rerun = subprocess.Popen(recmd, shell=False)
    rerun.wait()

    if rerun.returncode != 0:
        print('Indel realigner creator failed; Exiting program')
        sys.exit()

    #Base quality score recalibration
    bqs_path = '{0}/{1}.table'.format(out_path, os.path.splitext(os.path.basename(sam_path))[0])
     
    bqscmd = ['java', '-jar', gatk_path, '-T', 'BaseRecalibrator', '-R', ref_path,
        '-I', ral_path, '-o', bqs_path, '-nct', '1', '-cov', 'ReadGroupCovariate',
        '-knownSites', dbsnp_path]

    bqsrun = subprocess.Popen(bqscmd, shell=False)
    bqsrun.wait()

    if bqsrun.returncode != 0:
        print('Base quality score recalibrator failed; Exiting program')
        sys.exit()
    
    #Print Reads
    fbam_path = '{0}/{1}_final.bam'.format(out_path, os.path.splitext(os.path.basename(sam_path))[0])
    prcmd = ['java', '-jar', gatk_path, '-T', 'PrintReads', '-R', ref_path, '-I',
            ral_path, '-o', fbam_path, '-BQSR', bqs_path, '-nct', '1']


    prrun = subprocess.Popen(prcmd, shell=False)
    prrun.wait()

    if prrun.returncode != 0:
        print('Print reads failed; Exiting program')
        sys.exit()


    #Haplotype caller
    vcf_path = '{0}/variants.vcf'.format(out_path, os.path.splitext(os.path.basename(sam_path))[0])

    hcmd = ['java', '-jar', gatk_path, '-T', 'HaplotypeCaller', '-R', ref_path,
        '-I', fbam_path, '--dbsnp', dbsnp_path, '-o', vcf_path, '-nct', '1', 
        '-gt_mode', 'DISCOVERY']

    hrun = subprocess.Popen(hcmd, shell=False)
    hrun.wait()
    
    if hrun.returncode != 0:
        print('Haplotype caller failed; Exiting program')
        sys.exit()


    print('Variant call pipeline completed')
    print('VCF file can be found at {0}'.format(vcf_path))
    return

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='VariantCaller')
    parser.add_argument('-t', '--trimmomatic', dest='trim_path', type=str, help='Path to Trimmomatic')
    parser.add_argument('-b', '--bowtie', dest='bowtie_path', type=str, help='Path to Bowtie')
    parser.add_argument('-p', '--picard', dest='picard_path', type=str, help='Path to Picard')
    parser.add_argument('-g', '--gatk', dest='gatk_path', type=str, help='Path to GATK')
    parser.add_argument('-i', '--inputs', dest='input_path', nargs='+', type=str, help='Path to Paired end reads')
    parser.add_argument('-w', '--index', dest='index_path', type=str, help='Path to Reference bowtie index')
    parser.add_argument('-d', '--dbsnp', dest='dbsnp_path', type=str, help='Path to dbSNP vcf file')
    parser.add_argument('-r', '--refrence', dest='ref_path', type=str, help='Path to Reference file')
    parser.add_argument('-a', '--adapter', dest='adapter_path', type=str, help='Path to Adapter file')
    parser.add_argument('-o', '--outpath', dest='out_path', type=str, help='Path to Ouput directory')
    args = parser.parse_args()

main(args.trim_path, args.bowtie_path, args.picard_path, args.gatk_path,
    args.input_path, args.index_path, args.dbsnp_path, args.adapter_path,
    args.ref_path, args.out_path)
