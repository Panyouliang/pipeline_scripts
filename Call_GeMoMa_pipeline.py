#! /usr/bin/env python3
# -*- coding:utf-8 -*-


import sys,argparse,os
import subprocess,glob
import shutil


#======================================================================================================
version = "v1.0"
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''
======================================================================
This pipeline is used for call GeMoMa annotating pipeline.

Version: v1.0
Author: Panyouliang, panyouliang@genomics.cn
Date: 2023-08-12, yyyy-mm-dd
======================================================================''')

parser.add_argument('-v', '--version', action='version', version=version)
parser.add_argument('-ref', metavar='fasta', type=str, required=True, help='Please input the reference genome [.fa, .fasta].')
parser.add_argument('-homolog_list', metavar='gff', type=str, required=True, help='Please input the homolog_list.txt.')
parser.add_argument('-bam_dir', metavar='path', type=str, required=True, help='Please input the bamfile path.')
parser.add_argument('-threads', metavar='int', type=str, default='4', required=False, help='Please input the threads number, default=4, [OPTIONAL].')
parser.add_argument('-maxintron', metavar='int', type=str, default='100000',required=False, help='Please input the max intron [int], default=100000, [OPTIONAL].')
parser.add_argument('-strand', metavar='rna-strandness', type=str, default='FR_UNSTRANDED', required=False, help='Defines whether the reads are stranded, range={FR_UNSTRANDED, FR_FIRST_STRAND, FR_SECOND_STRAND}, default = FR_UNSTRANDED, [OPTIONAL]')
parser.add_argument('-MEM', metavar='int', type=str,default='20',required=False, help='Please assign the Memory size, default=20G, (OPTIONAL).')

args = parser.parse_args()
#======================================================================================================


GeMoMa='/tools/GeMoMa/GeMoMa-1.9.jar'
run_java = '/tools/JAVA/bin/java'
mmseq_path = '/tools/mmseqs2/bin/'



def check_dir_exists(directory_path):
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)
        print(f"The directory {directory_path} has been created.")
    else:
        print(f"The directory {directory_path} already exists.")
        try:
            for files in os.listdir(directory_path):
                file_path = os.path.join(directory_path, files)
                if os.path.isfile(file_path):
                    os.remove(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
        except Exception as e:
            print(f"An error occurred: {e}")


def list_subdirectories(base_directory):
    subdirectories = []
    for root, dirs, files in os.walk(base_directory):
        for dir_name in dirs:
            relative_path = os.path.relpath(os.path.join(root, dir_name), base_directory)
            subdirectories.append(relative_path)
    return subdirectories


def extractor():

    check_dir_exists('extractor_data')

    with open(args.homolog_list,'r') as f:
        for line in f:
            line = line.strip()
            sp,gff,genome = line.split()[0],line.split()[1],line.split()[2]
            Extract = [run_java,'-Xms20G','-Xmx400G','-jar',GeMoMa,'CLI','Extractor', 'a='+gff,'g='+genome,'p=true','c=true','r=true','outdir='+'extractor_data/'+sp]
            subprocess.run(Extract)

def Pipeline():
    GeMoMaPipeline = [run_java,'-Xms20G','-Xmx400G','-jar',GeMoMa,'CLI', 'GeMoMaPipeline','t='+args.ref,'threads='+args.threads,'tblastn=false','p=true','pc=true','o=true','r=MAPPED','d=DENOISE','m='+mmseq_path,'AnnotationFinalizer.n=false','AnnotationFinalizer.u=YES','AnnotationFinalizer.r=NO','GeMoMa.Score=ReAlign','ERE.v=STRICT','ERE.mil=20','DenoiseIntrons.m='+args.maxintron,'GeMoMa.m='+args.maxintron,'ERE.s='+args.strand,'outdir=output']

    homolog = list_subdirectories('./extractor_data')

    for subdir in homolog:

        homolog_unit = ['s=pre-extracted','i='+subdir,'c=extractor_data/'+subdir+'/cds-parts.fasta','a=extractor_data/'+subdir+'/assignment.tabular']
        GeMoMaPipeline += homolog_unit

    bamdir = glob.glob(args.bam_dir+'/*.bam')
    bamlist = []
    for bam in bamdir:
        bamlist.append('ERE.m='+bam)

    GeMoMaPipeline += bamlist

    run_feature = open('GeMoMaPipeline.sh','w')
    run_feature.write(' '.join(x for x in GeMoMaPipeline)+'\n')
    run_feature.close()
    os.chmod('GeMoMaPipeline.sh', 0o775)


def main():
    extractor()
    Pipeline()



if __name__ == '__main__':
    main()


