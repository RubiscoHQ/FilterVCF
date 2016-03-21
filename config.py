# wannovar file title, update on 2016-03-20
default_annovar_info_names = ['Func.refgene', 'Gene.refgene', 'GeneDetail.refgene', 'ExonicFunc.refgene',
                                'AAChange.refgene', '1000G_ALL', '1000G_AFR', '1000G_AMR', '1000G_EAS',
                                '1000G_EUR', '1000G_SAS', 'ExAC_Freq', 'ExAC_AFR', 'ExAC_AMR', 'ExAC_EAS',
                                'ExAC_FIN', 'ExAC_NFE', 'ExAC_OTH', 'ExAC_SAS', 'ESP6500si_ALL', 'ESP6500si_AA',
                                'ESP6500si_EA', 'CG46', 'NCI60', 'dbSNP', 'COSMIC_ID', 'COSMIC_DIS', 'ClinVar_SIG',
                                'ClinVar_DIS', 'ClinVar_STATUS', 'ClinVar_ID', 'ClinVar_DB', 'ClinVar_DBID',
                                'GWAS_DIS', 'GWAS_OR', 'GWAS_BETA', 'GWAS_PUBMED', 'GWAS_SNP', 'GWAS_P',
                                'SIFT_score', 'SIFT_pred', 'Polyphen2_HDIV_score', 'Polyphen2_HDIV_pred',
                                'Polyphen2_HVAR_score', 'Polyphen2_HVAR_pred', 'LRT_score', 'LRT_pred',
                                'MutationTaster_score', 'MutationTaster_pred', 'MutationAssessor_score',
                                'MutationAssessor_pred', 'FATHMM_score', 'FATHMM_pred', 'RadialSVM_score',
                                'RadialSVM_pred', 'LR_score', 'LR_pred', 'VEST3_score', 'CADD_raw', 'CADD_phred',
                                'GERP++_RS', 'phyloP46way_placental', 'phyloP100way_vertebrate', 'SiPhy_29way_logOdds']

default_annovar_ids = ['Chr', 'Start', 'End', 'Ref', 'Alt']

default_vcf_title = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']


def check_annovar_title(route):
    f = open(route)
    title = f.readline()
    title_list = title.strip().split('\t')
    dic = {}
    if title_list[:5] == default_annovar_ids \
       and title_list[5:69] == default_annovar_info_names \
       and title_list[72:81] == default_vcf_title:
        dic = {'ids': [i for i in range(5)],
               'infos': [i for i in range(5, 69)],
               'vcfs': [i for i in range(72, 81)],
               'samples': [i for i in range(81, len(title_list))]}
    else:
        dic = {'ids': [], 'infos': [], 'vcfs': [], 'samples': []}
        for i in title_list:
            if i in default_annovar_ids:
                dic['ids'].append(i)
            elif i in default_annovar_info_names:
                dic['infos'].append(i)
            elif i in default_vcf_title:
                dic['vcfs'].append(i)
            else:
                print 'Unknown column:', i
        print 'Annovar title is not default, please revise it or specify specific column types or revise config.py.'
        raise
    return dic


genotype_empty_symble = ['0/0', './.']


gene_name_column_in_wannovar = 'Gene.refgene'


info_na_char = '.'