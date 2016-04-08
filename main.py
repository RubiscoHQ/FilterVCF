import data
import filter
import config
import os
from sys import argv


def interpret_cmd(cmd_list):
    cmd_dict = {}
    flag = ''
    for cmd in cmd_list:
        if cmd[0] == '-' or flag == '':
            if cmd in ['-input', '-I']:
                flag = 'input'
            elif cmd in ['-sample_info', '-SI']:
                flag = 'sample_info'
            elif cmd in ['-sample', '-S']:
                flag = 'sample'
                if flag not in cmd_dict:
                    cmd_dict[flag] = []
            elif cmd in ['-gene', '-G']:
                flag = 'gene'
                if flag not in cmd_dict:
                    cmd_dict[flag] = []
            elif cmd in ['-region', '-R']:
                flag = 'region'
                if flag not in cmd_dict:
                    cmd_dict[flag] = []
            elif cmd in ['-column_filter', '-CF']:
                flag = 'column_filter'
                if flag not in cmd_dict:
                    cmd_dict[flag] = []
            elif cmd in ['-model', '-M']:
                flag = 'model'
                if flag not in cmd_dict:
                    cmd_dict[flag] = []
            elif cmd in ['-total_logic', '-TL']:
                flag = 'total_logic'
            elif cmd in ['-output', '-O']:
                flag = 'output'
            else:
                print 'Error: can not interpret command:', cmd
                raise

        elif flag == 'input':
            cmd_dict[flag] = cmd
            flag = ''
        elif flag == 'sample_info':
            cmd_dict[flag] = cmd
            flag = ''
        elif flag == 'sample':
            cmd_dict[flag].append(cmd)
        elif flag == 'gene':
            cmd_dict[flag].append(cmd)
        elif flag == 'region':
            cmd_dict[flag].append(cmd)
        elif flag == 'column_filter':
            cmd_dict[flag].append(cmd)
        elif flag == 'model':
            cmd_dict[flag].append(cmd)
        elif flag == 'total_logic':
            cmd_dict[flag] = cmd
            flag = ''
        elif flag == 'output':
            cmd_dict[flag] = cmd
            flag = ''
        else:
            print 'Error: can not interpret flag:', flag
            raise

    if 'column_filter' in cmd_dict:
        if len(cmd_dict['column_filter']) == 4 and 'total_logic' not in cmd_dict:
            cmd_dict['total_logic'] = 'ALL_TRUE'

    return cmd_dict


def load_files(froute, sample_file):
    ext = os.path.splitext(froute)[1]
    if ext == '.vcf':
        vcf = data.Vcf(froute)
        return vcf, {}
    elif ext in ['.annovar', '.txt']:

        annovar_title_dict = config.check_annovar_title(froute, sample_file)
        annovar = data.Annovar(froute, sample_columns=annovar_title_dict['samples'],
                               info_columns=annovar_title_dict['infos'],
                               vcf_columns=annovar_title_dict['vcfs'],
                               id_columns=annovar_title_dict['ids'])
        return annovar, annovar_title_dict
    else:
        print 'Error: Unknown file type:', ext, 'For file:', os.path.split(froute)[1]
        raise


def make_column_filter(cmd_list, total_logic, variants):
    if len(cmd_list) % 4 != 0:
        print 'Error: Commend length for column filter is not correct. ' \
              'Each --column_filter or -CF should have 4 parameters.' \
              'Your cmd_list:', cmd_list
        raise TypeError
    else:
        cmd_dic = {'name':[],'logic':[],'query_key':[],'na':[]}
        flag = 'name'
        for i in cmd_list:
            if flag == 'name':
                cmd_dic[flag].append(i)
                flag = 'logic'
            elif flag == 'logic':
                cmd_dic[flag].append(i)
                flag = 'query_key'
            elif flag == 'query_key':
                cmd_dic[flag].append(i)
                flag = 'na'
            elif flag == 'na':
                cmd_dic[flag].append(i)
                flag = 'name'

    #print cmd_dic
    remain_variants = []
    for variant in variants:
        if filter.combine_simple_filter(variant=variant, names=cmd_dic['name'], logics=cmd_dic['logic'],
                                        total_logic=total_logic, query_keys=cmd_dic['query_key'],
                                        na_remains=cmd_dic['na']):
            remain_variants.append(variant)
    return remain_variants


def apply_model(froute, variants, models):
    sample_group = data.SampleGroup(froute)
    cohort = filter.Cohort(sample_group.cases, sample_group.ctrls, variants, 'Model_cohort')
    remain_var = []
    for model in models:
        if model == 'Dom':
            remain_var += cohort.dominant_var
        elif model == 'ResHom':
            remain_var += cohort.recessive_hom_var
        elif model == 'ResComp':
            remain_var += cohort.recessive_compound_var
        else:
            print 'Error: Unknown disease model:', model
            raise
    return remain_var


print '#==> STEP 1: Phrasing command'
command_dict = interpret_cmd(argv[1:])
print '#   ', command_dict

print '#==> STEP 2: Load input files'

database, title_dict = load_files(command_dict['input'], command_dict['sample_info'])
variants = database.variants
print '#   ', len(variants), 'variants remained'

print '#==> STEP 3: Simple columns filter'
try:
    variants = make_column_filter(command_dict['column_filter'], command_dict['total_logic'], variants)
except KeyError:
    print '#    Skip simple columns filter'
print '#   ', len(variants), 'variants remained'


print '#==> STEP 4: Get sample(s), gene(s) and  region(s) variant.'
try:
    variants = filter.get_samples_variants(variants, command_dict['sample'])
    print '#   ', len(variants), 'variants remained by sample filter'
except KeyError:
    print '#    Skip sample filter for not provide commands.'
try:
    variants = filter.get_gene_variants(variants, command_dict['gene'])
    print '#   ', len(variants), 'variants remained by gene filter'
except KeyError:
    print '#    Skip gene filter for not provide commands.'
try:
    variants = filter.get_regions_variants(variants, command_dict['region'])
    print '#   ', len(variants), 'variants remained by region filter'
except KeyError:
    print '#    Skip region filter for not provide commands.'


print '#==> STEP 5: Apply disease model.'
try:
    variants = apply_model(command_dict['sample_info'], variants, command_dict['model'])
except KeyError:
    print '    Skip model filter'

print '#   ', len(variants), 'variants remained'

try:
    print '#==> writing variants into file:', command_dict['output']
    f = open(command_dict['output']+'.annovar.txt', 'w')
    f.write(database.no_otherinfo_title)
    for variant in variants:
        f.write(variant.return_annovar_line(title=database.title.strip().split('\t'),
                                            sample_columns=title_dict['samples'],
                                            info_columns=title_dict['infos'],
                                            vcf_columns=title_dict['vcfs']))
except KeyError:
    print database.title
    for variant in variants:
        print variant.return_annovar_line(title=database.title.strip().split('\t'),
                                          sample_columns=title_dict['samples'],
                                          info_columns=title_dict['infos'],
                                          vcf_columns=title_dict['vcfs'])

print '#==> End'



