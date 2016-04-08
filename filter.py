#coding=utf-8

import config
import re


def get_value(variant, name):
    name = name.lower()
    try:
        value = variant.infos[name]
    except KeyError:
        print "Error: variant information dict do not contain name:", name
        raise
    return value


def check_na(value, na_char):
    if value in na_char:
        return True
    else:
        return False


def search_words(pat, text):
    pattern = re.compile(r'%s' %pat,re.S)
    match = pattern.search(text)
    if match:
        return match.group()
    else:
        return None

    # search a string and return the FIRST matched word in a list.


def get_samples_gt_from_variant(sample_list, variant):
    rl = []
    for i in sample_list:
        gt = variant.samples[i.id.lower()]['GT']

        if gt in config.genotype_empty_symble:
            rl.append('no')
            continue
        gtl = gt.split('/')
        if len(gtl) != 2:
            rl.append('other')
        elif gtl[0] == gtl[1]:
            rl.append('hom')  # gt == '1/1', '2/2', ...
        else:
            rl.append('het')  # gt == '0/1', '1/2', ...
    return rl


def make_combination(lists, max_com=2):
    total_list = []
    len_lists = len(lists)
    if len_lists == max_com:
        return [lists]
    elif max_com == 1:
        for t in lists:
            total_list.append([t])
        return total_list
    elif max_com > len_lists:
        print 'Error max_com'
    else:
        for i in range(len_lists-max_com+1):
            settle = lists[i]
            unsettle = make_combination(lists[i+1:], max_com-1)
            for k in unsettle:
                total_list.append([settle]+k)
    return total_list


def simple_filter(variant, name, logic, query_key, na_char=config.info_na_char, na_remain=True):

    value = get_value(variant, name)

    if check_na(value, na_char):
        if na_remain == 'T':
            return True
        elif na_remain == 'F':
            return False

    elif isinstance(query_key, str):
        if logic == 'in':
            return value in query_key
        elif logic == '!in':
            return value not in query_key
        elif logic == 'include':
            return query_key in value
        elif logic == '!include':
            return query_key not in value
        elif logic == 'is':
            return query_key == value
        elif logic == '!is':
            return query_key != value
        else:
            print 'Error: Unknown logic term:', logic, '\nWith query key:', query_key
            raise

    elif isinstance(query_key, (int, float)):
        try:
            valuen = float(value)
        except ValueError:
            print 'Error for value:', value, type(value)
            raise
        if logic == '>':
            return valuen > query_key
        elif logic == '<':
            return valuen < query_key
        elif logic == '=':
            return valuen == query_key
        elif logic == '!=':
            return valuen != query_key
        elif logic == '<=':
            return valuen <= query_key
        elif logic == '>=':
            return valuen >= query_key
        else:
            print 'Error: Unknown logic term:', logic, '\nWith query key:', query_key
            raise


def combine_simple_filter(variant, names, logics, total_logic, query_keys, na_remains, na_char=config.info_na_char):
    judge_list = []
    for i in range(len(names)):
        try:
            qk = float(query_keys[i])
        except ValueError:
            qk = query_keys[i]
        judge_list.append(
            simple_filter(variant, names[i], logics[i], qk, na_char=na_char, na_remain=na_remains[i])
        )

    if total_logic == 'ALL_TRUE':
        return list(set(judge_list)) == [True]
    elif total_logic == 'NOT_ALL_TRUE':
        return False in judge_list
    elif total_logic == 'ALL_FALSE':
        return list(set(judge_list)) == [False]
    elif total_logic == 'NOT_ALL_FALSE':
        return True in judge_list
    elif search_words(pat='\d+_TRUE', text=total_logic) != None:
        return judge_list.count(True) == int(total_logic.split('_')[0])
    elif search_words(pat='\d+_FALSE', text=total_logic) != None:
        return judge_list.count(False) == int(total_logic.split('_')[0])
    else:
        print 'Error: Unknown total logic term:', total_logic
        raise
# input single variant, filter rules list, combine several single column filter and give a final result.


def get_gene_variants(variants, genes, gene_column=config.gene_name_column_in_wannovar):
    remain_var = []
    for var in variants:
        flag = False
        for gene in genes:
            if simple_filter(var, gene_column, 'include', gene, na_remain=False):
                flag = True
        if flag:
            remain_var.append(var)
    return remain_var
# input single gene, return gene's variant in a given database (file).


def get_samples_variants(variants, samples, empty_sample_GT=config.genotype_empty_symble):
    remain_var = []

    for var in variants:
        flag = False
        for sam in samples:
            try:
                if var.samples[sam]['GT'] not in empty_sample_GT:
                    flag = True
            except KeyError:
                print 'Error: Unknow sample id:', sam, ' when finding sample variants.'
        if flag:
            remain_var.append(var)

    return remain_var
# input samples, return sample's variant in a given database (file).


def get_regions_variants(variants, regions):
    remain_var = []
    if len(regions) % 3 != 0:
        print 'Error: Commend length for region filter is not correct. ' \
              'Each --region or -R should have 3 parameters.'
        raise

    for var in variants:
        flag = False
        rl = ['', '', '']
        for i in range(len(regions)):
            rl[i % 3] = regions[i]
            if (i % 3 == 2):
                if combine_simple_filter(var, names=['chrom', 'vcfpos', 'vcfpos'], logics=['is', '>=', '<='],
                                         total_logic='ALL_TRUE', query_keys=rl, na_remains=[False]*3):
                    flag = True
        if flag:
            remain_var.append(var)

    return remain_var


class Cohort():
    def __init__(self, case_id_list, ctrl_id_list, variants, cohort_name):
        self.cases = case_id_list
        self.ctrls = ctrl_id_list
        self.sample_list = case_id_list+ctrl_id_list
        self.variants = variants
        self.name = cohort_name
        self.other_var = []
        self.dominant_model()
        self.recessive_compound_model()
        self.recessive_hom_model()

    def dominant_model(self):
        self.dominant_var = []
        self.dominant_gene = {}
        for var in self.variants:
            cases_gt = get_samples_gt_from_variant(self.cases, var)
            ctrls_gt = get_samples_gt_from_variant(self.ctrls, var)
            if 'other' in cases_gt or 'other' in ctrls_gt:
                self.other_var.append(var)
                continue

            if 'no' not in cases_gt and list(set(ctrls_gt)) == ['no']:
                self.dominant_var.append(var)
                genes = var.infos[config.gene_name_column_in_wannovar].split(',')
                for gene in genes:
                    try:
                        self.dominant_gene[gene] += 1
                    except KeyError:
                        self.dominant_gene[gene] = 1
        return

    def recessive_hom_model(self):
        self.recessive_hom_var = []
        self.recessive_hom_gene = {}
        for var in self.variants:
            cases_gt = get_samples_gt_from_variant(self.cases, var)
            ctrls_gt = get_samples_gt_from_variant(self.ctrls, var)
            if 'other' in cases_gt or 'other' in ctrls_gt:
                self.other_var.append(var)
                continue
            if list(set(cases_gt)) == ['hom'] and 'hom' not in ctrls_gt:
                self.recessive_hom_var.append(var)
                genes = var.infos[config.gene_name_column_in_wannovar].split(',')
                for gene in genes:
                    try:
                        self.recessive_hom_gene[gene] += 1
                    except KeyError:
                        self.recessive_hom_gene[gene] = 1
        return

    def recessive_compound_model(self):
        self.recessive_compound_var = []
        self.recessive_compound_gene = {}
        ctrls_count = len(self.ctrls)
        gvd = {}
        for var in self.variants:
            cases_gt = get_samples_gt_from_variant(self.cases, var)
            ctrls_gt = get_samples_gt_from_variant(self.ctrls, var)
            if 'other' in cases_gt or 'other' in ctrls_gt:
                self.other_var.append(var)
                continue
            if ('no' in cases_gt) or (list(set(cases_gt)) == ['hom']) or ('hom' in ctrls_gt):
                continue

            genes = var.infos[config.gene_name_column_in_wannovar].split(',')
            for gene in genes:
                if not gvd. has_key(gene):
                    gvd[gene] = {}

                try:
                    gvd[gene][var] = {'cases_gt': cases_gt,
                                      'ctrls_gt': ctrls_gt}
                except KeyError:
                    gvd[gene] = {var: {'cases_gt': cases_gt,
                                       'ctrls_gt': ctrls_gt}}
        for gene in gvd:
            variants = gvd[gene].keys()
            if len(variants) < 2:
                continue
            var_combination = make_combination(variants, 2)  # make pairs of variation combination
            remain_var = []
            for pair in var_combination:
                for i in range(ctrls_count):
                    if gvd[gene][pair[0]]['ctrls_gt'][i] == 'no' or \
                       gvd[gene][pair[1]]['ctrls_gt'][i] == 'no':
                        remain_var += pair
            remain_var = list(set(remain_var))
            if remain_var != []:
                self.recessive_compound_var += remain_var
                self.recessive_compound_gene[gene] = len(remain_var)
        return
        # recessive compound het variants means at least tow mutations
        # occur in one gene that accord with the following rules:
        # 1. any cases should carry this variant but not all with gt=hom (all hom is included in recessive hom model);
        # 2. not any ctrls carry this variant with gt=hom;
        # 3. for pairs of variants in one genes, every ctrl do not carry these pairs at the same time


    def return_vcf(self, head, model, route=''):
        return

    def return_annovar(self):
        return


