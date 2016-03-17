import data
import config
import re


def simple_filter(variant, name, logic, query_key, na_char='.', na_remain=True):

    value = get_value(variant, name)

    if check_na(value, na_char):
        return na_remain

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
        valuen = float(value)
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


def get_value(variant, name):
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


def combine_simple_filter(variant, names, logics, total_logic, query_keys, na_remains, na_char='.'):
    judge_list = []
    for i in range(len(names)):
        judge_list.append(
            simple_filter(variant, names[i], logics[i], query_keys[i], na_char=na_char, na_remain=na_remains[i])
        )

    if total_logic == 'ALL TRUE':
        return list(set(judge_list)) == [True]
    elif total_logic == 'NOT ALL TRUE':
        return False in judge_list
    elif total_logic == 'ALL FALSE':
        return list(set(judge_list)) == [False]
    elif total_logic == 'NOT ALL FALSE':
        return True in judge_list
    elif search_words(pat='\d+ TRUE', text=judge_list) != None:
        return judge_list.count(True) == int(total_logic.split()[0])
    elif search_words(pat='\d+ FALSE', text=judge_list) != None:
        return judge_list.count(False) == int(total_logic.split()[0])
    else:
        print 'Error: Unknown combine logic term:', total_logic
        raise








test1 = data.Vcf('/Users/rubisco/Desktop/test.vcf')


test2 = data.Annovar('/Users/rubisco/Desktop/test.annovar.txt', sample_columns=[81],
                info_columns=config.default_annovar_infos, vcf_columns=config.default_annovar_vcfs)
#print test1.title
#print test2.title


for variant in test2.variants:
    print variant.infos['Gene.refgene']
    print variant.infos
    print combine_simple_filter(variant, ['GERP++_RS'], ['is'], 'ALL TRUE', ['.'], na_remains=[True])

    break


