#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import numpy as np
import scipy.stats
import os
import sys
import statsmodels.sandbox.stats.multicomp

assert sys.version_info >= (2, 7)


def parse_args():
    starting_parser = argparse.ArgumentParser(description="This tool counts statistics from linker file containing two tab-delimited columns without a header. The first column must contain file path, the second - group name.")
    starting_parser.add_argument("-g", "--groupdata", required=True,
                                 help="Text file without a header containing two tab-delimited columns: file path and group name")
    starting_parser.add_argument("-i", "--index", required=True,
                                 help="Column name to merge across tables")
    starting_parser.add_argument("-v", "--value", required=True,
                                 help="Column name to count statistics")
    starting_parser.add_argument("-a", "--alpha", type=float, default=0.05,
                                 help="(Optional) Alpha value to count statistics, 0.05 by default")
    starting_parser.add_argument("-n", "--no_filter", default=False, action='store_true',
                                 help="(Optional) Disables pre-filtration")
    starting_parser.add_argument("-o", "--output", required=True,
                                 help="Output directory")
    return starting_parser.parse_args()


def ends_with_slash(string):
    if string.endswith("/"):
        return string
    else:
        return str(string + "/")


def parse_namespace():
    namespace = parse_args()
    return namespace.groupdata, namespace.index, namespace.value, namespace.alpha, namespace.no_filter, ends_with_slash(namespace.output)


def is_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        print(exception)


def find_existing_files(paths_list):
    output_list = []
    for file_name in paths_list:
        if os.path.isfile(file_name):
            output_list.append(file_name)
        else:
            print("Not found: " + file_name)
    if len(output_list) == 0:
        print("No valid files to process! Exiting...")
        sys.exit(2)
    print("Queue files number:", len(output_list))
    return output_list


def file_to_list(file):
    file_wrapper = open(file, 'rU')
    output_list = sorted(filter(None, [file_row.replace('\r', '').replace('\n', '') for file_row in file_wrapper]))
    file_wrapper.close()
    return output_list


def append_dict(dictionary, key, value):
    if key in dictionary:
        dictionary[key].append(value)
    else:
        dictionary[key] = [value]
    return dictionary


def list2pairwise_tuples(flat_list):
    tuples = []
    for list_name1 in flat_list:
        for list_name2 in flat_list:
            if list_name1 != list_name2 and sorted([list_name1, list_name2]) not in [sorted(sublist) for sublist in tuples]:
                tuples.append([list_name1, list_name2])
    return tuples


def multi_core_queue(function_to_parallelize, queue):
    import multiprocessing
    pool = multiprocessing.Pool()
    output = pool.map(function_to_parallelize, queue)
    pool.close()
    pool.join()
    return output


def merge_dicts(*args):
    output_dict = {}
    for arg in args:
        output_dict.update(arg)
    return output_dict


def dict2pd_series(dictionary):
    output = pd.Series()
    for key in dictionary:
        output = output.set_value(key, dictionary[key])
    return output


def values2pvals_dict(list_1, list_2):
    output_dict = {}
    for method_name, method_function in zip(["t-test", "u-test", "wilcoxon", "h-test"], [scipy.stats.ttest_ind, scipy.stats.mannwhitneyu, scipy.stats.ranksums, scipy.stats.kruskal]):
        if np.count_nonzero(list_1) == 0 and np.count_nonzero(list_2) == 0:
            output_dict[method_name] = {"statistic": np.nan, "p-value": np.nan}
        else:
            output_dict[method_name] = {i: j for (i, j) in zip(["statistic", "p-value"], method_function(list_1, list_2))}
    return output_dict


def unpack_methods():
    output_dict = {}
    for input_pair_string in list(groupDictPairsDict):
        for statistical_test in list(groupDictPairsDict[input_pair_string]):
            try:
                output_dict[statistical_test].update({input_pair_string: groupDictPairsDict[input_pair_string][statistical_test]})
            except KeyError:
                output_dict[statistical_test] = {input_pair_string: groupDictPairsDict[input_pair_string][statistical_test]}
    return output_dict


def row2multitest(pvals_list, input_colnames_list):
    output_dict = {}
    alphac_sidak = 0
    alphac_bonf = 0
    for method in [  # Description: http://www.statsmodels.org/devel/generated/statsmodels.sandbox.stats.multicomp.multipletests.html#statsmodels.sandbox.stats.multicomp.multipletests
        "bonferroni",  # one-step correction
        "sidak",  # one-step correction
        "holm-sidak",  # step down method using Sidak adjustments
        "holm",  # step-down method using Bonferroni adjustments
        "simes-hochberg",  # step-up method  (independent)
        "hommel",  # closed method based on Simes tests (non-negative)
        "fdr_bh",  # Benjamini/Hochberg  (non-negative)
        "fdr_by",  # Benjamini/Yekutieli (negative)
        "fdr_tsbh",  # two stage fdr correction (non-negative)
        "fdr_tsbky",  # two stage fdr correction (non-negative)
    ]:
        rejects_list, pvals_corrected_list, alphac_sidak, alphac_bonf = statsmodels.sandbox.stats.multicomp.multipletests(pvals_list, method=method, alpha=globalAlpha)
        output_dict.update(zip([k + "_is_rejected_by_" + method for k in input_colnames_list], [str(i) for i in list(rejects_list)]))
        output_dict.update(zip([l + "_corrected_for_" + method for l in input_colnames_list], [float(j) for j in list(pvals_corrected_list)]))
    output_dict.update(zip(["raw_alpha", "alpha_corrected_for_Sidak_method", "alpha_corrected_for_Bonferroni_method"], [globalAlpha, alphac_sidak, alphac_bonf]))
    return output_dict


def mp_group_name2df(group_name):
    def read_table_slice(table_file_name):
        df = pd.read_table(table_file_name, sep='\t', header=0, engine='python')
        output_dataframe = df.loc[:, [indexColName, valueColName]]
        return output_dataframe.rename(columns={valueColName: table_file_name})
    group_dataframe = pd.DataFrame()
    group_files_list = groupDataDF.loc[groupDataDF['group_name'] == group_name, 'file_name'].values.tolist()
    if len(group_files_list) < 2:
        print("Excluded group due to the files number smaller than 2:", group_name)
    else:
        for file_name in group_files_list:
            try:
                if len(group_dataframe) != 0:
                    group_dataframe = pd.merge(group_dataframe, read_table_slice(file_name), on=indexColName, how='outer')
                else:
                    group_dataframe = read_table_slice(file_name)
            except FileNotFoundError:
                print("Not found file in group", group_name, ":", file_name)
        pivot_file_name = outputDir + '_'.join(["pivot", "by", valueColName, group_name + ".tsv"])
        group_dataframe.to_csv(pivot_file_name, sep='\t', index=False)
        print("Dumped the group '" + group_name + "' pivot:", pivot_file_name)
    return group_name, group_dataframe


def mp_index2pvals(index_string):
    values_lists_dict = {i: j for (i, j) in zip(groupPairList, [groupMergedDF.loc[index_string, values_columns_list].fillna(0).astype('float64').values.tolist() for values_columns_list in [valuesColumnsList_1, valuesColumnsList_2]])}
    base_metrics_dict = {indexColName: index_string}
    for group_name in values_lists_dict:
        values_list = [float(i) for i in values_lists_dict[group_name]]
        base_metrics_dict["total_values_" + group_name] = len(values_list)
        base_metrics_dict["non-zero_values_" + group_name] = np.count_nonzero(values_list)
        base_metrics_dict["non-zero_percentage_" + group_name] = float(base_metrics_dict["non-zero_values_" + group_name]) / float(base_metrics_dict["total_values_" + group_name])
        base_metrics_dict["sum_" + group_name] = np.sum(values_list)
        base_metrics_dict["mean_" + group_name] = np.mean(values_list)
        base_metrics_dict["median_" + group_name] = np.median(values_list)
        base_metrics_dict["sd_" + group_name] = np.std(values_list)
        base_metrics_dict["variance_" + group_name] = np.var(values_list)
    if base_metrics_dict["sum_" + groupPairList[0]] > base_metrics_dict["sum_" + groupPairList[1]]:
        base_metrics_dict["prevalent"] = groupPairList[0]
    elif base_metrics_dict["sum_" + groupPairList[0]] < base_metrics_dict["sum_" + groupPairList[1]]:
        base_metrics_dict["prevalent"] = groupPairList[1]
    else:
        base_metrics_dict["prevalent"] = 'equivalent'
    p_values_dict = values2pvals_dict(values_lists_dict[groupPairList[0]], values_lists_dict[groupPairList[1]])
    return {i: dict2pd_series(merge_dicts(base_metrics_dict, p_values_dict[i])) for i in p_values_dict}


def mp_adjust_method_pvals(index_string):
    adjusted_pvals_dict = {indexColName: index_string}
    adjusted_pvals_dict.update(row2multitest(valuesSubDF.loc[index_string, :].astype('float64').values.tolist(), list(valuesSubDF)))
    return adjusted_pvals_dict


if __name__ == '__main__':
    inputGroupDataFile, indexColName, valueColName, globalAlpha, noFilterBool, outputDir = parse_namespace()
    # inputGroupDataFile, indexColName, valueColName, globalAlpha, noFilterBool, outputDir = ['/data1/bio/projects/Danilova_Natalya/IGC_KOs_for_acetate.groupdata', 'reference_id', 'id_mapped_reads_per_million_sample_total_reads', 0.05, False, '/data1/bio/projects/Danilova_Natalya/test/']
    # inputGroupDataFile, indexColName, valueColName, globalAlpha, noFilterBool, outputDir = ['/data2/bio/Metagenomes/SampleData/CARD_HP_C_1_2_3.groupdata', 'reference_id', 'id_mapped_reads_per_million_sample_total_reads', 0.05, False, '/data2/bio/Metagenomes/ABR/CARD/pvals/']
    print("Creating directory:", outputDir)
    is_path_exists(outputDir)
    print("Parsing groups")
    groupDataDF = pd.read_table(inputGroupDataFile, sep='\t', header='infer', names=['file_name', 'group_name'], engine='python')
    groupsRawList = sorted(groupDataDF['group_name'].unique().tolist())
    print("Merging data from groups:", "'" + "', '".join(groupsRawList) + "'")

    groupsDict = {}
    for parsedGroupName, parsedGroupDF in multi_core_queue(mp_group_name2df, groupsRawList):
        if len(parsedGroupDF) > 0:
            groupsDict[parsedGroupName] = parsedGroupDF
    if len(groupsDict) == 0:
        raise ValueError("Nothing to parse!")
    elif len(groupsDict) == 1:
        raise ValueError("At least 2 groups are required to compare!")
    print("Successfully parsed groups:", "'" + "', '".join(list(groupsDict)) + "'")
    groupsPairs2DArray = list2pairwise_tuples(list(groupsDict))
    print("Pairs to count p-values:", "'" + "', '".join(list(["_vs_".join(i) for i in groupsPairs2DArray])) + "'")

    comparedIndexesMethodsDictsLists2DArray = []

    groupDictPairsDict = {}
    for groupPairList in groupsPairs2DArray:
        print("Comparing groups:", "'" + "' vs '".join(str(i) for i in groupPairList) + "'")
        groupPairDF_1, groupPairDF_2 = [groupsDict[i] for i in groupPairList]
        groupMergedDF = pd.merge(groupPairDF_1, groupPairDF_2, on=indexColName, how='outer').fillna(0).set_index(indexColName)
        valuesColumnsList_1, valuesColumnsList_2 = [[i for i in list(groupPairDF_1) if i != indexColName], [j for j in list(groupPairDF_2) if j != indexColName]]
        comparedIndexesMethodsDictsLists2DArray.append([groupPairList, multi_core_queue(mp_index2pvals, groupMergedDF.index.tolist())])

    def mp_assemble_pvals_tables_dict(input_list):
        group_pair_list, compared_indexes_methods_dicts_list = input_list
        # Assemble dataframe
        statistical_tests_dfs_dict = {i: pd.DataFrame() for i in values2pvals_dict([], [])}
        for compared_indexes_method_dict in compared_indexes_methods_dicts_list:
            for statistical_test in statistical_tests_dfs_dict:
                if len(statistical_tests_dfs_dict[statistical_test]) == 0:
                    statistical_tests_dfs_dict[statistical_test] = pd.DataFrame([compared_indexes_method_dict[statistical_test]])
                else:
                    statistical_tests_dfs_dict[statistical_test] = statistical_tests_dfs_dict[statistical_test].append(compared_indexes_method_dict[statistical_test], ignore_index=True)
        print("Completed combining p-values methods data")
        for statistical_test in statistical_tests_dfs_dict:
            statistical_tests_df = statistical_tests_dfs_dict[statistical_test]
            statistical_tests_df['is_rejected'] = statistical_tests_df['p-value'] < globalAlpha
            if not noFilterBool:
                statistical_tests_df = statistical_tests_df.loc[statistical_tests_df['p-value'] < globalAlpha]
            statistical_tests_df.to_csv(outputDir + "comparison_by_" + valueColName + "_" + "_vs_".join(str(i) for i in group_pair_list) + "_" + statistical_test + ".tsv", sep='\t', index=False)
            statistical_tests_dfs_dict[statistical_test] = statistical_tests_df
        global groupDictPairsDict
        groupDictPairsDict["_vs_".join(str(i) for i in group_pair_list)] = statistical_tests_dfs_dict


    multi_core_queue(mp_assemble_pvals_tables_dict, comparedIndexesMethodsDictsLists2DArray)

    methodsDictDFsDict = unpack_methods()
    print("Adjusting p-values")

    for methodsDictDFsDictKey in methodsDictDFsDict:
        valuesMethodDF = pd.DataFrame()
        methodsDFsDict = methodsDictDFsDict[methodsDictDFsDictKey]
        for methodsDFsDictKey in methodsDFsDict:
            methodDF = methodsDFsDict[methodsDFsDictKey].rename(columns={i: i + '_' + methodsDFsDictKey for i in list(methodsDFsDict[methodsDFsDictKey]) if i != indexColName})
            if len(valuesMethodDF) == 0:
                valuesMethodDF = methodDF
            else:
                valuesMethodDF = pd.merge(valuesMethodDF, methodDF, on=indexColName, how='outer')
        valuesSubDF = valuesMethodDF.loc[:, [i for i in list(valuesMethodDF) if i == indexColName or 'p-value' in i]].set_index(indexColName).fillna(1)
        valuesSubDFIndexList = valuesSubDF.index.tolist()
        if len(valuesSubDFIndexList) == 0:
            print("None values to analyze left for method", methodsDictDFsDictKey)
        else:
            adjustedValuesDictsList = multi_core_queue(mp_adjust_method_pvals, valuesSubDFIndexList)

            def mp_assemble_adjusted_tables_dict(adjusted_values_dict):
                adjusted_values_method_df = pd.DataFrame()
                if len(adjusted_values_method_df) == 0:
                    adjusted_values_method_df = pd.DataFrame([dict2pd_series(adjusted_values_dict)])
                else:
                    adjusted_values_method_df = adjusted_values_method_df.append(dict2pd_series(adjusted_values_dict), ignore_index=True)
                adjusted_total_df = pd.merge(valuesMethodDF.loc[:, [indexColName] + [i for i in list(valuesMethodDF) if any(j in i for j in["prevalent", "p-value", "is_rejected"])]], adjusted_values_method_df, on=indexColName, how='left')
                adjusted_total_df['non_zero_hypothesis_count'] = adjusted_total_df.loc[:, [i for i in list(adjusted_total_df) if i != indexColName]].applymap(lambda var: str(var) == 'True').sum(axis=1)
                adjusted_total_df.sort_values(by='non_zero_hypothesis_count', ascending=False).to_csv(outputDir + '_'.join(["adjusted_p-values_for_groups", *sorted([i for i in groupsDict]), "for_method", methodsDictDFsDictKey, "by", valueColName, "for_alpha", str(globalAlpha) + ".tsv"]), sep='\t', index=False)
                print("Adjusted p-values for", methodsDictDFsDictKey)

            multi_core_queue(mp_assemble_adjusted_tables_dict, adjustedValuesDictsList)

    print("Completed. Check the directory:", outputDir[:-1])
