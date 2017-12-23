#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys
import re
import subprocess
import pandas as pd
from groupdata2statistics import ends_with_slash, is_path_exists, file_to_list, find_existing_files, list2pairwise_tuples, dict2pd_series
from collections import Counter
from matplotlib import pyplot as plt


def parse_args():
    starting_parser = argparse.ArgumentParser(description="Given two tab-delimited tables: the first table contains aliases, path prefixes and suffixes, the second contains abstract sample names and corresponding group IDs")
    starting_parser.add_argument("-s", "--suffixes", required=True,
                                 help="Table without a header containing 3 columns: alias, path prefixes and suffixes for each sample")
    starting_parser.add_argument("-g", "--groupdata", required=True,
                                 help="Text file without a header containing two tab-delimited columns: file path and group name")
    starting_parser.add_argument("-f", "--filter", default=None,
                                 help="Text file containing words to filter one per line")
    starting_parser.add_argument("-t", "--top", type=int, default=10,
                                 help="(Optional) Number of top rows to process under own names, all other would be processed as 'other'")
    starting_parser.add_argument("-c", "--control", default=[], nargs='+',
                                 help="(Optional) First and second control to compose")
    starting_parser.add_argument("-m", "--method", default='u-test',
                                 help="(Optional) Method to compose, 'u-test' by default")
    starting_parser.add_argument("-r", "--fdr", default='fdr_bh',
                                 help="(Optional) Method to adjust multiple comparisons, 'fdr_bh' by default. See http://www.statsmodels.org/dev/generated/statsmodels.sandbox.stats.multicomp.multipletests.html")
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


def parse_namespace():
    namespace = parse_args()
    is_path_exists(namespace.output)
    return namespace.suffixes, namespace.groupdata, namespace.filter, namespace.top, namespace.control, namespace.method, namespace.fdr, namespace.index, namespace.value, namespace.alpha, namespace.no_filter, ends_with_slash(namespace.output)


def file_append(string, file_to_append):
    file = open(file_to_append, 'a+')
    file.write(string)
    file.close()


def external_route(input_direction_list, output_direction):
    process = subprocess.Popen(input_direction_list, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (output, error) = process.communicate()
    process.wait()
    if error:
        print(error)
    if not output_direction:
        return re.sub('[\r\n]', '', output.decode("utf-8"))
    else:
        file_append(output.decode("utf-8"), output_direction)


def external_route_wrapper(args):
    print(external_route(*args))


def counter_df2plot(df, plot_title, path):
    if not os.path.isdir(path):
        plot_title = path.split('/')[-1] + "_" + plot_title
        path = "/".join(path.split('/')[:-1])
        is_path_exists(path)
    path = ends_with_slash(path)
    is_path_exists(path)
    df.to_csv(path + plot_title + ".tsv", sep='\t', index=False)
    try:
        df = df.set_index("keyword")
        # Remove zero-only columns
        df = df.loc[:, [i for i in list(df) if df.loc[:, i].sum() != 0]]
    except KeyError:
        # The dataframe 'adjusted_rows' is not a real table, but a dump
        if "adjusted_rows" not in plot_title:
            print("Cannot make plots for table:", path + plot_title + ".tsv")
        return
    fig = plt.figure(figsize=(12, 5))
    rotated_pairs_df = df.apply(lambda x: 100 * x / float(x.sum())).transpose()
    ax1 = rotated_pairs_df.plot(kind='bar', stacked=True, label='Percentage', width=.99)
    plt.xticks(rotation=0)
    ax1.set_xlabel("Comparison pairs")
    ax1.set_ylabel("Percentage", color='black')
    ax1.tick_params('y', colors='black')
    ax1.set_ylim([0, 100])
    labels = []
    for j in rotated_pairs_df.columns:
        for i in rotated_pairs_df.index:
            if rotated_pairs_df.loc[i][j] > 0:
                labels.append(j)
            else:
                labels.append("")
    patches = ax1.patches
    for label, rect in zip(labels, patches):
        width = rect.get_width()
        if width > 0:
            ax1.text(rect.get_x() + width / 2., rect.get_y() + rect.get_height() / 2., label, ha='center', va='center', size='xx-small')
    for label, rect in zip(df.sum(), patches):
        ax1.text(rect.get_x() + rect.get_width() / 2., 101, str(int(label)), ha='center', va='bottom')
    ax1.legend(loc=1).set_visible(False)
    plt.title(plot_title, y=1.05)
    plt.savefig(path + "/plots/" + plot_title + ".png", dpi=600)
    plt.cla()
    plt.clf()
    plt.close(fig)


def controls2pairing_dict(input_list):
    output_dict = {}
    # Will process only 2 control groups
    control_groups_list = [i for i in controlWordsList if i in input_list][:2]
    test_groups_list = [i for i in input_list if i not in controlWordsList]
    if len(input_list) > 3:
        # Index-based control groups attachment
        output_dict["near"] = [control_groups_list] + [[i, test_groups_list[test_groups_list.index(i) + 1]] for i in test_groups_list if test_groups_list.index(i) != (len(test_groups_list) - 1)]
        if len(control_groups_list) > 0:
            output_dict["on_first_control"] = [[control_groups_list[0], j] for j in test_groups_list]
        if len(control_groups_list) > 1:
            output_dict["on_second_control"] = [[control_groups_list[1], j] for j in test_groups_list]
    else:
        output_dict["'free-for-all"] = list2pairwise_tuples(input_list)
    return output_dict


def pair_list2omni_indexes_list(input_list):
    input_string = "_vs_".join(str(i) for i in input_list)
    reversed_string = "_vs_".join(str(i) for i in input_list[::-1])
    if input_string in alias_pair_strings_list:
        return {input_string: alias_df.loc[alias_df['_'.join(["p-value", input_string, "is_rejected_by", adjustmentMethod])] == True][indexColName].values.tolist()}
    if reversed_string in alias_pair_strings_list:
        return {input_string: alias_df.loc[alias_df['_'.join(["p-value", reversed_string, "is_rejected_by", adjustmentMethod])] == True][indexColName].values.tolist()}
    else:
        raise ValueError("Not in the main dataframe: " + '_'.join(["p-value", input_string, "is_rejected_by", adjustmentMethod]))


def pair_list2left_indexes_list(input_list):
    input_string = "_vs_".join(str(i) for i in input_list)
    reversed_string = "_vs_".join(str(i) for i in input_list[::-1])
    if input_string in alias_pair_strings_list:
        return {input_string: alias_df.loc[(alias_df['_'.join(["p-value", input_string, "is_rejected_by", adjustmentMethod])] == True) & (alias_df["prevalent_" + input_string] == input_list[0])][indexColName].values.tolist()}
    if reversed_string in alias_pair_strings_list:
        return {input_string: alias_df.loc[(alias_df['_'.join(["p-value", reversed_string, "is_rejected_by", adjustmentMethod])] == True) & (alias_df["prevalent_" + reversed_string] == input_list[0])][indexColName].values.tolist()}
    else:
        raise ValueError("Not in the main dataframe: " + '_'.join(["p-value", input_string, "is_rejected_by", adjustmentMethod]))


def pair_list2right_indexes_list(input_list):
    input_string = "_vs_".join(str(i) for i in input_list)
    reversed_string = "_vs_".join(str(i) for i in input_list[::-1])
    if input_string in alias_pair_strings_list:
        return {input_string: alias_df.loc[(alias_df['_'.join(["p-value", input_string, "is_rejected_by", adjustmentMethod])] == True) & (alias_df["prevalent_" + input_string] == input_list[1])][indexColName].values.tolist()}
    if reversed_string in alias_pair_strings_list:
        return {input_string: alias_df.loc[(alias_df['_'.join(["p-value", reversed_string, "is_rejected_by", adjustmentMethod])] == True) & (alias_df["prevalent_" + reversed_string] == input_list[1])][indexColName].values.tolist()}
    else:
        raise ValueError("Not in the main dataframe: " + '_'.join(["p-value", input_string, "is_rejected_by", adjustmentMethod]))


def digest_words_list(input_list):
    filtered_list = ["_".join([j for j in re.findall('[A-Za-z]{4,}', i) if j.lower() not in filterWordsList]) for i in input_list]
    keyword_lists_dict = {'species': [j[0] for j in [re.findall('[A-Z]{1}[a-z]{4,}_[a-z]{4,}', i) for i in filtered_list] if len(j) > 0],
                          'genera': [j[0] for j in [re.findall('[A-Z]{1}[a-z]{4,}', i) for i in filtered_list] if len(j) > 0],
                          'common_words': [k for j in [re.findall('[a-z]{4,}', i.lower()) for i in input_list] for k in j]}
    for keyword_list_key in keyword_lists_dict:
        keyword_list = keyword_lists_dict[keyword_list_key]
        missing_values_number = len(input_list) - len(keyword_list)
        if missing_values_number > 0:
            keyword_list = keyword_list + ["unknown",] * missing_values_number
        keyword_lists_dict[keyword_list_key] = keyword_list
    return keyword_lists_dict


def list4count2df(input_list, input_string, only_top_bool):
    input_dict = Counter(input_list)
    if len(input_dict) > 0:
        df = pd.DataFrame.from_dict(Counter(input_list), orient='index').reset_index()
        df = df.rename(columns={i: j for i, j in zip(list(df), ["keyword", input_string])}).sort_values(by=input_string, ascending=False)
        if only_top_bool:
            values_not_fit = len(df) - keyWordTopRowsNumber
            if values_not_fit > 0:
                top_df = pd.DataFrame([dict2pd_series({"keyword": "other", input_string: df.tail(values_not_fit)[input_string].sum()})])
                top_df = top_df.append(df.head(keyWordTopRowsNumber), ignore_index=True)
                return top_df
    else:
        df = pd.DataFrame(columns=["keyword", input_string])
    return df


def pd_merge_dfs_list(input_list, index_string):
    ldf = pd.DataFrame()
    for rdf in input_list:
        if len(ldf) == 0:
            ldf = rdf.copy()
        else:
            ldf = pd.merge(ldf, rdf, on=index_string, how="outer")
    return ldf


def process_keywords_lists_dict(input_dict):
    base_ldf = pd.DataFrame()
    digest_dfs_dict = digest_words_list([])
    top_dfs_dict = digest_words_list([])
    for pair_string in input_dict:
        base_rdf = pd.DataFrame(input_dict[pair_string], columns=[pair_string]).sort_values(pair_string)
        if len(base_ldf) == 0:
            base_ldf = base_rdf
        else:
            base_ldf = pd.concat([base_ldf, base_rdf], axis=1)
        digest_dict = digest_words_list(input_dict[pair_string])
        for digest_method in digest_dict:
            digest_dfs_dict[digest_method].append(list4count2df(digest_dict[digest_method], pair_string, False))
            top_dfs_dict[digest_method].append(list4count2df(digest_dict[digest_method], pair_string, True))
    digest_dfs_dict = {i: pd_merge_dfs_list(digest_dfs_dict[i], "keyword").set_index('keyword').fillna(0).astype('int').reset_index() for i in digest_dfs_dict}
    top_dfs_dict = {i: pd_merge_dfs_list(top_dfs_dict[i], "keyword").set_index('keyword').fillna(0).astype('int').reset_index() for i in top_dfs_dict}
    return base_ldf, digest_dfs_dict, top_dfs_dict


def unpack_dicts():
    output_dict = {}
    output_dict.update({"omni_" + i: omni_keywords_dict[i] for i in omni_keywords_dict})
    output_dict.update({"left_prevalent_" + i: left_keywords_dict[i] for i in left_keywords_dict})
    output_dict.update({"right_prevalent_" + i: right_keywords_dict[i] for i in right_keywords_dict})
    return output_dict


def export_lists_dict(input_dict, output_prefix):
    base_ldf, digest_dfs_dict, top_dfs_dict = process_keywords_lists_dict(input_dict)
    output_dict = {"adjusted_rows": base_ldf}
    output_dict.update({"raw_counter_" + i: digest_dfs_dict[i] for i in digest_dfs_dict})
    output_dict.update({"top_" + str(keyWordTopRowsNumber) + "_counter_" + i: top_dfs_dict[i] for i in top_dfs_dict})
    for key in output_dict:
        counter_df2plot(output_dict[key], key, output_prefix)


def mp_dicts_export(key):
    is_path_exists(outputDir + "pvals/" + alias_string + "/counters/plots/")
    export_lists_dict(mainPairingDict[key], outputDir + "pvals/" + alias_string + "/counters/" + key)


def multi_core_queue(function_to_parallelize, queue):
    import multiprocessing
    pool = multiprocessing.Pool()
    pool.map(function_to_parallelize, queue)
    pool.close()
    pool.join()


if __name__ == '__main__':
    suffixesFileName, groupDataFileName, filterWordsFileName, keyWordTopRowsNumber, controlWordsList, composeMethod, adjustmentMethod, indexColName, valueColName, inputAlpha, noFilterBool, outputDir = parse_namespace()
    scriptDir = ends_with_slash(os.path.dirname(os.path.realpath(sys.argv[0])))
    os.chdir(scriptDir)
    suffixesDF = pd.read_table(suffixesFileName, sep='\t', header='infer', names=["alias", "prefix", "suffix"]).set_index("alias")
    groupDataRawDF = pd.read_table(groupDataFileName, sep='\t', header='infer', names=["sample_name", "group_name"])
    filterWordsList = [i.lower() for i in file_to_list(filterWordsFileName)]
    for alias_string in suffixesDF.index.values:
        prefix_string, suffix_string = suffixesDF.ix[alias_string].tolist()
        groupdata_processed_df = groupDataRawDF.copy()
        groupdata_processed_df["file_name"] = groupdata_processed_df.loc[:, "sample_name"].apply(lambda x: prefix_string + x + suffix_string)
        existing_files_list = find_existing_files(groupdata_processed_df.loc[:, "file_name"].tolist())
        if len(existing_files_list) == 0:
            print("Excluded alias:", alias_string)
            continue
        groupdata_processed_file_name = outputDir + alias_string + '.groupdata'
        groupdata_processed_df = groupdata_processed_df.loc[groupdata_processed_df["file_name"].isin(existing_files_list)]
        groupdata_processed_df.loc[:, ["file_name", "group_name"]].to_csv(groupdata_processed_file_name, sep='\t', header=False, index=False)
        cmd_list = ["python3", scriptDir + "groupdata2statistics.py", "-g", groupdata_processed_file_name, "-i", indexColName, "-v", valueColName, "-a", str(inputAlpha), "-o", outputDir + "pvals/" + alias_string]
        if noFilterBool:
            cmd_list = cmd_list[:-2] + ["-n"] + cmd_list[-2:]
        # P-values count
        external_route_wrapper((cmd_list, outputDir + alias_string + "_groupdata2statistics.log"))
        # Make plots
        pairing_dict = controls2pairing_dict(sorted(pd.unique(groupdata_processed_df['group_name']).tolist()))
        is_path_exists(outputDir + "pvals/" + alias_string + "/counters")
        alias_df = pd.read_table(subprocess.getoutput("ls -d " + outputDir + "pvals/" + alias_string + "/adjusted*" + composeMethod + "*").split('\n')[0], sep='\t', header=0)
        alias_pair_strings_list = [re.findall('p-value_(.*)_is_rejected_by_' + adjustmentMethod, i)[0] for i in list(alias_df) if '_is_rejected_by_' + adjustmentMethod in i]
        # Select data to visualize
        omni_keywords_dict = {i: {} for i in pairing_dict}
        left_keywords_dict = {i: {} for i in pairing_dict}
        right_keywords_dict = {i: {} for i in pairing_dict}
        for pairing_method in pairing_dict:
            for pairing_list in pairing_dict[pairing_method]:
                omni_keywords_dict[pairing_method].update(pair_list2omni_indexes_list(pairing_list))
                left_keywords_dict[pairing_method].update(pair_list2left_indexes_list(pairing_list))
                right_keywords_dict[pairing_method].update(pair_list2right_indexes_list(pairing_list))
        # Prepare main queue
        mainPairingDict = unpack_dicts()
        multi_core_queue(mp_dicts_export, mainPairingDict)
        print("Processed alias:", alias_string)
