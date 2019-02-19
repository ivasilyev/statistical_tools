#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import numpy as np
import scipy.stats
import os
import re
import multiprocessing
import statsmodels.sandbox.stats.multicomp
import statistics


def parse_args():
    starting_parser = argparse.ArgumentParser(description="This tool counts statistics from linker file containing two tab-delimited columns without a header. The first column must contain file path, the second - group name.")
    starting_parser.add_argument("-g", "--groupdata", required=True,
                                 help="Text file without a header containing two tab-delimited columns: file path and group name"
                                      "If suffix and prefix are provided, combines path as <prefix><sample_name><suffix>")
    starting_parser.add_argument("-p", "--prefix", default="",
                                 help="(Optional) File path prefix")
    starting_parser.add_argument("-s", "--suffix", default="",
                                 help="(Optional) File path suffix")
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
    return namespace.groupdata, namespace.prefix, namespace.suffix, namespace.index, namespace.value, namespace.alpha, namespace.no_filter, namespace.output


def multi_core_queue(function_to_parallelize, queue):
    pool = multiprocessing.Pool()
    output = pool.map(function_to_parallelize, queue)
    pool.close()
    pool.join()
    return output


class ComparingGroupParser:
    def __init__(self, name):
        self.name = name
        self.pivot_df = self.group_name2df()
        self.samples_list = list(self.pivot_df)
        self.width = len(self.samples_list)
    @staticmethod
    def read_table_slice(table_file_name):
        df = pd.read_table(table_file_name, sep='\t', header=0, engine='python')
        columns = [indexColName, valueColName]
        output_dataframe = df.loc[:, columns]
        if any(i not in list(df) for i in columns):
            raise ValueError("Not all columns \"{}\" are present in table columns \"{}\"".format(columns, list(df)))
        return output_dataframe.rename(columns={valueColName: table_file_name})
    def group_name2df(self):
        group_dataframe = pd.DataFrame()
        group_files_list = groupDataDF.loc[groupDataDF["group_name"] == self.name, "file_name"].values.tolist()
        if len(group_files_list) < 2:
            print("Excluded group due to the files number smaller than 2:", self.name)
        else:
            for file_name in group_files_list:
                try:
                    if len(group_dataframe) != 0:
                        group_dataframe = pd.merge(group_dataframe, self.read_table_slice(file_name), on=indexColName, how='outer')
                    else:
                        group_dataframe = self.read_table_slice(file_name)
                except FileNotFoundError:
                    print("Not found file in group", self.name, ":", file_name)
            pivot_file_name = outputDir + '_'.join(["pivot", "by", valueColName, self.name + ".tsv"])
            group_dataframe.to_csv(pivot_file_name, sep='\t', index=False)
            print("Dumped the group '" + self.name + "' pivot:", pivot_file_name)
        return group_dataframe.set_index(indexColName).fillna(0)
    @staticmethod
    def dict2pd_series(dictionary):
        output = pd.Series()
        for key in dictionary:
            output.at[key] = dictionary[key]
        return output
    @staticmethod
    def list2base_stats(input_list):
        functions_dict = {"total_values": lambda x: len(input_list),
                          "sum": sum,
                          "mean": statistics.mean,
                          "median": statistics.median,
                          "sd": statistics.stdev,
                          "variance": statistics.variance}
        out = {k: functions_dict[k](input_list) for k in functions_dict}
        out["non-zero_values"] = len([i for i in input_list if float(i) > 0])
        out["non-zero_values_percentage"] = float(out["non-zero_values"]) / float(out["total_values"])
        return out
    def __len__(self):
        return len(self.pivot_df)
    def get_index_list(self):
        return self.pivot_df.index.values.tolist()
    def get_values_series(self, index_string):
        try:
            input_dict = self.pivot_df.loc[index_string].to_dict()
            output_dict = {}
            for k in input_dict:
                try:
                    tmp = float(input_dict[k])
                    del tmp
                    output_dict.update({k: input_dict[k]})
                except ValueError:
                    output_dict.update({k: 0})
        except KeyError:
            output_dict = {i: 0 for i in list(self.pivot_df)}
        return pd.Series(self.dict2pd_series(output_dict), name=index_string)
    def get_values_list(self, index_string):
        return self.get_values_series(index_string).values.tolist()
    def mp_count_base_stats(self, index_string):
        return pd.Series(self.dict2pd_series(self.list2base_stats(self.get_values_list(index_string))), name=index_string)
    def create_base_stats_df(self):
        df = pd.concat(multi_core_queue(self.mp_count_base_stats, indexList), axis=1, sort=False).transpose()
        df.index.names = [indexColName]
        return df
    def recreate_pivot_df(self):
        df = pd.concat(multi_core_queue(self.get_values_series, indexList), axis=1, sort=False).transpose()
        df.index.names = [indexColName]
        return df


class GroupComparator:
    def __init__(self, group_pair_list):
        self.group_pair_list = group_pair_list
        self.group_names_list = [i.name for i in group_pair_list]
        self.name = "_vs_".join(self.group_names_list)
        self.functions_dict = self.create_functions_dict()
        self.dfs_dict = self.create_pvals_dfs_dict()
    @staticmethod
    def create_functions_dict():
        return {"t-test": scipy.stats.ttest_ind, "u-test": scipy.stats.mannwhitneyu, "wilcoxon": scipy.stats.ranksums, "h-test": scipy.stats.kruskal}
    @staticmethod
    def create_comparisons_list():
        """
        Given a flat list of n items, returns one-side comparison tuples with a 'n * (n - 1) / 2' number
        :return out:
        List of tuples
        """
        out = []
        for a in groupsParsedList:
            for b in groupsParsedList:
                if a is not b and not (a, b) in out and not (b, a) in out:
                    out.append((a, b))
        return out
    def mp_count_single_pvals(self, index_string):
        """
        The function for per-line p-values calculation
        :param index_string:
        Single value of 'indexColName' column
        :return output_dict:
        Dictionary <p-value counting method>: Pandas Series object supplied with 'statistic' and 'p-value' fields and 'name' attribute
        """
        group_1 = self.group_pair_list[0]
        group_2 = self.group_pair_list[-1]
        list_1 = group_1.get_values_list(index_string)
        list_2 = group_2.get_values_list(index_string)
        output_dict = {}
        for function_name in self.functions_dict:
            if np.count_nonzero(list_1) == 0 and np.count_nonzero(list_2) == 0:
                function_dict = {"statistic": np.nan, "p-value": np.nan, "is_rejected": False}
            else:
                function_dict = {i: j for (i, j) in zip(["statistic", "p-value"], self.functions_dict[function_name](list_1, list_2))}
                function_dict["is_rejected"] = function_dict["p-value"] < globalAlpha
            # True and False are converted into 1 and 0 otherwise
            function_dict["is_rejected"] = str(function_dict["is_rejected"])
            output_dict.update({function_name: pd.Series(ComparingGroupParser.dict2pd_series(function_dict), name=index_string)})
        return output_dict
    def create_pvals_dfs_dict(self):
        """
        The wrapper for 'mp_count_single_pvals' method
        :return output_dict:
        Dictionary <p-value counting method>: Pandas DataFrame object containing 'statistic' and 'p-value' columns and 'indexColName' index
        """
        pvals_dicts_list = multi_core_queue(self.mp_count_single_pvals, indexList)
        output_dict = {i: [] for i in self.functions_dict}
        for pvals_dict in pvals_dicts_list:
            for function_name in output_dict:
                output_dict[function_name].append(pvals_dict[function_name])
        for function_name in output_dict:
            df = pd.concat(output_dict[function_name], axis=1, sort=False).transpose()
            df.index.names = [indexColName]
            output_dict[function_name] = df
        return output_dict
    def mp_get_prevalent(self, series):
        d = series.to_dict()
        key_1, key_2 = d
        if d[key_1] > d[key_2]:
            output_dict = {"prevalent": key_1}
        elif d[key_1] < d[key_2]:
            output_dict = {"prevalent": key_2}
        else:
            output_dict = {"prevalent": np.nan}
        # Also np.nan == np.nan -> False
        return pd.Series(ComparingGroupParser.dict2pd_series(output_dict), name=series.name)
    def create_prevalents_df(self):
        group_1 = self.group_pair_list[0]
        group_2 = self.group_pair_list[-1]
        prevalent_function_name = "sum"
        df = pd.concat([i.loc[:, [prevalent_function_name]].rename(columns={c: j for c in list(i)}) for (i, j) in [(o.create_base_stats_df(), o.name) for o in [group_1, group_2]]], axis=1, sort=False)
        queue = [df.loc[i, :] for i in df.index.tolist()]
        output_df = pd.concat(multi_core_queue(self.mp_get_prevalent, queue), axis=1, sort=False).transpose()
        output_df.index.names = [indexColName]
        return output_df


class MultipleComparator:
    def __init__(self):
        self.single_pvals_functions_list = list(GroupComparator.create_functions_dict())
        self.single_pvals_dfs_dict = self.create_single_pvals_dfs_dict()
        self.multi_pvals_functions_list = self.create_multi_pvals_functions_list()
        self.multi_pvals_dfs_dict = self.create_multi_pvals_dfs_dict()
    def create_single_pvals_dfs_dict(self):
        """
        The function for per-method merging p-values tables.
        Given comparisons list, each object supplied with dictionary {<p-value counting method>: Pandas DataFrame object containing 'statistic' and 'p-value' columns and 'indexColName' index} in 'dfs_dict' attribute
        :return:
        Dictionary <p-value counting method>: Pandas DataFrame object containing '<group 1>_vs_<group 2>_p-value' columns renamed corresponding 'name' attribute of comparison object and 'indexColName' index
        """
        single_pvals_dfs_dict = {i: [] for i in self.single_pvals_functions_list}
        for single_pvals_method in single_pvals_dfs_dict:
            for comparison in groupComparisonsList:
                df = comparison.dfs_dict[single_pvals_method].loc[:, ["p-value"]]
                df.rename(columns={"p-value": comparison.name}, inplace=True)
                single_pvals_dfs_dict[single_pvals_method].append(df)
        # Replacing 'NA' with 1
        return {k: pd.concat(single_pvals_dfs_dict[k], axis=1, sort=False).fillna(1) for k in single_pvals_dfs_dict}
    @staticmethod
    def create_multi_pvals_functions_list():
        return [  # Description: http://www.statsmodels.org/devel/generated/statsmodels.sandbox.stats.multicomp.multipletests.html#statsmodels.sandbox.stats.multicomp.multipletests
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
            ]
    def mp_count_multi_pvals(self, input_dict):
        """
        The function for per-line p-values correction
        :param input_dict:
        Dictionary with 'index_string', 'pvals_list', 'colnames_list' keys
        :return:
        Pandas Series object containing '<group 1>_vs_<group 2>_is_rejected_by_<adjustment_method>' logic values columns,
        '<group 1>_vs_<group 2>_corrected_for_<adjustment_method>' numeric column, alpha values columns and 'name' attribute
        """
        output_dict = {}
        alphac_sidak = alphac_bonf = globalAlpha
        for function_name in self.multi_pvals_functions_list:
            rejects_list, pvals_corrected_list, alphac_sidak, alphac_bonf = statsmodels.sandbox.stats.multicomp.multipletests(pvals=input_dict["pvals_list"], method=function_name, alpha=globalAlpha)
            # True and False are converted into 1 and 0 otherwise
            output_dict.update({i + "_is_rejected_by_" + function_name: j for (i, j) in zip(input_dict["colnames_list"], [str(k) for k in rejects_list])})
            output_dict.update({i + "_corrected_for_" + function_name: j for (i, j) in zip(input_dict["colnames_list"], pvals_corrected_list)})
        output_dict.update({"raw_alpha": globalAlpha, "alpha_corrected_for_Sidak_method": alphac_sidak, "alpha_corrected_for_Bonferroni_method": alphac_bonf})
        return pd.Series(ComparingGroupParser.dict2pd_series(output_dict), name=input_dict["index_string"])
    def create_multi_pvals_dfs_dict(self):
        """
        The wrapper for 'mp_count_multi_pvals' method
        :return output_dict:
        Dictionary <p-value counting method>: Pandas DataFrame object with 'indexColName' index
        For n comparisons and dataframe shall contain:
            - n * m '<group 1>_vs_<group 2>_is_rejected_by_<adjustment_method>' columns;
            - n * m '<group 1>_vs_<group 2>_corrected_for_<adjustment_method>' columns;
            - 3 columns containing alpha values;
        where m = len(self.multi_pvals_functions_list)
        E.g. for 15 comparisons and 10 correction methods dataframe shall contain 303 columns
        """
        output_dict = {}
        for single_pvals_method in self.single_pvals_dfs_dict:
            df = self.single_pvals_dfs_dict[single_pvals_method]
            queue = [{"index_string": i, "pvals_list": df.loc[i].values.tolist(), "colnames_list": list(df)} for i in df.index.tolist()]
            output_dict.update({single_pvals_method: pd.concat(multi_core_queue(self.mp_count_multi_pvals, queue), axis=1, sort=False).transpose()})
        return output_dict


class Exporter:
    @staticmethod
    def prepare_string(s: str):
        return re.sub("@+", "_", s)
    @staticmethod
    def merge_pivot_df():
        output_list = []
        for group in groupsParsedList:
            df = group.recreate_pivot_df()
            output_list.append(df.rename(columns={i: "{}@{}".format(group.name, Exporter.prepare_string(i)) for i in list(df)}))
        out = pd.concat(output_list, axis=1, sort=False)
        out.index.names = [indexColName]
        return out
    @staticmethod
    def merge_base_stats_df():
        output_list = []
        for group in groupsParsedList:
            df = group.create_base_stats_df()
            output_list.append(df.rename(columns={i: "{}@{}".format(group.name, Exporter.prepare_string(i)) for i in list(df)}))
        out = pd.concat(output_list, axis=1, sort=False)
        out.index.names = [indexColName]
        return out
    @staticmethod
    def merge_single_pvals_df():
        output_list = []
        for comparison in groupComparisonsList:
            for single_pvals_method in comparison.dfs_dict:
                df = comparison.dfs_dict[single_pvals_method]
                output_list.append(df.rename(columns={i: "_".join([comparison.name, single_pvals_method, i]) for i in list(df)}))
        out = pd.concat(output_list, axis=1, sort=False)
        out.index.names = [indexColName]
        return out
    @staticmethod
    def merge_prevalents_df():
        output_list = []
        for comparison in groupComparisonsList:
            df = comparison.create_prevalents_df()
            output_list.append(df.rename(columns={i: "_".join([comparison.name, i]) for i in list(df)}))
        out = pd.concat(output_list, axis=1, sort=False)
        out.index.names = [indexColName]
        return out
    @staticmethod
    def merge_multi_pvals_df():
        output_list = []
        for single_pvals_method in correctionObject.multi_pvals_dfs_dict:
            df = correctionObject.multi_pvals_dfs_dict[single_pvals_method]
            output_list.append(df.rename(columns={i: i + "_for_" + single_pvals_method for i in list(df)}))
        out = pd.concat(output_list, axis=1, sort=False)
        out.index.names = [indexColName]
        return out
    def merge_everything(self):
        output_list = [i() for i in [self.merge_pivot_df, self.merge_base_stats_df, self.merge_single_pvals_df, self.merge_prevalents_df, self.merge_multi_pvals_df]]
        out = pd.concat(output_list, axis=1, sort=False)
        out.index.names = [indexColName]
        return out
    @staticmethod
    def _mp_count_rejected_null_hypothesis(series):
        output_dict = {"null_hypothesis_rejections_counter": len([i for i in series.values.tolist() if i == 'True'])}
        return pd.Series(ComparingGroupParser.dict2pd_series(output_dict), name=series.name)
    def concat_rejected_counter(self, input_df):
        queue = [input_df.loc[i, :] for i in input_df.index.tolist()]
        df = pd.concat(multi_core_queue(self._mp_count_rejected_null_hypothesis, queue), axis=1, sort=False).transpose()
        output_df = pd.concat([input_df, df], axis=1).sort_values(by="null_hypothesis_rejections_counter", ascending=False)
        return output_df


if __name__ == '__main__':
    inputGroupDataFileName, mainPrefix, mainSuffix, indexColName, valueColName, globalAlpha, noFilterBool, outputDir = parse_namespace()
    outputDir = ends_with_slash(outputDir)
    os.makedirs(outputDir, exist_ok=True)

    print("Parsing groups data")
    groupDataDF = pd.read_table(inputGroupDataFileName, sep='\t', header='infer', names=["sample_name", "group_name"], engine='python')
    groupDataDF["file_name"] = mainPrefix.strip() + groupDataDF["sample_name"] + mainSuffix.strip()
    filesToValidate = groupDataDF["file_name"].values.tolist()
    print("Validating files")
    groupDataDF = groupDataDF.loc[groupDataDF["file_name"].map(lambda x: os.path.isfile(x))]
    if len(groupDataDF) == 0:
        raise ValueError("No valid files found!")
    invalidFiles = [i for i in filesToValidate if i not in groupDataDF["file_name"].values.tolist()]
    if len(invalidFiles) > 0:
        print("Warning! The files have not been found: \n{}\n".format(("\n").join(invalidFiles)))

    groupsRawList = sorted(groupDataDF["group_name"].unique().tolist())
    groupsParsedList = [ComparingGroupParser(i) for i in groupsRawList]
    groupsParsedList = [i for i in groupsParsedList if len(i) > 0 and i.width > 0]
    print("Groups to process: " + "'" + "', '".join([i.name for i in groupsParsedList]) + "'")
    indexList = sorted(list(set([k for j in [i.get_index_list() for i in groupsParsedList] for k in j])))
    groupsTuplesList = GroupComparator.create_comparisons_list()
    groupComparisonsList = [GroupComparator(i) for i in groupsTuplesList]
    correctionObject = MultipleComparator()
    exporter = Exporter()
    exporterDF = exporter.merge_everything()
    outputDF = exporter.concat_rejected_counter(exporterDF).reset_index()
    outputDF.loc[:, [i for i in list(outputDF) if len(i.strip()) > 0]].to_csv(
        outputDir + '_'.join([i.name for i in groupsParsedList]) + "_total_dataframe.tsv", sep='\t', index=False)
    print("Completed")
