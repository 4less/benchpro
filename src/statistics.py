from src.taxon_profile import Profile, default_ranks, Rank, RankTemp
from src.binary_statistics import BinaryStatistics
import scipy as sc
import sys
import math
from loguru import logger
from collections import namedtuple
from skbio.diversity import alpha
import numpy as np


class AbundanceStatistics:
    def __init__(self, name):
        self.name = name
        self.pearson_correlation_intersection = 0
        self.pearson_correlation_union = 0
        self.spearman_correlation_intersection = 0
        self.spearman_correlation_union = 0
        self.bray_curtis_intersection = 0
        self.bray_curtis_union = 0
        self.l2_intersection = 0
        self.l2_union = 0
        self.l2_log_intersection = 0
        self.l2_log_union = 0

class MiscStatistics:
    def __init__(self, name: str):
        self.name = name
        self.statistics = dict()

def get_abundance_dict(profile, rank):
    abundance_dict = dict()
    for row in profile.Rows(rank):
        name = row.GetLineage().Get(rank)
        if name == "":
            continue
        if name not in abundance_dict:
            abundance_dict[name] = 0.0
        abundance_dict[name] += row.GetAbundance()

    return abundance_dict


def bray_curtis(a, b):
    total_a = sum(a)
    total_b = sum(b)
    shared_a = sum(va for va, vb in zip(a, b) if va > 0 and vb > 0)
    shared_b = sum(vb for va, vb in zip(a, b) if va > 0 and vb > 0)

    bc = (2*min(shared_a, shared_b) / (total_a + total_b))

    # print("Bray: {}".format(bc))
    return bc

def euclidian(a,b):
    sum = 0
    for a1, a2 in zip(a, b):
        sum += (a1 - a2) ** 2
    return math.sqrt(sum)

def get_rank_statistics(gold_profile: Profile, prediction_profile: Profile,
                        ranks=default_ranks, workbook_logger=None, metadata=None):
    MiscStats = namedtuple('MiscStats', ['name', 'metrics'])

    logger.info(f"Get rank statistics for {gold_profile.name} and {prediction_profile.name}")
    binary_statistics_rank = dict()
    abundance_statistics_rank = dict()
    miscellaneous_statistics_rank = dict()


    logger.info(f"Process ranks: {ranks}")
    for rank in ranks:
        logger.info(f"Process rank: {rank}")

        stats = BinaryStatistics(prediction_profile.name)
        abundance = AbundanceStatistics(prediction_profile.name)
        misc_stats = MiscStatistics(prediction_profile.name)

        prediction_set = set(row.GetLineage().Get(rank) for row in prediction_profile.Rows(rank))
        gold_set = set(row.GetLineage().Get(rank) for row in gold_profile.Rows(rank))


        shared = prediction_set.intersection(gold_set)
        pred_only = prediction_set.difference(gold_set)
        gold_only = gold_set.difference(prediction_set)
        pred_no_shared = prediction_set.difference(shared)
        gold_no_shared = gold_set.difference(shared)

        if len(shared) == 0:
            logger.warning("No shared taxa between gold profile and prediction profile")
            logger.warning("prediction_set: {}".format(prediction_set))
            logger.warning("gold_set: {}" .format(gold_set))

        gold_set = set(row.GetLineage().Get(rank) for row in gold_profile.Rows())
        stats.tp = len(shared)
        stats.fp = len(pred_no_shared)
        stats.fn = len(gold_no_shared)



        ##########################################################################################
        # Get abundance dicts for abundance statistics

        gold_dict = get_abundance_dict(gold_profile, rank)
        prediction_dict = get_abundance_dict(prediction_profile, rank)
        gold_vec = list(gold_dict.values())
        pred_vec = list(prediction_dict.values())
        gold_vec_shared = [gold_dict[key] for key in shared]
        pred_vec_shared = [prediction_dict[key] for key in shared]

        logger.info("Shared taxa {}".format(len(gold_vec_shared)))
        logger.info("Shared taxa: GoldStd Sum {}, Prediction Sum {}".format(sum(gold_vec_shared), sum(pred_vec_shared)))

        gold_vec_union = gold_vec_shared + [value for key, value in gold_dict.items() if key not in prediction_dict]
        pred_vec_union = pred_vec_shared + [0.0 for key in gold_dict if key not in prediction_dict]
        pred_vec_union += [value for key, value in prediction_dict.items() if key not in gold_dict]
        gold_vec_union += [0.0 for key in prediction_dict if key not in gold_dict]

        logger.info("Union taxa {}".format(len(gold_vec_union)))
        logger.info("Union taxa abundances {}".format(gold_vec_union))
        logger.info("Union taxa: GoldStd Sum {}, Prediction Sum {}".format(sum(gold_vec_union), sum(pred_vec_union)))


        if workbook_logger:
            dataset = metadata[prediction_profile.name]['Dataset'] if metadata else "unspecified"
            tool = metadata[prediction_profile.name]['Tool'] if metadata else "unspecified"
            workbook_logger.add_sample(gold_profile.name, rank.value.name, dataset, tool, gold_dict, prediction_dict)


        ##########################################################################################
        # Pearson correlation
        abundance.pearson_correlation_intersection = sc.stats.pearsonr(gold_vec_shared, pred_vec_shared).statistic if len(gold_vec_shared) > 1 else 0
        abundance.pearson_correlation_union = sc.stats.pearsonr(gold_vec_union, pred_vec_union).statistic if len(gold_vec_shared) > 1 else 0

        logger.info("Pearson Correlation Shared {} Union {}".format(abundance.pearson_correlation_intersection, abundance.pearson_correlation_union))

        ##########################################################################################
        # Spearman correlation
        abundance.spearman_correlation_intersection = sc.stats.spearmanr(gold_vec_shared, pred_vec_shared).statistic if len(gold_vec_shared) > 1 else 0
        abundance.spearman_correlation_union = sc.stats.spearmanr(gold_vec_union, pred_vec_union).statistic if len(gold_vec_shared) > 1 else 0

        logger.info("Spearman Correlation Shared {} Union {}".format(abundance.spearman_correlation_intersection, abundance.spearman_correlation_union))

        ##########################################################################################
        # Bray-Curtis
        abundance.bray_curtis_intersection = bray_curtis(gold_vec_shared, pred_vec_shared)
        abundance.bray_curtis_union = bray_curtis(gold_vec_union, pred_vec_union)

        logger.info("Bray-Curtis Shared {} Union {}".format(abundance.bray_curtis_intersection, abundance.bray_curtis_union))

        ##########################################################################################
        # Euclidean (L2)
        abundance.l2_intersection = euclidian(gold_vec_shared, pred_vec_shared)
        abundance.l2_union = euclidian(gold_vec_union, pred_vec_union)

        logger.info("L2 Shared {} Union {}".format(abundance.l2_intersection, abundance.l2_union))

        ##########################################################################################
        # Euclidean (L2) (Math Log10 dont forget +1)
        abundance.l2_log_intersection = euclidian([a + 1 for a in gold_vec_shared], [a + 1 for a in pred_vec_shared])
        abundance.l2_log_union = euclidian([a+1 for a in gold_vec_union], [a+1 for a in pred_vec_union])

        logger.info("L2 Log Shared {} Union {}".format(abundance.l2_log_intersection, abundance.l2_log_union))

        ##########################################################################################
        # Shannon diversity
        misc_stats.statistics["ShannonDiversity"] = alpha.shannon(np.array(pred_vec))
        misc_stats.statistics["ShannonDiversityGold"] = alpha.shannon(np.array(gold_vec))
        misc_stats.statistics["ShannonDiversityDiff"] = abs(misc_stats.statistics["ShannonDiversityGold"]- misc_stats.statistics["ShannonDiversity"])
        misc_stats.statistics["ShannonDiversityTP"] = alpha.shannon(np.array(pred_vec_shared))
        misc_stats.statistics["ShannonDiversityGoldTP"] = alpha.shannon(np.array(gold_vec_shared))
        misc_stats.statistics["ShannonDiversityDiffTP"] = abs(misc_stats.statistics["ShannonDiversityGoldTP"]- misc_stats.statistics["ShannonDiversityTP"])
        misc_stats.statistics["RichnessTP"] = len(pred_vec_shared)
        misc_stats.statistics["RichnessGoldTP"] = len(gold_vec_shared)
        misc_stats.statistics["RichnessTPDiff"] = abs(len(gold_vec) - len(pred_vec_shared))

        logger.info("Shannon diversity Prediction {} vs Gold {}".format(misc_stats.statistics["ShannonDiversity"], misc_stats.statistics["ShannonDiversityGold"]))
        logger.info("Shannon diversity Prediction {} vs Gold {}".format(misc_stats.statistics["ShannonDiversityTP"], misc_stats.statistics["ShannonDiversityGoldTP"]))
        ##########################################################################################



        abundance_statistics_rank[rank] = abundance
        binary_statistics_rank[rank] = stats
        miscellaneous_statistics_rank[rank] = misc_stats

        try:
            stats.Validate()
        except:
            logger.exception(f"{gold_profile.name}\t{prediction_profile.name}\t{rank}")
            exit(9)

    return binary_statistics_rank, abundance_statistics_rank, miscellaneous_statistics_rank


def print_binary_stats(stats_dict, tool: str = '', sample: str = '', dataset: str = '', sep: str = '\t', file_handle=sys.stdout, extra_metadata=None):

    for rank, stats in stats_dict.items():
        sample = stats.name
        if extra_metadata:
            dataset = extra_metadata[sample]['Dataset']
            tool = extra_metadata[sample]['Tool']

        extra_metadata = None
        file_handle.write('{}\t{}\t{}\t{}\t{}\t{}{}\n'.format(sample, dataset, tool, rank.value.name, 'TP', stats.tp, '\t' + extra_metadata if extra_metadata else ""))
        file_handle.write('{}\t{}\t{}\t{}\t{}\t{}{}\n'.format(sample, dataset, tool, rank.value.name, 'FP', stats.fp, '\t' + extra_metadata if extra_metadata else ""))
        file_handle.write('{}\t{}\t{}\t{}\t{}\t{}{}\n'.format(sample, dataset, tool, rank.value.name, 'FN', stats.fn, '\t' + extra_metadata if extra_metadata else ""))
        file_handle.write('{}\t{}\t{}\t{}\t{}\t{}{}\n'.format(sample, dataset, tool, rank.value.name, 'Sensitivity', stats.Sensitivity(), '\t' + extra_metadata if extra_metadata else ""))
        file_handle.write('{}\t{}\t{}\t{}\t{}\t{}{}\n'.format(sample, dataset, tool, rank.value.name, 'Precision', stats.Precision(), '\t' + extra_metadata if extra_metadata else ""))
        file_handle.write('{}\t{}\t{}\t{}\t{}\t{}{}\n'.format(sample, dataset, tool, rank.value.name, 'F1', stats.F1(), '\t' + extra_metadata if extra_metadata else ""))


def print_abundance_stats(stats_dict, tool: str = '', dataset: str = '', sep: str = '\t', file_handle=sys.stdout, extra_metadata=None):
    for rank, stats in stats_dict.items():

        sample = stats.name
        if extra_metadata:
            dataset = extra_metadata[sample]['Dataset']
            tool = extra_metadata[sample]['Tool']

        extra_metadata = None

        file_handle.write('{}\t{}\t{}\t{}\t{}\t{}{}\n'.format(sample, dataset, tool, rank.value.name, 'PearsonCorrelationIntersect',
                                                  stats.pearson_correlation_intersection, '\t' + extra_metadata if extra_metadata else ""))
        file_handle.write('{}\t{}\t{}\t{}\t{}\t{}{}\n'.format(sample, dataset, tool, rank.value.name, 'PearsonCorrelationUnion',
                                                              stats.pearson_correlation_union, '\t' + extra_metadata if extra_metadata else ""))
        file_handle.write('{}\t{}\t{}\t{}\t{}\t{}{}\n'.format(sample, dataset, tool, rank.value.name, 'SpearmanCorrelationIntersect',
                                                              stats.spearman_correlation_intersection, '\t' + extra_metadata if extra_metadata else ""))
        file_handle.write('{}\t{}\t{}\t{}\t{}\t{}{}\n'.format(sample, dataset, tool, rank.value.name, 'SpearmanCorrelationUnion',
                                                              stats.spearman_correlation_union, '\t' + extra_metadata if extra_metadata else ""))
        file_handle.write('{}\t{}\t{}\t{}\t{}\t{}{}\n'.format(sample, dataset, tool, rank.value.name, 'BrayCurtisIntersect',
                                                              stats.bray_curtis_intersection, '\t' + extra_metadata if extra_metadata else ""))
        file_handle.write('{}\t{}\t{}\t{}\t{}\t{}{}\n'.format(sample, dataset, tool, rank.value.name, 'BrayCurtisUnion',
                                                              stats.bray_curtis_union, '\t' + extra_metadata if extra_metadata else ""))
        file_handle.write('{}\t{}\t{}\t{}\t{}\t{}{}\n'.format(sample, dataset, tool, rank.value.name, 'L2Intersect',
                                                              stats.l2_intersection, '\t' + extra_metadata if extra_metadata else ""))
        file_handle.write('{}\t{}\t{}\t{}\t{}\t{}{}\n'.format(sample, dataset, tool, rank.value.name, 'L2Union',
                                                              stats.l2_union, '\t' + extra_metadata if extra_metadata else ""))


def print_misc_stats(stats_dict, tool: str = '', dataset: str = '', sep: str = '\t', file_handle=sys.stdout, extra_metadata=None):
    for rank, stats in stats_dict.items():

        sample = stats.name
        if extra_metadata:
            dataset = extra_metadata[sample]['Dataset']
            tool = extra_metadata[sample]['Tool']

        extra_metadata = None

        for metric, value in stats.statistics.items():
            file_handle.write('{}\t{}\t{}\t{}\t{}\t{}{}\n'.format(sample, dataset, tool, rank.value.name, metric,
                                                                  value, '\t' + extra_metadata if extra_metadata else ""))



def write_statistics(output: str, all_statistics_dict, all_abundance_statistics_dict, all_misc_statistics_dict, header=True, extra_metadata=None):
    with open(output, 'w') as out:
        if header:
            out.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                "sample", "dataset", "tool", "rank", "metric", "value"
            ))
        for tool, statistics_dict_list in all_statistics_dict.items():
            abundance_dict_list = all_abundance_statistics_dict[tool]
            misc_dict_list = all_misc_statistics_dict[tool]
            dataset_num = 0
            for statistics_dict, abundance_dict, misc_dict in zip(statistics_dict_list, abundance_dict_list, misc_dict_list):
                dataset = str(dataset_num)
                print_binary_stats(statistics_dict, tool, dataset, file_handle=out, extra_metadata=extra_metadata)
                print_abundance_stats(abundance_dict, tool, dataset, file_handle=out, extra_metadata=extra_metadata)
                print_misc_stats(misc_dict, tool, dataset, file_handle=out, extra_metadata=extra_metadata)
                # print_binary_stats(statistics_dict, tool, dataset, extra_metadata=extra_metadata)
                # print_abundance_stats(abundance_dict, tool, dataset, extra_metadata=extra_metadata)
                dataset_num += 1
                