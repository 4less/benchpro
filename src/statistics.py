from src.profile import Profile, default_ranks
from src.binary_statistics import BinaryStatistics
from src.profile import RankTemp
from src.profile import Rank
import scipy as sc
import sys
import math



class AbundanceStatistics:
    def __init__(self):
        self.pearson_correlation_intersection = 0
        self.pearson_correlation_union = 0
        self.bray_curtis = 0
        self.l2 = 0


def get_abundance_dict(profile, rank):
    abundance_dict = dict()
    for row in profile.Rows():
        name = row.GetLineage().Get(rank)
        if name not in abundance_dict:
            abundance_dict[name] = 0.0
        abundance_dict[name] += row.GetAbundance()

    return abundance_dict


def bray_curtis(a, b):
    total_a = sum(a)
    total_b = sum(b)
    shared_a = sum(va for va, vb in zip(a, b) if va > 0 and vb > 0)
    shared_b = sum(vb for va, vb in zip(a, b) if va > 0 and vb > 0)
    #
    # print("Bray_Curtis")
    # print("gold: {}".format(a))
    # print("pred: {}".format(b))
    # print("total_gold: {}, total_pred: {}\tshared_gold: {}, pred_gold: {}".format(total_a, total_b, shared_a, shared_b))
    # print("min@ {}".format(min(shared_a, shared_b)))
    # print("(total_a + total_b): {}".format((total_a + total_b)))
    bc = (2*min(shared_a, shared_b) / (total_a + total_b))

    # print("Bray: {}".format(bc))
    return bc

def euclidian(a,b):
    sum = 0
    for a1, a2 in zip(a, b):
        sum += (a1 - a2) ** 2
    return math.sqrt(sum)

def get_rank_statistics(gold_profile: Profile, prediction_profile: Profile,
                        ranks=default_ranks, labeled_profile_output=None):
    rank_dict = dict()
    abundance_stats = dict()

    for rank in ranks:
        stats = BinaryStatistics()
        abundance = AbundanceStatistics()
        # print(rank.name)


        prediction_set = set(row.GetLineage().Get(rank) for row in prediction_profile.Rows())
        gold_set = set(row.GetLineage().Get(rank) for row in gold_profile.Rows())

        shared = prediction_set.intersection(gold_set)
        pred_only = prediction_set.difference(gold_set)
        gold_only = gold_set.difference(prediction_set)

        if len(shared) == 0:
            print("prediction_set: {}" .format(prediction_set))
            print("gold_set: {}" .format(gold_set))
            input()

        gold_set = set(row.GetLineage().Get(rank) for row in gold_profile.Rows())
        stats.tp = len(shared)
        stats.fp = len(pred_only)
        stats.fn = len(gold_only)



        ##########################################################################################
        # Pearson correlation on shared taxa

        gold_dict = get_abundance_dict(gold_profile, rank)
        prediction_dict = get_abundance_dict(prediction_profile, rank)
        gold_vec = [gold_dict[key] for key in shared]
        pred_vec = [prediction_dict[key] for key in shared]

        if labeled_profile_output:
            print(shared)
            print(pred_only)
            print(gold_only)
            for label in shared:
                labeled_profile_output.write("{}\t{}\t{}\t{}\t{}\n".format(prediction_profile.name, "TP", "Prediction", label, prediction_dict[label]))
                labeled_profile_output.write("{}\t{}\t{}\t{}\t{}\n".format(gold_profile.name, "TP", "Gold", label, gold_dict[label]))
            for label in pred_only:
                labeled_profile_output.write("{}\t{}\t{}\t{}\t{}\n".format(prediction_profile.name, "FP", "Prediction", label, prediction_dict[label]))
            for label in gold_only:
                labeled_profile_output.write("{}\t{}\t{}\t{}\t{}\n".format(gold_profile.name, "FN", "Gold", label, gold_dict[label]))


        # print(gold_vec)
        # print(pred_vec)
        if len(shared) > 1:
            # print(sc.stats.pearsonr(gold_vec, pred_vec))
            abundance.pearson_correlation_intersection = sc.stats.pearsonr(gold_vec, pred_vec).statistic
        else:
            abundance.pearson_correlation_intersection = 0

        ##########################################################################################
        # Pearson correlation on union taxa

        gold_vec += [gold_dict[key] for key in gold_only]
        gold_vec += [0.0 for _ in pred_only]
        pred_vec += [0.0 for _ in gold_only]
        pred_vec += [prediction_dict[key] for key in pred_only]

        # print(gold_vec)
        # print(pred_vec)

        if len(gold_vec) > 1:
            # print(sc.stats.pearsonr(gold_vec, pred_vec))
            abundance.pearson_correlation_union = sc.stats.pearsonr(gold_vec, pred_vec).statistic
        else:
            abundance.pearson_correlation_union = 0


        ##########################################################################################
        # Bray-Curtis
        abundance.bray_curtis = bray_curtis(gold_vec, pred_vec)
        ##########################################################################################
        # Euclidean (L2)
        abundance.l2 = euclidian(gold_vec, pred_vec)
        ##########################################################################################

        # print("False negatives: {}".format(gold_set.difference(prediction_set)))
        # print("False positives: {}".format(prediction_set.difference(gold_set)))

        abundance_stats[rank] = abundance
        rank_dict[rank] = stats

        try:
            stats.Validate()
        except:
            print(f"{gold_profile.name}\t{prediction_profile.name}\t{rank}")
            exit(9)

    return rank_dict, abundance_stats


def print_binary_stats(stats_dict, name: str = '', dataset: str = '', sep: str = '\t', file_handle=sys.stdout):
    prefix = '{}{}'.format(name, sep) if name else ''
    dataset_prefix = '{}{}'.format(dataset, sep) if name else ''

    for rank, stats in stats_dict.items():
        rank = '{}{}'.format(rank.value.name, sep) if rank else ''

        file_handle.write('{}{}{}{}\t{}\n'.format(prefix, dataset_prefix, rank, 'TP', stats.tp))
        file_handle.write('{}{}{}{}\t{}\n'.format(prefix, dataset_prefix, rank, 'FP', stats.fp))
        file_handle.write('{}{}{}{}\t{}\n'.format(prefix, dataset_prefix, rank, 'FN', stats.fn))
        file_handle.write('{}{}{}{}\t{}\n'.format(prefix, dataset_prefix, rank, 'Sensitivity', stats.Sensitivity()))
        file_handle.write('{}{}{}{}\t{}\n'.format(prefix, dataset_prefix, rank, 'Precision', stats.Precision()))
        file_handle.write('{}{}{}{}\t{}\n'.format(prefix, dataset_prefix, rank, 'F1', stats.F1()))


def print_abundance_stats(stats_dict, name: str = '', dataset: str = '', sep: str = '\t', file_handle=sys.stdout):
    prefix = '{}{}'.format(name, sep) if name else ''
    dataset_prefix = '{}{}'.format(dataset, sep) if name else ''

    for rank, stats in stats_dict.items():
        rank = '{}{}'.format(rank.value.name, sep) if rank else ''

        file_handle.write('{}{}{}{}\t{}\n'.format(prefix, dataset_prefix, rank, 'PearsonCorrelationIntersect',
                                                  stats.pearson_correlation_intersection))
        file_handle.write('{}{}{}{}\t{}\n'.format(prefix, dataset_prefix, rank, 'PearsonCorrelationUnion',
                                                  stats.pearson_correlation_union))
        file_handle.write('{}{}{}{}\t{}\n'.format(prefix, dataset_prefix, rank, 'BrayCurtis',
                                                  stats.bray_curtis))
        file_handle.write('{}{}{}{}\t{}\n'.format(prefix, dataset_prefix, rank, 'L2',
                                                  stats.l2))
