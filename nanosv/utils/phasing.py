from utils import parse_breakpoints as breakpoint
from utils import parse_reads as read
from utils import parse_bam as bam
from utils import create_vcf as c_vcf
import NanoSV
import random
matrix = []

def make_matrix(sv_id, windows):
    """
    Create matrix of positions and reads with their ref/alt/- as variable. Parse phasing result and select
    best result
    :param sv_id:
    :param windows:
    :return best result with the purity, phasing score and snps used:
    """
    global matrix
    x = 0
    scores = [0,0,0,0]
    bp = -1
    for window in windows:
        matrix = []
        if x == 0:
            chr = breakpoint.structural_variants[sv_id].chr
        else:
            chr = breakpoint.structural_variants[sv_id].chr2
        bin_start = int(window[0] / NanoSV.opts_variant_bin_size)
        bin_end = int(window[1] / NanoSV.opts_variant_bin_size)
        if bin_start < 0:
            bin_start = 0
        if bin_end > c_vcf.vcf.contigs[chr][1]:
            bin_end = int(c_vcf.vcf.contigs[chr][1] / NanoSV.opts_variant_bin_size)
        sv_reads = get_sv_reads(sv_id, x)
        for bin in range(bin_start, bin_end+1):
            for variant_position in bam.variants[chr][bin]:
                if int(window[0]) <= variant_position <= int(window[1]):
                    matrix.append([])
                    matrix_fill_position(breakpoint.structural_variants[sv_id].ref_qname_clips[x], variant_position, sv_id, chr, bin, x)
                    matrix_fill_position(sv_reads, variant_position, sv_id, chr, bin, x)
        matrix = list(map(list, zip(*matrix)))
        if len(matrix) == 0:
            x += 1
            continue
        phasing_result= clustering(matrix, list(range(len(breakpoint.structural_variants[sv_id].ref_qname_clips[x]),len(matrix))), x)
        if sum(phasing_result[:1]) > sum(scores[:1]):
            scores = phasing_result
            scores.append(len(matrix[0]))
            bp += 1
        x += 1
        deleted = 0
        for segment in range(len(matrix)):
            if len(set(matrix[segment-deleted])) <= 1 and matrix[segment-deleted][0] == "-":
                del matrix[segment-deleted]
                deleted += 1
    return scores


def get_sv_reads(sv_id, x):
    """
    gets reads supporting the sv. Used later to evaluate clustering
    :param sv_id:
    :param x:
    :return list with sv_reads:
    """
    sv_reads = []
    for bp_id in breakpoint.structural_variants[sv_id].breakpoints:
        if x == 0:
            sv_reads.append(read.breakpoints[bp_id].segment_1['id'])
        else:
            sv_reads.append(read.breakpoints[bp_id].segment_2['id'])
    return sv_reads


def matrix_fill_position(sv_or_ref, variant_position, sv_id, chr, bin, x):
    """
    fill position in matrix with correct value for either ref or alt base, or a - with unknown haplotype
    :param sv_or_ref:
    :param variant_position:
    :param sv_id:
    :param chr:
    :param bin:
    :param x:
    """
    global matrix
    for segment_id in sv_or_ref:
        if segment_id[2] in bam.variants[chr][bin][variant_position].segments:
            matrix[len(matrix) - 1].append(bam.variants[chr][bin][variant_position].segments[segment_id[2]][0])
        else:
            if breakpoint.structural_variants[sv_id].flag1 == 0 and x == 0:
                if bam.segments[segment_id[0]][segment_id[1]][segment_id[2]].pos > variant_position:
                    matrix[len(matrix) - 1].append('-')
                else:
                    matrix[len(matrix) - 1].append('=')
            elif breakpoint.structural_variants[sv_id].flag1 == 16 and x == 0:
                if bam.segments[segment_id[0]][segment_id[1]][segment_id[2]].end < variant_position:
                    matrix[len(matrix) - 1].append('-')
                else:
                    matrix[len(matrix) - 1].append('=')
            elif breakpoint.structural_variants[sv_id].flag2 == 0 and x == 1:
                if bam.segments[segment_id[0]][segment_id[1]][segment_id[2]].end < variant_position:
                    matrix[len(matrix) - 1].append('-')
                else:
                    matrix[len(matrix) - 1].append('=')
            elif breakpoint.structural_variants[sv_id].flag2 == 16 and x == 1:
                if bam.segments[segment_id[0]][segment_id[1]][segment_id[2]].pos > variant_position:
                    matrix[len(matrix) - 1].append('-')
                else:
                    matrix[len(matrix) - 1].append('=')


def clustering(matrix, sv_reads, bp_id):
    """
    creates clustering matrix by calling make_clustering_matrix(). Clusters most similar reads and sends
    result to judge_clustering to be evaluated. Result comes back and is returned.
    :param matrix:
    :param sv_reads:
    :param bp_id:
    :return phasing result:
    """
    clustering_matrix = make_clustering_matrix(matrix)
    while len(clustering_matrix) > 2:
        keys = []
        for x in clustering_matrix:
            keys.append(x)
        highest_score = 0
        readA = 0
        readB = 0
        for i in keys:
            for key, value in clustering_matrix[i].items():
                if value >= highest_score:
                    highest_score = value
                    readA = key
                    readB = i
        if highest_score < NanoSV.opts_clustering_cutoff:
            break
        merged_name = str(readA) + "," + str(readB)
        merged_dict = {}
        for j in keys:
            if j == readA or j == readB:
                continue
            sum_of_scores = 0
            for read in [readA, readB]:
                if max(map(int, read.split(","))) >= max(map(int, j.split(","))):
                    sum_of_scores += clustering_matrix[str(read)][str(j)]
                else:
                    sum_of_scores += clustering_matrix[str(j)][str(read)]
            merged_dict[str(j)] = sum_of_scores / 2
        del_list = []
        for item, merged_value in merged_dict.items():
            if max(map(int, item.split(","))) <= max(map(int, readB.split(","))):
                continue
            clustering_matrix[str(item)][str(merged_name)] = merged_value
            del_list.append(item)
        del clustering_matrix[str(readA)]
        del clustering_matrix[str(readB)]
        for read in clustering_matrix:
            if readA in clustering_matrix[read]:
                del clustering_matrix[read][readA]
            if readB in clustering_matrix[read]:
                del clustering_matrix[read][readB]
        for item in del_list:
            del merged_dict[str(item)]
        clustering_matrix[merged_name] = merged_dict
    breakpoint_result = judge_clustering(clustering_matrix, sv_reads, len(matrix), bp_id)
    return breakpoint_result


def judge_clustering(clustering_matrix, sv_reads, total_reads, x):
    """
    Judges clustering result and calculates purity and phasing scores. Also calls randomise() 100000 times to get a
    random result. With these random results a p-value can be calculated.
    :param clustering_matrix:
    :param sv_reads:
    :param total_reads:
    :param x:
    :return list with purity score, phasing score and p-value:
    """
    purity_chart = []
    phasing_chart = []
    clusters = []
    for key in clustering_matrix:
        clusters.append(key.split(","))
    longest_clusters = []
    for length_cluster in [len(clusters), len(clusters)-1]:
        long_cluster = []
        for j in range(length_cluster):
            if len(long_cluster) <= len(clusters[j]):
                long_cluster = clusters[j]
        del clusters[clusters.index(long_cluster)]
        longest_clusters.append(long_cluster)
    amounts_per_cluster = []
    clusternr = 0
    for cluster in longest_clusters:
        amounts = [0, 0]
        clusternr += 1
        for read in cluster:
            if int(read) in sv_reads:
                amounts[0] += 1
            else:
                amounts[1] += 1
        amounts_per_cluster.append(amounts)
    pur_percentage_per_cluster = []
    phasing_percentage_per_cluster = []
    for i in range(len(amounts_per_cluster)):
        purity_percentages = []
        purity_percentages.append(amounts_per_cluster[i][0] / (amounts_per_cluster[i][0] + amounts_per_cluster[i][1]) * 100)
        purity_percentages.append(amounts_per_cluster[i][1] / (amounts_per_cluster[i][0] + amounts_per_cluster[i][1]) * 100)
        pur_percentage_per_cluster.append(purity_percentages)
        phasing_percentages = []
        phasing_percentages.append(amounts_per_cluster[i][0] / len(sv_reads) * 100)
        try:
            phasing_percentages.append(amounts_per_cluster[i][1] / (total_reads - len(sv_reads)) * 100)
        except ZeroDivisionError:
            phasing_percentages.append(0)
        phasing_percentage_per_cluster.append(phasing_percentages)


    pur_sv_score = (pur_percentage_per_cluster[0][0] - pur_percentage_per_cluster[1][0])
    pur_ref_score = (pur_percentage_per_cluster[0][1] - pur_percentage_per_cluster[1][1])
    if pur_sv_score < 0:
        pur_sv_score = pur_sv_score * -1
    if pur_ref_score < 0:
        pur_ref_score = pur_ref_score * -1
    phasing_sv_score = (phasing_percentage_per_cluster[0][0] + phasing_percentage_per_cluster[1][0])
    phasing_ref_score = (phasing_percentage_per_cluster[0][1] + phasing_percentage_per_cluster[1][1])

    purity_score = (pur_sv_score + pur_ref_score) / 2
    phasing_score = (phasing_sv_score + phasing_ref_score) / 2
    random_purities = []
    for i in range(100000):
        random_purities.append(randomise(longest_clusters, sv_reads))
    teller = 1
    for value in random_purities:
        if float(value) >= purity_score:
            teller += 1
    pvalue = teller/len(random_purities)
    purity_chart.append(int(purity_score))
    phasing_chart.append(int(phasing_score))
    return [int(purity_score), int(phasing_score), pvalue]


def randomise(longest_clusters, sv_reads):
    """
    Creates a random result and calculates scores to be used with a calculation of the p-value.
    :param longest_clusters:
    :param sv_reads:
    :return purity score:
    """
    random_options = []
    for ref in range(len(matrix) - len(sv_reads)):
        random_options.append(1)
    for alt in range(len(sv_reads)):
        random_options.append(2)
    random.shuffle(random_options)
    new_clusters = [[], [], []]
    for clusternr in range(0,3):
        if clusternr == 0:
            for zero in range(len(matrix) - (len(longest_clusters[0])+len(longest_clusters[1]))):
                new_clusters[clusternr].append(random_options[-1])
                del random_options[-1]
        else:
            for one_or_two in range(len(longest_clusters[clusternr-1])):
                new_clusters[clusternr].append(random_options[-1])
                del random_options[-1]

    pur_percentage_per_cluster = []
    amounts_per_cluster = [[new_clusters[1].count(2),new_clusters[1].count(1)], [new_clusters[2].count(2), new_clusters[2].count(1)]]
    for i in range(len(amounts_per_cluster)):
        purity_percentages = []
        purity_percentages.append(amounts_per_cluster[i][0] / (amounts_per_cluster[i][0] + amounts_per_cluster[i][1]) * 100)
        purity_percentages.append(amounts_per_cluster[i][1] / (amounts_per_cluster[i][0] + amounts_per_cluster[i][1]) * 100)
        pur_percentage_per_cluster.append(purity_percentages)
    pur_sv_score = (pur_percentage_per_cluster[0][0] - pur_percentage_per_cluster[1][0])
    pur_ref_score = (pur_percentage_per_cluster[0][1] - pur_percentage_per_cluster[1][1])
    if pur_sv_score < 0:
        pur_sv_score = pur_sv_score * -1
    if pur_ref_score < 0:
        pur_ref_score = pur_ref_score * -1

    purity_score = (pur_sv_score + pur_ref_score) / 2
    return purity_score


def make_clustering_matrix(matrix):
    """
    calculates similarity between reads to be put in the clustering matrix
    :param matrix:
    :return filled clustering matrix:
    """
    clustering_matrix = {}
    for i in range(len(matrix)):
        clustering_matrix[str(i)] = {}
        for j in range(i + 1):
            mutations_in_common = 0
            if j == i:
                clustering_matrix[str(i)][str(j)] = 0
            else:
                amount_positions = len(matrix[i])
                for pos in range(len(matrix[i])):
                    if matrix[i][pos] == '-' and matrix[j][pos] == '-':
                        amount_positions -= 1
                        continue
                    if matrix[i][pos] == matrix[j][pos]:
                        mutations_in_common += 1
                if amount_positions == 0:
                    clustering_matrix[str(i)][str(j)] = 0
                else:
                    clustering_matrix[str(i)][str(j)] = mutations_in_common / amount_positions
    return clustering_matrix
