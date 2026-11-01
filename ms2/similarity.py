"""Square-root cosine with greedy, closest-first fragment matching."""

import numba as nb
import numpy as np

from .protocol import INTENSITY_POWER, PRECURSOR_MARGIN_DA, TOLERANCE_DA, TOLERANCE_PPM


@nb.njit(cache=True)
def clean_peaks(mz, intensity, start, stop, precursor):
    """Remove peaks above the cutoff and normalise the retained intensities."""
    masses = np.empty(stop - start)
    values = np.empty(stop - start)
    maximum = 0.0
    for i in range(start, stop):
        if mz[i] <= precursor - PRECURSOR_MARGIN_DA and intensity[i] > maximum:
            maximum = intensity[i]
    if maximum <= 0:
        return masses[:0], values[:0]
    count = 0
    for i in range(start, stop):
        if mz[i] > precursor - PRECURSOR_MARGIN_DA:
            continue
        if intensity[i] / maximum < 0:
            continue
        masses[count] = mz[i]
        values[count] = intensity[i] / maximum
        count += 1
    return masses[:count], values[:count]


@nb.njit(cache=True)
def match_peaks(first, second):
    """Match each peak at most once, resolving candidates by mass difference."""
    candidate_first = []
    candidate_second = []
    differences = []
    lower = 0
    for i in range(len(first)):
        tolerance = max(TOLERANCE_DA, first[i] * TOLERANCE_PPM * 1e-6)
        while lower < len(second) and second[lower] < first[i] - tolerance:
            lower += 1
        j = lower
        while j < len(second) and second[j] <= first[i] + tolerance:
            candidate_first.append(i)
            candidate_second.append(j)
            differences.append(abs(second[j] - first[i]))
            j += 1
    if len(differences) == 0:
        return np.empty(0, np.int64), np.empty(0, np.int64)
    distance_array = np.empty(len(differences))
    for k in range(len(differences)):
        distance_array[k] = differences[k]
    # Keep NumPy/Numba's default argsort, including its treatment of ties.
    order = np.argsort(distance_array)
    used_first = np.zeros(len(first), np.bool_)
    used_second = np.zeros(len(second), np.bool_)
    selected_first = []
    selected_second = []
    for k in order:
        i = candidate_first[k]
        j = candidate_second[k]
        if not used_first[i] and not used_second[j]:
            used_first[i] = True
            used_second[j] = True
            selected_first.append(i)
            selected_second.append(j)
    result_first = np.empty(len(selected_first), np.int64)
    result_second = np.empty(len(selected_second), np.int64)
    for k in range(len(selected_first)):
        result_first[k] = selected_first[k]
        result_second[k] = selected_second[k]
    return result_first, result_second


@nb.njit(cache=True)
def score_pairs(mz, intensity, offsets, precursors, queries, references):
    """Score directed spectrum pairs; an empty retained spectrum gives NaN."""
    scores = np.empty(len(queries))
    for pair in range(len(queries)):
        query = queries[pair]
        reference = references[pair]
        first, v1 = clean_peaks(
            mz, intensity, offsets[query], offsets[query + 1], precursors[query]
        )
        second, v2 = clean_peaks(
            mz, intensity, offsets[reference], offsets[reference + 1], precursors[reference]
        )
        if len(first) == 0 or len(second) == 0:
            scores[pair] = np.nan
            continue
        w1 = v1 ** INTENSITY_POWER
        w2 = v2 ** INTENSITY_POWER
        matched_first, matched_second = match_peaks(first, second)
        numerator = 0.0
        for k in range(len(matched_first)):
            numerator += w1[matched_first[k]] * w2[matched_second[k]]
        scores[pair] = numerator / np.sqrt((w1 * w1).sum() * (w2 * w2).sum())
    return scores


def maximum_scores(mz, intensity, offsets, precursors, queries, reference_sets):
    """Return each maximum and its first attaining reference in input order."""
    if len(queries) != len(reference_sets):
        raise ValueError("Every query needs one reference set")
    if not len(queries):
        return np.empty(0), np.empty(0, dtype=np.int64)
    lengths = np.asarray([len(refs) for refs in reference_sets], dtype=np.int64)
    if np.any(lengths == 0):
        raise ValueError("A query has no eligible reference spectra")
    first = np.repeat(np.asarray(queries, dtype=np.int64), lengths)
    second = np.concatenate(reference_sets).astype(np.int64, copy=False)
    owner = np.repeat(np.arange(len(queries)), lengths)
    scores = score_pairs(mz, intensity, offsets, precursors, first, second)
    scores = np.where(np.isnan(scores), -1.0, scores)
    maxima = np.full(len(queries), -1.0)
    np.maximum.at(maxima, owner, scores)
    best = np.full(len(queries), -1, dtype=np.int64)
    for group, reference in zip(owner[scores == maxima[owner]], second[scores == maxima[owner]]):
        if best[group] < 0:
            best[group] = reference
    return maxima, best
