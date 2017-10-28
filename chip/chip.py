import numpy as np
import argparse
import matplotlib.pyplot as plt


def read_mappings( fname ):
    """Reads mappings from text file"""
    forward = []
    reverse = []
    with open(fname) as f:
        for line in f:
            parts = line.rstrip().split(' ')
            if len(parts) < 3:
                continue

            if parts[2] == '+': #forward sring
                forward.append((int(parts[0]), int(parts[1])))
            else:
                reverse.append((int(parts[0]), int(parts[1])))

    return forward, reverse


def smooth_array( data, k):
    """Smooths an array based with an average on k, which specifies the size of the filter 2k + 1"""
    padd_data = np.concatenate((np.zeros(k), data, np.zeros(k)))
    filter = np.ones(2 * k + 1) / (2 * k + 1)
    result = np.zeros(data.shape)
    for i in range(data.shape[0]):
        window = padd_data[i: i + (2 * k) + 1]
        avg = sum(window * filter)
        result[i] = avg
    return result


def smooth_array_2(data, k):
    """The same as smooth array but faster"""
    k = 2 * k + 1
    bound = k / 2
    return np.convolve(data, np.ones((k,)) / k)[bound:-bound]

def get_arg_maximas(data, threshold=None):
    """Returns and array with the indexes of the maximas, if threshold, removes all peaks
    under the threshold value"""
    peaks = []
    for i in range(1, data.shape[0] - 1):
        prev = data[i - 1]
        next = data[i + 1]
        if data[i] >= prev and data[i] > next:
            peaks.append(i)

    # if threshold
    if threshold is not None:
        return  np.array([ peak for peak in peaks if data[peak] >= threshold])

    # otherwise
    return np.array(peaks)

def fstore(data, file):
    """stores data in file"""
    with open(file, 'w') as f:
        for line in data:
            f.write(' '.join([str(cell) for cell in line]) + '\n')

def plot(forward_data, reverse_data, peaks):
    """Plots both strands and shadow areas for enriched regions"""
    x = range(forward_data.shape[0])
    y = forward_data
    y2 = reverse_data * -1.

    fig, ax = plt.subplots()

    # lines for strands
    plt.fill(x, y, '-', linewidth=2,
                     label='Forward strand', color='b')

    plt.fill(x, y2, '-', linewidth=2,
                     label='Reverse strand', color='r')

    # fills for enriched regions

    for init, end in peaks:
        plt.axvspan(init, end, color='g', alpha=0.5)
        #plt.text(init, y2[end],'Enriched area', fontsize=10)


    plt.yticks([])
    plt.show()

def merge_pairs(a, b=None):
    """This function merge pairs, recursively,  of data when overlaping occurs, this means, when the end of first pair
    is in between the init of the second pair and the end, this way, a bigger pair is created"""
    if b is None: # first condition
        b = a
        a = [b.pop(0)]
        return merge_pairs(a, b)
    elif b == []: # end condition
        return a
    else: # working condition
        if a[-1][1] > b[0][0]: ## merge # if end of a in  between init and end of next pair
            a[-1][1] = b[0][1]
            b.pop(0)
        else: ## do not merge, #otherwise, just keep this pair
            a.append(b.pop(0))
        return merge_pairs(a, b)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file",
                        required=True,
                        type=str,
                        default=None,
                        help="Text file containing the mappings")


    parser.add_argument("-o", "--output",
                        required=True,
                        type=str,
                        default=None,
                        help="Output file")

    args = parser.parse_args()

    f_maps, r_maps = read_mappings(args.file)

    forward_strand = np.zeros(1000) # 1000 base pairs
    reverse_strand = np.zeros(1000)

    # a)
    for init, end in f_maps:
        forward_strand[init:end] += 1

    for init, end in r_maps:
        reverse_strand[init:end] += 1

    ## CHECK thiS
    #reverse_strand = reverse_strand[::-1]
    # b)
    k = 4
    # naive
    #forward_smooth = smooth_array(forward_strand, k)
    #reverse_smooth = smooth_array(reverse_strand, k)

    # optimized
    forward_smooth = smooth_array_2(forward_strand, k)
    reverse_smooth = smooth_array_2(reverse_strand, k)

    # c and d
    #peaks = get_arg_maximas(forward_smooth)
    peaks_forward = get_arg_maximas(forward_smooth, threshold=100)
    peaks_reverse = get_arg_maximas(reverse_smooth, threshold=100)
    # e
    peak_pairs = []
    from_ = 120 #100
    to_ = 200
    for peak in peaks_forward:
        candidates = peaks_reverse[(peaks_reverse >= peak + from_) & (peaks_reverse < peak + to_)]

        if len(candidates) > 0:
            best_one = np.argmax(reverse_smooth[candidates])
            peak_pairs.append([peak, candidates[best_one]])


    ## CHECK merging pairs
    #peak_pairs = merge_pairs(peak_pairs)
    print peak_pairs
    fstore(peak_pairs, args.output)

    plot(forward_smooth, reverse_smooth, peak_pairs)

