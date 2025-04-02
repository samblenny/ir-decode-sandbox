#!/usr/bin/python3
"""
Calculate clustered histograms for IR remote control pulse lengths from logic
analyzer capture files. The point of this is to help me identify suitable
timing thresholds to use for writing a state machine to do feature extraction.
"""
import csv
import sys

NEC_FILES = """
captures/insignia-power.csv
captures/insignia-1-repeating.csv
captures/insignia-1.csv
captures/insignia-11-fast.csv
captures/insignia-11-slow.csv
captures/insignia-power.csv
""".split()

SAMSUNG_FILES = """
captures/samsung-1-repeating.csv
captures/samsung-1-slow.csv
captures/samsung-11-slow.csv
captures/samsung-11.csv
""".split()

SIRC_FILES = """
captures/sony-BD-1-fast.csv
captures/sony-BD-1-repeating.csv
captures/sony-BD-11-fast.csv
captures/sony-BD-power.csv
""".split()


# Timing Constants (units are all µs, HIGH/LOW are for active LOW receiver)
#
# These timing constants were measured on a Vishay TSOP18638 which has their
# AG6 gain control and noise filtering. Supposedly AGC2 receivers are more
# suitable for Sony SIRC codes in particular. But, whatever. This works. Just
# be aware that my timing constants may be a bit off compared to a TSOP38238.
#
SIRC_HEADER = 2383     # Initial LOW pulse, no subsequent leader pulse
SIRC_HI = 590          # SIRC HIGH pulses are all the same width
SIRC_LO_0 = 609        # SIRC 0 bit LOW time
SIRC_LO_1 = 1208       # SIRC 1 bit LOW time
SIRC_GAP_20 = 14800    # SIRC gap after end of 20-bit code before repeat code
SIRC_PERIOD = 45000    # SIRC start-of-code to start-of-next-repeat is 45ms
NEC_HEADER = 8984      # Initial LOW pulse before HIGH leader pulse
NEC_LEADER = 4500      # NEC HIGH pulse after header (not part of data)
NEC_REPEAT_HI = 2246   # NEC HIGH pulse for held down key (NEC repeat code)
NEC_LO = 586           # NEC LOW pulses (except header)
NEC_HI_0 = 512         # NEC 0 bit HIGH time
NEC_HI_1 = 1681        # NEC 1 bit HIGH time
NEC_GAP_FIRST = 43550  # NEC gap before first repeat code
NEC_GAP_REST = 96222   # NEC gap before the rest of the repeat codes
SAMSUNG_HEADER = 4528  # Initial LOW pulse before HIGH leader pulse
SAMSUNG_LEADER = 4509  # Samsung HIGH pulse after header (not part of data)
SAMSUNG_LO = 578       # Samsung LOW pulses (except header)
SAMSUNG_HI_0 = 532     # Samsung 0 bit HIGH time
SAMSUNG_HI_1 = 1683    # Samsung 1 bit HIGH time
SAMSUNG_GAP = 46930    # Samsung HIGH gap between repeats for held down key

# This is used for clustering histogram bins (unit is µs)
TOLERANCE = 100


class StateMachine:

    def __init__(self, debug):
        self.debug = debug
        self.prev = None
        self.histogram = {}
        self.file_count = 0
        self.total_edges = 0
        self.current_file = ""

    def add_edge(self, time, pin_state):
        self.total_edges += 1
        if self.prev is None:
            print(f"--- {self.current_file}")
            self.prev = (pin_state, time)
        else:
            µs = round((time - self.prev[1]) * 1e6)
            # For pin state here, use the state that the pin was in before the
            # edge transition because the time delta measures the length of the
            # pin's previous state.
            prev_pin_state = self.prev[0]
            k = (prev_pin_state, µs)
            self.histogram[k] = self.histogram.get(k, 0) + 1
            self.prev = (pin_state, time)
            if self.debug:
                print(f"{k[0]} {k[1]:9d}")

    def start_new_file(self, filename):
        self.current_file = filename
        self.file_count += 1
        self.prev = None

    def combine_bins_into_clusters(self):
        """Combine raw histogram bins into clusters of similar pulse lengths"""
        bins = []
        current = None
        prev_µs = None
        prev_ps = None
        for key in sorted(self.histogram):
            (pin_state, µs) = key
            n = self.histogram[key]
            gap_µs = 0
            # Avoid a discontinuity when crossing from the longest LOW pulse
            # to the shortest HIGH pulse (remember we're sorting by tuples)
            if prev_ps is None:
                prev_ps = pin_state
                current = []
            elif prev_ps != pin_state:
                prev_ps = pin_state
                prev_µs = None
                if current and len(current) > 0:
                    bins.append(current)
                    current = []
            # Calculate the gap in µs between this bin and the previous bin
            if prev_µs is None:
                prev_µs = µs
            else:
                gap_µs = µs - prev_µs
                prev_µs = µs
            # Start new cluster if gap is too big
            if current and len(current) > 0 and gap_µs > TOLERANCE:
                bins.append(current)
                current = []
            # Add this bin's info to current cluster
            current.append((pin_state, µs, n))
        if current and len(current) > 0:
            bins.append(current)
        clustered_bins = bins
        if self.debug:
            print("Clustered Histogram Bins:")
            print("\n".join([str(n) for n in bins]))
        # Calculate average µs for each cluster
        averaged = {}
        frequency_total = 0
        for list_ in clustered_bins:
            total = 0
            n = 0
            for (pin_state, µs, frequency) in list_:
                total += µs * frequency
                n += frequency
            avg = round(total / n)
            averaged[(pin_state, avg)] = n
            frequency_total += n
        print("Cluster Averaged Histogram:")
        for key in averaged:
            (pin_state, µs) = key
            frequency = averaged[key]
            print(f" {pin_state} {µs:7d} {frequency:4d}")
        print(f"TOTAL EDGES: {self.total_edges}")
        print(f"TOTAL GAPS: {frequency_total}")
        print(f"TOTAL FILES: {self.file_count}")

    def print_raw_histogram(self):
        if not self.debug:
            return
        print("Raw Histogram:")
        prev_µs = None
        prev_ps = None
        for key in sorted(self.histogram):
            (pin_state, µs) = key
            n = self.histogram[key]
            gap_µs = 0
            if prev_ps is None:
                prev_ps = pin_state
            elif prev_ps != pin_state:
                prev_ps = pin_state
                prev_µs = None
            if prev_µs is None:
                prev_µs = µs
            else:
                gap_µs = µs - prev_µs
                prev_µs = µs
            print(f" {pin_state} {µs:7d} {n:3d}  {gap_µs}")


def parse_data(filename):
    events = []
    with open(filename) as f:
        reader = csv.reader(f)
        next(reader)  # Skip header
        for row in reader:
            events.append((float(row[0]), int(row[1])))  # (time, pin_state)
    return events

def main(file_list, debug):
    s = StateMachine(debug)
    for f in file_list:
        events = parse_data(f)
        s.start_new_file(f)
        for (time, pin_state) in events:
            s.add_edge(time, pin_state)
    s.print_raw_histogram()
    s.combine_bins_into_clusters()
    print()


print_help = True
args = sys.argv[1:]
debug = 'debug' in args
if 'nec' in args:
    main(NEC_FILES, debug)
    print_help = False
if 'sirc' in args:
    main(SIRC_FILES, debug)
    print_help = False
if 'samsung' in args:
    main(SAMSUNG_FILES, debug)
    print_help = False
if print_help:
    print("Usage: python3 histogram.py (nec|sirc|samsung) [debug]")
