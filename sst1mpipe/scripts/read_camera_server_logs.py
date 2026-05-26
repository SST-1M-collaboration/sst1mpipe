import argparse
import re


def main():

    parser = argparse.ArgumentParser(description="Extract information on logs of the camera server")
    parser.add_argument('file', type=str, nargs='?', help='Path to log file')
    args = parser.parse_args()

    pattern = re.compile(
          r'^(INFO|DEBUG|WARN|ERROR)\s+' # LOG LEVEL
          r'(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2})\s+' # DATE
          r'(\w+)\s+'
          r'([\w.]+):\s+\[([0-9A-F]+):([0-9A-F]+)\]\s+'
          r'(-?\d+)\s+' # RXRING drop
          r'(-?\d+)\s+' # PACKET ERRS gaps+
          r'(-?\d+)\s+' # PACKET ERRS gaps-
          r'(-?\d+)\s+' # UDP INGEST pckts
          r'(-?\d+)\s+' # UDP INGEST MBs
          r'(-?\d+(?:\.\d+)?)\s+' # CLOCK drift
          r'(-?\d+)\s+'  # TRIGGER RATES part
          r'(-?\d+)\s+'  # TRIGGER RATES camera
          r'(-?\d+)\s+'  # TRIGGER RATES array
          r'(-?\d+)\s+'  # TRIGGER RATES muon
          r'(-?\d+)\s+'  # TRIGGER RATES hill
          r'([A-Z]?)\s+' # MERGE OFFSET L
          r'(-?\d+)\s+'  # MERGE OFFSET min
          r'(-?\d+)\s+'  # MERGE OFFSET max
          r'([A-Z]?)\s+' # MERGE OFFSET E
          r'(-?\d+)\s+'  # MERGE OFFSET mean
          r'([A-Z]?)\s+' # MERGE OFFSET P
          r'(-?\d+)\s+'  # BUFFERS/QUEUES ARRQ
          r'(-?\d+)\s+'  # BUFFERS/QUEUES BUNQ
          r'(\d+%)\s+'   # BUFFERS/QUEUES UARN
          r'(\d+%)\s+'   # BUFFERS/QUEUES PARN
          r'(-?\d+)\s+'  # DAQ RATES obj
          r'(-?\d+)\s+'  # DAQ RATES MB/s

    )

    n_packets = 0
    n_array = 0
    n_daq = 0
    n_camera = 0
    n_part = 0
    n_gaps_p = 0
    n_gaps_m = 0
    merge_late = False
    merge_early = False
    merge_p = False


    with open(args.file) as f:
        lines = f.readlines()
        for line in lines:



            m = pattern.search(line)
            if m:
                # level = str(m.group(1))
                # date = m.group(2)
                # process = m.group(3)
                # version = m.group(4)
                # id_1, id_2 = m.group(5), m.group(6)

                # rxring = int(m.group(7))
                gaps_p = int(m.group(8))
                gaps_m = int(m.group(9))
                udp_pckts = int(m.group(10))
                # udp_mb = m.group(11)
                # clock_drift = m.group(12)

                rate_part = int(m.group(13))
                rate_camera = int(m.group(14))
                rate_array = int(m.group(15))
                # rate_muon = int(m.group(16))
                # rate_hill = int(m.group(17))

                if m.group(18) == 'L':
                    merge_late = True

                # merge_min = int(m.group(19))
                # merge_max = int(m.group(20))
                if m.group(21) == 'E':
                    merge_early = True

                # merge_mean = int(m.group(22))
                if m.group(23) == 'P':
                    merge_p = True

                rate_daq = int(m.group(28))

                n_packets += udp_pckts
                n_array += rate_array
                n_camera += rate_camera
                n_daq += rate_daq
                n_part += rate_part
                n_gaps_p += gaps_p
                n_gaps_m += gaps_m

    print("Total number of UDP packets: ", n_packets)
    print("Total number of camera events: ", n_camera)
    print("Total number of parts events: ", n_part)
    print("Total number of arrays events: ", n_array)
    print("Total number of gaps+ : ", n_gaps_p)
    print("Total number of gaps- : ", n_gaps_m)
    print("Total number of DAQ obj : ", n_daq)
    print(f"Merge late {merge_late}")
    print(f"Merge early {merge_early}")
    print(f"Merge P {merge_p}")
