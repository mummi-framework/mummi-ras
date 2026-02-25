#!/usr/bin/env python3
# Authors: Loic Pottier <pottier1@llnl.gov>

import argparse
import csv
from datetime import datetime
import logging
import os
import time
from typing import List

import numpy as np

logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', level=logging.INFO)

def read_csv(csv_file: str) -> List[List[str]]:
    if not os.path.isfile(csv_file):
        return []
    result = []
    i = 0
    with open(csv_file, newline='') as f:
        for row in csv.reader(f, delimiter=','):
            if i == 0 or row == []:
                i += 1
                continue
            result.append(row)
    return result

def timestamp_fmt(short=False):
    if short:
        return "%H%M%S"
    return "%Y%m%d-%H%M%S"

def cleanup_x_data(x_data):
    """
    If all timestamps are all during one day, we remove the day.
    Example of one timestamp: 20240728-101201
    """
    if len(set([s.split('-')[0] for s in x_data])) == 1:
        return [s.split('-')[-1] for s in x_data]
    else:
        return x_data

def select_data(data, convert_ts=False):
    """
    Select the correct data to form X and Y for the plot
    We select number of createsims and number of CG for Y
    """
    x_data, y_data, y_data1, y_data2 = [], {}, [], []
    for row in data:
        #  One row: 2136900,1000,20240728-123853,1,1,0,0,0,0,0,0,0,0,0,4111,0,0,0,0,0,0,0
        if convert_ts:
            x = datetime.strptime(row[2], timestamp_fmt(short=False))
        else:
            x = row[2]
        x_data.append(x)
        y_data1.append(int(row[6]))
        y_data2.append(int(row[14]))

    y_data["createsims"] = np.array(y_data1)
    y_data["cg"] = np.array(y_data2)
    return x_data, y_data


def plot_stackedbar(x_data, y_data, outfile = "mummi-gantt"):
    try:
        import matplotlib.pyplot as plt
        from matplotlib.dates import DateFormatter, DayLocator, HourLocator, drange
    except Exception as e:
        logging.error(f"Did not find Matplolib")
        return

    fig, ax = plt.subplots()

    bottom = np.zeros(len(x_data))
    x_data = np.array(x_data)

    i = 0
    colors = ["darkorange", "darkgreen"]
    for job, count in y_data.items():
        ax.bar(x_data, count, align='center', label=job, width=0.1, bottom=bottom, color=colors[i])
        bottom += count
        i += 1

    ax.xaxis.set_major_formatter(DateFormatter(timestamp_fmt()))
    ax.fmt_xdata = DateFormatter('%Y-%m-%d %H:%M:%S')
    fig.autofmt_xdate()

    ax.set_xlabel('Timestamp')
    ax.set_ylabel('Number of createsims / CG')
    ax.legend(prop={'size': 12})

    ax.autoscale_view()
    ax.legend(fancybox=True, framealpha=0)
    #plt.show()
    timestr = time.strftime("%Y%m%d-%H%M%S")
    plt.savefig(f"{outfile}-{timestr}.png", bbox_inches = "tight", transparent = True, dpi = 600)

def main(args):
    data = read_csv(args.file)
    x_data, y_data = select_data(data, convert_ts=True)
    # x_data = cleanup_x_data(x_data)
    plot_stackedbar(x_data, y_data)

def parse_args() -> dict:
    parser = argparse.ArgumentParser(description="Plot the Gantt Diagram for a MuMMI run")
    parser.add_argument(
        '-f',
        '--file',
        help="CSV status file (status_XX.csv where XX is a timestamp)",
        required=True
    )

    args = parser.parse_args()

    return args

if __name__ == "__main__":
    args = parse_args()
    main(args)

