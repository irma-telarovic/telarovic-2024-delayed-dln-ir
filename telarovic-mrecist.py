#!/usr/bin/python

import csv
import math

# The first day from which tumor volumes are taken into account when 
# calculating the best response and the best average response.
FIRST_DAY = 7 

INPUT_FILE_NAME = 'input-example.csv'

# Converts string values to numbers. Returns -1 if unsuccessful.
def parseValue(strVal):
    if strVal == '':
        return -1
    try:
        return float(strVal)
    except:
        return -1

# Reads mouse IDs and tumor volumes from the .csv file.
# Expects the mouse ID in the first row.
# Subseqent rows contain tumor volumes, one row for each day of measurement.
def readInput():
    csvFile = open(INPUT_FILE_NAME)
    reader = csv.reader(csvFile)
    ids = []
    volumes = []
    first = True
    for row in reader:
        if first:
            first = False
            ids = row
            continue
        volumes.append(list(map(parseValue, row)))
    return ids, volumes

# The response call returned here is defined in the Statistical analysis
# section of the publication.
def calcResponseCall(bestResponse, bestAverageResponse):
    if bestResponse < -0.95 and bestAverageResponse < -0.4:
        return 'mCR'
    elif bestResponse < -0.5 and bestAverageResponse < -0.2:
        return 'mPR'
    elif bestResponse < 0.35 and bestAverageResponse < 0.3:
        return 'mSD'
    else:
        return 'mPD'

def main():
    mouseIDs, tumorVolumes = readInput()
    numDays = len(tumorVolumes)
    numMice = len(tumorVolumes[0])
    print('Number of mice: %d' % numMice)

    bestResponses = [math.inf] * numMice
    bestAvgResponses = [math.inf] * numMice
    sumDeltaV = [0] * numMice

    # Looping over all mice to calculate the best response and the best average
    # response.
    # Formulas are defined in the Statistical analysis section of the
    # publication.
    for mouseIndex in range(numMice):
        startingTumorVolume = tumorVolumes[0][mouseIndex]
        seenValues = 0
        for day in range(1, numDays):
            volume = tumorVolumes[day][mouseIndex]
            if volume < 0:
                continue
            seenValues += 1
            deltaV = (volume - startingTumorVolume) / startingTumorVolume
            sumDeltaV[mouseIndex] += deltaV
            if day < FIRST_DAY:
                continue
            bestResponses[mouseIndex] = min(bestResponses[mouseIndex], deltaV)
            avgDeltaV = sumDeltaV[mouseIndex] / seenValues
            bestAvgResponses[mouseIndex] = min(bestAvgResponses[mouseIndex], avgDeltaV)

    # Output.
    print("Mouse ID\tBest response\tBest avg response\tmRECIST")
    for mouseIndex in range(numMice):
        bestResponse = bestResponses[mouseIndex]
        bestAvgResponse = bestAvgResponses[mouseIndex]
        print(
            "{}\t\t{:.2f}\t\t{:.2f}\t\t\t{}".format(
            mouseIDs[mouseIndex],
            bestResponse * 100,
            bestAvgResponse * 100,
            calcResponseCall(bestResponse, bestAvgResponse)))

main()