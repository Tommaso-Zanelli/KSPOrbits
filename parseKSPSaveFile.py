#!/usr/bin/env python3
# import cProfile, pstats, io
# import pickle
import pyksp as pk


inputFile = r"C:\Users\t.zanelli\Downloads\TEMPDELETEME\all-good.sfs"
outputFile = r"C:\Users\t.zanelli\Downloads\TEMPDELETEME\all-good-check.sfs"

# Create timer object
tt = pk.timer()


def main():
    # Load file as dictionary
    tt.tic()
    dataDict = pk.sfsFileToDict(inputFile)
    #fileId = open(r"C:\Users\t.zanelli\Downloads\TEMPDELETEME\all-good.pkl", 'wb')
    #pickle.dump(dataDict, fileId)
    #fileId.close()
    tt.toc("Time elapsed loading file (s):")

    # Output file for checking functionality
    tt.tic()
    pk.dictTosfsFile(dataDict, outputFile)
    tt.toc("Time elapsed writing file (s):")

# pr = cProfile.Profile()
# pr.enable()
main()
# pr.disable()

# s = io.StringIO()
# ps = pstats.Stats(pr, stream=s).sort_stats('cumulative')
# ps.print_stats()  # top 20 lines
# print(s.getvalue())