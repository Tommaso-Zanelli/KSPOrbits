#!/usr/bin/env python3
import cProfile, pstats, io
import pickle

# Create timer object
import time
class timer:
    def __init__(self):
        self.t1 = 0.0
        self.t2 = 0.0
        self.dt = 0.0
    def tic(self):
        self.t1 = time.time()
    def toc(self, outStr = None):
        self.t2 = time.time()
        self.dt = self.t2 - self.t1
        if outStr != None:
            print(outStr, self.dt)
        return self.dt
tt = timer()

def remIndent(line):
    stripped = line.lstrip('\t')
    return stripped, len(line) - len(stripped)

def parseLines(lines, verbose = False):
    numLines = len(lines)
    blockLevel = 0
    if verbose:
        listsOfPropVals = 0
    
    # Dictionary accumulators initialization
    dataDict = {}
    currentDict = dataDict
    previousDict = []
    
    # Loop over each line
    for ii, line in enumerate(lines):
        # Remove indentation and store its level
        lineStrip, indentLevel = remIndent(line)
        
        # If line starts a sub-block
        if lineStrip == '{':
            blockLevel += 1
            
            # If a block with that name is already present in current block
            if blockName in currentDict:

                if not isinstance(currentDict[blockName], list):
                    currentDict[blockName] = [currentDict[blockName].copy()]
                currentDict[blockName].append({})
                previousDict.append(currentDict)
                currentDict = currentDict[blockName][-1]
            
            # If the block is currently unique
            else:
                currentDict[blockName] = {}
                previousDict.append(currentDict)
                currentDict = currentDict[blockName]
        
        # If line ends a sub-block
        elif lineStrip == '}':
            blockLevel -= 1
            currentDict = previousDict.pop()
        
        # For every other line
        else:
            # Verbose checks
            if verbose:
                # Check nesting level
                if blockLevel != len(previousDict):
                    print(f"WARNING: dictionary nesting level ({len(previousDict):d}) does not match block nesting level ({blockLevel:d})")
                # Check indentation level
                if blockLevel != indentLevel:
                    print(f"WARNING: indentation level ({indentLevel:d}) does not match block nesting level ({blockLevel:d})")
            
            # Check if line is block name
            if ii + 1 < numLines and lines[ii + 1].lstrip('\t') == '{':
                blockName = lineStrip

            # Check if line is property/value pair
            elif '=' in lineStrip:
                
                # Split property and value
                propVal = lineStrip.split('=', 1)
                prop = propVal[0]
                if prop[-1] == ' ':
                    prop = prop[:-1]
                val  = propVal[1]
                if val[0] == ' ':
                    val = val[1:]
                
                # If a property with that name is already present in current block
                if prop in currentDict:
                    if not isinstance(currentDict[prop], list):
                        currentDict[prop] = [currentDict[prop]]
                        if verbose:
                            listsOfPropVals += 1
                    currentDict[prop].append(val)
                
                # If the property is currently unique
                else:
                    currentDict[prop] = val

            # Check if line is trailing empty line
            elif len(lineStrip) == 0:
                pass

            # Line not recognized
            else:
                if verbose:
                    print(f"WARNING: Unclassifiable line n. {ii + 1:d} of {numLines:d}: \"{line:s}\"")

    if verbose:
        if listsOfPropVals > 0:
            print(f"WARNING: number of lists of property/value pairs with same property name: {listsOfPropVals:d}")

    return dataDict

def dictToLines(dataDict, verbose = False):
    lines = []
    blockLevel = 0
    
    # Dictionary accumulators initialization
    currentDict = dataDict
    previousDict = []
    previousIndex = []
    previousListIndex = []
    previousWasList = []
    
    iterate = True
    inList = False
    ii = 0
    jj = 0
    dictLen  = len(dataDict.keys())
    dictKeys = list(dataDict.keys())
    while iterate:

        # Check if the current dictionary is over
        if ii >= dictLen:
            
            if blockLevel > 0:
                currentDict = previousDict.pop()
                dictLen = len(currentDict.keys())
                dictKeys = list(currentDict.keys())
                jj = previousListIndex.pop()
                inList = previousWasList.pop()
                ii = previousIndex.pop()
                if inList:
                    jj += 1
                else:
                    ii += 1
                    jj = 0
                blockLevel -= 1
                lines.append(f"{'\t'*(blockLevel):s}{'}':s}")
                if verbose:
                    print(lines[-1])
                
            else:
                iterate = False
        
        else:
            if isinstance(currentDict[dictKeys[ii]], dict):
                lines.append(f"{'\t'*(blockLevel):s}{dictKeys[ii]:s}")
                if verbose:
                    print(lines[-1])
                lines.append(f"{'\t'*(blockLevel):s}{'{':s}")
                if verbose:
                    print(lines[-1])
                previousDict.append(currentDict)
                previousIndex.append(ii)
                previousListIndex.append(jj)
                previousWasList.append(False)
                currentDict = currentDict[dictKeys[ii]]
                ii = 0
                blockLevel += 1
                dictLen = len(currentDict.keys())
                dictKeys = list(currentDict.keys())
            elif isinstance(currentDict[dictKeys[ii]], list):
                if jj < len(currentDict[dictKeys[ii]]):
                    if isinstance(currentDict[dictKeys[ii]][jj], dict):
                        lines.append(f"{'\t'*(blockLevel):s}{dictKeys[ii]:s}")
                        if verbose:
                            print(lines[-1])
                        lines.append(f"{'\t'*(blockLevel):s}{'{':s}")
                        if verbose:
                            print(lines[-1])
                        previousDict.append(currentDict)
                        previousIndex.append(ii)
                        previousListIndex.append(jj)
                        previousWasList.append(True)
                        currentDict = currentDict[dictKeys[ii]][jj]
                        ii = 0
                        blockLevel += 1
                        dictLen = len(currentDict.keys())
                        dictKeys = list(currentDict.keys())
                    else:
                        lines.append(f"{'\t'*(blockLevel):s}{dictKeys[ii]:s} = {currentDict[dictKeys[ii]][jj]:s}")
                        if verbose:
                            print(lines[-1])
                        jj += 1
                else:
                    ii += 1
                    jj = 0
            else:
                lines.append(f"{'\t'*(blockLevel):s}{dictKeys[ii]:s} = {currentDict[dictKeys[ii]]:s}")
                if verbose:
                    print(lines[-1])
                ii += 1
                jj = 0

    lines.append('')
    return lines

def sfsFileToDict(fileName, verbose = False):
    fileId = open(fileName, 'rb')
    lines  = fileId.read().decode('utf-8').split("\r\n")
    fileId.close()
    return parseLines(lines, verbose)

def dictTosfsFile(dataDict, fileName, verbose = False):
    lines  = dictToLines(dataDict, verbose)
    fileId = open(fileName, 'wb')
    fileId.write('\r\n'.join(lines).encode('utf-8'))
    fileId.close()

def main():
    # Load file as dictionary
    tt.tic()
    dataDict = sfsFileToDict(r"C:\Users\t.zanelli\Downloads\TEMPDELETEME\all-good.sfs")
    fileId = open(r"C:\Users\t.zanelli\Downloads\TEMPDELETEME\all-good.pkl", 'wb')
    pickle.dump(dataDict, fileId)
    fileId.close()
    tt.toc("Time elapsed loading file (s):")
    
    # Output file for checking functionality
    tt.tic()
    dictTosfsFile(dataDict, r"C:\Users\t.zanelli\Downloads\TEMPDELETEME\all-good-check.sfs")
    tt.toc("Time elapsed writing file (s):")

# pr = cProfile.Profile()
# pr.enable()
main()
# pr.disable()

# s = io.StringIO()
# ps = pstats.Stats(pr, stream=s).sort_stats('cumulative')
# ps.print_stats()  # top 20 lines
# print(s.getvalue())