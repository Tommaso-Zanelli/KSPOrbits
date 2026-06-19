#!/bin/bash

# =============================================================
# Variables

# Define seatch path
dirPath="./OctaveScripts"

# Build full lists of files
allmFilesList=$(echo "$(find ${dirPath} -type f -name "*.m")")
allmatFilesList=$(echo "$(find ${dirPath} -type f -name "*.mat")")

# =============================================================
# Functions

# Check if element is already in the list
isNotInList()
{
    checkElement=$1
    list=$2
    notInList=True
    for element in ${list}; do
        if [[ "${element}" == "${checkElement}" ]]; then
            notInList=False
        fi
    done
    echo ${notInList}
}

# Add entry to a list if not present already
addToListIfNotInAlready()
{
    newElement=$1
    oldList=$2
    if [[ "$(isNotInList ${newElement} "${oldList}")" == "True" ]]; then
        echo ${oldList} ${newElement}
    else
        echo ${oldList}
    fi
}

# Check file dependencies: for each function check if it is called in the file
checkDependencies()
{
    inputFile=$1
    funList=$2
    depFileList=""
    # Loop over all functions
    for fun in ${funList}; do
        funcName=$(echo $fun | cut -f 1 -d'@')
        fileName=$(echo $fun | cut -f 2 -d'@')
        if [[ ! -z $(grep -w ${funcName} ${inputFile}) ]]; then
            depFileList=$(addToListIfNotInAlready "${fileName}" "${depFileList}")
        fi
    done
    echo "${depFileList}"
}

# Check recursively for file dependencies
checkRecursive()
{
    inputFile=$1
    funList=$2
    oldList=""
    newList="$(checkDependencies ${inputFile} "${funList}")"
    while [[ ! "${oldList}" == "${newList}" ]]; do
        oldList="${newList}"
        for element in ${oldList}; do
            addList=$(checkDependencies ${element} "${funList}")
            for newElement in ${addList}; do
                if [[ "$(isNotInList ${newElement} "${newList}")" == "True" ]]; then
                    newList="${newList} ${newElement}"
                fi
            done
        done
    done
    echo ${oldList}
}

# =============================================================
# Main

# Get full list of standalone matlab functions
echo "Getting full list of standalone matlab functions..."
time {
    stdAlFunList="";
    searchRegex="function[ ]{1,}[^=]{1,}=[ ]{1,}([^\(]{1,})\([^\)]*\)"
    while read line; do 
        fileName=$(basename $(echo $line | cut -f 1 -d':') | cut -f 1 -d '.')
        funcName=$(echo $line | cut -f 2 -d':')
        if [[ "${fileName}" == "${funcName}" ]]; then
            stdAlFunList="${stdAlFunList} ${funcName}@$(echo $line | cut -f 1 -d':')"
        fi
    done < <(grep -rioE "${searchRegex}" ${dirPath} | sed -r "s/${searchRegex}/\1/g")
}
echo; echo

# Get list of essential .m files
echo "Getting list of essential .m files..."
time essentialmFiles=$(checkRecursive "./orbits.py" "${stdAlFunList}")
echo; echo

# Get list of non-essential .m files
echo "Getting list of non-essential .m files..."
time {
    uselessmFilesList=""
    for file in ${allmFilesList}; do
        if [[ "$(isNotInList ${file} "${essentialmFiles}")" == "True" ]]; then
            uselessmFilesList="${uselessmFilesList} ${file}"
        fi    
    done
}
echo; echo

# Get list of essential and non-essential .mat files
echo "Getting list of essential and non-essential .mat files..."
time {
    essentialmatFiles=""
    uselessmatFilesList=""
    for matfile in ${allmatFilesList}; do
        isNotInUse=True
        for mfile in ${essentialmFiles}; do
            if [[ ! -z $(grep -w "${matfile}" ${mfile}) ]]; then
                isNotInUse=False
            fi
        done
        if [[ "${isNotInUse}" == "True" ]]; then
            uselessmatFilesList="${uselessmatFilesList} ${matfile}"
        else
            essentialmatFiles="${essentialmatFiles} ${matfile}"
        fi
    done
}
echo; echo

for mFileToAdd in ${essentialmFiles}; do
    echo "git add ${mFileToAdd}"
done
for matFileToAdd in ${essentialmatFiles}; do
    echo "git add ${matFileToAdd}"
done
echo
for mFileToDel in ${uselessmFilesList}; do
    echo "rm -rfv ${mFileToDel}"
done
for matFileToDel in ${uselessmatFilesList}; do
    echo "rm -rfv ${matFileToDel}"
done
echo