#!/usr/bin/env python

'''
Made by Ethan Gyllenhaal (egyllenhaal@unm.edu)
Last updated 21 November 2020

Script for converting haploid SNAPP input to diploid input.
Quick and dirty script, uses internal values for inputs and outputs instead of command line arguments.
'''

# Generic function for reading tab delimited file in as a table and making a list of names
# Input is the name of the text/tsv file
def getData(name):
    import csv
    data=[]
    namelist=[]
    # opens the file fed into function
    with open(name) as text:
        # makes csv reader for file
        database = csv.reader(text, delimiter=' ')
        # for each line, append the data (in this case 0/1 entries) to one array and the names to another
        for line in database:    
            data.append(list(line[2]))
            namelist.append(line[0])
    # returns two lists
    return data, namelist


# Main function for driving all operations and running the basic math to make "diploid" individuals    
def main():
    import csv
    # gets the data and names using getData function
    data, namelist = getData('snapp_dip.txt')
    # sets a length variable for the data, i.e. number of loci
    length = len(data[0])
    # initializes arrays for names and output data array
    array = []
    names = []
    # for each entry, get the name and add the "diploid" value for a given locus to the output array
    for i in range(0,len(data)//2):
        # a number for finding the data entries
        num = 2*i
        # temporary array for storing newly converted diploid data
        temp = []
        # appends the sample name to the name list
        names.append(namelist[num])
        # for each locus, add the two values together (i.e. 0/1 --> 0/1/2 format)
        for j in range(0,length):
            temp.append(int(data[num][j])+int(data[num+1][j]))
        # append the temporary array to the output array
        array.append(temp)
    output = open('snapp_out.txt', 'w')
    # puts together output, with some formatting
    for row in range(0,len(array)):
        output.write(names[row]+' '+''.join(str(e) for e in array[row])+'\n')
            
if __name__ == '__main__':
    main()
