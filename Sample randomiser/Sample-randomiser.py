"""
This script takes an input number of reagents (e.g. 1-6x valency) and dosages
of each (e.g. 8 dosages, serial four-fold dilution of 115 nM scFv) and randomises
their position on the plate to generate an Excel/ CSV file. The script is run
repeatedly to generate new plate layouts, and these layouts are manually moved
to a sensible location for access by the data analysis script, and renamed
(numbered).
"""

from contextlib import contextmanager
import os
import sys
import csv
import random

@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)

def get_input():
    # This function usually gets user input and unpacks it into variables for the
    # rest of the script. For the sake of clarity, I've hard-coded the user
    # inputs that were given to generate the layouts used in Fryer et al. 2021.
    while True:
        output_folder = ''
        # output_folder = str(input("Where would you like the output files from "
        #                          "this script to be written to? This will "
        #                          "default to the script's run directory if "
        #                          "left blank."))
        if output_folder == '':
            output_folder = os.getcwd()
        if os.path.exists(output_folder):
            print("Writing output files to " + output_folder + ".")
            break
        else:
            print("Sorry, I can't find the directory " + output_folder + ","
                    "please try again.")
            continue
    while True:
        try:
            # reagent_num = int(input("How many different reagents do you "
            #                        "have? These will form the columns of "
            #                        "your plates (1-12 per plate)."))
            reagent_num = 6
        except TypeError:
            print("That's not an integer, try again!")
            continue
        break
    while True:
        try:
            # dosage_num = int(input("How many different dosages of each "
            #                      "reagent do you have? These will form the"
            #                      "rows of your plates (A-H per plate)."))
            dosage_num = 8
        except TypeError:
            print("That's not an integer, try again!")
            continue
        break
    while True:
        lookup_folder = ''
        # lookup_folder = str(input("If you like, you can enter the PATH "
        #                         "for a folder containing"
        #                         ".csv files which convert the reagent "
        #                         "and dosage numbers into actual values."
        #                         "These files should be named 'reagents.csv' "
        #                         "and 'dosages.csv', and should be simple lists "
        #                         "of a single conversion value per row. "
        #                         "If conversion is not desired, simply hit "
        #                         "Enter, or to use the scripts' directory "
        #                         "type '.': "))
        lookup_folder = '.'
        if lookup_folder == '':
            lookup_reagents = ''
            lookup_dosages = ''
            break
        if lookup_folder == '.':
            lookup_folder = os.getcwd()
        if os.path.exists(lookup_folder):
            if os.path.isfile(os.path.join(lookup_folder, "reagents.csv")):
                lookup_reagents = os.path.join(lookup_folder,
                                                 "reagents.csv")
            else:
                print("I can't find the file 'reagents.csv' in " +
                      lookup_folder + ", please try again.")
                continue
            if os.path.isfile(os.path.join(lookup_folder, "dosages.csv")):
                lookup_dosages = os.path.join(lookup_folder,
                                                 "dosages.csv")
            else:
                print("I can't find the file 'dosages.csv' in " +
                      lookup_folder + ", please try again.")
                continue
            break
        else:
            print("I can't find the PATH: " + lookup_folder + ", please try "
                                                              "again")
            continue
    return reagent_num, dosage_num, lookup_folder, lookup_reagents,\
           lookup_dosages, output_folder

def make_lists(nreagents, ndosages, lfolder, lreagents, ldosages):
    # This function constructs two default lists of numbers, which are
    # converted into text lists if lookup tables have been provided.
    reagents_list = []
    dosages_list = []
    for i in range(nreagents):
        reagents_list.append(i)
    for i in range(ndosages):
        dosages_list.append(i)
    if lfolder != '':
        with cd(lfolder):
            with open(lreagents, newline='') as c:
                reagents_csv = csv.reader(c, delimiter = ',', quotechar ='"')
                for rown, value in enumerate(reagents_csv):
                    reagents_list[rown] = str(value[0])
            with open(ldosages, newline='') as c:
                dosages_csv = csv.reader(c, delimiter = ',', quotechar ='"')
                for rown, value in enumerate(dosages_csv):
                    dosages_list[rown] = str(value[0])
    return reagents_list, dosages_list



def randomise_lists_to_arrays(rlist, dlist, out_folder):
    # This is the core of the script, which randomises the list of reagents
    # first, then constructs an array of randomised dosage orders for each
    # reagent. It doesn't yet handle more than 12 constructs sensibly
    # (e.g. by splitting into a list of lists (= plates) of size 12).
    csv_array = []
    random.shuffle(rlist)
    for rnum, reagent in enumerate(rlist):
        csv_array.append([reagent])
        random.shuffle(dlist)
        for dosage in dlist:
            csv_array[rnum].append(dosage)
    print(csv_array)
    # Then transpose the list of lists for easy csv writing
    csv_array_t = [list(i) for i in zip(*csv_array)]
    print(csv_array_t)
    with cd(out_folder):
        with open('shuffle_output.csv', 'w', newline='') as csv_file:
            writer = csv.writer(csv_file, delimiter = ',', quotechar = '"')
            writer.writerows(csv_array_t)


def main():
    nreagents, ndosages, lfolder, lreagents, ldosages, ofolder = get_input()
    rlist, dlist = make_lists(nreagents, ndosages, lfolder, lreagents, ldosages)
    randomise_lists_to_arrays(rlist, dlist, ofolder)

if __name__ == "__main__":
    main()